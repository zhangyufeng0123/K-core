#include<bits/stdc++.h>

using namespace std;

const string filename = "data1.txt";

struct Vertex{
    int val;    //这个点的编号
    vector<int> connectV;   //与这个点连接的点
    int weight; //该点的权重
    //int degrees; //该点的度数
};

//the whole graph
struct Graph{
    //vector<int> vertex; //标注顺序下是哪个点的位置
    vector<Vertex> nodes;
    map<int, int> location; //用来存储一个点在第几列
    map<int, int> weights;
};

//将string类型转换为int类型
int stringToInt(string str){
    int res = 0;
    for(int i = 0; i < str.size(); i++){
        if(str[i] < '0' || str[i] > '9')    continue;
        res = res * 10 + int(str[i] - '0');
    }
    return res;
}

//判断输入的数字是否合法
bool judgeNum(string str){
    for(int i = 0; i < str.size(); i++){
        if(str[i] < '0' || str[i] > '9')    return false;
    }
    return true;
}

//读取数据
bool readData(Graph &G_out, Graph &G_in){
    ifstream infile;
    infile.open(filename);
    if(!infile){
        cout << "打开文件出错！" << endl;
    }
    string lineStr;
    int row = 0;    //记录现在是第几行

    //存储出度
    while(getline(infile, lineStr)){
        //int col = 0;    //记录现在是第几列
        row++;
        //存成二维表结构
        stringstream ss(lineStr);
        string str;
        vector<string> lineArray;
        //按照空格分开
        while(getline(ss, str, ' ')){
            lineArray.push_back(str);
        }

        //if(lineArray.size() == 0)   continue;   //只有一个结点，说明没有出度，可以直接忽略
        //将数据放到vector容器中
        Vertex temp;
        for(int i = 0; i < lineArray.size(); i++){
            /*
            if(!judgeNum(lineArray[i])){//判断输入的数据是否报错
                cout << "第" << row << "行，第" << i << "个数据出现了问题" << endl;
                return false;
            }
             */
            int num = stringToInt(lineArray[i]);
            //int tmpWeight;
            //当i=0时代表输入num的是这个节点，其余的点是这个点一步就可以到达的点
            if(i == 0){
                temp.val = num;
            }else if(i == 1){
                temp.weight = num;
                //G_out.weight[temp.val] = num;//赋值权重
                //tmpWeight = num;
            }else{
                temp.connectV.push_back(num);
            }
        }
        G_out.nodes.push_back(temp);
        //G_out.vertex.push_back(temp.val);
        G_out.location[temp.val] = G_out.location.size();
        G_out.weights[temp.val] = temp.weight;
        //G_out.nodes.back().degrees = G_out.nodes.back().connectV.size();
    }

    //存储入度
    Vertex tmp;
    //初始化入度
    for(int i = 0; i < row; i++){
        tmp.val = G_out.nodes[i].val;
        tmp.weight = G_out.weights[tmp.val];
        G_in.nodes.push_back(tmp);
        G_in.location[tmp.val] = i;
        G_in.weights[tmp.val] = tmp.weight;
        //G_in.vertex.push_back(tmp.val);
    }

    for(int i = 0; i < row; i++){
        for(int j = 0; j < G_out.nodes[i].connectV.size(); j++){
            int point = G_out.nodes[i].connectV[j];
            int locate = G_in.location[point];
            G_in.nodes[locate].connectV.push_back(G_out.nodes[i].val);
        }
    }
}

//二分查找
int binarySearch(vector<int> &nums, int target){
    int left = 0, right = nums.size() - 1;
    while(right >= left){
        int mid = (left + right) >> 1;
        if(nums[mid] == target){
            return mid;
        }else if(nums[mid] > target){
            right = mid - 1;
        }else{
            left = mid + 1;
        }
    }
    return -1;
}

//根据点的权重进行升序
static bool cmp_weight(Vertex &a, Vertex&b){
    return a.weight < b.weight;
}

//对点的出度进行升序
static bool cmp_out(Vertex &a, Vertex &b){
    return a.connectV.size() < b.connectV.size();
}

//对点的入度进行升序
static bool cmp_in(Vertex &a, Vertex &b){
    return a.connectV.size() < b.connectV.size();
}

void DeleteSingleNode(Graph &G, int nodeVal, int targetNode){
    int point = G.location[targetNode];
    int i = 0;
    for(; i < G.nodes[point].connectV.size(); i++){
        if(G.nodes[point].connectV[i] == nodeVal)   break;
    }

    swap(G.nodes[point].connectV[i], G.nodes[point].connectV.back());
    G.nodes[point].connectV.pop_back();
}


//根据点集获取图
void pointsChangeToGraph(Graph G_out, Graph G_in, Graph &G_outr, Graph &G_inr, vector<int> points){
    map<int, bool> judge;
    //初始化judge
    for(int i = 0; i < G_out.nodes.size(); i++){
        judge[G_out.nodes[i].val] = false;
    }

    //将目标点集所拥有的点置为true
    for(int i = 0; i < points.size(); i++){
        judge[points[i]] = true;
    }

    for(int i = 0; i < points.size(); i++){
        Vertex tmp;
        tmp.val = points[i];
        int locate = G_out.location[tmp.val];
        tmp.weight = G_out.nodes[locate].weight;
        vector<int> connectV_tmp;
        for(int j = 0; j < G_out.nodes[locate].connectV.size(); j++){
            if(judge[G_out.nodes[locate].connectV[j]]){
                connectV_tmp.push_back(G_out.nodes[locate].connectV[j]);
            }
        }
        G_outr.nodes.push_back(tmp);
        G_outr.location[tmp.val] = i;
        G_outr.weights[tmp.val] = tmp.weight;
        connectV_tmp.clear();
        for(int j = 0; j < G_in.nodes[locate].connectV.size(); j++){
            if(judge[G_in.nodes[locate].connectV[j]]){
                connectV_tmp.push_back(G_in.nodes[locate].connectV[j]);
            }
        }
        G_inr.nodes.push_back(tmp);
        G_inr.location[tmp.val] = i;
        G_inr.weights[tmp.val] = tmp.weight;
    }
}

/*
 * 删除图2中
 * 先找到第一张图要删的点的位置
 * 将这个点的位置与最后一个点的位置交换
 */
void DeleteNode(Graph &G1, Graph &G2, int targetNode, vector<int> &nodes){
    int point = G1.location[targetNode];

    //删除所有这个点所有出边或入边
    for(int i = 0; i < G1.nodes[point].connectV.size(); i++){
        int temppoint = G1.nodes[point].connectV[i];
        nodes.push_back(temppoint);
        DeleteSingleNode(G2, targetNode, temppoint);
    }

    //获取nodes最后一个点的信息
    int lastNodeVal = G1.nodes.back().val;
    int lastNodeLocate = G1.location[lastNodeVal];
    swap(G1.nodes[point], G1.nodes[lastNodeLocate]);
    G1.location[lastNodeVal] = point;
    G1.nodes.pop_back();
    G1.location[targetNode] = -1;
}

//Algorithm 1:QueryDcore
//k是入度，l是出度
bool QueryDcore(Graph G_out, Graph G_in, int k, int l, int q, Graph &G_outq, Graph &G_inq){

    int flag;
    vector<int> useless;
    while(1){
        flag = 1;
        //先将点按照出度进行排序
        sort(G_out.nodes.begin(), G_out.nodes.end(), cmp_out);
        for(int i = 0; i < G_out.nodes.size(); i++){
            G_out.location[G_out.nodes[i].val] = i;
        }
        for(int i = 0; i < G_out.nodes.size(); i++){
            if(G_out.nodes[i].connectV.size() < l){
                //delete node
                int e = G_out.nodes[i].val;
                DeleteNode(G_out, G_in, e, useless);
                DeleteNode(G_in, G_out, e, useless);
                flag = 0;
            }
        }

        //将点按照入度进行排序
        sort(G_in.nodes.begin(), G_in.nodes.end(), cmp_in);
        for(int i = 0; i < G_in.nodes.size(); i++){
            G_in.location[G_in.nodes[i].val] = i;
        }
        for(int i = 0; i < G_in.nodes.size(); i++){
            if(G_in.nodes[i].connectV.size() < k){
                //delete node
                int e = G_in.nodes[i].val;
                DeleteNode(G_out, G_in, e, useless);
                DeleteNode(G_in, G_out, e, useless);
                flag = 0;
            }
        }
        if(flag)    break;
    }

    if(G_out.location[q] == -1) return false;
    vector<int> arrVal;
    queue<int> que;
    map<int, int> stamp;
    for(int i = 0; i < G_out.nodes.size(); i++){
        stamp[G_out.nodes[i].val] = 0;
    }
    stamp[q] = 1;
    que.push(q);
    arrVal.push_back(q);
    while(!que.empty()){
        int tmp = que.front();
        que.pop();
        //在出度的图中找相连的点
        int locate = G_out.location[tmp];
        for(int i = 0; i < G_out.nodes[locate].connectV.size(); i++){
            int nodeVal = G_out.nodes[locate].connectV[i];
            if(stamp[nodeVal] == 0){
                que.push(nodeVal);
                arrVal.push_back(nodeVal);
                stamp[nodeVal] = 1;
            }
        }

        //在入度的图中找出相连的点
        locate = G_in.location[tmp];
        for(int i = 0; i < G_in.nodes[locate].connectV.size(); i++){
            int nodeVal = G_in.nodes[locate].connectV[i];
            if(stamp[nodeVal] == 0){
                que.push(nodeVal);
                arrVal.push_back(nodeVal);
                stamp[nodeVal] = 1;
            }
        }
    }

    for(int i = 0; i < arrVal.size(); i++){
        Vertex tmp;
        tmp.val = arrVal[i];
        int locate = G_out.location[tmp.val];
        for(int j = 0; j < G_out.nodes[locate].connectV.size(); j++){
            int t = G_out.nodes[locate].connectV[j];
            if(stamp[t]){
                tmp.connectV.push_back(t);
            }
        }
        G_outq.nodes.push_back(tmp);
        G_outq.location[tmp.val] = i;
    }

    for(int i = 0; i < arrVal.size(); i++){
        Vertex tmp;
        tmp.val = arrVal[i];
        int locate = G_in.location[tmp.val];
        for(int j = 0; j < G_in.nodes[locate].connectV.size(); j++){
            int t = G_in.nodes[locate].connectV[j];
            if(stamp[t]){
                tmp.connectV.push_back(t);
            }
        }
        G_inq.nodes.push_back(tmp);
        G_inq.location[tmp.val] = i;
    }
    return true;
}

//k是入度，l是出度
void Peel(Graph G_out, Graph G_in, int k, int l, int q, Graph &G_outr, Graph &G_inr){
    Graph G_outq, G_inq;
    //获得极大子图
    bool flag = QueryDcore(G_out, G_in, k, l, q, G_outq, G_inq);
    if(!flag){
        cout << "目标点不能构成k-core" << endl;
        return;
    }
    //对极大子图中的点进行排序
    sort(G_outq.nodes.begin(), G_outq.nodes.end(), cmp_weight);
    queue<int> Q;
    vector<int> D, L;
    //用一个hash来记录L中的每个点的位置
    //map<int, int> locate;
    //L.push_back(G_outq.nodes.front().val);

    //初始化L
    for(int i = 0; i < G_outq.nodes.size(); i++){
        L.push_back(G_outq.nodes[i].val);
        //locate[L.back()] = i;
    }

    map<int, bool> judge;
    for(int i = 0; i < G_outq.nodes.size(); i++){
        judge[G_outq.nodes[i].val] = false;
        //judge_tmp[G_outq.nodes[i].val] = false;//标记还有哪些点还剩余
    }
    //获取删减之后的顶点
    while (!L.empty()){
        int tmp = L.front();
        D.clear();
        Q.push(tmp);
        judge[tmp] = true;
        while(!Q.empty()){
            int pointVal = Q.front();
            judge[pointVal] = false;
            Q.pop();
            if(pointVal == q){
//                int res = G_out.nodes[G_out.location[L.back()]].weight;
//                for(int i = 0; i < D.size(); i++){
//                    res = max(res, G_out.nodes[G_out.location[D[i]]].weight);
//                }
//                return res;
                while(!L.empty()){
                    int tmpe = L.back();
                    D.push_back(tmpe);
                    L.pop_back();
                }
                break;
            }
            vector<int> nodes;//记录哪些与pointVal相连的边被删除
            DeleteNode(G_outq, G_inq, pointVal, nodes);
            nodes.clear();
            DeleteNode(G_inq, G_outq, pointVal, nodes);
            for(int i = 0; i < nodes.size(); i++){
                if(G_outq.nodes[G_outq.location[nodes[i]]].connectV.size() < l && !judge[G_outq.location[nodes[i]]]){
                    Q.push(nodes[i]);
                    judge[G_outq.location[nodes[i]]] = true;
                    continue;
                }
                if(G_inq.nodes[G_inq.location[nodes[i]]].connectV.size() < k && !judge[G_inq.location[nodes[i]]]){
                    Q.push(nodes[i]);
                    judge[G_inq.location[nodes[i]]] = true;
                }
            }
            int location = binarySearch(L, pointVal);
            L.erase(L.begin() + location);
            D.push_back(pointVal);
        }
    }
    pointsChangeToGraph(G_out, G_in, G_outr, G_inr, D);
//    return -1;
//    for(int i = 0; i < D.size(); i++){
//        judge_tmp[D[i]] = true;
//    }
//    for(int i = 0; i < D.size(); i++){
//        Vertex tmp;
//        tmp.val = D[i];
//        int locate = G_out.location[tmp.val];
//        for(int j = 0; j < G_out.nodes[locate].connectV.size(); j++){
//            int t = G_out.nodes[locate].connectV[j];
//            if(judge_tmp[t]){
//                tmp.connectV.push_back(t);
//            }
//        }
//        G_outr.nodes.push_back(tmp);
//        G_outr.location[tmp.val] = i;
//    }
//    for(int i = 0; i < D.size(); i++){
//        Vertex tmp;
//        tmp.val = D[i];
//        int locate = G_in.location[tmp.val];
//        for(int j = 0; j < G_in.nodes[locate].connectV.size(); j++){
//            int t = G_in.nodes[locate].connectV[j];
//            if(judge_tmp[t]){
//                tmp.connectV.push_back(t);
//            }
//        }
//        G_inr.nodes.push_back(tmp);
//        G_inr.location[tmp.val] = i;
//    }
}

//G_out存储的是每个点能一步到达的点，G_in存储的是每个点能被一步到达的点
void batchPeelingAlgorithm(Graph G_out, Graph G_in, int k, int l, int q, Graph &G_outr, Graph &G_inr){
    Graph G_outq, G_inq;
    //获取存在q的极大子图
    bool flag = QueryDcore(G_out, G_in, k, l, q, G_outq, G_inq);//Line1
    if(!flag){
        cout << "目标点不能构成k-core" << endl;
        return;
    }

    //对极大子图中的点根据权重进行排序
    sort(G_outq.nodes.begin(), G_outq.nodes.end(), cmp_weight);
    vector<int> L, H;//L存储图中的所有点的val
    for(int i = 0; i < G_outq.nodes.size(); i++){
        L.push_back(G_outq.nodes[i].val);   //Line3
    }
    H.assign(L.begin(), L.begin() + (L.size() / 2));//Line4
    while(!H.empty()){
        vector<int> L0(L);  //Line6
        L.assign(L.begin() + (L.size() / 2), L.end());  //Line7
        Graph G_outq_tmp, G_inq_tmp;
        flag = QueryDcore(G_outq, G_inq, k, l, q, G_outq_tmp, G_inq_tmp);//Line8
        if(!flag){
            H.assign(H.begin(), H.begin() + (H.size() / 2));
            L.assign(L0.begin(), L0.end());
        }else{
            vector<int> L_star;
            sort(G_outq_tmp.nodes.begin(), G_outq_tmp.nodes.end(), cmp_weight);
            for(int i = 0; i < G_outq_tmp.nodes.size(); i++){
                L_star.push_back(G_outq_tmp.nodes[i].val);
            }
            H.assign(L_star.begin(), L_star.begin() + (L_star.size() / 2));
            L.assign(L_star.begin(), L_star.end());
        }
    }
    pointsChangeToGraph(G_outq, G_inq, G_outr, G_inr, L);
}

int main(){
    Graph G_out, G_in;
    bool j = readData(G_out, G_in);
    //cout << 1 << endl;

    //Graph temp = QueryDcore(G_out, G_in, 1, 1, 1);
    Graph G_outr, G_inr;
    //bool t = QueryDcore(G_out, G_in, 1, 1, 8, G_outq, G_inq);
    cout << 1 << endl;

    return 0;
}