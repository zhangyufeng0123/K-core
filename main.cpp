#include<bits/stdc++.h>

using namespace std;

const string filename = "data1.txt";

struct Vertex{
    int val;    //这个点的编号
    vector<int> connectV;   //与这个点连接的点
    //int degrees; //该点的度数
};

//the whole graph
struct Graph{
    //vector<int> vertex; //标注顺序下是哪个点的位置
    vector<Vertex> nodes;
    map<int, int> location; //用来存储一个点在第几列
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

            //当i=0时代表输入num的是这个节点，其余的点是这个点一步就可以到达的点
            if(i == 0){
                temp.val = num;
            }else{
                temp.connectV.push_back(num);
            }
        }
        G_out.nodes.push_back(temp);
        //G_out.vertex.push_back(temp.val);
        G_out.location[temp.val] = G_out.location.size();
        //G_out.nodes.back().degrees = G_out.nodes.back().connectV.size();
    }

    //存储入度
    Vertex tmp;
    //初始化入度
    for(int i = 0; i < row; i++){
        tmp.val = G_out.nodes[i].val;
        G_in.nodes.push_back(tmp);
        G_in.location[tmp.val] = i;
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

//对点的出度进行排序
static bool cmp_out(Vertex &a, Vertex &b){
    return a.connectV.size() < b.connectV.size();
}

//对点的入度进行排序
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

/*
 * 删除图2中
 * 先找到第一张图要删的点的位置
 * 将这个点的位置与最后一个点的位置交换
 */
void DeleteNode(Graph &G1, Graph &G2, int targetNode){
    int point = G1.location[targetNode];

    //删除所有这个点所有出边或入边
    for(int i = 0; i < G1.nodes[point].connectV.size(); i++){
        int temppoint = G1.nodes[point].connectV[i];
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
//k是出度，l是入度
bool QueryDcore(Graph G_out, Graph G_in, int k, int l, int q, Graph &G_outq, Graph &G_inq){

    int flag;
    while(1){
        flag = 1;
        //先将点按照出度进行排序
        sort(G_out.nodes.begin(), G_out.nodes.end(), cmp_out);
        for(int i = 0; i < G_out.nodes.size(); i++){
            G_out.location[G_out.nodes[i].val] = i;
        }
        for(int i = 0; i < G_out.nodes.size(); i++){
            if(G_out.nodes[i].connectV.size() < k){
                //delete node
                int e = G_out.nodes[i].val;
                DeleteNode(G_out, G_in, e);
                DeleteNode(G_in, G_out, e);
                flag = 0;
            }
        }

        //将点按照入度进行排序
        sort(G_in.nodes.begin(), G_in.nodes.end(), cmp_in);
        for(int i = 0; i < G_in.nodes.size(); i++){
            G_in.location[G_in.nodes[i].val] = i;
        }
        for(int i = 0; i < G_in.nodes.size(); i++){
            if(G_in.nodes[i].connectV.size() < l){
                //delete node
                int e = G_in.nodes[i].val;
                DeleteNode(G_out, G_in, e);
                DeleteNode(G_in, G_out, e);
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

void Peel(Graph G_out, Graph G_in, int k, int l, int q){
    Graph G_outq, G_inq;
    bool flag = QueryDcore(G_out, G_in, k, l, q, G_outq, G_inq);
    if(!flag){
        cout << "目标点不能构成k-core" << endl;
        return;
    }
    
}

int main(){
    Graph G_out, G_in;
    bool j = readData(G_out, G_in);
    //cout << 1 << endl;

    //Graph temp = QueryDcore(G_out, G_in, 1, 1, 1);
    Graph G_outq, G_inq;
    bool t = QueryDcore(G_out, G_in, 1, 1, 8, G_outq, G_inq);
    cout << 1 << endl;

    return 0;
}