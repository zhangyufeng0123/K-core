#include<bits/stdc++.h>

using namespace std;

const string filename = "data.txt";

int Kkmax, Llmax;

typedef pair<vector<int>, vector<int>> PVV;
typedef pair<int, int> PII;

struct Vertex {
    int val;    //这个点的编号
    vector<int> connectV;   //与这个点连接的点
    int weight; //该点的权重
    //int degrees; //该点的度数
};

//the whole graph
struct Graph {
    //vector<int> vertex; //标注顺序下是哪个点的位置
    vector<Vertex> nodes;
    map<int, int> location; //用来存储一个点在第几列
    map<int, int> weights;
};

//将string类型转换为int类型
int stringToInt(string str) {
    int res = 0;
    for (int i = 0; i < str.size(); i++) {
        if (str[i] < '0' || str[i] > '9') continue;
        res = res * 10 + int(str[i] - '0');
    }
    return res;
}

//根绝权值进行排序
static bool cmp_initWeight(PII &a, PII &b) {
    if (a.second == b.second) return a.first < b.first;
    return a.second < b.second;
}

//判断输入的数字是否合法
bool judgeNum(string str) {
    for (int i = 0; i < str.size(); i++) {
        if (str[i] < '0' || str[i] > '9') return false;
    }
    return true;
}

//读取数据2
bool readDataTmp(Graph &G_out, Graph &G_in){
    ifstream infile;
    infile.open(filename);
    if(!infile){
        cout << "打开文件错误!" << endl;
        return false;
    }
    int points, edges;
    infile >> points >> edges;

    vector<PII> weights;
    vector<Vertex> initTmp(points);
    G_out.nodes.assign(initTmp.begin(), initTmp.end());
    G_in.nodes.assign(initTmp.begin(), initTmp.end());
    for(int i = 0; i < points; i++){
        //初始化权值begin
        PII tmp;
        tmp.first = i + 1;
        tmp.second = 0;
        weights.push_back(tmp);
        //初始化权值end

        //初始化图begin
        G_out.nodes[i].val = i + 1;
        G_out.location[i + 1] = i;
        G_in.nodes[i].val = i + 1;
        G_in.location[i + 1] = i;
        //初始化图end
    }

    //输入边
    for(int i = 0; i < edges; i++){
        int a, b;
        infile >> a >> b;
        weights[a - 1].second--;
        weights[b - 1].second++;
        G_out.nodes[a - 1].connectV.push_back(b);
        G_in.nodes[b - 1].connectV.push_back(a);
    }

    sort(weights.begin(), weights.end(), cmp_initWeight);
    Kkmax = 0, Llmax = 0;
    for(int i = 0; i < weights.size(); i++){
        Kkmax = max(Kkmax, int(G_in.nodes[weights[i].first - 1].connectV.size()));
        Llmax = max(Llmax, int(G_out.nodes[weights[i].first - 1].connectV.size()));

        G_out.nodes[weights[i].first - 1].weight = i;
        G_in.nodes[weights[i].first - 1].weight = i;
    }

    //fill the weights
    for(int i = 0; i < points; i++){
        G_out.weights[i + 1] = G_out.nodes[i].weight;
        G_in.weights[i + 1] = G_in.nodes[i].weight;
    }
    return true;
}

//读取数据
bool readData(Graph &G_out, Graph &G_in) {
    ifstream infile;
    infile.open(filename);
    if (!infile) {
        cout << "打开文件出错！" << endl;
    }
    string lineStr;
    int row = 0;    //记录现在是第几行
    bool flag = false;   //判断是不是第一行
    vector<PII> weights;    //存储代替权值的数值

    while (getline(infile, lineStr)) {
        if (!flag) {
            //存成二维表结构
            stringstream ss(lineStr);
            string str;
            //按照空格分开
            vector<string> lineArray;
            while (getline(ss, str, ' ')) {
                lineArray.push_back(str);
            }

            vector<Vertex> initTmp(stringToInt(lineArray[1]));
            G_out.nodes.assign(initTmp.begin(), initTmp.end());
            G_in.nodes.assign(initTmp.begin(), initTmp.end());
            for (int i = 0; i < G_out.nodes.size(); i++) {
                PII tmp;
                tmp.first = i + 1;
                tmp.second = 0;
                weights.push_back(tmp);

                G_out.nodes[i].val = i + 1;
                G_out.location[i + 1] = i;
                G_in.nodes[i].val = i + 1;
                G_in.location[i + 1] = i;
            }
            flag = true;
        } else {
            //存成二维表结构
            stringstream ss(lineStr);
            string str;
            //按照空格分开
            vector<string> lineArray;
            while (getline(ss, str, ' ')) {
                lineArray.push_back(str);
            }
            for (int i = 0; i < lineArray.size(); i++) {
                int a = stringToInt(lineArray[0]), b = stringToInt(lineArray[1]);
                weights[a - 1].second--;
                weights[b - 1].second++;
                G_out.nodes[a].connectV.push_back(b);
                G_in.nodes[b].connectV.push_back(a);
            }
        }
    }

    //计算各点的权值
    sort(weights.begin(), weights.end(), cmp_initWeight);
    Kkmax = 0;
    Llmax = 0;
    for (int i = 0; i < weights.size(); i++) {
        Kkmax = max(Kkmax, int(G_in.nodes[weights[i].first - 1].connectV.size()));
        Llmax = max(Llmax, int(G_out.nodes[weights[i].first - 1].connectV.size()));

        G_out.nodes[weights[i].first - 1].weight = i;
        G_in.nodes[weights[i].first - 1].weight = i;
    }

    //存储出度
//    while (getline(infile, lineStr)) {
//        //int col = 0;    //记录现在是第几列
//        row++;
//        //存成二维表结构
//        stringstream ss(lineStr);
//        string str;
//        vector<string> lineArray;
//        //按照空格分开
//        while (getline(ss, str, ' ')) {
//            lineArray.push_back(str);
//        }
//
//        //if(lineArray.size() == 0)   continue;   //只有一个结点，说明没有出度，可以直接忽略
//        //将数据放到vector容器中
//        Vertex temp;
//        for (int i = 0; i < lineArray.size(); i++) {
//            /*
//            if(!judgeNum(lineArray[i])){//判断输入的数据是否报错
//                cout << "第" << row << "行，第" << i << "个数据出现了问题" << endl;
//                return false;
//            }
//             */
//            int num = stringToInt(lineArray[i]);
//            //int tmpWeight;
//            //当i=0时代表输入num的是这个节点，其余的点是这个点一步就可以到达的点
//            if (i == 0) {
//                temp.val = num;
//            } else if (i == 1) {
//                temp.weight = num;
//                //G_out.weight[temp.val] = num;//赋值权重
//                //tmpWeight = num;
//            } else {
//                temp.connectV.push_back(num);
//            }
//        }
//        G_out.nodes.push_back(temp);
//        //G_out.vertex.push_back(temp.val);
//        G_out.location[temp.val] = G_out.location.size();
//        G_out.weights[temp.val] = temp.weight;
//        //G_out.nodes.back().degrees = G_out.nodes.back().connectV.size();
//    }
//
//    //存储入度
//    Vertex tmp;
//    //初始化入度
//    for (int i = 0; i < row; i++) {
//        tmp.val = G_out.nodes[i].val;
//        tmp.weight = G_out.weights[tmp.val];
//        G_in.nodes.push_back(tmp);
//        G_in.location[tmp.val] = i;
//        G_in.weights[tmp.val] = tmp.weight;
//        //G_in.vertex.push_back(tmp.val);
//    }
//
//    for (int i = 0; i < row; i++) {
//        for (int j = 0; j < G_out.nodes[i].connectV.size(); j++) {
//            int point = G_out.nodes[i].connectV[j];
//            int locate = G_in.location[point];
//            G_in.nodes[locate].connectV.push_back(G_out.nodes[i].val);
//        }
//    }
}

//二分查找
int binarySearch(vector<int> &nums, int target) {
    int left = 0, right = nums.size() - 1;
    while (right >= left) {
        int mid = (left + right) >> 1;
        if (nums[mid] == target) {
            return mid;
        } else if (nums[mid] > target) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }
    return -1;
}

//根据点的权重进行升序
static bool cmp_weight(Vertex &a, Vertex &b) {
    return a.weight < b.weight;
}

//根据lmax（v, k)排序
static bool cmp_lmax(PII &a, PII &b) {
    if (a.second == b.second) return a.first < b.first;
    return a.second < b.second;
}

//对点的出度进行升序
static bool cmp_out(Vertex &a, Vertex &b) {
    return a.connectV.size() < b.connectV.size();
}

//对点的入度进行升序
static bool cmp_in(Vertex &a, Vertex &b) {
    return a.connectV.size() < b.connectV.size();
}

void DeleteSingleNode(Graph &G, int nodeVal, int targetNode) {
    int point = G.location[targetNode];
    int i = 0;
    for (; i < G.nodes[point].connectV.size(); i++) {
        if (G.nodes[point].connectV[i] == nodeVal) break;
    }

    swap(G.nodes[point].connectV[i], G.nodes[point].connectV.back());
    G.nodes[point].connectV.pop_back();
}

/*
 * 根据点集获取图
 * G_out,G_in是完整的图
 * G_outr,G_inr是得到后的图
 */
void pointsChangeToGraph(Graph G_out, Graph G_in, Graph &G_outr, Graph &G_inr, vector<int> points) {
    map<int, bool> judge;
    //初始化judge
    for (int i = 0; i < G_out.nodes.size(); i++) {
        judge[G_out.nodes[i].val] = false;
    }

    //将目标点集所拥有的点置为true
    for (int i = 0; i < points.size(); i++) {
        judge[points[i]] = true;
    }

    for (int i = 0; i < points.size(); i++) {
        Vertex tmp;
        tmp.val = points[i];
        int locate = G_out.location[tmp.val];
        tmp.weight = G_out.nodes[locate].weight;
        vector<int> connectV_tmp;
        for (int j = 0; j < G_out.nodes[locate].connectV.size(); j++) {
            if (judge[G_out.nodes[locate].connectV[j]]) {
                connectV_tmp.push_back(G_out.nodes[locate].connectV[j]);
            }
        }
        G_outr.nodes.push_back(tmp);
        G_outr.location[tmp.val] = i;
        G_outr.weights[tmp.val] = tmp.weight;
        connectV_tmp.clear();
        for (int j = 0; j < G_in.nodes[locate].connectV.size(); j++) {
            if (judge[G_in.nodes[locate].connectV[j]]) {
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
void DeleteNode(Graph &G1, Graph &G2, int targetNode, vector<int> &nodes) {
    int point = G1.location[targetNode];

    //删除所有这个点所有出边或入边
    for (int i = 0; i < G1.nodes[point].connectV.size(); i++) {
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
bool QueryDcore(Graph G_out, Graph G_in, int k, int l, int q, Graph &G_outq, Graph &G_inq) {
    int flag;
    vector<int> useless;//存储与被删的点的关联点
    while (1) {
        flag = 1;

        for (int i = 0; i < G_out.nodes.size(); i++) {
            if (G_out.nodes[i].connectV.size() < l) {
                //delete node
                int e = G_out.nodes[i].val;
                DeleteNode(G_out, G_in, e, useless);
                DeleteNode(G_in, G_out, e, useless);
                flag = 0;
            }
        }

        for (int i = 0; i < G_in.nodes.size(); i++) {
            if (G_in.nodes[i].connectV.size() < k) {
                //delete node
                int e = G_in.nodes[i].val;
                DeleteNode(G_out, G_in, e, useless);
                DeleteNode(G_in, G_out, e, useless);
                flag = 0;
            }
        }
        if (flag) break;
    }

    if (G_out.location[q] == -1) return false;
    vector<int> arrVal;
    queue<int> que;
    map<int, int> stamp;
    for (int i = 0; i < G_out.nodes.size(); i++) {
        stamp[G_out.nodes[i].val] = 0;
    }
    stamp[q] = 1;
    que.push(q);
    arrVal.push_back(q);
    while (!que.empty()) {
        int tmp = que.front();
        que.pop();
        //在出度的图中找相连的点
        int locate = G_out.location[tmp];
        for (int i = 0; i < G_out.nodes[locate].connectV.size(); i++) {
            int nodeVal = G_out.nodes[locate].connectV[i];
            if (stamp[nodeVal] == 0) {
                que.push(nodeVal);
                arrVal.push_back(nodeVal);
                stamp[nodeVal] = 1;
            }
        }

        //在入度的图中找出相连的点
        locate = G_in.location[tmp];
        for (int i = 0; i < G_in.nodes[locate].connectV.size(); i++) {
            int nodeVal = G_in.nodes[locate].connectV[i];
            if (stamp[nodeVal] == 0) {
                que.push(nodeVal);
                arrVal.push_back(nodeVal);
                stamp[nodeVal] = 1;
            }
        }
    }

    for (int i = 0; i < arrVal.size(); i++) {
        Vertex tmp;
        tmp.val = arrVal[i];
        int locate = G_out.location[tmp.val];
        tmp.weight = G_out.weights[tmp.val];
        for (int j = 0; j < G_out.nodes[locate].connectV.size(); j++) {
            int t = G_out.nodes[locate].connectV[j];
            if (stamp[t]) {
                tmp.connectV.push_back(t);
            }
        }
        G_outq.nodes.push_back(tmp);
        G_outq.location[tmp.val] = i;
        G_outq.weights[tmp.val] = tmp.weight;
    }

    for (int i = 0; i < arrVal.size(); i++) {
        Vertex tmp;
        tmp.val = arrVal[i];
        tmp.weight = G_in.weights[tmp.val];
        int locate = G_in.location[tmp.val];
        for (int j = 0; j < G_in.nodes[locate].connectV.size(); j++) {
            int t = G_in.nodes[locate].connectV[j];
            if (stamp[t]) {
                tmp.connectV.push_back(t);
            }
        }
        G_inq.nodes.push_back(tmp);
        G_inq.location[tmp.val] = i;
        G_inq.weights[tmp.val] = tmp.weight;
    }

    return true;
}

//k是入度，l是出度
void Peel(Graph G_out, Graph G_in, int k, int l, int q, Graph &G_outr, Graph &G_inr) {
    Graph G_outq, G_inq;
    //获得极大子图
    bool flag = QueryDcore(G_out, G_in, k, l, q, G_outq, G_inq);
    if (!flag) {
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
    for (int i = 0; i < G_outq.nodes.size(); i++) {
        L.push_back(G_outq.nodes[i].val);
        //locate[L.back()] = i;
    }

    map<int, bool> judge;
    for (int i = 0; i < G_outq.nodes.size(); i++) {
        judge[G_outq.nodes[i].val] = false;
        //judge_tmp[G_outq.nodes[i].val] = false;//标记还有哪些点还剩余
    }
    //获取删减之后的顶点
    while (!L.empty()) {
        int tmp = L.front();
        D.clear();
        Q.push(tmp);
        judge[tmp] = true;
        while (!Q.empty()) {
            int pointVal = Q.front();
            judge[pointVal] = false;
            Q.pop();
            if (pointVal == q) {
//                int res = G_out.nodes[G_out.location[L.back()]].weight;
//                for(int i = 0; i < D.size(); i++){
//                    res = max(res, G_out.nodes[G_out.location[D[i]]].weight);
//                }
//                return res;
                while (!L.empty()) {
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
            for (int i = 0; i < nodes.size(); i++) {
                if (G_outq.nodes[G_outq.location[nodes[i]]].connectV.size() < l && !judge[G_outq.location[nodes[i]]]) {
                    Q.push(nodes[i]);
                    judge[G_outq.location[nodes[i]]] = true;
                    continue;
                }
                if (G_inq.nodes[G_inq.location[nodes[i]]].connectV.size() < k && !judge[G_inq.location[nodes[i]]]) {
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
}

//G_out存储的是每个点能一步到达的点，G_in存储的是每个点能被一步到达的点
void batchPeelingAlgorithm(Graph G_out, Graph G_in, int k, int l, int q, Graph &G_outr, Graph &G_inr) {
    Graph G_outq, G_inq;
    //获取存在q的极大子图
    bool flag = QueryDcore(G_out, G_in, k, l, q, G_outq, G_inq);//Line1
    if (!flag) {
        cout << "目标点不能构成k-core" << endl;
        return;
    }

    //对极大子图中的点根据权重进行排序
    sort(G_outq.nodes.begin(), G_outq.nodes.end(), cmp_weight);
    vector<int> L, H, S;//L存储图中的所有点的val
    for (int i = 0; i < G_outq.nodes.size(); i++) {
        L.push_back(G_outq.nodes[i].val);   //Line3
        if(G_outq.nodes[i].weight < G_outq.nodes[G_outq.location[q]].weight)
            S.push_back(G_outq.nodes[i].val);
    }
    H.assign(S.begin(), S.end());
    //H.assign(L.begin(), L.begin() + (L.size() / 2));//Line4
    while(!S.empty() || !H.empty()){
        vector<int> L0(L);
        L.assign(L.begin() + H.size(), L.end());
        //将L的点集转化为图

        Graph G_Outq, G_Inq;
        pointsChangeToGraph(G_outq, G_inq, G_Outq, G_Inq, L);

        //转化成的图为G_Outq和G_Inq
        //将转换好的图复制到G_outq和G_inq begin

        swap(G_Outq.location, G_outq.location);
        swap(G_Outq.weights, G_outq.weights);
        swap(G_Inq.location, G_Inq.location);
        swap(G_Inq.weights, G_Inq.weights);
        G_outq.nodes.assign(G_Outq.nodes.begin(), G_Outq.nodes.end());
        G_inq.nodes.assign(G_Inq.nodes.begin(), G_Inq.nodes.end());
        //end

        //update L with the degree constraint
        Graph G_outq_tmp, G_inq_tmp;
        flag = QueryDcore(G_outq, G_inq, k, l, q, G_outq_tmp, G_inq_tmp);//Line8

        //更新图begin
        swap(G_outq.location, G_outq_tmp.location);
        swap(G_outq.weights, G_outq_tmp.weights);
        swap(G_inq.location, G_inq_tmp.location);
        swap(G_inq.weights, G_inq_tmp.weights);
        G_outq.nodes.assign(G_outq_tmp.nodes.begin(), G_outq_tmp.nodes.end());
        G_inq.nodes.assign(G_inq_tmp.nodes.begin(), G_inq_tmp.nodes.end());
        //更新图end

        sort(G_outq.nodes.begin(), G_outq.nodes.end(), cmp_weight);

        vector<int> LStar;
        S.clear();
        for(int i = 0; i < G_outq.nodes.size(); i++){
            LStar.push_back(G_outq.nodes[i].val);
            if(G_outq.nodes[i].weight < G_outq.nodes[G_outq.location[q]].weight)
                S.push_back(G_outq.nodes[i].val);
        }
        if(!flag){
            H.assign(S.begin(), S.begin() + int(S.size() / 2));
            L.assign(L0.begin(), L0.end());
        }else{
            H.assign(H.begin(), H.begin() + int(H.size() / 2));
            L.assign(LStar.begin(), LStar.end());
        }
    }
    pointsChangeToGraph(G_outq, G_inq, G_outr, G_inr, L);
}

//index的子函数by zyf
void OutCoreDecom1(Graph G_out, Graph G_in, int k, map<int, int> &Lmax) {
    int ltmp = 1;
    vector<int> useless;//存储与被删的点的关联点
    while (!G_out.nodes.empty()) {
        vector<int> points;//记录每一轮中的被删的点
        int flag = 1;
        while (1) {
            flag = 1;
            //入度部分
            //sort(G_in.nodes.begin(), G_in.nodes.end(), cmp_out);
            for (int i = 0; i < G_in.nodes.size(); i++) {
                G_in.location[G_in.nodes[i].val] = i;
            }
            for (int i = 0; i < G_in.nodes.size(); i++) {
                if (G_in.nodes[i].connectV.size() < k) {
                    //delete node
                    int e = G_in.nodes[i].val;
                    DeleteNode(G_out, G_in, e, useless);
                    DeleteNode(G_in, G_out, e, useless);
                    points.push_back(e);
                    flag = 0;
                } else break;
            }

            //出度部分
            //sort(G_out.nodes.begin(), G_out.nodes.end(), cmp_out);
            for (int i = 0; i < G_out.nodes.size(); i++) {
                G_out.location[G_out.nodes[i].val] = i;
            }
            for (int i = 0; i < G_out.nodes.size(); i++) {
                if (G_out.nodes[i].connectV.size() < ltmp) {
                    int e = G_out.nodes[i].val;
                    DeleteNode(G_out, G_in, e, useless);
                    DeleteNode(G_in, G_out, e, useless);
                    points.push_back(e);
                    flag = 0;
                } else break;
            }

            if (flag = 1) {
                break;
            }
        }

        for (int i = 0; i < points.size(); i++) {
            Lmax[points[i]] = ltmp - 1;
        }
        points.clear();
        ltmp++;
    }
}

//index的子函数2 by zyf
void OutCoreDecom2(Graph G_out, Graph G_in, int k, map<int, int> &Lmax) {
    int ltmp = 1;
    vector<int> useless;
    map<int, bool> flag;
    for (int i = 0; i < G_out.nodes.size(); i++) {
        flag[G_out.nodes[i].val] = false;
    }
    while (!G_out.nodes.empty()) {
        vector<int> points;
        queue<int> Q;
        for (int i = 0; i < G_out.nodes.size(); i++) {
            if (G_out.nodes[i].connectV.size() == ltmp) {
                Q.push(G_out.nodes[i].val);
                flag[G_out.nodes[i].val] = true;
            }
        }
        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();
            flag[u] = false;
            Lmax[u] = ltmp;

            int locate = G_in.location[u];
            //入度part
            for (int i = 0; i < G_in.nodes[locate].connectV.size(); i++) {
                int v = G_in.nodes[locate].connectV[i];
                int vtmp = G_out.location[v];
                if (G_out.nodes[vtmp].connectV.size() <= ltmp + 1 && !flag[v]) {
                    Q.push(v);
                    flag[v] = true;
                }
            }

            //出度part
            locate = G_out.location[u];
            for (int i = 0; i < G_out.nodes[locate].connectV.size(); i++) {
                int v = G_out.nodes[locate].connectV[i];
                int vtmp = G_in.location[v];
                if (G_in.nodes[vtmp].connectV.size() <= ltmp + 1 && !flag[v]) {
                    Q.push(v);
                    flag[v] = true;
                }
            }

            DeleteNode(G_out, G_in, u, useless);
            DeleteNode(G_in, G_out, u, useless);
        }
        ltmp++;
    }
}

void QueryDcore(Graph G_out, Graph G_in, int k, int l, Graph &G_outq, Graph &G_inq) {
    int flag;
    vector<int> useless;
    while (1) {
        flag = 1;
//        sort(G_out.nodes.begin(), G_out.nodes.end(), cmp_out);
//        for (int i = 0; i < G_out.nodes.size(); i++) {
//            G_out.location[G_out.nodes[i].val] = i;
//        }
        for (int i = 0; i < G_out.nodes.size(); i++) {
            if (G_out.nodes[i].connectV.size() < l) {
                //delete node
                int e = G_out.nodes[i].val;
                DeleteNode(G_out, G_in, e, useless);
                DeleteNode(G_in, G_out, e, useless);
                flag = 0;
            }
        }

        //将点按照入度进行排序
//        sort(G_in.nodes.begin(), G_in.nodes.end(), cmp_in);
//        for (int i = 0; i < G_in.nodes.size(); i++) {
//            G_in.location[G_in.nodes[i].val] = i;
//        }
        for (int i = 0; i < G_in.nodes.size(); i++) {
            if (G_in.nodes[i].connectV.size() < k) {
                //delete node
                int e = G_in.nodes[i].val;
                DeleteNode(G_out, G_in, e, useless);
                DeleteNode(G_in, G_out, e, useless);
                flag = 0;
            }
        }
        if (flag) break;
    }

    for (int i = 0; i < G_out.nodes.size(); i++) {
        G_outq.nodes.push_back(G_out.nodes[i]);
        G_outq.location[G_out.nodes[i].val] = i;
        G_outq.weights[G_out.nodes[i].val] = G_out.nodes[i].weight;
    }
    for (int i = 0; i < G_in.nodes.size(); i++) {
        G_inq.nodes.push_back(G_in.nodes[i]);
        G_inq.location[G_in.nodes[i].val] = i;
        G_inq.weights[G_in.nodes[i].val] = G_in.nodes[i].weight;
    }
}

/*
 * 建立所有index索引的数组
 * I为三维数组，第一维为k值，第二维为点，第三维为该点下的所有点以及与该点对应的l值
 * location存储每个k值下每个点的位置
 */
void IndexConstructionAlgorithm(Graph G_out, Graph G_in, vector<vector<vector<PII>>> &I, vector<map<int, int>> &location) {
    int k = 1, kmax = 0;
    for (int i = 0; i < G_in.nodes.size(); i++) {
        kmax = max(kmax, int(G_in.nodes[i].connectV.size()));
    }
    //vector<PVV> I;
    //vector<map<int, int>> location;  //存储每个k值下每个点的位置
    while (k <= kmax) {
        Graph G_outr, G_inr;
        QueryDcore(G_out, G_in, k, 1, G_outr, G_inr);
        map<int, int> Lmax;
        OutCoreDecom2(G_outr, G_inr, k, Lmax);
        vector<int> G_points;
        map<int, int> PLTmp;
        for (int i = 0; i < G_outr.nodes.size(); i++) {
            G_points.push_back(G_outr.nodes[i].val);
            PLTmp[G_outr.nodes[i].val] = i;
        }
        location.push_back(PLTmp);
        vector<PII> P;
        vector<vector<PII>> Ptmp;
        PII tmp;
        for (int i = 0; i < G_outr.nodes.size(); i++) {//u
            int LmaxTmp = Lmax[G_outr.nodes[i].val];
            for (int j = 0; j < G_outr.nodes[i].connectV.size(); j++) {//v
                if (Lmax[G_outr.nodes[i].connectV[j]] >= 1) {
                    tmp.first = G_outr.nodes[i].connectV[j];
                    //.first.push_back(G_outr.nodes[i].connectV[j]);
                    if (Lmax[G_outr.nodes[i].connectV[j]] > LmaxTmp) {
                        tmp.second = LmaxTmp;
                        //P.second.push_back(LmaxTmp);
                    } else {
                        tmp.second = Lmax[G_outr.nodes[i].connectV[j]];
                        //P.second.push_back(Lmax[G_outr.nodes[i].connectV[j]]);
                    }
                    P.push_back(tmp);
                }
            }
            sort(P.begin(), P.end(), cmp_lmax);
            Ptmp.push_back(P);

        }
        I.push_back(Ptmp);
    }
}

/*
 * 查询前需要搞定I的所有值
 * 需要提供一个index，index应包含所有的有k值情况下的点
 * 需要提供q点、l值
 * location提供的是每个点的位置
 */
vector<int> IndexQueryAlgorithm(vector<vector<PII>> I, int q, int l, map<int, int> location) {
    queue<int> que;
    vector<int> res;    //存储满足条件的所有点
    res.push_back(q);
    que.push(q);
    vector<bool> visited(I.size(), false);
    visited[location[q]] = true;
    while (!que.empty()) {
        int u = que.front(), locate = location[u];
        que.pop();
        for (int i = 0; i < I[locate].size(); i++) {
            if (I[locate][i].second >= l) {
                res.push_back(I[locate][i].first);//增加满足条件的点
                if (!visited[location[I[locate][i].first]]) {
                    que.push(I[locate][i].first);
                    visited[location[I[locate][i].first]] = true;
                }
            } else break;
        }
    }
    sort(res.begin(), res.end());
    vector<int> ret;
    ret.push_back(res.front());
    for (int i = 1; i < res.size(); i++) {
        if (res[i] != res[i - 1]) {
            ret.push_back(res[i]);
        }
    }
    return ret;
}

/*
 *查询前置
 *在此函数中获取查询所需要的条件
*/
void BehindQuery(Graph G_out, Graph G_in) {
    vector<vector<vector<PII>>> I;
    vector<map<int, int>> location;
    IndexConstructionAlgorithm(G_out, G_in, I, location);
    /*手动输入查询的信息begin*/
    int k, l, q;
    cout << "请输入你所要查询的k值、l值和q值：" << endl;
    cin >> k;
    cout << "输入k值：" << endl;
    if(k <= 0){
        cout << "输入的k值有误!" << endl;
        return;
    }
    cout << "输入l值：" << endl;
    cin >> l;
    cout << "输入q值:" << endl;
    cin >> q;
    if(G_out.location.find(q) == G_out.location.end()){
        cout << "输入的q值有误" << endl;
        return;
    }
    /*手动输入查询的信息end*/
    vector<vector<PII>> Itmp;
    Itmp = I[k - 1];
    vector<int> points = IndexQueryAlgorithm(Itmp, q, l, location[k - 1]);
}


int main() {
    //验证basic peel是否能成功使用
    Graph G_out, G_in;
    readDataTmp(G_out, G_in);

    //测试Query Community
    Graph G_outq, G_inq;
    bool flag = QueryDcore(G_out, G_in, 2, 3, 1, G_outq, G_inq);
    if(!flag){
        cout << "q不满足" << endl;
    }
    //get edges in the new graph
//    for(int i = 0; i < G_outq.nodes.size(); i++){
//        for(int j = 0; j < G_outq.nodes[i].connectV.size(); j++){
//            cout << G_outq.nodes[i].val << ' ' << G_outq.nodes[i].connectV[j] << endl;
//        }
//    }

    return 0;
}