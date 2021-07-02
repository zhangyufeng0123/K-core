#include<bits/stdc++.h>

using namespace std;

const string filename = "chess.txt"; //change

int KLmax = 1;

typedef long long ll;
typedef pair<int, int> PII;
typedef vector<PII> VPII;

struct Vertex {
    int val;    //点的编号
    vector<int> connectV;   //存储与该点连接的点的编号
    int weight; //点的权重
};

struct Graph {
    vector<Vertex> nodes;   //存储图中的所有点的信息
    vector<int> location;   //存储图中的所有点的位置
    vector<int> weights;    //存储图中的所有点的权重
};

//用计数排序对cmp_weight进行优化
void sort_weight(vector<Vertex> &num) {
    int wmax = num.front().weight, wmin = num.front().weight;
    for (int i = 1; i < num.size(); i++) {
        wmax = max(wmax, num[i].weight);
        wmin = min(wmin, num[i].weight);
    }
    vector<Vertex> tmp(wmax - wmin + 1);
    for (int i = 0; i < num.size(); i++) {
        tmp[num[i].weight - wmin] = num[i];
    }
    int k = 0;
    for (int i = 0; i <= (wmax - wmin); i++) {
        if (tmp[i].weight != 0) {
            num[k++] = tmp[i];
        }
    }
}

static bool cmp_weight(Vertex &a, Vertex &b) {
    return a.weight > b.weight;
}

static bool cmp_val(PII &a, PII &b) {
    if (a.first == b.first) return a.second < b.second;
    return a.first < b.first;
}

//用计数排序对sort进行优化
void sort_lmax(vector<PII> &num) {
    int lmax = num.front().second, lmin = num.front().second;
    //确认其区间
    for (int i = 1; i < num.size(); i++) {
        lmax = max(lmax, num[i].second);
        lmin = min(lmin, num[i].second);
    }
    vector<vector<PII>> tmp(lmax - lmin + 1);
    for (int i = 0; i < num.size(); i++) {
        tmp[num[i].second - lmin].push_back(num[i]);
    }
    int k = 0;
    for (int i = (lmax - lmin); i >= 0; i--) {
        while (!tmp[i].empty()) {
            num[k++] = tmp[i].back();
            tmp[i].pop_back();
        }
    }
}

//二分查找
int binarySearch(vector<PII> nums, int targetWeight) {
    int left = 0, right = nums.size() - 1;
    while (right >= left) {
        int mid = (left + right) >> 1;
        if (nums[mid].second == targetWeight) {
            return mid;
        } else if (nums[mid].second > targetWeight) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }
    return -1;
}

/*
 * 对初始化的权重进行排序，不需要考虑耗时，所以使用sort
 * */
static bool cmp_initWeight(PII &a, PII &b) {
    if (a.second == b.second) return a.first < b.first;
    return a.second < b.second;
}

/*
 * 读取数据
 * */
bool ReadData(Graph &G_out, Graph &G_in) {
    ifstream infile;
    infile.open(filename);
    if (!infile) {
        cout << "打开文件错误！" << endl;
        return false;
    }
    int points, edges;
    infile >> points >> edges;
    VPII weights;

    //初始化nodes的长度
    vector<Vertex> initTmp(points + 1);
    G_out.nodes.assign(initTmp.begin(), initTmp.end());
    G_in.nodes.assign(initTmp.begin(), initTmp.end());

    //初始化location和weights
    vector<int> lw(points + 1);
    G_in.location.assign(lw.begin(), lw.end());
    G_out.location.assign(lw.begin(), lw.end());
    G_in.weights.assign(lw.begin(), lw.end());
    G_out.weights.assign(lw.begin(), lw.end());

    for (int i = 0; i < points + 1; i++) {
        //初始化权值begin
        PII tmp;
        tmp.first = i;
        tmp.second = 0;
        weights.push_back(tmp);
        //初始化权值end

        //初始化图 begin
        G_out.nodes[i].val = i;
        G_out.location[i] = i;
        G_in.nodes[i].val = i;
        G_in.location[i] = i;
        //初始化图 end
    }

    //输入边 begin
    for (int i = 0; i < edges; i++) {
        int a, b;
        infile >> a >> b;
        if (a == b)continue;
        weights[a].second--;
        weights[b].second--;
        G_out.nodes[a].connectV.push_back(b);
        G_in.nodes[b].connectV.push_back(a);
    }
    //输入边 end

    int k = 1;
    sort(weights.begin(), weights.end(), cmp_initWeight);

    for (auto a : weights) {
        G_out.nodes[a.first].weight = k;
        G_in.nodes[a.first].weight = k++;
    }

    //fill the weights
    for (int i = 0; i < points; i++) {
        G_out.weights[i] = G_out.nodes[i].weight;
    }
    for (int i = 0; i < points; i++) {
        G_in.weights[i] = G_in.nodes[i].weight;
    }
    return true;
}

/*
 * 获取图中的点中最大的编号
 * */
int GetMaxNodeVal(vector<Vertex> nodes) {
    int maxVal = 0;
    for (auto node : nodes) {
        maxVal = max(maxVal, node.val);
    }
    return maxVal;
}

/*
 * 删除单个点
 * */
void DeleteSingleNode(Graph &G, int nodeVal, int targetNode) {
    int point = G.location[targetNode];
    if (point == -1) return;
    int i = 0;
    for (; i < G.nodes[point].connectV.size(); i++) {
        if (G.nodes[point].connectV[i] == nodeVal) break;
    }
    G.nodes[point].connectV[i] = G.nodes[point].connectV.back();
    G.nodes[point].connectV.pop_back();
}

/*
 * 1.删除G2中的targetNode的点
 * 2.1 先找到第一张图要删点的位置
 * 2.2 将这个点与最后一个点进行位置交换
 * nodes记录删除了哪些点
 * */
void DeleteNodes(Graph &G1, Graph &G2, int targetNode, vector<int> &nodes) {
    int pointLocation = G1.location[targetNode];

    //删除图2与targetNode相连的点
    for (int i = 0; i < G1.nodes[pointLocation].connectV.size(); i++) {
        int tempPoint = G1.nodes[pointLocation].connectV[i];
        nodes.push_back(tempPoint);
        DeleteSingleNode(G2, targetNode, tempPoint);
    }

    int lastNodeVal = G1.nodes.back().val;
    G1.nodes[pointLocation] = G1.nodes.back();
    G1.nodes.pop_back();
    G1.location[lastNodeVal] = pointLocation;
    G1.location[targetNode] = -1;
}

/*
 * 将A的入度出度图转到B的入度出度图
 * */
void GraphAToB(Graph G_out, Graph G_in, Graph &G_outr, Graph &G_inr) {
    G_outr.nodes = move(G_out.nodes);
    G_inr.nodes = move(G_in.nodes);
    G_outr.location = move(G_out.location);
    G_inr.location = move(G_in.location);
    G_outr.weights = move(G_out.weights);
    G_inr.weights = move(G_in.weights);
}

/*
 * 不包含q的(k, l)-core查询
 * */
bool QueryCommunity(Graph G_out, Graph G_in, int k, int l, Graph &G_outr, Graph &G_inr) {
    //寻找图中的最大的编号
    int maxVal = GetMaxNodeVal(G_out.nodes);
    vector<bool> judge(maxVal + 1, false);
    queue<int> Q;

    for (auto node : G_out.nodes) {
        if (node.connectV.size() < l && !judge[node.val]) {
            Q.push(node.val);
            judge[node.val] = true;
        }
    }

    for (auto node : G_in.nodes) {
        if (node.connectV.size() < k && !judge[node.val]) {
            Q.push(node.val);
            judge[node.val] = true;
        }
    }

    vector<int> useless;
    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();

        int uOutLocation = G_out.location[u], uInLocation = G_in.location[u];
        for (int i = 0; i < G_out.nodes[uOutLocation].connectV.size(); i++) {
            int point = G_out.nodes[uOutLocation].connectV[i];  //获取该点的编号
            if ((G_out.nodes[G_out.location[point]].connectV.size() < l ||
                 G_in.nodes[G_in.location[point]].connectV.size() < k + 1) && !judge[point]) {
                Q.push(point);
                judge[point] = true;
            }
        }

        for (int i = 0; i < G_in.nodes[uInLocation].connectV.size(); i++) {
            int point = G_in.nodes[uInLocation].connectV[i];    //获取该点编号
            if ((G_out.nodes[G_out.location[point]].connectV.size() < l + 1 ||
                 G_in.nodes[G_in.location[point]].connectV.size() < k) && !judge[point]) {
                Q.push(point);
                judge[point] = true;
            }
        }
        DeleteNodes(G_out, G_in, u, useless);
        DeleteNodes(G_in, G_out, u, useless);
    }
    GraphAToB(G_out, G_in, G_outr, G_inr);
}

/*
 * 包含q的(k,l)-core查询
 * */
bool QueryCommunity(Graph G_out, Graph G_in, int k, int l, int q, Graph &G_outq, Graph &G_inq) {
    //寻找图中的最大的编号
    int maxVal = GetMaxNodeVal(G_out.nodes);
    vector<bool> judge(maxVal + 1, false);
    queue<int> Q;

    for (auto node : G_out.nodes) {
        if (node.connectV.size() < l && !judge[node.val]) {
            Q.push(node.val);
            judge[node.val] = true;
        }
    }

    for (auto node : G_in.nodes) {
        if (node.connectV.size() < k && !judge[node.val]) {
            Q.push(node.val);
            judge[node.val] = true;
        }
    }

    vector<int> useless;
    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();

        int uOutLocation = G_out.location[u], uInLocation = G_in.location[u];
        for (int i = 0; i < G_out.nodes[uOutLocation].connectV.size(); i++) {
            int point = G_out.nodes[uOutLocation].connectV[i];  //获取该点的编号
            if ((G_out.nodes[G_out.location[point]].connectV.size() < l ||
                 G_in.nodes[G_in.location[point]].connectV.size() < k + 1) && !judge[point]) {
                Q.push(point);
                judge[point] = true;
            }
        }

        for (int i = 0; i < G_in.nodes[uInLocation].connectV.size(); i++) {
            int point = G_in.nodes[uInLocation].connectV[i];    //获取该点编号
            if ((G_out.nodes[G_out.location[point]].connectV.size() < l + 1 ||
                 G_in.nodes[G_in.location[point]].connectV.size() < k) && !judge[point]) {
                Q.push(point);
                judge[point] = true;
            }
        }

        DeleteNodes(G_out, G_in, u, useless);
        DeleteNodes(G_in, G_out, u, useless);
    }
    if (G_out.location[q] == -1) return false;
    vector<int> arrVal;
    for (auto a : judge) {
        a = false;
    }
    judge[q] = true;
    Q.push(q);
    arrVal.push_back(q);
    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();
        //在出度图中找相连的点
        for (auto a : G_out.nodes[G_out.location[u]].connectV) {
            if (!judge[a]) {
                judge[a] = true;
                Q.push(a);
                arrVal.push_back(a);
            }
        }

        //在入度图中找相连的点
        for (auto a : G_in.nodes[G_in.location[u]].connectV) {
            if (!judge[a]) {
                judge[a] = true;
                Q.push(a);
                arrVal.push_back(a);
            }
        }
    }

    maxVal = 0;
    for (auto a : arrVal) {
        maxVal = max(maxVal, a);
    }
    vector<Vertex> initmp1(arrVal.size());
    G_outq.nodes.assign(initmp1.begin(), initmp1.end());
    G_inq.nodes.assign(initmp1.begin(), initmp1.end());
    vector<int> initmp2(maxVal + 1);
    G_outq.location.assign(initmp2.begin(), initmp2.end());
    G_outq.weights.assign(initmp2.begin(), initmp2.end());
    G_inq.location.assign(initmp2.begin(), initmp2.end());
    G_inq.weights.assign(initmp2.begin(), initmp2.end());

    k = 0;
    for (auto a : arrVal) {
        Vertex tmp;
        tmp.val = a;
        tmp.weight = G_out.weights[a];
        for (auto t : G_out.nodes[G_out.location[a]].connectV) {
            if (judge[t]) tmp.connectV.push_back(t);
        }

        G_outq.nodes[k] = tmp;
        G_outq.location[a] = k++;
        G_outq.weights[a] = tmp.weight;
    }

    k = 0;
    for (auto a : arrVal) {
        Vertex tmp;
        tmp.val = a;
        tmp.weight = G_in.weights[a];
        for (auto t : G_in.nodes[G_in.location[a]].connectV) {
            if (judge[t]) tmp.connectV.push_back(t);
        }

        G_inq.nodes[k] = tmp;
        G_inq.location[a] = k++;
        G_inq.weights[a] = tmp.weight;
    }

    return true;
}

/*
 * 构建index索引的子函数
 * */
void OutCoreDecom(Graph G_out, Graph G_in, int k, unordered_map<int, int> &Lmax) {
    int ltmp = 1;
    vector<int> useless;
    unordered_map<int, bool> flag;
    for (auto a : G_out.nodes) {
        flag[a.val] = false;
    }

    while (!G_out.nodes.empty()) {
        vector<int> points;
        queue<int> Q;
        for (int i = 0; i < G_out.nodes.size(); i++)
            if (G_out.nodes[i].connectV.size() == ltmp) {
                Q.push(G_out.nodes[i].val);
                flag[G_out.nodes[i].val] = true;
            }

        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();
            flag[u] = false;
            Lmax[u] = ltmp;

            //入度part
            for (auto v : G_in.nodes[G_in.location[u]].connectV) {
                if (G_out.nodes[G_out.location[v]].connectV.size() <= ltmp + 1 && !flag[v]) {
                    Q.push(v);
                    flag[v] = true;
                }
            }

            //出度part
            for (auto v : G_out.nodes[G_out.location[u]].connectV) {
                if (G_in.nodes[G_in.location[v]].connectV.size() <= k && !flag[v]) {
                    Q.push(v);
                    flag[v] = true;
                }
            }

            DeleteNodes(G_out, G_in, u, useless);
            DeleteNodes(G_in, G_out, u, useless);
        }
        ltmp++;
    }
}

/*
 * 建立所有index索引的数组
 * I分为三维，第一维为k，第二维为点的编号，第三维为(与第二维相连接的点, lmax)
 * location存储每个k值下每个点的位置，(point,location)
 * */
void IndexConstructAlgorithm(Graph G_out, Graph G_in, vector<vector<vector<PII>>> &I, vector<vector<int>> &location) {
    int k = 1;
    while (k <= KLmax) {
        Graph G_outr, G_inr;
        QueryCommunity(G_out, G_in, k, 1, G_outr, G_inr);
        unordered_map<int, int> Lmax;
        OutCoreDecom(G_outr, G_inr, k, Lmax);
        int maxVal = GetMaxNodeVal(G_outr.nodes);
        vector<int> locationtmp(maxVal + 1, -1);
        for (int i = 0; i < G_outr.nodes.size(); i++) {
            locationtmp[G_outr.nodes[i].val] = i;
        }
        location.push_back(locationtmp);
        VPII P;
        vector<vector<PII>> PLtmp;
        PII tmp;
        unordered_map<int, bool> judgePoint;
        for (auto node : G_outr.nodes) {
            judgePoint[node.val] = false;
        }
        vector<int> points;
        for (int i = 0; i < G_outr.nodes.size(); i++) {
            P.clear();

            for (int j = 0; j < points.size(); j++) {
                judgePoint[points[j]] = false;
            }
            points.clear();
            int point = G_outr.nodes[i].val;
            int LmaxTmp = Lmax[point];

            for (int j = 0; j < G_outr.nodes[i].connectV.size(); j++) {
                if (Lmax[G_outr.nodes[i].connectV[j]] >= 1 && !judgePoint[G_outr.nodes[i].connectV[j]]) {
                    judgePoint[G_outr.nodes[i].connectV[j]] = true;
                    points.push_back(G_outr.nodes[i].connectV[j]);
                    tmp.first = G_outr.nodes[i].connectV[j];
                    if (Lmax[G_outr.nodes[i].connectV[j]] > LmaxTmp) {
                        tmp.second = LmaxTmp;
                    } else {
                        tmp.second = Lmax[G_outr.nodes[i].connectV[j]];
                    }
                }
                P.push_back(tmp);
            }

            int locate = G_inr.location[point];
            for (int j = 0; j < G_inr.nodes[locate].connectV.size(); j++) {
                if (Lmax[G_inr.nodes[locate].connectV[j]] >= 1 && !judgePoint[G_inr.nodes[locate].connectV[j]]) {
                    tmp.first = G_inr.nodes[locate].connectV[j];
                    if (Lmax[G_inr.nodes[locate].connectV[j]] > LmaxTmp) {
                        tmp.second = LmaxTmp;
                    } else {
                        tmp.second = Lmax[G_inr.nodes[locate].connectV[j]];
                    }
                    P.push_back(tmp);
                }
            }

            sort_lmax(P);
            PLtmp.push_back(P);
        }
        I.push_back(PLtmp);
        k++;
    }
}

/*
 * 计算整个图的kmax，作为数据所需
 * */
int cal_kmax(Graph G_out, Graph G_in) {
    int k = 1;
    while (1) {
        Graph G_outr, G_inr;
        QueryCommunity(G_out, G_in, k, 0, G_outr, G_inr);
        if (G_outr.nodes.empty()) {
            return k - 1;
        }
        k++;
        GraphAToB(G_outr, G_inr, G_out, G_in);
    }
}

/*
 * 计算整个图的lmax，作为数据所需
 * */
int cal_lmax(Graph G_out, Graph G_in) {
    int l = 1;
    while (1) {
        Graph G_outr, G_inr;
        QueryCommunity(G_out, G_in, 0, l, G_outr, G_inr);
        if (G_outr.nodes.empty()) {
            return l - 1;
        }
        l++;
        GraphAToB(G_outr, G_inr, G_out, G_in);
    }
}

/*
 * 计算整个图的klmax，用来满足乘以系数之后的kl能有一个满足的图
 * */
void cal_KLmax(Graph G_out, Graph G_in) {
    while (1) {
        Graph G_outr, G_inr;
        QueryCommunity(G_out, G_in, KLmax, KLmax, G_outr, G_inr);
        if (G_outr.nodes.empty()) {
            KLmax--;
            return;
        }
        KLmax++;
    }
}

/*
 * 点集转换为图集
 * */
void PointsChangeToGraph(Graph G_out, Graph G_in, Graph &G_outr, Graph &G_inr, vector<int> points) {
    int maxVal = GetMaxNodeVal(G_out.nodes);
    vector<bool> judge(maxVal + 1, false);
    for (int i = 0; i < points.size(); i++) {
        judge[points[i]] = true;
    }
    G_outr.location.assign(maxVal + 1, -1);
    G_outr.weights.assign(maxVal + 1, -1);
    G_inr.location.assign(maxVal + 1, -1);
    G_inr.weights.assign(maxVal + 1, -1);

    vector<Vertex> initmp(points.size());
    G_inr.nodes.assign(initmp.begin(), initmp.end());
    G_outr.nodes.assign(initmp.begin(), initmp.end());

    int k = 0;
    for (auto point : points) {
        Vertex tmp;
        tmp.val = point;
        int locate = G_out.location[point];
        tmp.weight = G_out.weights[point];

        for (auto node : G_out.nodes[locate].connectV) {
            if (judge[node])
                tmp.connectV.push_back(node);
        }

        G_outr.nodes[k] = tmp;
        G_outr.location[point] = k;
        G_outr.weights[point] = tmp.weight;

        tmp.connectV.clear();
        locate = G_in.location[point];
        for (auto node : G_in.nodes[locate].connectV) {
            if (judge[node])
                tmp.connectV.push_back(node);
        }
        G_inr.nodes[k] = tmp;
        G_inr.location[point] = k++;
        G_inr.weights[point] = tmp.weight;
    }
}

void PointsChangeToGraph(Graph G_out, Graph G_in, Graph &G_outr, Graph &G_inr, vector<PII> points) {
    int maxVal = GetMaxNodeVal(G_out.nodes);
    vector<bool> judge(maxVal + 1, false);
    for (int i = 0; i < points.size(); i++) {
        judge[points[i].first] = true;
    }

    G_outr.location.assign(maxVal + 1, -1);
    G_outr.weights.assign(maxVal + 1, -1);
    G_inr.location.assign(maxVal + 1, -1);
    G_inr.weights.assign(maxVal + 1, -1);

    vector<Vertex> initmp(points.size());
    G_inr.nodes.assign(initmp.begin(), initmp.end());
    G_outr.nodes.assign(initmp.begin(), initmp.end());

    int k = 0;
    for (auto point : points) {
        Vertex tmp;
        tmp.val = point.first;
        int locate = G_out.location[point.first];
        tmp.weight = G_out.weights[point.first];

        for (auto node : G_out.nodes[locate].connectV) {
            if (judge[node])
                tmp.connectV.push_back(node);
        }

        G_outr.nodes[k] = tmp;
        G_outr.location[point.first] = k;
        G_outr.weights[point.first] = tmp.weight;

        tmp.connectV.clear();
        locate = G_in.location[point.first];
        for (auto node : G_in.nodes[locate].connectV) {
            if (judge[node])
                tmp.connectV.push_back(node);
        }
        G_inr.nodes[k] = tmp;
        G_inr.location[point.first] = k++;
        G_inr.weights[point.first] = tmp.weight;
    }
}

/*
 * baseline
 * */
void BasicPeelingAlgorithm(Graph G_out, Graph G_in, int k, int l, int q, Graph &G_outr, Graph &G_inr) {
    //备份G_out, G_in
    Graph G_outTmp, G_inTmp;

    //对极大子图中的点进行排序
    sort_weight(G_out.nodes);

    //更新G_out中的location
    for (int i = 0; i < G_out.nodes.size(); i++) {
        G_out.location[G_out.nodes[i].val] = i;
    }

    GraphAToB(G_out, G_in, G_outTmp, G_inTmp);

    queue<int> Q;
    vector<int> D;
    vector<PII> L;
    for (auto node : G_out.nodes) {
        PII tmpe;
        tmpe.first = node.val;
        tmpe.second = node.weight;
        L.push_back(tmpe);
    }

    int maxVal = GetMaxNodeVal(G_out.nodes);
    vector<int> judge(maxVal + 1, 0);

    while (!L.empty()) {
        int tmp = L.front().first;
//        for (int i = 0; i < D.size(); i++) {
//            judge[D[i]] = 0;
//        }
        D.clear();
        Q.push(tmp);
        judge[tmp] = true;
        while (!Q.empty()) {
            int pointVal = Q.front();
            judge[pointVal] = 0;
            Q.pop();

            if (pointVal == q) {
                while (!L.empty()) {
                    int tmpe = L.back().first;
                    D.push_back(tmpe);
                    L.pop_back();
                }
                break;
            }

            vector<int> nodes;
            DeleteNodes(G_out, G_in, pointVal, nodes);
            DeleteNodes(G_in, G_out, pointVal, nodes);

            for (int i = 0; i < nodes.size(); i++) {
                if (G_out.nodes[G_out.location[nodes[i]]].connectV.size() < l && !judge[nodes[i]]) {
                    Q.push(nodes[i]);
                    judge[nodes[i]] = 1;
                    continue;
                }
                if (G_in.nodes[G_in.location[nodes[i]]].connectV.size() < k && !judge[nodes[i]]) {
                    Q.push(nodes[i]);
                    judge[nodes[i]] = 1;
                }
            }
            int targetWeight = G_out.weights[pointVal];
            int location = binarySearch(L, targetWeight);
            L.erase(L.begin() + location);
            D.push_back(pointVal);
        }
    }
    PointsChangeToGraph(G_outTmp, G_inTmp, G_outr, G_inr, D);
}

bool JudgeQBelongToLH(Graph G_out, Graph G_in, int k, int l, int q, int wpiovt) {
    for (auto a : G_out.nodes[G_out.location[q]].connectV) {
        if (G_out.weights[a] > wpiovt) l--;
        if (!l) break;
    }
    if (l) return false;
    for (auto a : G_in.nodes[G_in.location[q]].connectV) {
        if (G_in.weights[a] > wpiovt) k--;
        if (!k) break;
    }
    if (k) return false;
    return true;
}

void BatchPeelingAlgorithm(Graph G_out, Graph G_in, int k, int l, int q, Graph &G_outr, Graph &G_inr) {
    Graph G_outTmp, G_inTmp;
    sort_weight(G_out.nodes);//根据权重对G_out.nodes进行升序排序
    for (int i = 0; i < G_out.nodes.size(); i++) {
        G_out.location[G_out.nodes[i].val] = i;
    }

    GraphAToB(G_out, G_in, G_outTmp, G_inTmp);//备份图
    vector<PII> L(G_out.nodes.size()), H, S;
    for (int i = 0; i < G_out.nodes.size(); i++) {
        PII tmp;
        tmp.first = G_out.nodes[i].val;
        tmp.second = G_out.nodes[i].weight;
        L[i] = tmp;
    }
    int tmpWeight = G_out.weights[q];
    for (auto l : L) {
        if (l.second < tmpWeight)
            S.push_back(l);
        else break;
    }

    int maxVal = GetMaxNodeVal(G_out.nodes);
    vector<int> judgePoint(maxVal + 1, 0);

    H.assign(S.begin(), S.begin() + int((S.size()) / 2));
    while (!S.empty() || !H.empty()) {
        bool flag = JudgeQBelongToLH(G_out, G_in, k, l, q, H.empty() ? 0 : H.back().second);
        if (!flag) {
            H.assign(H.begin(), H.begin() + int((H.size()) / 2));
            continue;
        }
        vector<PII> L0(L);
        L.assign(L.begin() + H.size(), L.end());

        Graph G_OUT, G_IN;
        PointsChangeToGraph(G_outTmp, G_inTmp, G_OUT, G_IN, L);
        flag = QueryCommunity(G_OUT, G_IN, k, l, q, G_out, G_in);
        if (!flag) {
            H.assign(H.begin(), H.begin() + int((H.size()) / 2));
            //L.assign(L0.begin(), L0.end());
            L = move(L0);
        } else {
            //GraphAToB(G_out, G_in, G_outTmp, G_inTmp);
            sort_weight(G_out.nodes);
            vector<PII> LStar(G_out.nodes.size());
            for (int i = 0; i < G_out.nodes.size(); i++) {
                PII tmp;
                tmp.first = G_out.nodes[i].val;
                tmp.second = G_out.nodes[i].weight;
                LStar[i] = tmp;
            }

            int u = LStar.front().first;
            vector<int> D;
            queue<int> Q;
            Q.push(u);
            while (!Q.empty()) {
                int utmp = Q.front();
                judgePoint[utmp] = 1;
                Q.pop();
                if (utmp == q) {
                    while (!LStar.empty()) {
                        D.push_back(LStar.back().first);
                        LStar.pop_back();
                    }
                    PointsChangeToGraph(G_outTmp, G_inTmp, G_outr, G_inr, D);
                    return;
                }

                vector<int> nodes;
                DeleteNodes(G_out, G_in, utmp, nodes);
                DeleteNodes(G_in, G_out, utmp, nodes);

                for (int i = 0; i < nodes.size(); i++) {
                    if (G_out.nodes[G_out.location[nodes[i]]].connectV.size() < l && !judgePoint[nodes[i]]) {
                        Q.push(nodes[i]);
                        judgePoint[nodes[i]] = 1;
                        continue;
                    }
                    if (G_in.nodes[G_in.location[nodes[i]]].connectV.size() < k && !judgePoint[nodes[i]]) {
                        Q.push(nodes[i]);
                        judgePoint[nodes[i]] = 1;
                    }
                }

                int targetWeight = G_out.weights[utmp];
                int location = binarySearch(L, targetWeight);
                LStar.erase(LStar.begin() + location);
                D.push_back(utmp);
            }

            S.clear();
            for (int i = 0; i < LStar.size(); i++) {
                if (LStar[i].second < tmpWeight)
                    S.push_back(LStar[i]);
                else break;
            }
            H.assign(S.begin(), S.begin() + int((S.size()) / 2));
            //L.assign(LStar.begin(), LStar.end());
            L = move(LStar);
        }
    }
}

vector<int> IndexQueryAlgorithm(vector<vector<PII>> I, int q, int l, vector<int> location) {
    queue<int> que;
    vector<int> res;
    res.push_back(q);
    que.push(q);
    vector<int> visited(location.size(), 0);
    visited[q] = 1;
    while (!que.empty()) {
        int u = que.front(), locate = location[u];
        que.pop();
        for (int i = 0; i < I[locate].size(); i++) {
            if (I[locate][i].second >= l) {
                if (!visited[location[I[locate][i].first]]) {
                    res.push_back(I[locate][i].first);
                    que.push(I[locate][i].first);
                    visited[location[I[locate][i].first]] = true;
                }
            } else break;
        }
    }

    return res;
}

void solve() {
    double c[5] = {0.1, 0.3, 0.5, 0.7, 0.9};    //c是实验所需要的系数

    //统计实验的时间变量
    clock_t start_QC, start_BasicP, start_BatchP, start_Ind;
    clock_t end_QC, end_BasicP, end_BatchP, end_Ind;
    clock_t QCBasicP = 0, QCBatchP = 0, IndBatchP = 0;
    Graph G_out, G_in;
    //获取数据
    ReadData(G_out, G_in);

    //输出kmax和lmax，用来debug
    cout << "Kmax = " << cal_kmax(G_out, G_in) << endl;
    cout << "Lmax = " << cal_lmax(G_out, G_in) << endl;
    cal_KLmax(G_out, G_in);
    cout << "KLmax = " << KLmax << endl;
    //Index 预处理 begin
    vector<vector<vector<PII>>> I; //I分为三维，第一维为k，第二维为点的编号，第三维为(与第二维相连接的点, lmax)
    vector<vector<int>> location;
    IndexConstructAlgorithm(G_out, G_in, I, location);
    //Index 预处理 end

    for (int i = 0; i < 5; i++) {
        int k = KLmax * c[i];
        if (KLmax * c[i] > k) k++;
        int l = k;
        //k = 3, l = 3;
        int weight_BasicP = G_out.nodes.size(), weight_BatchP = G_out.nodes.size(), weight_Ind = G_out.nodes.size();
        QCBasicP = 0, QCBatchP = 0, IndBatchP = 0;
        cout << k << ' ' << l << ":" << endl;
        int q;
        //获取G_outq, G_inq的值
        Graph G_outq, G_inq;
        Graph G_otmp, G_itmp;
        QueryCommunity(G_out, G_in, k, l, G_otmp, G_itmp);
        if (G_otmp.nodes.empty()) {
            cout << i << endl;
            cout << "不存在(" << k << "," << l << ")-core" << endl;
            break;
        }
        int e = rand() % G_otmp.nodes.size();
        //q = 3;
        q = G_otmp.nodes[e].val;
        cout << "Finish the random and q = " << q << endl;

        start_QC = clock();
        bool flag = QueryCommunity(G_out, G_in, k, l, q, G_outq, G_inq);
        end_QC = clock();
        QCBatchP += (end_QC - start_QC);
        QCBasicP += (end_QC - start_QC);
//
//        for(int i = 0; i < G_outq.nodes.size(); i++){
//            cout << i << ' ' << G_outq.nodes[i].val << ' ' << G_outq.nodes[i].weight << endl;
//        }

        //Baseline Part
        Graph G_outBasicP, G_inBasicP;
        start_BasicP = clock();
        BasicPeelingAlgorithm(G_outq, G_inq, k, l, q, G_outBasicP, G_inBasicP);
        end_BasicP = clock();

        //打印生成的图中的所有点
        cout << G_outBasicP.nodes.size() << endl;
        sort(G_outBasicP.nodes.begin(), G_outBasicP.nodes.end(), cmp_weight);
        for (auto node : G_outBasicP.nodes) {
            cout << node.val << ' ';
        }
        cout << endl;
        //打印结束

        //BatchPeeling Part
        Graph G_outBatchP, G_inBatchP;
        start_BatchP = clock();
        BatchPeelingAlgorithm(G_outq, G_inq, k, l, q, G_outBatchP, G_inBatchP);
        end_BatchP = clock();

        //打印生成的图中的所有点
        cout << G_outBatchP.nodes.size() << endl;
        sort(G_outBatchP.nodes.begin(), G_outBatchP.nodes.end(), cmp_weight);
        for (auto node : G_outBatchP.nodes) {
            cout << node.val << ' ';
        }
        cout << endl;

        Graph G_outIndBatch, G_inIndBatch;
        start_Ind = clock();
        vector<vector<PII>> Itmp;
        Itmp = I[k - 1];
        vector<int> points = IndexQueryAlgorithm(Itmp, q, l, location[k - 1]);
        PointsChangeToGraph(G_out, G_in, G_outIndBatch, G_inIndBatch, points);
        Graph G_outr, G_inr;
        BatchPeelingAlgorithm(G_outIndBatch, G_inIndBatch, k, l, q, G_outr, G_inr);
        end_Ind = clock();

        //打印生成的图中的所有点
        cout << G_outr.nodes.size() << endl;
        sort(G_outr.nodes.begin(), G_outr.nodes.end(), cmp_weight);
        for (auto node : G_outr.nodes) {
            cout << node.val << ' ';
        }
        cout << endl;
    }
}

int main() {
    srand(time(NULL));
    solve();
    return 0;
}
