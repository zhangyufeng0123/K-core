#include<bits/stdc++.h>

using namespace std;

const string filename = "data.txt";

struct Vertex{
    int val;    //这个点的编号
    vector<int> connectV;   //与这个点连接的点
    int degrees; //该点的度数
};

//the whole graph
struct Graph{
    vector<int> vertex; //标注顺序下是哪个点的位置
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
        G_out.vertex.push_back(temp.val);
        G_out.location[temp.val] = G_out.location.size();
        G_out.nodes.back().degrees = G_out.nodes.back().connectV.size();
    }

    //存储入度
    Vertex tmp;
    //初始化入度
    for(int i = 0; i < row; i++){
        tmp.val = G_out.nodes[i].val;
        G_in.nodes.push_back(tmp);
        G_in.location[tmp.val] = i;
        G_in.vertex.push_back(tmp.val);
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
    return a.degrees < b.degrees;
}

//对点的入度进行排序
static bool cmp_in(Vertex &a, Vertex &b){
    return a.degrees < b.degrees;
}

//Algorithm 1:QueryDcore
//k是出度，l是入度
Graph QueryDcore(Graph G_out, Graph G_in, int k, int l){
    //先对图中的点的度数排序

    int flag = 1;
    while(1){
        flag = 1;

        if(flag)    break;
    }
}

int main(){
    Graph G_out, G_in;
    bool j = readData(G_out, G_in);
    //cout << 1 << endl;

    return 0;
}