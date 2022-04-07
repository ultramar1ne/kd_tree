/*
todo:   1. dim on varinace?
        2. median partition
 */

#include <iostream>
#include<algorithm>
#include<queue>
using namespace std;

int N=129;   // N nodes
const int k = 3; // k-d

inline int murmurhash(int u )
{
    size_t v = u * 3935559000370003845ul + 2691343689449507681ul;
    v ^= v >> 21;
    v ^= v << 37;
    v ^= v >>  4;
    v *= 4768777513237032717ul;
    v ^= v << 20;
    v ^= v >> 41;
    v ^= v <<  5;
    return v;
}

struct kd_node{
    float Xs[k];
    kd_node* left = nullptr;
    kd_node* right = nullptr;
};

void creat_kds(kd_node* kd_node_list){
    for(int i = 0; i<N; i++){
        int rand = murmurhash(i);
        for(int j = 0; j<k; j++){
            kd_node_list[i].Xs[j]=(rand%997 * murmurhash(j)%997)%997;
        }
    }
}

//todo:
float select_pivot(kd_node* kds, int n, int dim, int sample_size = 5){
    if(n<sample_size){
        return kds[0].Xs[dim]/2.0 + kds[1].Xs[dim]/2.0;
    }
    float avg =0;
    for(int i = 0; i<sample_size; i++){
        int rand = murmurhash(i)%n;
        avg += (kds[rand].Xs[dim])/float(sample_size);
    }
    return avg;
}



void build_kd_tree(kd_node* left_kds, int left_n, kd_node* right_kds, int right_n, int height, kd_node* parent){
    int dim = height %k;//todo: select dim on variance?

    if (left_n==1){
        parent->left=&left_kds[0];
        left_n=0;
    }else if(left_n>=2){
        std::sort(left_kds,left_kds+left_n, [&dim](kd_node x, kd_node y){return x.Xs[dim] < y.Xs[dim];});
        parent->left=&left_kds[left_n/2];
        build_kd_tree(left_kds,left_n/2,left_kds+left_n/2+1,left_n-left_n/2-1,height+1,parent->left);
    }

    if (right_n==1){
        parent->right=&right_kds[0];
        right_n=0;
    }else if(false){
        float avg = select_pivot(right_kds,right_n,dim);
        auto mid = std::partition(right_kds,right_kds+right_n,[&dim,&avg](kd_node x){return x.Xs[dim]<=avg;});
        parent->right=mid;
        int index_r = distance(right_kds,mid);
        build_kd_tree(right_kds,index_r,right_kds+index_r+1,right_n-index_r-1,height+1,parent->right);
    }else if(right_n>=2){
        std::sort(right_kds,right_kds+right_n, [&dim](kd_node x, kd_node y){return x.Xs[dim] < y.Xs[dim];});
        parent->right=&right_kds[right_n/2];
        build_kd_tree(right_kds,right_n/2,right_kds+right_n/2+1,right_n-right_n/2-1,height+1,parent->right);
    }

    if(left_n == 0 and right_n == 0){
        return;
    }
}

//copied from lc
void levelOrder(kd_node* root) {
    if (root == nullptr)
        return ;
    kd_node* curNode = nullptr;
    vector<float*> perLevelNodes;
    vector<vector<int>> result;
    int count = 0;
    std::queue<kd_node*> nodeQueue;
    nodeQueue.push(root);
    int height = 0;
    while (!nodeQueue.empty()) // 每轮循环，都访问一层节点
    {
        perLevelNodes.clear();
        count = nodeQueue.size(); // 当前层有count个节点
        for (int i = 0; i < count; ++i)
        { // 循环弹出count次节点进行访问，也就是访问该层所有节点；同时将下一层所有节点压入栈
            curNode = nodeQueue.front();
            nodeQueue.pop();
            perLevelNodes.push_back(curNode->Xs);
            if (curNode->left != nullptr) // 左子节点非空压入栈
                nodeQueue.push(curNode->left);
            if (curNode->right != nullptr) // 右子节点非空压入栈
                nodeQueue.push(curNode->right);
        }
        cout<<endl<<"level :" <<height++ <<endl;

        for(int i = 0; i< perLevelNodes.size(); i++){
            for(int j = 0; j<k;j++){
                cout<<perLevelNodes[i][j]<<" ";
            }
            cout<<endl;
        }
    }
    return;
}

int main(){
    kd_node* kds=new kd_node[N];
    creat_kds(kds);

    //build
    std::sort(kds,kds+N, [](kd_node x, kd_node y){return x.Xs[0] < y.Xs[0];});
    build_kd_tree(kds,N/2,kds+N/2+1,N-N/2-1,1,&kds[N/2]);
    cout<<"built, nice"<<endl;

    //level order
    levelOrder(&kds[N/2]);
    return 0;
}
