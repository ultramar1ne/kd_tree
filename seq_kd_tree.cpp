#include <iostream>
#include<algorithm>
#include<queue>
#include<cmath>
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

    bool isVisited = false;
    kd_node* parent_ptr = nullptr;

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

//todo: appro median causes imbalance?
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
        left_kds[0].parent_ptr=parent;
        left_n=0;
    }else if(left_n>=2){
        std::sort(left_kds,left_kds+left_n, [&dim](kd_node x, kd_node y){return x.Xs[dim] < y.Xs[dim];});
        parent->left=&left_kds[left_n/2];
        left_kds[left_n/2].parent_ptr=parent;
        build_kd_tree(left_kds,left_n/2,left_kds+left_n/2+1,left_n-left_n/2-1,height+1,parent->left);
    }

    if (right_n==1){
        parent->right=&right_kds[0];
        right_kds[0].parent_ptr=parent;
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
        right_kds[right_n/2].parent_ptr=parent;
        build_kd_tree(right_kds,right_n/2,right_kds+right_n/2+1,right_n-right_n/2-1,height+1,parent->right);
    }

    if(left_n == 0 and right_n == 0){
        return;
    }
}

float euclid_dis(float* X, float* Y){
    float res = 0;
    for(int i = 0; i<k; i++){
        res+=(X[i]-Y[i])*(X[i]-Y[i]);
    }
    return sqrt(res);
}

//todo: after reaching a leaf node, we may need to trace back if
const bool NeedTrackBack = true;
bool nn_query_down(float* X, kd_node* root, int* height, float* cur_min_dis, kd_node** cur_node){
    int dim = *height % k;
    if(*cur_min_dis>= euclid_dis(X,root->Xs)){
        *cur_min_dis= euclid_dis(X,root->Xs);
    }
    *cur_node = root;

    if(X[dim]>root->Xs[dim]){
        if(root->right){
            *height+=1;
            return nn_query_down(X, root->right,height, cur_min_dis, cur_node);
        }
    }else{
        if(root->left){
            *height+=1;
            return nn_query_down(X, root->left,height, cur_min_dis, cur_node);
        }
    }

    // if the hyper plant "intersection"s the hyper ball :
    if (root->parent_ptr && ( abs(X[abs(dim-1)%k] - root->parent_ptr->Xs[abs(dim-1)%k]) > *cur_min_dis ) ) {
        cout<<"Need Trace Back!"<<abs( X[abs(dim-1)%k] - root->parent_ptr->Xs[abs(dim-1)%k] )<< endl;
        return NeedTrackBack;
    } else {
        return false;
    }
}

float* nn_query(float* X, kd_node* root ){
    kd_node* cur_node = nullptr;
    int height = 0; float cur_min_dis = 999999;
    nn_query_down(X,root,&height,&cur_min_dis,&cur_node);
    return cur_node->Xs;
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
    cout<<"built, nice"<<endl<<"Level Traverse:";

    //level order
    levelOrder(&kds[N/2]);

    //1-NN on  k-d tree
    cout<<endl<<"NN query"<<endl;
    auto * X = new float [k]; for(int i = 0; i<k; i++){ X[i] = murmurhash(i)%666;}

    float* res = nn_query(X,&kds[N/2]);
    for(int i = 0; i<k; i++){ cout<<"dim "<<i<<": "<<res[i]<<" "<<endl;}
    cout<< "euclid dis: "<< euclid_dis(X,res)<<endl;

    return 0;
}
