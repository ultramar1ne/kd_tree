#include <iostream>
#include<algorithm>
#include<queue>
#include<cmath>
#include <thread>         // std::this_thread::sleep_for
#include <queue>          // std::priority_queue


using namespace std;

int N = 12345;   // N nodes
const int k = 5; // k-d

inline int murmurhash(int u) {
    size_t v = u * 3935559000370003845ul + 2691343689449507681ul;
    v ^= v >> 21;
    v ^= v << 37;
    v ^= v >> 4;
    v *= 4768777513237032717ul;
    v ^= v << 20;
    v ^= v >> 41;
    v ^= v << 5;
    return v;
}

struct kd_node {
    float Xs[k];
    string info;
    bool isVisited = false;

    kd_node *parent_ptr = nullptr;
    kd_node *left = nullptr;
    kd_node *right = nullptr;
};

struct dis_info_pair{
    float distance;
    string info;
    bool operator <(dis_info_pair a) const  {  return distance < a.distance; }
    bool operator >(dis_info_pair a) const  {  return distance > a.distance; }
};

void creat_kds(kd_node *kd_node_list) {
    for (int i = 0; i < N; i++) {
        int rand = murmurhash(i);
        for (int j = 0; j < k; j++) {
            kd_node_list[i].Xs[j] = (rand % 997 * murmurhash(j) % 997) % 997;
        }
    }
}

//todo: appro median causes imbalance?
float select_pivot(kd_node *kds, int n, int dim, int sample_size = 5) {
    if (n < sample_size) {
        return kds[0].Xs[dim] / 2.0 + kds[1].Xs[dim] / 2.0;
    }
    float avg = 0;
    for (int i = 0; i < sample_size; i++) {
        int rand = murmurhash(i) % n;
        avg += (kds[rand].Xs[dim]) / float(sample_size);
    }
    return avg;
}

void build_kd_tree(kd_node *left_kds, int left_n, kd_node *right_kds, int right_n, int height, kd_node *parent) {
    int dim = height % k;//todo: select dim on variance?

    if (left_n == 1) {
        parent->left = &left_kds[0];
        left_kds[0].parent_ptr = parent;
        left_n = 0;
    } else if (left_n >= 2) {
        std::sort(left_kds, left_kds + left_n, [&dim](kd_node x, kd_node y) { return x.Xs[dim] < y.Xs[dim]; });
        parent->left = &left_kds[left_n / 2];
        left_kds[left_n / 2].parent_ptr = parent;
        build_kd_tree(left_kds, left_n / 2, left_kds + left_n / 2 + 1, left_n - left_n / 2 - 1, height + 1,
                      parent->left);
    }

    if (right_n == 1) {
        parent->right = &right_kds[0];
        right_kds[0].parent_ptr = parent;
        right_n = 0;
    } else if (right_n >= 2) {
        std::sort(right_kds, right_kds + right_n, [&dim](kd_node x, kd_node y) { return x.Xs[dim] < y.Xs[dim]; });
        parent->right = &right_kds[right_n / 2];
        right_kds[right_n / 2].parent_ptr = parent;
        build_kd_tree(right_kds, right_n / 2, right_kds + right_n / 2 + 1, right_n - right_n / 2 - 1, height + 1,
                      parent->right);
    }

    if (left_n == 0 and right_n == 0) {
        return;
    }
}

float euclid_dis(float *X, float *Y) {
    float res = 0;
    for (int i = 0; i < k; i++) {
        res += (X[i] - Y[i]) * (X[i] - Y[i]);
    }
    return sqrt(res);
}


void nn_query_down(float *X, kd_node *root, int height, float *cur_min_dis) {
    if (!root){
        return;
    }
    int dim = height % k;
    if (*cur_min_dis > euclid_dis(X, root->Xs)) {
        *cur_min_dis = euclid_dis(X, root->Xs);
    }
    //root->isVisited = true;
    bool goLeft = false;

    if (X[dim] > root->Xs[dim]) {
        nn_query_down(X, root->right, height + 1, cur_min_dis);
    } else {
        goLeft = true;
        nn_query_down(X, root->left, height + 1, cur_min_dis);
    }

    if( abs(X[dim] - root->Xs[dim]) < *cur_min_dis ){
        if (goLeft){
            nn_query_down(X, root->right, height + 1, cur_min_dis);
        }else{
            nn_query_down(X, root->left, height + 1, cur_min_dis);
        }
    }
}

float bruteforce_nn(float *X, kd_node *kds) {
    float *resX = new float[k];
    float res = 99999;
    for (int i = 0; i < N; i++) {
        if (res > euclid_dis(X, kds[i].Xs)) {
            res = euclid_dis(X, kds[i].Xs);
            resX = kds[i].Xs;
        }
    }
    return res;
}

float* bruteforce_nn(float *X, kd_node *kds, int fixedSize) {
    float *resX = new float[fixedSize+1];
    for(int i = 0; i<fixedSize+1; i++){
        resX[i] = 9999999;
    }
    for (int i = 0; i < N; i++) {
        resX[fixedSize] = euclid_dis(X, kds[i].Xs);
        sort(resX,resX+fixedSize+1);
    }
    return resX;
}


float nn_query(float *X, kd_node *root) {
    int height = 0;
    float cur_min_dis = 999999; //todo: C++ type.inf?
    nn_query_down(X, root, height, &cur_min_dis);
    return cur_min_dis;
}

void nn_query_down(float *X, kd_node *root, int height, int fixedSize, std::priority_queue<dis_info_pair>* node_pairs);

std::priority_queue<dis_info_pair> nn_query(float *X, kd_node *root, int fixedSize) {
    int height = 0;
    float cur_min_dis = 999999; //todo: C++ inf?
    std::priority_queue<dis_info_pair> heap;
    nn_query_down(X, root, height, fixedSize, &heap);
    return heap;
}
void nn_query_down(float *X, kd_node *root, int height, int fixedSize, std::priority_queue<dis_info_pair> *node_heap) {
    if (!root){
        return;
    }
    int dim = height % k;
    dis_info_pair temp{euclid_dis(X, root->Xs), root->info};
    bool goLeft = false;

    if(node_heap->size()<fixedSize){
        node_heap->push(temp);
    }  else if(temp < node_heap->top()){
        node_heap->pop();
        node_heap->push(temp);
    }


    if (X[dim] > root->Xs[dim]) {
        nn_query_down(X, root->right, height + 1,fixedSize,node_heap);
    } else {
        goLeft = true;
        nn_query_down(X, root->left, height + 1,fixedSize,node_heap);
    }

    if( abs(X[dim] - root->Xs[dim]) < node_heap->top().distance ){
        if (goLeft){
            nn_query_down(X, root->right, height + 1,fixedSize,node_heap);
        }else{
            nn_query_down(X, root->left, height + 1,fixedSize,node_heap);
        }
    }
}


//copied from lc
void levelOrder(kd_node *root) {
    if (root == nullptr)
        return;
    kd_node *curNode = nullptr;
    vector<float *> perLevelNodes;
    vector<vector<int>> result;
    int count = 0;
    std::queue<kd_node *> nodeQueue;
    nodeQueue.push(root);
    int height = 0;
    while (!nodeQueue.empty()) // 每轮循环，都访问一层节点
    {
        perLevelNodes.clear();
        count = nodeQueue.size(); // 当前层有count个节点
        for (int i = 0; i < count; ++i) { // 循环弹出count次节点进行访问，也就是访问该层所有节点；同时将下一层所有节点压入栈
            curNode = nodeQueue.front();
            nodeQueue.pop();
            perLevelNodes.push_back(curNode->Xs);
            if (curNode->left != nullptr) // 左子节点非空压入栈
                nodeQueue.push(curNode->left);
            if (curNode->right != nullptr) // 右子节点非空压入栈
                nodeQueue.push(curNode->right);
        }
        cout << endl << "level :" << height++ << endl;

        for (int i = 0; i < perLevelNodes.size(); i++) {
            for (int j = 0; j < k; j++) {
                cout << perLevelNodes[i][j] << " ";
            }
            cout << endl;
        }
    }
    return;
}

int main() {
    kd_node *kds = new kd_node[N];
    creat_kds(kds);

    //build
    std::sort(kds, kds + N, [](kd_node x, kd_node y) { return x.Xs[0] < y.Xs[0]; });
    build_kd_tree(kds, N / 2, kds + N / 2 + 1, N - N / 2 - 1, 1, &kds[N / 2]);
    cout << "built, nice" << endl;

    //level order
    //levelOrder(&kds[N/2]);

    //1-NN on  k-d tree
    auto *X = new float[k];
    X[0] = 152;
    X[1] = 779;
    X[2] = 100;
    X[3] = 330;

    for(int i = 0; i<1234; i++){
        X[0] = (murmurhash(i)* murmurhash(i)/213)%997;
        X[1] = (murmurhash(i)* murmurhash(i)/333)%997;
        X[2] = (murmurhash(i)* murmurhash(i)/444)%997;
        X[3] = (murmurhash(i)* murmurhash(i)/555)%997;
        if(nn_query(X,&kds[N/2]) != bruteforce_nn(X, kds)){
            cout<<"Wrong Answer!"<<bruteforce_nn(X, kds)<<"  "  <<nn_query(X,&kds[N/2])<< endl  << X[0]<<" "<<X[1]<<" "<<X[2]<<" "<<X[3];
        }
    }
    cout<<"passed NN test"<<endl;

    for(int i = 0; i<5678; i++){
        X[0] = (murmurhash(i)* murmurhash(i)/213)%997;
        X[1] = (murmurhash(i)* murmurhash(i)/333)%997;
        X[2] = (murmurhash(i)* murmurhash(i)/444)%997;
        X[3] = (murmurhash(i)* murmurhash(i)/555)%997;

        std::priority_queue<dis_info_pair> A = nn_query(X,&kds[N/2],5);
        float* B = bruteforce_nn(X,kds,5);
        float* AA = new float[5];
        for(int i = 4; i>=0; i--){
            float a = A.top().distance;
            A.pop();
            AA[i] = a;
            if (AA[i]!=B[i]){
                cout<<"failed knn"<<endl;
                for(int j = 0; j<5; j++){
                    cout<<AA[j]<<" "<<B[j]<<endl;
                }
            }
        }

    }

    cout<<"passed kNN test"<<endl;
    return 0;
}
