

# current result    Aug 3

Qun Lou 

## seq Knn
### 1. build kd_tree from points list
Now: use std::nth_element in every layer
Todo: 

  	sort the points list on different dimension in advance?  
  	select dim on k-dim tree according to "variance" instead of round-robin?    

### 2. k-nearest-neighbor search

```cpp
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
```
I didn't tell the difference between "leaf" nodes and "parent" node. I am afraid that the current realization is too naïve to be parallelized.

### 3. correctness&speed test

0. my_knn
1. CGAL
2. brute force

####   	dataset:

4d points: https://archive.ics.uci.edu/ml/datasets/iris

3d points: https://archive.ics.uci.edu/ml/datasets/haberman's+survival

manual generated:  randomly distributed or highly repeated points 

#### 	results example:

N = 1023  K-D = 4  K-NN = 7

```
cgal time 0.002647
my_knn time 5.2e-05
brute force time 0.000224
KNN correctness: correct
```

### parallel try
I simply add "cilk_spawn" and get a much slower result now, even without lock to guarantee correct ...