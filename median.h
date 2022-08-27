#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <cmath>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <algorithm> 
using namespace std;


float* ext; 
float* original_A;
bool isInitialized;

int leq_threshold;
int seq_quicksort_threshold; 
int seq_scan_threshold; 
int default_sample_size; 

inline size_t murmurhash(size_t u )
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



template <typename T>
// Using avg. instead of real median pivot
T select_pivot(T* A, size_t n, int sample_size = default_sample_size){
	T pivot = 0;
	for(size_t i = 0; i<sample_size ; i++){ 
		pivot +=  A[(murmurhash(i))%n]/sample_size;  
	}
	return pivot;
}



template <typename T>
// todo: delayed sequence??
void set_flag(T* A, size_t n, T pivot, int* flag){
	cilk_for (int i = 0; i < n; i++){
		flag[i] = (A[i]<=pivot) ? 0 : 1;
	}
	return ;
}


template <typename T>   //decrease&conquer
void prefix_sum(int* flag, size_t n, T* Out){
	if (n<seq_scan_threshold){   //todo: ????? it seems non parallel is the fastest?
		Out[0]=flag[0];
		for( int i = 1; i<n; i++){
			Out[i]=flag[i]+Out[i-1];
		}
		return ;
	}

	int* B = new int[size_t(n/2)];
	int* C = new int[n];

	cilk_for (size_t i = 0; i <size_t(n/2); i++){
		B[i] = flag[2*i] + flag[2*i+1];
	}
	
	prefix_sum(B,size_t(n/2),C); 

	cilk_for (size_t i = 0; i<n; i++){
		if (i%2==0){
			Out[i] =  C[(i-1)/2]+flag[i];
		}else{
			Out[i] = C[i/2];
		}
	}
//todo: How to avoding malloc B C?   hwo to GC?
	return;
}

template <typename T>
size_t para_partition(T* A, size_t n, T pivot) {
	int* flag = new int[n];
	set_flag(A,n,pivot,flag);  //todo:  delay sequence??

	size_t* ps = new size_t[n];
	prefix_sum(flag,n,ps);
	int offset = (A - original_A) ;

	int num_gt = ps[n-1];  
	int num_leq = n- num_gt;
	if ( num_gt < leq_threshold || num_leq < leq_threshold){
		return num_leq;
	}

	if(ps[0] == 0){
		ext[offset] = A[0];
	}else{
		ext[offset+num_leq] = A[0];
	}


	//done:  cilk for here
	cilk_for (size_t i = 1; i<n; i++){  
		if (ps[i]!=ps[i-1]){
			ext[offset+ps[i-1]+num_leq] = A[i];
		}else{
			ext[offset+ i-ps[i]] = A[i];
		}
	}

	cilk_for(int i = 0; i<n; i++){
		A[i]=ext[offset+i];
	}
	return num_leq;
}

template <typename T>
void median(T* A, int n, int k) {
	// using global var instead of allocating space recursively 
	if (!isInitialized){
		original_A = A;
		ext = new T[n];
		isInitialized = true;
	}
    
    int n_leq =  para_partition(A,n,select_pivot(A,n));
    if(n_leq < k ){
        if(n_leq < seq_quicksort_threshold){
            nth_element(A+n_leq,A+k,A+n);
            return;
        }
	    median(A+n_leq,n-n_leq,k-n_leq);
    }else if (n_leq > k){
        if(n_leq < seq_quicksort_threshold){
            nth_element(A,A+k,A+n_leq);
            return;
        }
	    median(A,n_leq,k);
    }else{
        return;
    }
}
