#ifndef MATRIX_H:q
#define MATRIX_H

#include "vector.h"

using namespace std;

template<typename TP, int N>
class Matrix
{
    public:
        
        matrix<TP> data;

        Matrix(void);
        Matrix(Matrix m);

        Matrix operator*(Matrix m);
        Vector operator*(Vector v);

};

LD matinfnorm(int n,LD **A){
	int i,j;
	LD r,s;
	r=0;
	for(i=0;i<n;i++){
		s=0;
		for(j=0;j<n;j++) s+=fabs(A[i][j]);
		if(s>r) r=s;
	}
	return r;
}

void matdot(int n, LD *u, LD **A, LD *rhs){
	int i,j;
	for(i=0;i<n;i++) rhs[i]=0.0L;
	for(i=0;i<n;i++) for(j=0;j<n;j++) rhs[i]+=A[i][j]*u[j];
}

void matmult(int n, LD **A, LD (*B)(int,int)){
    int i,j,k;
    LD **X=(LD**)malloc(n*sizeof(LD*));
    for(i=0;i<n;i++) X[i]=(LD*)malloc(n*sizeof(LD)); 
    for(i=0;i<n;i++) for(j=0;j<n;j++){
        X[i][j]=0.0L;
        for(k=0;k<n;k++) X[i][j]+=A[j][k]*B(k,j);
    }
}

void transpose(int n, LD **A, LD **rhs){
    int i,j;
    for(i=0;i<n;i++) for(j=0;j<n;j++) rhs[i][j]=0.0L;
    for(i=0;i<n;i++) for(j=0;j<n;j++) rhs[i][j]=A[j][i];
}

void matprintf(int n, LD **A){
	int i,m;
	if(n>5) m=5;
	else m=n;
	for(i=0;i<m;i++) vecprintf(n,A[i]);
    if(n!=m) printf("\n .\n. \n.");
}

#endif
