#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef long double LD;

LD vecdot(int n, LD *u, LD *v){
	int i;
    LD x=0.0L;
	for(i=0;i<n;i++) x+=u[i]*v[i];
    return x;
}

void vecsub(int n, LD *u, LD *v, LD t, LD *rhs){
	int i;
    for(i=0;i<n;i++) rhs[i]=0.0L;
	for(i=0;i<n;i++) rhs[i]=u[i]-t*v[i];
}

void vecadd(int n, LD *u, LD *v, LD t, LD *rhs){
	int i;
    for(i=0;i<n;i++) rhs[i]=0.0L;
	for(i=0;i<n;i++) rhs[i]=u[i]+t*v[i];
}

LD vecnorm(int n, LD *v){
	return sqrt(vecdot(n,v,v));
}

void vecprintf(int n, LD *u){
    int i,m;
	if(n<3) m=n;
	else m=3;
    for(i=0;i<m;i++) printf("%15.15Le ",u[i]);
	if(n!=m) printf("...");
    printf("\n");
}
