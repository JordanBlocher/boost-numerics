#ifndef HESSENBERG_H
#define HESSENBERG_H

#include "deps.h"

mat hessenberg(mat A, int n, int m){
	unsigned int i,j,k;
    mat D(n,m),C(n,m),E(n,m),B(n,m),U(n,m),Uconj(n,m),temp(n,m);
    identity_matrix<val> I(n);
    zero_matrix<val> zeros(n,n);
    vec e1(m), v(m);
	val beta,alpha;
    for(i=0;i<e1.size();i++) e1(i)=0.0; e1(0)=1.0;
    for(k=0;k<A.size2()-2;k++){
        B=project(A,range(0,k+1),range(0,k+1));
        D.clear();D.resize(n-k-1,k+1);
        project(D,range(0,n-k-1),range(k,k+1))=project(A,range(k+1,n),range(k,k+1));
        e1.resize(D.size1());
        C=project(A,range(0,k+1),range(k+1,m));
        E=project(A,range(k+1,n),range(k+1,m));
        beta=-1.0*(D(0,k)/abs(D(0,k)))*norm_2(vec(column(D,k)));
        alpha=sqrt(2)/norm_2(vec(column(D,k))-beta*e1);
        v=alpha*(vec(column(D,k))-beta*e1);
        I.resize(B.size1());
        U.clear();Uconj.clear();
        project(U,range(0,k+1),range(0,k+1))=I;
        project(Uconj,range(0,k+1),range(0,k+1))=I;
        I.resize(v.size());
        mat vmat=I-outer_prod(v,conj(v));
        project(U,range(k+1,n),range(k+1,m))=I-outer_prod(v,conj(v));
        project(Uconj,range(k+1,n),range(k+1,m))=conj(I-outer_prod(v,conj(v)));
        temp=prod(A,conj(trans(U)));A=prod(U,temp);
    }
    return A;
}

#endif
