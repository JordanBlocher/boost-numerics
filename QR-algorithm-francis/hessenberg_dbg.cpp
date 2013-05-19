#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/xpressive/detail/core/matcher/action_matcher.hpp>
#include "./util/print.h"

using namespace boost::numeric::ublas;

typedef double val;

typedef boost::numeric::ublas::matrix<double> mat;
typedef boost::numeric::ublas::vector<double> vec;

int main(){
	unsigned int i,j,k;
    int m,n;
	scanf("%d %d",&n,&m);
    mat A(n,m),D(n,m),C(n,m),E(n,m),B(n,m),U(n,m),Uconj(n,m),temp(n,m);
    identity_matrix<val> I(n);
    zero_matrix<val> zeros(n,n);
    vec e1(m), v(m);
	val beta,alpha;
    printf("Reading an %d x %d matrix...\n",n,m);
    for(i=0;i<A.size1();i++) for(j=0;j<A.size2();j++) scanf("%lf", &A(i,j));
    for(i=0;i<e1.size();i++) e1(i)=0.0; e1(0)=1.0;
    printf("A= \n");matprintf(A);
    for(k=0;k<A.size2()-1;k++){
        printf("Iteration %d\n\n",k);
        B=project(A,range(0,k+1),range(0,k+1));
        printf("B= \n");matprintf(B);
        D.clear();D.resize(n-k-1,k+1);
        project(D,range(0,n-k-1),range(k,k+1))=project(A,range(k+1,n),range(k,k+1));
        e1.resize(D.size1());
        printf("D= \n");matprintf(D);
        C=project(A,range(0,k+1),range(k+1,m));
        printf("C= \n");matprintf(C);
        E=project(A,range(k+1,n),range(k+1,m));
        printf("E= \n");matprintf(E);
        beta=-1.0*(D(0,k)/fabs(D(0,k)))*norm_2(vec(column(D,k)));
        printf("beta= %lf\n",beta);
        alpha=sqrt(2)/norm_2(vec(column(D,k))-beta*e1);
        printf("alpha= %lf\n",alpha);
        v=alpha*(vec(column(D,k))-beta*e1);
        printf("v= \n");vecprintf(v);
        I.resize(B.size1());
        printf("I= \n");matprintf(I);
        U.clear();Uconj.clear();
        project(U,range(0,k+1),range(0,k+1))=I;project(Uconj,range(0,k+1),range(0,k+1))=I;
        I.resize(v.size());
        printf("I= \n");matprintf(I);
        mat vmat=I-outer_prod(v,conj(v));
        printf("vmat= \n");matprintf(vmat);
        project(U,range(k+1,n),range(k+1,m))=I-outer_prod(v,conj(v));
        project(Uconj,range(k+1,n),range(k+1,m))=conj(I-outer_prod(v,conj(v)));
        printf("U= \n");matprintf(U);
        printf("Uconj= \n");matprintf(Uconj);
        temp=prod(A,conj(trans(U)));A=prod(U,temp);
        printf("UAU*= \n");matprintf(A);
    }
    exit(0);
}
