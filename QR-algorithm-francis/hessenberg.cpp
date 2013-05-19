#include "./util/deps.h"

mat hessenberg(mat A, int n, int m){
	unsigned int i,j,k;
    mat D(n,m),C(n,m),E(n,m),B(n,m),U(n,m),Uconj(n,m),temp(n,m);
    identity_matrix<val> I(n);
    zero_matrix<val> zeros(n,n);
    vec e1(m), v(m);
	val beta,alpha;
    e1.clear();e1(0)=1.0;
    printf("\n Reducing....\n");
    for(k=0;k<A.size2()-2;k++){
        printf(" Iteration %d\n",k+1);
        B=project(A,range(0,k+1),range(0,k+1));
        D.clear();D.resize(n-k-1,k+1);
        project(D,range(0,n-k-1),range(k,k+1))=project(A,range(k+1,n),range(k,k+1));
        printf(" d= %4.2s"," ");vecprintf(vec(column(D,k)));
        e1.resize(D.size1());
        C=project(A,range(0,k+1),range(k+1,m));
        E=project(A,range(k+1,n),range(k+1,m));
        beta=-1.0*(D(0,k)/abs(D(0,k)))*norm_2(vec(column(D,k)));
        printf(" beta=  %15.8lf\n",beta);
        alpha=sqrt(2)/norm_2(vec(column(D,k))-beta*e1);
        printf(" alpha= %15.8lf\n",alpha);
        v=alpha*(vec(column(D,k))-beta*e1);
        printf(" v= %4.2s"," ");vecprintf(v);
        I.resize(B.size1());
        U.clear();Uconj.clear();
        project(U,range(0,k+1),range(0,k+1))=I;
        project(Uconj,range(0,k+1),range(0,k+1))=I;
        I.resize(v.size());
        mat vmat=I-outer_prod(v,conj(v));
        printf(" U= \n");matprintf(vmat);
        project(U,range(k+1,n),range(k+1,m))=I-outer_prod(v,conj(v));
        project(Uconj,range(k+1,n),range(k+1,m))=conj(I-outer_prod(v,conj(v)));
        temp=prod(A,conj(trans(U)));A=prod(U,temp);
        printf(" UAU*= \n");matprintf(A);printf("\n");
    }printf("\n");
    return A;
}

int main(){
    int m,n,i,j;
	scanf("%d %d",&n,&m);
    mat A(n,m), hessenberg_form(n,m);
    printf(" Reading an %d x %d matrix...\n",n,m);
    for(i=0;i<A.size1();i++) for(j=0;j<A.size2();j++) scanf("%lf", &A(i,j));
    printf(" A= \n");matprintf(A);
    hessenberg_form = hessenberg(A,n,m);
    printf(" Matrix reduced to Hessenberg form:\n");matprintf(hessenberg_form);
    exit(0);
}
