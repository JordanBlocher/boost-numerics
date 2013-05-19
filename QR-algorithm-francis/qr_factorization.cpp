#include "./util/deps.h"

#define N 20

mat triangle(mat A){
    int i,j;
    mat R(A.size1(),A.size2());
    R.clear();
    for(i=0;i<R.size1();i++) for(j=i;j<R.size2();j++)
        R(i,j)= j<i?0:A(i,j);
    return R;
}

mat QR_factorization(mat A, int n, int m){
    int i,j;
    mat D(n,m),Q(n,m),R(n,m),Qhat(n,m),Rhat(n,m),H(n,m),temp(n,m);
    identity_matrix<val> I(n);
    vec tau(n),v(n);
    printf(" Factorizing..\n");
    for(i=0;i<N+1;i++){
        D.clear();
	    lapack::geqrf(A,tau);
        //printf(" After call to lapack:\n A%d=\n",i+1);matprintf(A);
        //printf("tau=");vecprintf(tau);
        R=triangle(A);
        for(j=0;j<m;j++) if(abs(R(j,j))>0.0) D(j,j)=R(j,j)/abs(R(j,j)); else D(j,j)=1.0;
        temp=prod(conj(D),R);
        R=temp;
        //printf("Rhat%d=\n",i+1);matprintf(R);
        Q=I;
	    lapack::ormqr('L','N',A,tau,Q,lapack::optimal_workspace());
        temp=prod(Q,D);
        Q=temp;
        //printf("Qhat%d=\n",i+1);matprintf(Q);
        A=prod(R,Q);
        if(i==1||i==9){
            printf(" A%d=\n",i+1);matprintf(A);
            printf("    .\n    .\n    .\n");
        }
    }
        printf(" A%d=\n",N);matprintf(A);
    return A;
}

int main(){
    int m,n,i,j;
	scanf("%d %d",&n,&m);
    mat A(n,m),T(n,m);
    vec eigenv(n);
    printf(" Reading an %d x %d matrix...\n",n,m);
    for(i=0;i<A.size1();i++) for(j=0;j<A.size2();j++) scanf("%lf", &A(i,j));
    printf(" A= \n");matprintf(A);
    T=QR_factorization(A,n,m);
    exit(0);
}

