#include "./util/hessenberg.h"

mat triangle(mat A){
    int i,j;
    mat R(A.size1(),A.size2());
    R.clear();
    for(i=0;i<R.size1();i++) for(j=i;j<R.size2();j++)
        R(i,j)= j<i?0:A(i,j);
    return R;
}

vec QR_shifted(mat input, int n, int m){
    int i=0,j,k,M=m;
    val z,lambda;
    mat A(n,m),Q(n,m),R(n,m),temp(n,m);
    identity_matrix<val> I(n);
    vec tau(n),eigenv(n);
    A=hessenberg(input,n,m);
    printf(" Hessenberg form:\n");matprintf(A);
    printf(" Factorizing..\n");
    for(k=0;k<M;k++){
        I.resize(n);tau.resize(n);
        temp.resize(n,m);R.resize(n,m);Q.resize(n,m);
        while(norm_2(vec(row(A,n-1)))>abs(A(n-1,m-1))){i++;
            z=A(n-1,m-1);
            temp=A-z*I;
            A=temp;
//            lapack::geqrf(A,tau);
            R=triangle(A);
            Q=I;
//            lapack::ormqr('L','N',A,tau,Q,lapack::optimal_workspace());
            A=prod(R,Q)+z*I;
        }
        m-=1;n-=1;
        eigenv(k)=A(n,m);
        if(k<M-1){ printf(" A%d=\n",i+1);matprintf(A);}
        if(k<M-2) printf(" Deflating...\n");
        A.resize(m,n);
    }
    return eigenv;
}

int main(){
    double inr,ini;
    int m,n,i,j;
	scanf("%d %d",&n,&m);
    mat A(n,m),R(n,m);
    vec eigenv(n);
    printf(" Reading an %d x %d matrix...\n",n,m);
    for(i=0;i<A.size1();i++) for(j=0;j<A.size2();j++){ scanf("%lf", &A(i,j));}
    printf(" A= \n");matprintf(A);
    eigenv=QR_shifted(A,n,m);
    printf("Eigenvalues from Shifted QR-Algorithm:\n");vecprintf(eigenv);
    exit(0);
}

