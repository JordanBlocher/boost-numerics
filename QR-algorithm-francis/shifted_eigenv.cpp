#include "./util/hessenberg.h"

#define N 2

mat triangle(mat A){
    int i,j;
    mat R(A.size1(),A.size2());
    R.clear();
    for(i=0;i<R.size1();i++) for(j=i;j<R.size2();j++)
        R(i,j)= j<i?0:A(i,j);
    return R;
}

vec eigenvector(mat A, val lambda){
    int i;
    vec x(A.size2()),y(A.size2());
    mat temp(A.size2(),A.size2());
    identity_matrix<val> I(A.size2());
    x(0)=1.0;
    temp=A-lambda*I;
    for(i=0;i<10;i++){//while(fabs((norm_2(y)/norm_2(x))-eigenvalue)>1e-8){
        lapack::gesv(temp,x);
        x=y/norm_2(y);
    }
    printf("Eigenvalue= %lf",norm_2(y)/norm_2(x));
    return x;
}

mat deflate(mat A, int n, int m){
    val lambda,i,j,k,alpha,beta;
    mat U(n,m),temp(n,m);
    vec x(n),e1(n),v(n),b(n);
    identity_matrix<val> I(n);
    lambda=A(n-1,n-1);
    printf(" lambda= %lf\n",lambda);
    //eigenv(k)=lambda;
    if(fabs(A(0,0))<1e-8) beta=1.0; else beta=A(0,0)/fabs(A(0,0));
    printf(" beta=  %15.8lf\n",beta);
    alpha=sqrt(2)/norm_2(x-beta*e1);
    printf(" alpha= %15.8lf\n",alpha);
    v=alpha*(x-beta*e1);
    printf(" v= %4.2s"," ");vecprintf(v);
    U=I-outer_prod(v,conj(v));
    printf(" U= \n");matprintf(U);
    temp=prod(A,conj(trans(U)));A=prod(U,temp);
    printf(" UAU*= \n");matprintf(A);printf("\n");
    return A; 
}

mat QR_shifted(mat input, int n, int m){
    int i=0,j,k;
    val z,alpha,beta,lambda;
    mat A(n,m),D(n,m),Q(n,m),R(n,m),temp(n,m),U(n,m),Uconj(n,m);
    identity_matrix<val> I(n);
    vec tau(n),eigenv(n),x(n),b(n),e1(n),v(n);
    e1.clear();e1(0)=1.0;
    A=hessenberg(input,n,m);
    printf(" Factorizing..\n");
    for(k=0;k<A.size2()-1;k++){
        I.resize(n);tau.resize(n);
        temp.resize(n,m);D.resize(n,m);R.resize(n,m);Q.resize(n,m);
        while(norm_2(vec(row(A,n-1)))>fabs(A(n-1,m-1))){i++;
            z=A(n-1,m-1);
            temp=A-z*I;
            A=temp;
            lapack::geqrf(A,tau);
            R=triangle(A);
            for(j=0;j<m;j++) if(R(j,j)>0.0||R(j,j)<0.0) D(j,j)=R(j,j)/fabs(R(j,j)); else D(j,j)=1.0;
            temp=prod(conj(D),R);
            R=temp;
            Q=I;
            lapack::ormqr('L','N',A,tau,Q,lapack::optimal_workspace());
            temp=prod(Q,D);
            Q=temp;
            A=prod(R,Q)+z*I;
        }
        m-=1;n-=1;
        eigenv(k)=A(n,m);
        printf(" A%d=\n",i);matprintf(A);
        //temp=deflate(A,A.size1(),A.size2());
        //A=temp;
        A.resize(m,n);
        printf(" A deflated%d=\n",i);matprintf(A);
    }
    return A;
}

int main(){
    int m,n,i,j;
	scanf("%d %d",&n,&m);
    mat A(n,m),R(n,m),test(n,m);
    vec eigenv(n);
    printf(" Reading an %d x %d matrix...\n",n,m);
    for(i=0;i<A.size1();i++) for(j=0;j<A.size2();j++) scanf("%lf", &A(i,j));
    printf(" A= \n");matprintf(A);
    A=QR_shifted(A,n,m);
    //for(i=0;i<R.size1();i++) eigenv(i)=R(i,i);
    //printf("Eigenvalues from QR-Algorithm:\n");vecprintf(eigenv);
    exit(0);
}

