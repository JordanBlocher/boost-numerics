#include "deps.h"
#include "print.h"

#define NU 0.00001
#define N 400
#define T 2.0
#define K 100
#define r 0.5

val f(val x){
    return 3.0*sin(4.0*M_PI*x);
}

main(int argc, char* argv[]){ 

    int k,n,i,j;

    std::string op=argv[1], prob=argv[3];
    std::istringstream tin(argv[2]);
    val t; tin>>t;

	val dx=1.0/K, dt=(val)T/N, xk, tn, tmp;
    val R=NU*dt/dx;    // R = NU*dt/dx
    val alpha, beta;
    
    static banded_matrix<val> U(K+2, K+2, 2, 4);
    vector<fortran_int_t> p(K+2);    
    
    static vec u(K+2), w(K+2), z(K+2), v(K+2);

    if(prob=="a"){
    // Tridiagonal Matrix U
    for(i=0; i<U.size1(); i++){
            U(i,i)=1.0-R-2*r; 
            k=std::max(i-1,0);
            U(k,k+1)=r;
            U(k+1,k)=R+r;
    }
    // Boundary Conditions
    U(0,0)+=r; U(K+1,K+1)+=(R+r); 
   
    w(0)=1.0; w(K+1)=-1.0; z(0)=r; z(K+1)=-(R+r); // Sherman-Morrison
    }
   
    // Test Boundary Conditions
    if(op=="test"){ 
        mat wzt = outer_prod(w, z);
        mat B = U - wzt; 
        printf("Matrix Q1\n"); matprintf(B);
        printf("Matrix B\n"); matprintf(U);
        printf("Matrix wz^T\n"); matprintf(wzt);
    }
   
    // Initial conditions
    for(k=0;k<u.size();k++){
        xk=k*dx;
        u(k)=f(xk);
    }
   
    lapack::gbtrf(U, p); // LU-decompostion
    lapack::gbtrs(U, p, w); //B^-1w
    alpha=1.0/(1.0-inner_prod(z, w)); // alpha = 1/(1-z^T(b^-1w))

    for(n=0;n<=N;n++){
        tn=n*dt;
        lapack::gbtrs(U, p, u); // B^-1u
        beta=alpha*inner_prod(z, u); // beta = alpha*z^T(B^-1u)
        u+=beta*w;
        if(op=="pipe") plotu(u, tn, K);
        if(op=="approx" && tn==t) printu(u, tn, K);
    }
    return 0;
}
