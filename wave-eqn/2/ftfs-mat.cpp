#include "deps.h"
#include "print.h"

val f(val x){
    val sinx=1.0;
    for(int i=0;i<40;i++) sinx*=sin(M_PI*x);
    return sinx;
}

main(int argc, char* argv[]){ 
  
    int k,n,i,j;
    static int N, K, T;
    static int a=-1.0;

    std::string op=argv[1], prob=argv[3];
    std::istringstream tin(argv[2]);
	val t; tin>>t;
    if(prob=="a"){ N=25; K=5; T=1.0;}
    else if(prob=="b"){ N=2500; K=100; T=20.0;}


    val dx=1.0/K, dt=(val)T/N, xk, tn;
    val R=a*dt/dx;    // r = a*dt/dx
    
    static banded_matrix<val> U(K+2, K+2, 2, 4);
    vector<fortran_int_t> p(K+2);
    
    // Initialize matrix
    for(i=0; i<U.size1(); i++){
            U(i,i)=1.0-R; 
            k=std::max(i-1,0);
            U(k,k+1)=1.0*R;
    }
    // Boundary Conditions
    U(0,0)+=R;

    static vec u(K+2), w(K+2), z(K+2), v(K+2);

    // Initial conditions
    for(k=0;k<u.size();k++){
        xk=k*dx;
        u(k)=f(xk);
    }

    val alpha, beta; // Sherman-Morrison 
    w(0)=1.0; w(K+1)=-1.0; z(0)=R; z(K+1)=0.0;
    
    // Test Boundary Conditions
    if(op=="test"){
        mat wzt = outer_prod(w, z);
        mat B = U - wzt; 
        printf("Matrix Q1\n"); matprintf(B);
        printf("Matrix B\n"); matprintf(U);
        printf("Matrix wz^T\n"); matprintf(wzt);
    }
    
    lapack::gbtrf(U, p); // LU-decompostion
    lapack::gbtrs(U, p, w); //B^-1w
    alpha=1.0/(1.0-inner_prod(z, w)); // alpha = 1/(1-z^T(B^-1w)

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
