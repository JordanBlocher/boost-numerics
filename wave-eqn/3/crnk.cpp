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

    if(prob=="a"){ N=25; K=20; T=1.0;}
    else if(prob=="b"){ N=2500; K=100; T=20.0;}

	val dx=1.0/K, dt=(val)T/N, xk, tn, tmp;
    val R=a*dt/dx;    // R = a*dt/dx
    
    static banded_matrix<val> U(K, K, 2, 4);
    vector<fortran_int_t> p(K);    
    
    // Tridiagonal Matrix U
    for(i=0; i<U.size1(); i++){
            U(i,i)=1.0; 
            k=std::max(i-1,0);
            U(k,k+1)=R/4.0;
            U(k+1,k)=-R/4.0;
    }
    // Boundary Conditions
    U(0,0)+=R/4.0; U(K-1,K-1)-=R/4.0; 
    
    static vec u(K), w(K), z(K), r(K);
  
    // Initial conditions
    for(k=0;k<u.size();k++){
        xk=k*dx;
        u(k)=f(xk);
    }
    if(op=="energy") {printf("t=0.0\n"); printenergy(u);}

    val alpha, beta; // Sherman-Morrison 
    w(0)=1.0; w(K-1)=-1.0; z(0)=R/4.0; z(K-1)=R/4.0;
    
    // Test Boundary Conditions
    if(op=="test"){ 
        mat wzt = outer_prod(w, z);
        mat B = U - wzt; 
        printf("Matrix U\n"); matprintf(U);
        printf("Matrix B\n"); matprintf(B);
        printf("Matrix wz^T\n"); matprintf(wzt);
    }

    lapack::gbtrf(U, p); // LU-decompostion
    lapack::gbtrs(U, p, w); //B^-1w
    alpha=1.0/(1.0-inner_prod(z,w)); // alpha = 1/(1-z^T(b^-1w))
    
    if(op=="test") printf("alpha=%f\n",alpha);

    #define slice(k) u(((k)+K)%K)

    for(n=0;n<=N;n++){
        for(k=0;k<U.size1();k++) r(k)=(R/4.0)*(slice(k-1)-slice(k+1))+u(k);
        tn=n*dt;
        lapack::gbtrs(U, p, r); // B^-1r^n
        beta=alpha*inner_prod(z, trans(r)); // beta = alpha*z^T(B^-1r^n)
        if(op=="test"){ printf("beta=%f\n",beta);  printavg(u);}
        u=r+beta*w;
        if(op=="pipe") plotu(u, tn, K);
        if(op=="approx" && tn==t)printu(u, tn, K);
    }
    
    if(op=="energy"){ printf("t=%d\n",T); printenergy(u);}
    
    return 0;
}
