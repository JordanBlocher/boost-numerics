#include "deps.h"
#include "output.h"
#include "print.h"

#define K 10
#define T 0.25
#define N 1
#define NU 0.00005
#define NMOD 200 

// Note: 
// Test program for gbsv used for parameters, 
// followed from svn boost: https://svn.boost.org/svn/boost/sandbox/numeric_bindings/libs/numeric/bindings/lapack/test/ublas_gbsv.cpp

val f(val x){
    return x*x*x*(-48.0*x*x + 112.0*x - 64.0);
}

main(int argc, char* argv[]){ 
    int k,n,i,j;
    std::string op=argv[1]; // command line argument
	val dx=1.0/K, dt=T/N, dx4=dx*dx*dx*dx, xk, tn; // discretization variables
    val rho=NU*dt/dx4;    // rho = nu*dt/dx^4
    static vec u(K-1), v(K+1);
    
    // We allocate 2 lower & 4 upper diagonal, according to the example
    static banded_matrix<val> U(K-1, K-1, 2, 4);
    vector<fortran_int_t> p(K-1);
    
    // Initialize matrix
    for(i=0; i<U.size1(); i++){
            U(i,i)=1.0+6.0*rho; 
            k=std::max(i-1,1);
            U(k,k-1)=U(k-1,k)=-4.0*rho;
            U(k,k+1)=U(k+1,k)=-4.0*rho;
            k=std::max(i-2,2);
            U(k,k-2)=U(k-2,k)=1.0*rho;
            U(k,k+2)=U(k+2,k)=1.0*rho;
        }
    // Boundary Conditions
    U(0,0)-=1.0*rho;
    U(K-2,K-2)-=1.0*rho;
    if(op=="matrix"){ printf("Pentadiagonal Matrix\n"); matprintf(U);}

    // Initial conditions
    for(k=0;k<=K;k++){
        xk=k*dx;
        v(k)=f(xk);
    }
    u=subrange(v,1,K);
    //printf("Original Vector\n"); vecprintf(u,dx);

    lapack::gbtrf(U, p); // LU-decompostion
    for(n=0;n<=N;n++){
        tn=n*dt;
        lapack::gbtrs(U, p, u); // Solve
        if(op=="plot0") plot0(u, tn, K-1, N);
        if(op=="plot1") plot1(u, tn, K-1, N);
        if(op=="plot3d" && (n-1)%NMOD==0) plot3d(u, tn, K-1, N);
        if(op=="approx") output(tn, 0.5, u(K/2-1), K-1, N);
    }
    return 0;
}
