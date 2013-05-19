#include "deps.h"
#include "output.h"
#include "print.h"

#define MU 0.1
#define NU 0.00005
#define K 2000
#define T 0.25
#define N 10000
#define NMOD 2000

// Function for the initial condition
val f(val x){
	return -48*pow(x,5) + 112*pow(x,4) - 64*pow(x,3);
}

main(int argc, char* argv[]){ 

    std::string op=argv[1]; // command line argument

    // Variable decloation
    static vec v(K+1), u(K-1), w(K-1);
    val dt = T/N;
    val dx = 1.0/K;
    val gam = MU*dt/(dx*dx);
    val rho = NU*dt/(dx*dx*dx*dx);
	val ar = dt/(2*dx);
	
    int i,j,k,n;
    val xk, tn;
    val un, up;

   static banded_matrix<val, row_major> A(K-1, K-1, 2, 4);

    // Initialize matrix
    for(i=0; i< A.size1(); i++){
            A(i,i)=1.0+6.0*rho-2*gam; 
            k=std::max(i-1,1);
            A(k,k+1)=A(k+1,k)=-4.0*rho+gam;
            k=std::max(i-2,2);
            A(k,k+2)=A(k+2,k)=A(k,k-2)=A(k-2,k)=1.0*rho;
        }
   // Boundary Conditions
   A(0,0)-=1.0*rho;
   A(K-2,K-2)-=1.0*rho;
   if(op=="matrix"){ printf("Pentadiagonal Matrix\n"); matprintf(A);}

    //Creating a permutation storage vector for the lapack functions
    vector<fortran_int_t> p(K+1);
	lapack::gbtrf(A, p); // LU-decompostion
	
    //Initialize u
	for (k=0;k<=K;k++){
        xk=k*dx;
        v(k)=f(xk);
    }
    u=subrange(v,1,K);
    //if(op=="approx") {printf("Initial Vector\n");vecprintf(u,dx);}

    for(n=0;n<=N;n++){
        tn=n*dt;
        for (k=1;k<K-2;k++){
            up=u(k-1); un=u(k+1);
            w(k)=u(k)*(1 - ar*(un-up));
        }   
        w(0)=u(0)*(1 + ar*u(1));
		w(K-2)=u(K-2)*(1 - ar*u(K-3));
        u=w;
		lapack::gbtrs(A, p, u); // Solve
        if(op=="approx") output(tn, 0.5, u(K/2-1), K-1, N);
        if(op=="plot0" && n%NMOD==0) plot0(u, tn, K-1, N);
        if(op=="plot1" && n%NMOD==0) plot1(u, tn, K-1, N);
        if(op=="plot3d" && n%NMOD==0) plot3d(u, tn, K-1, N);
	}
	return 0;
}
