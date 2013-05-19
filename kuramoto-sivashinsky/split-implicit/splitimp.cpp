#include "deps.h"
#include "output.h"
#include "print.h"

#define MU 0.1
#define NU 0.00005
#define K 1000
#define T 0.25
#define N 100000
#define NMOD1 1000
#define NMOD2 2000

// Function for the initial condition
val f(val x){
	return x*x*x*(-48*x*x + 112*x - 64);
}

int main(int argc, char* argv[]){ 
    std::string op=argv[1]; // command line argument
	// Variable declaration
	static vec v(K-1), u(K+1);
	val dt = T/N;
	val dx = 1.0/K;
	val gam = MU*dt/(dx*dx);
	val rho = NU*dt/(dx*dx*dx*dx);
	val ar = dt/(2*dx);
	int i,j,k;

	//Initialize U
	for (j=0;j<K-1;j++) v(j)=f((j+1)*dx);

    //printf("Initial Vector\n");vecprintf(v,dx);
    //printf("set terminal x11 noraise\nset yrange [-4:4]\nset style data lines\n\n");

    static banded_matrix<val, row_major> A(K-1, K-1, 2, 4);

	//Creating a permutation storage vector for the lapack functions

//	printf("tn=%3g , %3g\n",0.0,v(K/2-1));
	int n;
	for(n=0;n<=N;n++){
		val tn=n*dt;
 
    	// Initialize matrix
    	for(i=0; i<A.size1(); i++){
            A(i,i)=1.0-2*gam+6*rho; 
          	if(i+1<A.size1())A(i,i+1)=gam-4*rho+ar*v(i);
			if(i-1>=0)A(i,i-1)=gam-4*rho-ar*v(i);
            if(i+2<A.size1())A(i,i+2)=A(i+2,i)=rho;
			if(i-2>=0)A(i,i-2)=A(i-2,i)=rho;
		}
   	 // Boundary Conditions
    	A(0,0)=A(K-2,K-2)=1.0-2*gam+5*rho;
    	//printf("Pentadiagonal Matrix\n"); matprintf(A);

		vector<fortran_int_t> p(K-1);
	
		lapack::gbtrf(A, p); // LU-decompostion
		//matprintf(A);

		lapack::gbtrs(A, p, v); // Solve
		//printf("tn=%24.15e %24.15e\n",tn,v(K/2-1));
            
        if(op=="plot0" && n%NMOD1==0) plot0(v, tn, K-1, N);
        if(op=="plot1"  && n%NMOD1==0) plot1(v, tn, K-1, N);
        if(op=="plot3d" && n%NMOD2==0) plot3d(v, tn, K-1, N);
        if(op=="approx") output(tn, 0.5, v(K/2-1), K-1, N);
	}
	
	//Reassigns the data from the lapack vector v to the spacial vector u
	for (i=1;i<K;i++) u(i)=v(i-1);
	u(0)=u(K)=0.0;
    //Prints the desired record
	//printf("At x=0.5 and t=0.25 the heat is %24.15e\n",u(K/2));
	return 0;
}
