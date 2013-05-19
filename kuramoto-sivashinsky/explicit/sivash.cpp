#include "deps.h"
#include "output.h"
#include "print.h"

#define MU 0.1
#define NU 0.00005
#define K 10
#define T 0.25
#define N 100000000

val f(val x){
	return -48*pow(x,5) + 112*pow(x,4) - 64*pow(x,3);
}

val delta0(vec u, int i){
	if (i==0||i==K) return 0.0;
	return u(i+1) - u(i-1); 
}

val delta2(vec u, int i){
	if (i==0||i==K) return 0.0;
	return u(i+1) - 2*u(i) + u(i-1);
}

val delta4(vec u, int i){
	if (i==0||i==K) return 0.0;
	if (i==1) return u(3) - 4*u(2) + 5*u(1);
	if (i==K-1) return 5*u(K-1) - 4*u(K-2) + u(K-3);
	return u(i+2) - 4*u(i+1) + 6*u(i) - 4*u(i-1) + u(i-2);
}

int main(){
	static vec u(K+1), un(K+1);
	val dt = T/N;
	val dx = 1.0/K;
	val gam = MU*dt/(dx*dx);
	val rho = NU*dt/(dx*dx*dx*dx);
	val ar = dt/(2*dx);
	int j,n=1;

	//Initialize U
	for (j=0;j<K+1;j++) u(j)=f(j*dx);
    //printf("Initial vector\n"); vecprintf(u,K);

	val tn = n*dt;
	while (tn<=T||n<N){
		for (j=0;j<K+1;j++) un(j) = u(j) - gam*delta2(u,j) - rho*delta4(u,j) - ar*u(j)*delta0(u,j);
		for (j=0;j<K+1;j++) u(j) = un(j);
		n+=1;
		tn = n*dt;
	}
	printf("At t = 0.25 and space 0.5 the 'heat' is %24.15e\n", u(K/2)); 
	return 0;
}
