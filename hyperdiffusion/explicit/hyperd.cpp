#include "deps.h"
#include "output.h"
#include "print.h"

#define K 10
#define T 0.25
#define N 5
#define NU 0.00005
#define NMOD 1 

val f(val x){
    return x*x*x*(-48.0*x*x + 112.0*x - 64.0);
}
val delta4(vec u, int i){
	if (i==0||i==K) return 0.0;
	else if (i==1) return u(3) - 4*u(2) + 5*u(1);
	else if (i==K-1) return 5*u(K-1) - 4*u(K-2) + u(K-3);
	else return u(i+2) - 4*u(i+1) + 6*u(i) - 4*u(i-1) + u(i-2);
}

main(int argc, char* argv[]){ 
    static vec u(K+1);
    int k,n;
    std::string op=argv[1]; // command line argument
	val dx=1.0/K, dt=T/N, dx4=dx*dx*dx*dx, xk, tn; // discretization variables

    // Initialize
    for(k=0;k<=K;k++){
        xk=k*dx;
        u(k)=f(xk);
    }
    if(op=="plot") printf("set terminal x11 noraise\nset yrange [-5:5]\nset style data lines\n\n");

    val rho=NU*dt/dx4;    // rho = nu*dt/dx^4
    vec temp(K+1);
    temp=u;

    for(n=0;n<=N;n++){
		tn=n*dt;
        u=temp;
        for(k=0;k<=K;k++){   
            xk=k*dx;
            temp(k) = u(k) - rho*delta4(u,k);    
            //printu(u,tn,xk,K);
        }
        if(op=="plot0") plot0(u, tn, K-1, N);
        if(op=="plot1") plot1(u, tn, K-1, N);
        if(op=="plot3d") plot3d(u, tn, K-1, N);
        if(op=="approx") output(tn, 0.5, u(K/2), K-1, N);
    }

    return 0;
}
