#include "../deps.h"
#include "../print.h"
#include "../output.h"

#define M 10
#define T 1.0
#define N 50
#define NU 0.166666666666666

val f(val x){
	return sin(2*M_PI*x);
}
val exact_sol(val x, val t){
    return sin(2*M_PI*x)*exp(-NU*4.0*M_PI*M_PI*t);
}
void init(vec &u, vec &sol){
    int k;
	val dx=1.0/M, xk;
	for(k=0;k<=M;k++){
        xk=k*dx;
        u(k)=f(xk);
        sol(k)=exact_sol(xk,0);
    }
}
 
main(int argc, char* argv[]){
    static vec u(M+1), utp1(M+1), utp2(M+1);
    static vec sol(M+1);
	int k,n;
    std::string op=argv[1];
    init(u,sol);
	val dx=1.0/M, tn, dt=T/N, dx2=dx*dx, up, un;
    val r=NU*dt/dx2;    // r = nu*dt/dx^2
    for(n=1;n<=N;n++){
        utp2=u; // u^n-1
        tn=n*dt;
        u(0)=0.0;  // a(tn)
        for(k=1;k<M;k++){
           up=u(k-1);
           un=u(k+1);
           if(n==1) u(k)=exact_sol(k*dx,tn);
           else u(k)=utp1(k)+2*r*(un-2*u(k)+up); // u_k = 2*r*(u_k+1 - 2*u_k + u_k-1)
           sol(k)=exact_sol(k*dx,tn);
        }
        utp1=utp2;
        u(M)= 0.0;  // b(tn)
        output(tn, dx, u, sol, op);
    }
    return 0;
}

