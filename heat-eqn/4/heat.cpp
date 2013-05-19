#include "../deps.h"
#include "output.h"

#define T 1
#define NU 0.1

val f(val x){
	return 0.0;
}
val F(val x, val t){
    return sin(2*M_PI*x)*sin(4*M_PI*t);
}
val a(val t){
	return 0.0;
}
val b(val t){
    return 0.0;
}

int main(int argc, char* argv[]){ 
    std::istringstream mss(argv[2]);
    std::istringstream nss(argv[3]);
    static double M, N;
    mss>>M; nss>>N;
    static vec u(M+1);
	int k,n;
    bool erf=true;
    std::string op=argv[1];
	val dx=1.0/M, xk, tn, un, up;
	for(k=0;k<=M;k++){
        xk=k*dx;
        u(k)=f(xk);
    }
    val dt=(val)T/N, dx2=dx*dx;
    for(n=1;n<=2*N;n++){
        tn=n*dt;
        u(0)=a(tn);
        for(k=1;k<M;k++){
           un=u(k+1);
           up=u(k-1);
           u(k)+=NU*dt/dx2*(un-2*u(k)+up)+dt*F(k*dx,tn);
        }
        output(tn, dx, u, op);
        u(M)=b(tn);
        }
    return 0;
}

