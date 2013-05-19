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

	val dx=1.0/K, dt=(val)T/N, xk, tn, tmp; // discretization variables
    val R=a*dt/dx;    // R = a*dt/dx
    
    static vec u(K+2), v(K+2);
    
    // Initial conditions
    for(k=0;k<u.size();k++){
        xk=k*dx;
        u(k)=f(xk);
    }
    if(op=="energy") {printf("t=0.0\n"); printenergy(u);}

    for(n=0;n<=N;n++){
        tn=n*dt;
        for(k=1;k<=K;k++){
            v(k)=u(k)-R*(u(k+1)-u(k));
        }
        v(0)=v(K); v(K+1)=v(1);
        u=v;
        if(op=="pipe") plotu(u, tn, K);
        if(op=="approx" && t==tn) printu(u, tn, K);
    }

    if(op=="energy"){ printf("t=%d\n",T); printenergy(u);}

    return 0;
}
