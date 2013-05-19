#include "deps.h"
#include "print.h"

#define NU 0.00001
#define N 400
#define T 2
#define K 100
#define r 0.5

val f(val x){
    return 3.0*sin(4.0*M_PI*x);
}


main(int argc, char* argv[]){ 
    
    int k,n,i,j;

    std::string op=argv[1], prob=argv[3]; // command line argument
    std::istringstream tin(argv[2]);
    val t; tin>>t;

	val dx=1.0/K, dt=(val)T/N, xk, tn; // discretization variables
    val R=NU*dt/dx;    // R = NU*dt/dx
    
    static vec u(K+2), v(K+2);
    
    // Initial conditions
    for(k=0;k<u.size();k++){
        xk=k*dx;
        u(k)=f(xk);
    }
    
    if(prob=="a"){
    for(n=0;n<=N;n++){
        tn=n*dt;
        for(k=1;k<=K;k++){
            v(k)=u(k)-R*(u(k)-u(k-1))+r*(u(k-1)-2*u(k)+u(k+1));
        }
        v(0)=v(K+1)=0;
        u=v;
        if(op=="pipe") plotu(u, tn, K);
        if(op=="approx" && t==tn) printu(u, tn, K);
    }
    }
    if(prob=="b"){
    for(n=0;n<=N;n++){
        tn=n*dt;
        for(k=1;k<=K;k++){
            v(k)=u(k)-0.5*R*(u(k)-u(k-1))+0.5*R*(u(k-1)-\
            2*u(k)+u(k+1))+r*(u(k-1)-2*u(k)+u(k+1));
        }
        v(K+1)=v(0)=0;
        u=v;
        if(op=="pipe") plotu(u, tn, K);
        if(op=="approx" && t==tn) printu(u, tn, K);
    }
    }
    return 0;
}
