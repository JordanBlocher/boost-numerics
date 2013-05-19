#include "deps.h"
#include "print.h"

#define K 100
#define N 1000
#define T 0.2
#define L 10.0
#define NMOD 10
#define i cval(0,1)

cval f(val x){
    return 2*exp(-i*x*M_PI/5.0)-i*cos(3*x*M_PI/5.0);
}

int main(int argc, char* argv[]){
    
    std::string op=argv[1]; // Command line arguments
    
    int k,n;

    val dx=L/K, dt=T/N;
    static cvec u(K+2), w(K+2), unm1(K+2);

    // Initial condition
    for(k=0;k<K;k++) {
        val xk=k*dx;
        u[k]=f(xk);
    }
    // Boundary conditions
    u[K]=u[0]; u[K+1]=u[1];

    if(op=="energy") {printf("t=0.0\n"); printenergy(u);}
    cval rho=-i*2*dt/dx/dx;
    if(op=="test") printf("# rho = (%g,%g)\n",rho.real(),rho.imag());
    if(op=="pipe") plotu(u);

    for(k=1;k<=K;k++){
        w[k]=rho*(u[k+1]-2*u[k]+u[k-1])
            -i*2*dt*(2*u[k]*conj(u[k])*u[k]);
    }
    w[0]=w[K]; w[K+1]=w[1];
    for(k=0;k<=K+1;k++){
        unm1[k]=u[k]+w[k]/2;
    }
    for(n=0;n<=N;n++){
        for(k=1;k<=K;k++){
            w[k]=rho*(u[k+1]-2*u[k]+u[k-1])
                -i*4*dt*(u[k]*conj(u[k])*u[k]);
        }
        w[0]=w[K]; w[K+1]=w[1];
        for(k=0;k<=K+1;k++){
            cval tmp=u[k];
            u[k]=unm1[k]+w[k];
            unm1[k]=tmp;
            if(op=="approx" && n*dt==0.2 && k*dx==5.0)
                printf("NLS Explicit\nt=0.2, x=5\n u(k)=\
                %g+%gi\nN=%d, K=%d\n\n", \
                u(k).real(), u(k).imag(),N,K);
        }
        if(op=="pipe" && !((n+1)%NMOD)){ 
            printf("set title \"tn=%-15.5f\"\n",(n+1)*dt);
            plotu(u);
        }
        if(op=="plot" && n*dt==0.2){ 
            printf("#tn=%g\n\n",n*dt);
            printu(u);
        }
    }
    if(op=="energy") {printf("t=0.2\n"); printenergy(u);}

    return 0;
}
