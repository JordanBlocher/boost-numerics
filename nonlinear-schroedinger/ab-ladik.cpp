#include "deps.h"
#include "print.h"

#define K 100
#define N 1000
#define T 0.2
#define L 10.0
#define NMOD 10
#define i cval(0,1)

cval f(val x){
    if(x>L/2) x-=L;
    return 2*exp(-i*x*M_PI/5.0)-i*cos(3*x*M_PI/5.0);
}

int main(){
    int k;

    val dx=L/K, dt=T/N;
    static cvec u(K+2), w(K+2), z(K+2);
    static banded_matrix<val> U(K, K, 2, 4);
    vector<fortran_int_t> p(K);    
    cval e=2*i*dx*dx/dt;
    
    // Tridiagonal Matrix U
    for(i=0; i<U.size1(); i++){
            U(i,i)=-(2.0+e); 
            k=std::max(i-1,0);
            U(k,k+1)=1.0;
            U(k+1,k)=1.0;
    }
    // Boundary Conditions
    U(0,0)+=1.0; U(K-1,K-1)-=1.0; 
    
    // Initial condition
    for(k=0;k<K;k++) {
        val xk=k*dx;
        u[k]=f(xk);
    }
    // Boundary conditions
    u[K]=u[0]; u[K+1]=u[1];

    val alpha, beta; // Sherman-Morrison 
    w(0)=1.0; w(K-1)=-1.0; z(0)=R/4.0; z(K-1)=R/4.0;
    
    // Test Boundary Conditions
    if(op=="test"){ 
        mat wzt = outer_prod(w, z);
        mat B = U - wzt; 
        printf("Matrix U\n"); matprintf(U);
        printf("Matrix B\n"); matprintf(B);
        printf("Matrix wz^T\n"); matprintf(wzt);
    }

    lapack::gbtrf(U, p); // LU-decompostion
    lapack::gbtrs(U, p, w); //B^-1w
    alpha=1.0/(1.0-inner_prod(z,w)); // alpha = 1/(1-z^T(b^-1w))
    
    if(op=="test") printf("alpha=%f\n",alpha);

    #define slice(k) u(((k)+K)%K)

    for(n=0;n<=N;n++){
        for(k=0;k<U.size1();k++) r(k)=(R/4.0)*(slice(k-1)-slice(k+1))+u(k);
        tn=n*dt;
        lapack::gbtrs(U, p, r); // B^-1r^n
        beta=alpha*inner_prod(z, trans(r)); // beta = alpha*z^T(B^-1r^n)
        if(op=="test"){ printf("beta=%f\n",beta);  printavg(u);}
        u=r+beta*w;
        if(op=="pipe") plotu(u, tn, K);
        if(op=="approx" && tn==t)printu(u, tn, K);
    }
 

    printenergy(u);
    printf("# rho = (%g,%g)\n",rho.real(),rho.imag());
    plotu(u);
    int n;
    for(k=1;k<=K;k++){
        w[k]=rho*(u[k+1]-2*u[k]+u[k-1])
            -i*2*dt*(2*u[k]*conj(u[k])*u[k]);
    }
    w[0]=w[K]; w[K+1]=w[1];
    for(k=0;k<=K+1;k++){
        unm1[k]=u[k]+w[k]/2;
    }
    for(n=1;n<N;n++){
        for(k=1;k<=K;k++){
            w[k]=rho*(u[k+1]-2*u[k]+u[k-1])
                -i*4*dt*(u[k]*conj(u[k])*u[k]);
        }
        w[0]=w[K]; w[K+1]=w[1];
        for(k=0;k<=K+1;k++){
            cval tmp=u[k];
            u[k]=unm1[k]+w[k];
            unm1[k]=tmp;
        }
        if(!((n+1)%NMOD)){
            printf("set title \"tn=%-15.5f\"\n",(n+1)*dt);
            plotu(u);
        }
    }
    printenergy(u);

        return 0;
}
