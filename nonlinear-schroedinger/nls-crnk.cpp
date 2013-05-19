/*
 * Implicit Crank Nicolson Scheme 
 *
 * Source:
 * Thiab R. Taha, Mark J. Ablowitz, 
 * Analytical and Numerical Aspects of Certain
 * Nonlinear Evolution Equations II:
 * Numerical, Nonlinear Schr¨odinger Equation,
 * J. Comput. Phys., 55 (1984), pp. 203–230.
 */


#include "deps.h"
#include "print.h"

#define K 20
#define N 200
#define T 0.2
#define L 10.0
#define NMOD 10
#define i cval(0,1)

cval f(val x){
    return 2*exp(-i*x*M_PI/5.0)-i*cos(3*x*M_PI/5.0);
}

int main(int argc, char* argv[]){

    std::string op=argv[1]; // Command line arguments

    int k,n,j;
    val dx=L/K, dt=T/N, xk, tn;
    
    static banded_matrix<cval> B(K+1, K+1, 2, 4);
    vector<fortran_int_t> p(K+1);    
    
    cval lambda=(dt/dx/dx/2,0);
    if(op=="test")
        printf("# lambda = (%g,%g)\n",lambda.real(),lambda.imag());
    
    // Tridiagonal Matrix B
    for(j=0; j<B.size1(); j++){
        B(j,j)=i+lambda; 
        k=std::max(j-1,0);
        B(k,k+1)=-lambda/2.0;
        B(k+1,k)=-lambda/2.0;
    }
    // Boundary Conditions
    B(0,0)-=lambda/2.0; B(K,K)-=lambda/2.0; 
    
    static cvec u(K+1), w(K+1), z(K+1), v(K+1), r(K+1);

    // Initial condition
    for(k=0;k<u.size();k++) {
        xk=k*dx;
        u(k)=f(xk);
    }
    // Boundary conditions
    if(op=="pipe") plotu(u);

    cval alpha, beta; // Sherman-Morrison 
    w(0)=(-1.0,-1.0); w(K)=(-1.0,-1.0); 
    z(0)=(lambda/2.0,0); z(K)=(lambda/2.0,0);
    
    // Test Boundary Conditions
    if(op=="test"){ 
        cmat wzt = outer_prod(w, z);
        cmat A = B + wzt; 
        printf("Matrix A\n"); matprintf(A);
        printf("Matrix B\n"); matprintf(B);
        printf("Matrix wz^T\n"); matprintf(wzt);
    }

    lapack::gbtrf(B, p); // LU-decompostion
    lapack::gbtrs(B, p, w); //B^(-1)w
    alpha=1.0/(1.0-(z(0)*w(0)+z(K)*w(K))); // alpha = 1/(1-z^T(B^(-1)w))
    
    if(op=="test") printf("alpha=%g %g\n",alpha.real(), alpha.imag());

    if(op=="energy") {printf("t=0.0\n"); printenergy(u);}
    #define S(k) u(((k)+K)%K)
    
    v=u; // Initial iteration
    lapack::gbtrs(B, p, v); // B^(-1)r^n
    beta=alpha*inner_prod(z,v); //beta=alpha*z^T(B^(-1)r^n)
    if(op=="test") printf("#beta=%g %g\n",beta.real(), beta.imag());
    v+=beta*w;

    // Solve
    for(n=0;n<=N;n++){
        tn=n*dt;
        while(norm_inf(u-v)< 0.000000008){
            for(k=0;k<u.size();k++){ 
                r(k)=(lambda/2.0)*(S(j+1)+S(j-1))+(i-lambda)*u(k)
                    +dt*real(v(k)*conj(v(k)))*v(k)
                    +dt*real(u(k)*conj(u(k)))*u(k);
                if(op=="approx" && tn==0.2 && k*dx==5.0)
                    printf("NLS Crank Nicolson Method\nt=0.2, x=5\n u(k)=\
                    %g+%gi\nN=%d, K=%d\n\n",\
                    u(k).real(), u(k).imag(),N,K);
            }
            lapack::gbtrs(B, p, r); // B^(-1)r^n
            beta=alpha*inner_prod(z,r); //beta=alpha*z^T(B^(-1)r^n)
            if(op=="test") printf("#beta=%g %g\n",beta.real(), beta.imag());
            v=r+beta*w;
        }
        u=v;
        if(op=="pipe" && !((n+1)%NMOD)){ 
            printf("set title \"tn=%-15.5f\"\n",(n+1)*dt);
            plotu(u);
        }
        if(op=="plot" && tn==0.2){ 
            printf("#tn=%g\n\n",tn);
            printu(u);
        }
    }
    if(op=="energy") {printf("t=0.2\n"); printenergy(u);}

        return 0;
}
