#include "../deps.h"
#include "../print.h"
#include "../outputd.h"
#include "../output.h"

#define T 1.0

val f(val x){
	return sin(4*M_PI*x);
}
void init(vec &u, int m){
    int k;
	val dx=1.0/m, xk;
	for(k=0;k<=m;k++){
        xk=k*dx;
        u(k)=f(xk);
    }
}
 
main(int argc, char* argv[]){
    std::istringstream mss(argv[2]);
    std::istringstream nuss(argv[3]);
    std::istringstream nss(argv[4]);
    bool small=atoi(argv[5]);
    static double M, NU, N;
    mss>>M; nuss>>NU; nss>>N;
    static vec u(M+1);
    static vec sol(0);
	int k,n;
    val alpha;
    std::string op=argv[1];
    init(u,M);
    if(!small) alpha=1.0;
    else if(small && NU==1.0 && M==40) alpha=0.02; 
    else alpha=1.0; 
	val dx=1.0/(alpha*NU*M), tn, dt=T/N, dx2=dx*dx, up, un;
    val r=NU*dt/dx2;    // r = nu*dt/dx^2
    if(op=="approx") printf("Steps taken: N=%f with dx=%f, dt=%f\n\n",N,dx,dt);
    for(n=1;n<=N;n++){
        tn=n*dt;
        u(0)= 0.0;  // a(tn)
        for(k=1;k<M;k++){
           up=u(k-1);
           un=u(k+1);
           // u_k = dt/dx*(u_k+1 - u_k-1) + r*(u_k+1 - 2*u_k + u_k-1)
           u(k)+=-dt/dx*(un-up)+r*(un-2*u(k)+up);
        }
        u(M)= 0.0;  // b(tn)
        if(!small) outputd(tn, dx, u, op);
        else output(tn, dx, u, sol, op);
    }
    return 0;
}

