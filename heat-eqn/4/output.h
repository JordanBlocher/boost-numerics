#ifndef OUTPUT_H
#define OUTPUT_H

#include "../deps.h"
#include "../print.h"

void output(val tn, val dx, vec u, std::string op){
   vec p1(u.size()), p2(u.size()), e1(u.size()), e2(u.size());
   if(op=="approx" || op=="all" || op=="plot"){
         if(tn==0.1){
            p1=u;
        }
        if(tn==0.9){
            p2=u;
        }
        if(tn==2.0){
            if(op=="approx" || op=="all"){
                printf("%2s %2.8s 0.1 %10.8s 0.9 %11.8s 2.0\n%s %12.8s %15.8s %15.8s\n\n","time ="," "," "," ","x+dx","u(x+dx)","u(x+dx)","u(x+dx)");
                printf("Approximation\n");
                vecprintf2(p1,p2,u,dx);
            }
            else if(op=="plot"){
                    printf("\n#tn = %24.15e\n", 0.1);
                    vecprintf(p1,dx);
                    printf("\n#tn = %24.15e\n", 0.9);
                    vecprintf(p2,dx);
                    printf("\n#tn = %24.15e\n", 2.0);
                    vecprintf(u,dx);
           }
        }
   }
    if(op=="plot3d"){
        plotf(u,tn,dx);
    }
}
#endif
