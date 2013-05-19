#ifndef OUTPUT_H
#define OUTPUT_H

#include "deps.h"
#include "print.h"

void output(val tn, val dx, vec u, vec sol, std::string op){
   vec p1(u.size()), p2(u.size()), e1(u.size()), e2(u.size()), s1(u.size()), s2(u.size());
   if(op=="approx" || op=="all" || op=="plot" || op=="app+erf"){
         if(tn==0.06){
            p1=u;
        }
        if(tn==0.1){
            p2=u;
        }
        if(tn==0.9){
            if(op=="approx" || op=="all" || op=="app+erf"){
                printf("%2s %2.8s 0.06 %10.8s 0.1 %11.8s 0.9\n%s %12.8s %15.8s %15.8s\n\n","time ="," "," "," ","x+dx","u(x+dx)","u(x+dx)","u(x+dx)");
                printf("Approximation\n");
                vecprintf2(p1,p2,u,dx);
            }
            else if(op=="plot"){
                    printf("\n#tn = %24.15e\n", 0.06);
                    vecprintf(p1,dx);
                    printf("\n#tn = %24.15e\n", 0.1);
                    vecprintf(p2,dx);
                    printf("\n#tn = %24.15e\n", 0.9);
                    vecprintf(u,dx);
           }
        }
   }
   if(op=="sol" || op=="all" || op=="plotsol"){
        if(tn==0.06){
            s1=sol;
        }
        if(tn==0.1){
            s2=sol;
        }
        if(tn==0.9){
            if(op=="sol" || op=="all"){
                printf("Exact Solution\n");
                vecprintf2(s1,s2,sol,dx);
            }
            else if(op=="plotsol"){
                printf("\n#tn = %24.15e\n", 0.06);
                vecprintf(s1,dx);
                printf("\n#tn = %24.15e\n", 0.1);
                vecprintf(s2,dx);
                printf("\n#tn = %24.15e\n", 0.9);
                vecprintf(sol,dx);
            }
        }
    }
    if(tn==50.0 && op=="approx50"){
        printf("%2s %2.8s 50.0\n%s %10.8s\n","time ="," ","x+dx","u(x+dx)");
        vecprintf3(u,dx);
    }
    if(tn==50.0 && op=="plot50"){
        printf("\n#tn = %24.15e\n", 50.0);
        vecprintf(u,dx);
    }
    if(tn==50.0 && op=="plot50erf"){
        vecprintferr(u-sol,50.0);
    }
    if(op=="erf" || op=="ploterf" || op=="all" || op=="app+erf"){
        if(tn==0.06){
            e1=sol-u;
        }
        if(tn==0.1){
            e2=sol-u;
        }
        if(tn==0.9){
            if(op=="erf" || op=="all" || op=="app+erf"){
                printf("Error\n");
                vecprintf2err(e1,e2,sol-u,dx);
            }
            else if(op=="ploterf"){
                //printf("\n#tn = %24.15e\n", 0.06);
                vecprintferr(e1,0.06);
                //printf("\n#tn = %24.15e\n", 0.1);
                vecprintferr(e2,0.1);
                //printf("\n#tn = %24.15e\n", 0.9);
                vecprintferr(sol-u,0.9);
            }

        }
    }
    if(op=="plot3d"){
        plotf(u,tn,dx);
    }
    if(op=="plotsol3d"){
        plotf(sol,tn,dx);
    }
}
#endif
