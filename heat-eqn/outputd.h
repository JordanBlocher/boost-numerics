#ifndef OUTPUTSMALL_H
#define OUTPUTSMALL_H

#include "deps.h"
#include "print.h"

void outputd(val tn, val dx, vec u, std::string op){
   vec p1(u.size()), p2(u.size());
   if(op=="plot"){
         if(tn==0.06){
            p1=u;
        }
        if(tn==0.1){
            p2=u;
        }
        if(tn==0.9){
            printf("\n#tn = %24.15e\n", 0.06);
            vecprintf(p1,dx);
            printf("\n#tn = %24.15e\n", 0.1);
            vecprintf(p2,dx);
            printf("\n#tn = %24.15e\n", 0.9);
            vecprintf(u,dx);
        }
   }
   else if((op=="plotd")){
         if(tn==0.001){
            p1=u;
        }
        if(tn==0.002){
            p2=u;
        }
        if(tn==0.003){
            printf("\n#tn = %24.15e\n", 0.001);
            vecprintf(p1,dx);
            printf("\n#tn = %24.15e\n", 0.002);
            vecprintf(p2,dx);
            printf("\n#tn = %24.15e\n", 0.003);
            vecprintf(u,dx);
        }
   }
   else if(op=="plot3d" && tn<0.003) plotf(u,tn,dx);
}
#endif
