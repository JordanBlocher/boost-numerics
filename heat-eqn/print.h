#ifndef PRINT_H
#define PRINT_H

#include "deps.h"

void plotf(vec u, val tn, val dx){
    int i,m;
	if(u.size()>5) m=5;
	else m=u.size();
    for(i=0;i<u.size();i++) printf("%f %f %f\n",tn, i*dx, u(i));
    printf("\n");
}
void vecprintf(vec u, val dx){
    int i,m;
    for(i=0;i<u.size();i++) printf("%24.15e %24.15e\n",i*dx, u(i));
    printf("\n");
}
void vecprintferr(vec u, val tn){
    int i,m;
    //for(i=1;i<u.size()-1;i++) if( fabs(u(i))== norm_inf(u)) printf("%24.15e %24.15e\n",i*dx, fabs(u(i)));
    printf("%24.15e %24.15e\n",tn, norm_inf(u));
    printf("\n");
}
void vecprintf2(vec u, vec v, vec y, val dx){
    int i,m;
    for(i=0;i<u.size();i++) printf("%2.2f %15.8f %15.8f %15.8f\n",i*dx, u(i), v(i), y(i));
    printf("\n");
}
void vecprintf2err(vec u, vec v, vec y, val dx){
    int i,m;
    for(i=0;i<u.size();i++) printf("%2.2f %15.8f %15.8f %15.8f\n",i*dx, fabs(u(i)), fabs(v(i)), fabs(y(i)));
    printf("\n");
}
void vecprintf3(vec u, val dx){
    int i,m;
    for(i=0;i<u.size();i++) printf("%2.2f %19.8e\n",i*dx, u(i));
    printf("\n");
}
#endif
