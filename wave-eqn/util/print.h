#ifndef PRINT_H
#define PRINT_H

#include "deps.h"

void printenergy(vec u){
    val energy=inner_prod(u,u);
    printf("energy=%25.14e\n",energy/u.size());
}
// Vector print
void vecprintf(vec u, val dx){
    int i,m;
    for(i=0;i<u.size();i++) printf("%24.15e %24.15e\n",i*dx, u(i));
    printf("\n");
}
void printu(vec u, val tn, int K){
    int i;
    val dx=1.0/K;
    printf("\n#tn = %24.15e\n", tn);
    for(i=0;i<u.size();i++) printf("%24.15e %24.15e\n",i*dx, u(i));
    printf("\n");
}
void printavg(vec u){
    val avg = accumulate( u.begin(), u.end(), 0.0 )/ u.size();
    printf("Average value = %f\n", avg);
}
void plotu(vec u, val tn, int K){
    int k;
    val dx=1.0/K;
    if(tn==0.0) printf("set terminal x11 noraise\nset style data lines\n\n");
	printf("plot '-'\n");
    for(k=0;k<u.size();k++){
        val xk=k*dx;
        printf("%g %g\n",xk,u[k]);
    }
	printf("e\n");
}
// Matrix print
void printU(mat U, val tn, val xk){
    int i,m;
    printf("%2s t %10.5s x %9.5s u(x,t)\n\n"," "," "," ");
    printf("%2.5f %4s %2.5f %4s %2.5f\n",tn, " ", xk, " ", U(xk,tn));
    printf("\n");
}
void matprintrow(vec u){
    int i,m;
  	if(u.size()>10) m=10;
	else m=u.size();
	for(i=0;i<m;i++) printf("%12.8f",u(i));
    if(m<u.size()) printf("  ...");
    printf("\n");
}
void matprintf(mat A){
	int i,m;
	if(A.size1()>10) m=10;
	else m=A.size1();
	for(i=0;i<m;i++){ matrix_row<mat> row(A,i); matprintrow(row);}
    if(m<A.size1()) printf("  .\n  .\n  .\n");
    printf("\n");
}


#endif
