#ifndef PRINT_H
#define PRINT_H

#include "deps.h"


// Vector print
void printenergy(cvec u){
    val energy=real(inner_prod(u,conj(u)));
    printf("energy=%25.14e\n",energy/u.size());
}
void printu(cvec u){
    int i;
    val dx=1.0/u.size();
    for(i=0;i<u.size();i++) printf("%g %g %g\n",i*dx, u(i).real(), u(i).imag());
    printf("\n");
}
void plotu(cvec u){
    int i;
    val dx=1.0/u.size();
    printf("set terminal x11 noraise\nset style data lines\n\n");
	printf("set yrange [-5:5]\n");
    printf("set zrange [-5:5]\n");
    printf("splot '-'\n");
    for(i=0;i<u.size();i++) printf("%g %g %g\n",i*dx, u(i).real(), u(i).imag());
	printf("e\n");
    fflush(stdout);
}
// Matrix print
void matprintrow(cvec u){
    int i,m;
  	if(u.size()>10) m=10;
	else m=u.size();
    for(i=0;i<m;i++) printf("%8.8g+%1.8gi",u(i).real(), u(i).imag());
    if(m<u.size()) printf("  ...");
    printf("\n");
}
void matprintf(cmat A){
	int i,m;
	if(A.size1()>10) m=10;
	else m=A.size1();
	for(i=0;i<m;i++){ matrix_row<cmat> row(A,i); matprintrow(row);}
    if(m<A.size1()) printf("  .\n  .\n  .\n");
    printf("\n");
}

#endif
