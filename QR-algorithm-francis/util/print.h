#ifndef PRINT_H
#define PRINT_H

#include "deps.h"

void vecprintf(vec u){
    int i,m;
	if(u.size()>5) m=5;
	else m=u.size();
    for(i=0;i<m;i++){ 
        if(typeid(u[i])==typeid(std::complex<double>)){ 
            if(abs(u[i].imag())<1e-8) printf("%15.8lf %4s",u[i].real()," ");
            else if(u[i].imag()<1e-8) printf("%10.2lf-%2.2lfi ",u[i].real(),fabs(u[i].imag())); 
            else printf("%10.2lf+%2.2lfi ",u[i].real(),u[i].imag()); 
        }
        else printf("%15.8lf ",u[i]);
        }
	if(u.size()!=m) printf("...");
    printf("\n");
}
void matprintf(mat A){
	int i,m;
	if(A.size1()>5) m=5;
	else m=A.size1();
	for(i=0;i<m;i++) vecprintf(vec(row(A,i)));
    if(A.size1()!=m) printf("\n .\n. \n.");
}

#endif
