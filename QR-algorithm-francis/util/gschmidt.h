#include "deps.h"

mat gschmidt(mat A, int n, int m){
	int i,j,k;
	double r;
    mat At(m,n), Qt(m,n), IdmI(m,n);
    At=conj(trans(A));
    vec temp(n);
	for(j=0;j<n;j++) {
		for(i=0;i<m;i++) temp(i)=At(j,i);
		for(k=0;k<j;k++){
			r=inner_prod(vec(row(Qt,k)),vec(row(At,j)));
			for(i=0;i<m;i++) temp(i)-=r*Qt(k,i);
		}
		r=norm_2(temp);
		for(i=0;i<m;i++) Qt(j,i)=temp(i)/r;
	}
    // Check that Qt Q=In
	for(i=0;i<n;i++) for(j=0;j<n;j++) {
		IdmI(i,j)=inner_prod(vec(row(Qt,i)),vec(row(Qt,j)));
		if(i==j) IdmI(i,j)-=1.0;
	}
    if(norm_inf(IdmI)>10e-8){ 
        printf("Gram-Schmidt failed to orthogonalize the matrix A\n");
        printf("|I-Q*Q|=%lf\n",(float)norm_inf(IdmI));
    }
	
	return Qt;
}
