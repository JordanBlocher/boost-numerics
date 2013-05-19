#include "deps.h"

#define e 2.7

mat householder(mat A, int n, int m){
	int i,j,k;
	val alpha;
    mat Q(n,m);
    identity_matrix I(n);
    vec u(n),v(n),e1(n);
    e1.clear(); e1(0)=1.0;
    for(i=0;i<A.size2();i++){
        alpha=-pow(e,sqrt(-1))*norm_2(vec(column(A,i));
        u=vec(column(A,i))+alpha*e1;
        v=(1/norm_2(u))*u;
        Q=I-(1+((inner_prod(trans(conj(x)),v))/inner_prod(trans(conj(v)),x))*outer_product(v,conj(v));
    }
    return Q;
}
