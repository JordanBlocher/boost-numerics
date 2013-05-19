#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <typeinfo>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp> 


#include <boost/numeric/bindings/lapack/geqrf.hpp> 
#include <boost/numeric/bindings/lapack/ormqr.hpp> 
#include <boost/numeric/bindings/lapack/gesv.hpp> 
#include <boost/numeric/bindings/lapack/gels.hpp> 



#include "print.h"

using namespace boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

typedef std::complex<double> val;

typedef boost::numeric::ublas::matrix<val,column_major> mat;
typedef boost::numeric::ublas::vector<val> vec;


