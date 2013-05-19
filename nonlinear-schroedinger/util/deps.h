#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <complex>
#include <typeinfo>

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/bindings/bandwidth.hpp>

#include <boost/numeric/bindings/lapack/driver/gbsv.hpp> 
#include <boost/numeric/bindings/lapack/computational/gbtrf.hpp> 
#include <boost/numeric/bindings/lapack/computational/gbtrs.hpp> 
#include <boost/numeric/bindings/lapack/computational/getri.hpp> 
#include <boost/numeric/bindings/lapack/computational/getrf.hpp> 
#include <boost/numeric/bindings/lapack/driver/gtsv.hpp> 

#include "print.h"

using namespace boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

typedef std::complex<double> cval;
typedef double val;

typedef boost::numeric::ublas::matrix<val,column_major> mat;
typedef boost::numeric::ublas::matrix<cval,column_major> cmat;
typedef boost::numeric::ublas::vector<val> vec;
typedef boost::numeric::ublas::vector<cval> cvec;

