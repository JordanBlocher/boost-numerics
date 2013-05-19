//#ifndef DEPS_H
//#define DEPS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <typeinfo>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp> 

#include <boost/numeric/ublas/matrix.hpp> 
#include <boost/numeric/ublas/vector.hpp> 

#include <string.h>
#include <sstream>

using namespace boost::numeric::ublas;

typedef double val;

typedef boost::numeric::ublas::matrix<val,column_major> mat;
typedef boost::numeric::ublas::vector<val> vec;

//#endif

