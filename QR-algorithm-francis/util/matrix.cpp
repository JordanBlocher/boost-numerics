#include "../basic/matrix.h"

Matrix::Matrix(void)
{
    e11=e12=e13=e21=e22=e23=e31=e32=e33=0;
}

Matrix::Matrix(double r1c1, double r1c2, double r1c3, double r2c1, double r2c2, double r2c3, double r3c1, double r3c2, double r3c3)
{
    e11= r1c1;
    e12= r1c2;
    e13= r1c3;
    e21= r2c1;
    e22= r2c2;
    e23= r2c3;
    e31= r3c1;
    e32= r3c2;
    e33= r3c3;
}

Matrix Matrix::operator *(Matrix m)
{
    return Matrix( m.e11*e11 + m.e12*e21 + m.e13*e31,
                m.e11*e12 + m.e12*e22 + m.e13*e32,
                m.e11*e13 + m.e12*e23 + m.e13*e33,
                m.e21*e11 + m.e22*e21 + m.e23*e31,
                m.e21*e12 + m.e22*e22 + m.e23*e32,
                m.e21*e13 + m.e22*e23 + m.e23*e33,
                m.e31*e11 + m.e32*e21 + m.e33*e31,
                m.e31*e12 + m.e32*e22 + m.e33*e32,
                m.e31*e13 + m.e32*e23 + m.e33*e33);
}

Vector Matrix::operator*(Vector u)
{
    return Vector(e11*u.X() + e21*u.Y() + e31*u.Z(),
                e12*u.X() + e22*u.Y() + e23*u.Z(),
                e13*u.X() + e32*u.Y() + e33*u.Z());
}


