/**
Copyright © 2019 Alexey A. Shcherbakov. All rights reserved.

This file is part of sphereml.

sphereml is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

sphereml is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with sphereml. If not, see <https://www.gnu.org/licenses/>.
**/

#ifndef _MATRIX_H
#define _MATRIX_H

#include <iostream>
#include <complex>

#define _USE_MATH_DEFINES

using namespace std;

typedef complex<double> Complex;

const Complex j_(0.,1.);

class Vector;
class Matrix;

class Vector {
public:
     unsigned Nrow;
     Complex* Data;

     Vector(unsigned Mown = 1);
     Vector(const Vector&);
     ~Vector();

     Vector& operator = (const Vector&);
     Complex& operator() (unsigned i) const {
          if (i<Nrow) return Data[i];
          else {cout<<"vector: out of boundaries "<<i<<" "<<Nrow<<endl; return Data[0];}
     }

     Vector& operator += (const Vector&);
     Vector& operator -= (const Vector&);

     Vector operator * (const Complex);
     Vector& operator *= (const Complex tc) {(*this) = (*this)*tc; return *this;}

     Complex operator * (const Vector&) const;
     Vector operator * (const Matrix& B);
     Vector& operator *= (const Matrix& B) {*this = (*this)*B; return *this;}
     Vector& operator /= (const Matrix&);
     Vector& operator %= (const Matrix&);
     Vector operator + (const Vector& B) {Vector C(*this); C += B; return C;}
     Vector operator - (const Vector& B) {Vector C(*this); C -= B; return C;}
     Vector operator / (const Matrix& B) {Vector C(*this); C /= B; return C;}
     Vector operator % (const Matrix& B) {Vector C(*this); C %= B; return C;}

     double normF(unsigned int nn = 0) const;

     double cmp(const Vector&);

     friend Vector operator * (Complex tc, const Vector &V) {Vector VV(V); VV *= tc; return VV;};
};

class Matrix {
public:
     unsigned Nrow, Ncol;
     Complex* Data;

     Matrix(unsigned Mown = 1, unsigned Mext = 0);
     Matrix(const Matrix&);
     ~Matrix();

     Matrix& operator = (const Matrix&);
     Complex& operator () (unsigned i, unsigned j) const {
          if (i<Nrow && j<Ncol) return Data[i*Ncol + j];
          else {cout<<"matrix: out of boundaries "<<i<<" "<<j<<" "<<Nrow<<" "<<Ncol<<endl; return Data[0];}
     }

     double cmp(const Matrix&);
     double normF(int n1=-1) const;

     void eye();

     Matrix& operator += (const Matrix&);
     Matrix& operator -= (const Matrix&);
     Vector operator * (const Vector&) const;
     Matrix operator * (const Matrix&) const;
     Matrix& operator *= (const Matrix& B) {*this = (*this)*B; return *this;}
     Matrix& operator /= (const Matrix&);
     Matrix& operator %= (const Matrix&);
     Matrix operator + (const Matrix& B) {Matrix C(*this); C += B; return C;}
     Matrix operator - (const Matrix& B) {Matrix C(*this); C -= B; return C;}
     Matrix operator / (const Matrix& B) {Matrix C(*this); C /= B; return C;}
     Matrix operator % (const Matrix& B)     {Matrix C(*this); C %= B; return C;}

     Matrix transp() const;
     Matrix mconj() const;

     Complex trace() const;
};

#endif
