
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

#include <iostream>
#include <memory.h>

#include "./matrix.h"
#include <memory.h>

const double nz_ = 1.e-15;

Vector::Vector(unsigned Mown) {
     Nrow = Mown;
     Data = new Complex[Nrow];
}

Vector::Vector(const Vector& V) {
     Nrow = V.Nrow; 
     Data = new Complex[Nrow];
     memcpy(Data,V.Data,Nrow*sizeof(Complex));
}

Vector::~Vector() {delete[] Data;}

Vector& Vector::operator = (const Vector& M) {
     if (Nrow != M.Nrow) {if (Nrow) delete[] Data; Nrow = M.Nrow; Data = new Complex[Nrow];}
     memcpy(Data, M.Data, Nrow*sizeof(Complex));
     return *this;
}

Vector& Vector::operator += (const Vector& M) {
     Complex *D1 = Data;
     Complex *D2 = M.Data;
     Complex *D1end = Data + Nrow;
     if ((*this).Nrow != M.Nrow) {cout<<"Vector operator += : wrong vector length\n"; return *this;}
     while (D1 < D1end) *D1++ += *D2++;
     return *this;
}

Vector& Vector::operator -= (const Vector& M) {
     Complex *D1 = Data;
     Complex *D2 = M.Data;
     Complex *D1end = Data + Nrow;
     if ((*this).Nrow != M.Nrow) {cout<<"Vector operator -= : wrong vector length\n"; return *this;}
     while (D1 < D1end) *D1++ -= *D2++;
     return *this;
}

Vector Vector::operator * (const Complex tc) {
     Vector V(*this);
     for (unsigned int i=0; i<Nrow; i++) V.Data[i] *= tc;
     return V;
}

Complex Vector::operator * (const Vector& B) const {
     if (Nrow != B.Nrow) {cout<<"Vector operator * Vector : wrong vector length\n"; return 0.;}
     Complex tc = 0.;
#ifdef _OPENMP
     {
          if (Nrow < 1.e6) {
               double ta, taa, tb, tbb, *va, *vb, *vc, *vd;
               va = reinterpret_cast <double*>(Data); vb = reinterpret_cast <double*>(B.Data);
               vc = reinterpret_cast <double*>(&tc); vd=va+2*Nrow; 
               while (va<vd) {ta = *va++; taa = *va++; tb = *vb++; tbb = *vb++; *vc += ta*tb+taa*tbb; *(vc+1) += ta*tbb-taa*tb;} 
          }
          else {
               double tvr = 0., tvi = 0.;
#pragma omp parallel default(shared) reduction(+:tvr,tvi)
               {
                    int nt = omp_get_thread_num(), np = omp_get_num_threads(), nn = Nrow/np;
                    double ta, taa, tb, tbb, *va, *vb, *vd;
                    va = reinterpret_cast<double*>(Data) + 2*nt*nn; vb = reinterpret_cast<double*>(B.Data) + 2*nt*nn;
                    vd = (nt == np-1) ? (reinterpret_cast<double*>(Data) + 2*Nrow) : va + 2*nn;
                    while (va<vd) {
                         ta = *va++; taa = *va++;
                         tb = *vb++; tbb = *vb++;
                         tvr += ta*tb + taa*tbb;
                         tvi += ta*tbb - taa*tb;
                    }
               }
               tc = Complex(tvr,tvi);
          }
     }
#else
     {
          double ta, taa, tb, tbb, *va, *vb, *vc, *vd;
          va = reinterpret_cast <double*>(Data); vb = reinterpret_cast <double*>(B.Data);
          vc = reinterpret_cast <double*>(&tc); vd=va+2*Nrow; 
          while (va<vd) {ta = *va++; taa = *va++; tb = *vb++; tbb = *vb++; *vc += ta*tb+taa*tbb; *(vc+1) += ta*tbb-taa*tb;} 
     }
#endif
     return tc;
}

Vector Vector::operator * (const Matrix& B) {
     Vector C(B.Ncol);
     if (Nrow != B.Nrow) {cout<<"Vector operator * Matrix : wrong vector length\n"; return C;}
     double *va, *vb, *vc, *vaa, *vbb, *vcc, ta, taa, tb, tbb;
     va = reinterpret_cast <double*>(Data); vb = reinterpret_cast <double*>(B.Data); 
     vc = reinterpret_cast <double*>(C.Data); vbb = vb + 2*B.Nrow*B.Ncol;
     for (vcc=vc+2*B.Ncol; vc<vcc; vc+=2) {
          *vc = *(vc+1) = 0.; 
          for (vaa=vb; vaa<vbb; vaa += 2*B.Ncol) 
          {ta = *va++; taa = *va++; tb = *vaa; tbb = *(vaa+1); *vc += ta*tb-taa*tbb; *(vc+1) += ta*tbb+taa*tb;}
          va -= 2*Nrow; vb += 2;
     }
     return C;
}

Vector& Vector::operator /= (const Matrix& A) {
     if (A.Ncol != A.Nrow) {cout<<"Vector operator /= : matrix is not square\n"; return *this;}
     if (A.Nrow != Nrow) {cout<<"Vector operator /= : wrong matrix dimension\n"; return *this;}
     unsigned int i, ii;
     Matrix AA(Nrow); 
     for (i=0; i<Nrow; i++) for (ii=0; ii<Nrow; ii++) AA.Data[i*Nrow+ii] = A.Data[ii*Nrow+i]; 
     *this %= AA; return *this;
}

Vector& Vector::operator %= (const Matrix& A) {
     if (A.Ncol != A.Nrow) {cout<<"Vector operator %= : matrix is not square\n"; return *this;}
     if (A.Nrow != Nrow) {cout<<"Vector operator %= : wrong matrix dimension\n"; return *this;}

     double MaxMod, NewMod;
     const double NearZero = 1.e-15;
     Matrix R(A);
     unsigned int i, ii, j, jj/*, t = 2*Nrow*/;
     double *vva, *vvb, *vvc, *vv, va, vb, vc, vd, tv;
     Complex *a = R.Data, *b = Data, *aa, *bb, *aaa = R.Data + Nrow*Nrow - 1, ra; 

     for(i=0; i<Nrow-1; i++,a+=Nrow+1,b++) {
          MaxMod = abs(*a); j = i; vva = reinterpret_cast <double*>(a+Nrow);
          for (ii=i+1; ii<Nrow; ii++,vva+=2*Nrow) {
               va = *vva; vb = *(vva+1); NewMod = sqrt(va*va+vb*vb);
               if (MaxMod < NewMod) {MaxMod = NewMod; j = ii;}
          }
          if (MaxMod <= NearZero) {cout<<"Vector operator %= : matrix is zero\n"; return *this;}
          jj = 2*(Nrow-i);
          if (j > i) {
               vva = reinterpret_cast <double*>(a); vvc = reinterpret_cast <double*>(a+(j-i)*Nrow); vv = vvc+jj;
               while(vvc < vv) {tv = *vva; *vva++ = *vvc; *vvc++ = tv;}
               vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(b+(j-i)); 
               tv = *vva; *vva++ = *vvc; *vvc++ = tv; tv = *vva; *vva++ = *vvc; *vvc++ = tv;
          }
          jj -= 2; 
          //ra = -1./(*a);
          va = a->real(); vb = a->imag(); tv = va*va+vb*vb; ra = Complex(-va/tv,vb/tv);
          for (aa=a+Nrow,bb=b+1; aa<aaa; aa+=Nrow,bb++) {
               va = ra.real(); vb = ra.imag(); 
               tv = aa->real(); vd = aa->imag(); vc = va*tv-vb*vd; vd = va*vd+vb*tv;
               vva = reinterpret_cast <double*>(a+1); vvc = reinterpret_cast <double*>(aa+1); vv = vvc+jj;        
               while(vvc < vv) {va = *vva++; vb = *vva++; *vvc++ += va*vc-vb*vd; *vvc++ += va*vd+vb*vc;}
               vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(bb);         
               va = *vva++; vb = *vva; *vvc++ += va*vc-vb*vd; *vvc += va*vd+vb*vc;
          }
     }
     for (i=Nrow-1; i>=0; i--,a-=Nrow+1,b--) {
          //ra = 1./(*a); vc = ra.real(); vd = ra.imag(); 
          va = a->real(); vb = a->imag(); tv = va*va+vb*vb; vc = va/tv; vd = -vb/tv;
          vvc = reinterpret_cast <double*>(b);      
          va = *vvc; vb = *(vvc+1); *vvc++ = va*vc-vb*vd; *vvc = va*vd+vb*vc;
          vvb = reinterpret_cast <double*>(R.Data+i);    
          for (bb=Data; bb<b; vvb+=2*Nrow,bb++) {
               vc = *vvb; vd = *(vvb+1); 
               vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(bb);        
               va = *vva++; vb = *vva; *vvc++ -= va*vc-vb*vd; *vvc -= va*vd+vb*vc;
          }
     }
     return *this;
}

double Vector::normF(unsigned int nn) const {
     if ( nn >= Nrow ) return 0.;
     double tv = 0., tvr, tvi, *pv;
     pv = reinterpret_cast<double*>(this->Data) + 2*nn;
     for (unsigned int i=nn; i<Nrow; i++) {tvr = *pv++; tvi = *pv++; tv += tvr*tvr + tvi*tvi;}
     return sqrt(tv);
}

double Vector::cmp(const Vector &V) {
     if (Nrow != V.Nrow) return -1.;
     double tv = 0., tvv;
     for (unsigned int i=0; i<Nrow; i++) if (tv < (tvv = abs(V.Data[i] - Data[i]))) tv = tvv;
     return tv;
}

//////////////////////////////////////////////////////////////////////////

Matrix::Matrix(unsigned Mown, unsigned Mext) {
     Nrow = Mown;
     Ncol = (Mext) ? Mext : Mown;
     Data = new Complex[Nrow*Ncol];
}

Matrix::Matrix(const Matrix& M) {
     Nrow = M.Nrow; Ncol = M.Ncol;
     Data = new Complex[Nrow*Ncol];
     memcpy(Data, M.Data, Nrow*Ncol*sizeof(Complex));
}

Matrix::~Matrix() {delete [] Data;}

double Matrix::normF(int n1) const {
     double tv = 0.;
     //for (int i=n1; i<Nrow*Ncol; i++) tv += abs(Data[i]*Data[i]);
     if (n1 > -1) for (int i=0; i<n1; i++) tv += abs(Data[i]*Data[i]);
     else for (unsigned int i=0; i<Nrow*Ncol; i++) tv += abs(Data[i]*Data[i]);
     return sqrt(tv);//
}

double Matrix::cmp(const Matrix &M) {
     if ((Nrow != M.Nrow) || (Ncol != M.Ncol)) return -1.;
     double tv = 0., tvv;
     for (unsigned int i=0; i<Nrow*Ncol; i++) if (tv < (tvv = abs(M.Data[i] - Data[i]))) tv = tvv;
     return tv;
}

void Matrix::eye(void) {
     int N = (Nrow < Ncol) ? Nrow : Ncol;
     memset(Data,0,Nrow*Ncol*sizeof(Complex));
     for (int i=0; i<N; ++i) Data[i*Ncol+i] = 1.;
}

Matrix& Matrix::operator = (const Matrix& M) {
     if ((Nrow != M.Nrow)||(Ncol != M.Ncol))
          {if (Nrow*Ncol) delete[] Data; Nrow = M.Nrow; Ncol = M.Ncol; Data = new Complex[Nrow*Ncol];}
     memcpy(Data, M.Data, Nrow*Ncol*sizeof(Complex));
     return *this;
}

Matrix Matrix::transp() const {
     unsigned int i, j; Matrix M(Ncol,Nrow);
     for (i=0; i<Nrow; i++) for (j=0; j<Ncol; j++) M.Data[j*Nrow+i] = Data[i*Ncol+j];
     return M;
}

Matrix Matrix::mconj() const {
     unsigned int i, j; Matrix M(Ncol,Nrow);
     for (i=0; i<Nrow; i++) for (j=0; j<Ncol; j++) M.Data[j*Nrow+i] = conj(Data[i*Ncol+j]);
     return M;
}

Complex Matrix::trace() const {
     if (Nrow != Ncol) return 0.;
     Complex tc = 0.;
     for (unsigned int i=0; i<Nrow; i++) tc += Data[i*Nrow+i];
     return tc;
}

Matrix& Matrix::operator += (const Matrix& M) {
  Complex *D1 = Data;
  Complex *D2 = M.Data;
  Complex *D1end = Data + Nrow*Ncol;
     if ((*this).Nrow != M.Nrow) {cout<<"Matrix operator += : wrong dimension\n"; return *this;}
  if ((*this).Ncol != M.Ncol) {cout<<"Matrix operator += : wrong dimension\n"; return *this;}
  while (D1 < D1end) *D1++ += *D2++;
  return *this;
}

Matrix& Matrix::operator -= (const Matrix& M) {
  Complex *D1 = Data;
  Complex *D2 = M.Data;
  Complex *D1end = Data + Nrow*Ncol;
  if ((*this).Nrow != M.Nrow) {cout<<"Matrix operator -= : wrong dimension\n"; return *this;}
  if ((*this).Ncol != M.Ncol) {cout<<"Matrix operator -= : wrong dimension\n"; return *this;}
  while (D1 < D1end) *D1++ -= *D2++;
  return *this;
}

Vector Matrix::operator * (const Vector& B) const {
     Vector C(Nrow);
     //if (Nrow != Ncol) {cout<<"Matrix operator * Vector : wrong dimension\n"; return C;}
     if (Ncol != B.Nrow) {cout<<"Matrix operator * Vector : wrong dimension\n"; return C;}

     for (unsigned int i=0; i<Nrow; ++i) {C.Data[i] = 0.; for (unsigned int j=0; j<Ncol; ++j) C.Data[i] += Data[i*Ncol+j]*B(j);}
     return C;

     unsigned tt = 2*Nrow;
     double *va, *vb, *vc, *vbb, *vcc, ta, taa, tb, tbb;
     va = reinterpret_cast <double*>(Data); vb = reinterpret_cast <double*>(B.Data); 
     vc = reinterpret_cast <double*>(C.Data); 
     for (vcc = vc+2*Nrow; vc<vcc; vc+=2) {
          *vc = *(vc+1) = 0.; 
          for (vbb=vb+2*Ncol; vb<vbb;) 
               {ta = *va++; taa = *va++; tb = *vb++; tbb = *vb++; *vc += ta*tb-taa*tbb; *(vc+1) += ta*tbb+taa*tb;}
          vb -= tt;
     }
     return C;
}

Matrix Matrix::operator * (const Matrix& B) const {
     Matrix C(Nrow, B.Ncol);
     if (Ncol != B.Nrow) {cout<<"Matrix operator * Matrix : wrong dimension\n"; return C;}; // exception ???
//if ((Nrow < 0x200)&&(Nrow%2 == 0)&&(Ncol%2 == 0)&&(B.Ncol%2 == 0))
     {
          unsigned tt = 2*Nrow*C.Ncol, tc = 2*C.Ncol, t = 2*Ncol*C.Ncol;
          double *va, *vb, *vc, *vaa, *vbb, *vcc, tva, tvaa, tvb, tvbb;
          va = reinterpret_cast <double*>(Data); vb = reinterpret_cast <double*>(B.Data); vc = reinterpret_cast <double*>(C.Data); 
          for (vcc = vc+tt; vc<vcc;) *vc++ = *vc++ = 0.; vcc = vc-tt+tc;
          for (vbb=va+2*Nrow*Ncol; va<vbb; vcc+=tc, vb-=t) for (vaa=va+2*Ncol; va<vaa;) {
               tva = *va++; tvaa = *va++; vc = vcc-tc; 
               while (vc < vcc) {
                    tvb = *vb++; tvbb = *vb++; *vc++ += tva*tvb-tvaa*tvbb; *vc++ += tva*tvbb+tvaa*tvb;
                    //*vc++ += tva*(tvb = *vb++)-tvaa*(tvbb = *(++vb)); *vc++ += tva*tvbb+tvaa*tvb;
                    //tvb = tva*(tv = *vb++); tvbb = tvaa*(tvv = *vb++); *vc++ += tvb-tvbb; *vc++ += (tv+tva)*(tvv+tvaa)-tvb-tvbb;
               }
          }
     }
     return C;
     //else
     {  // Strassen
          int i, ii, n1 = Nrow/2, n2 = Ncol/2, n3 = B.Ncol/2; 
          Matrix S(n1,n2), SS(n2,n3), M1(n1,n3), M2(n1,n3), M3(n1,n3), M4(n1,n3), M5(n1,n3);
          for (ii=0; ii<n2; ii++) {
               for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[i*Ncol+ii]-Data[(i+n1)*Ncol+ii];
               for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[(ii+n2)*B.Ncol+n3+i]-B.Data[ii*B.Ncol+n3+i];
          }     
          M4 = S*SS;
          for (ii=0; ii<n2; ii++) {
               for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[i*Ncol+n2+ii];
               for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[(ii+n2)*B.Ncol+i];
          }     
          M3 = S*SS;
          for (ii=0; ii<n2; ii++) {
               for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[i*Ncol+ii];
               for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[ii*B.Ncol+i];
          }     
          M2 = S*SS;
          for (ii=0; ii<n2; ii++) {
               for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[(i+n1)*Ncol+ii]+Data[(i+n1)*Ncol+n2+ii];
               for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[ii*B.Ncol+n3+i]-B.Data[ii*B.Ncol+i];
          }     
          M5 = S*SS;
          for (ii=0; ii<n2; ii++) {
               for (i=0; i<n1; i++) S.Data[i*n2+ii] -= Data[i*Ncol+ii];
               for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[(ii+n2)*B.Ncol+n3+i]-SS.Data[ii*n3+i];
          }     
          M1 = S*SS+M2; M4 += M1; M1 += M5;
          for (i=0; i<n1; i++) for (ii=0; ii<n3; ii++) {
               C.Data[i*C.Ncol+ii] = M2.Data[i*n3+ii]+M3.Data[i*n3+ii];
               C.Data[(i+n1)*C.Ncol+n3+ii] = M4.Data[i*n3+ii]+M5.Data[i*n3+ii];
          }
          M5 = SS;
          for (ii=0; ii<n2; ii++) {
               for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[i*Ncol+n2+ii]-S.Data[i*n2+ii];
               for (i=0; i<n3; i++) SS.Data[ii*n3+i] = B.Data[(ii+n2)*B.Ncol+n3+i];
          }     
          M2 = S*SS;
          for (ii=0; ii<n2; ii++) {
               for (i=0; i<n1; i++) S.Data[i*n2+ii] = Data[(i+n1)*Ncol+n2+ii];
               for (i=0; i<n3; i++) SS.Data[ii*n3+i] = M5.Data[ii*n3+i]-B.Data[(ii+n2)*B.Ncol+i];
          }     
          M3 = S*SS;
          for (i=0; i<n1; i++) for (ii=0; ii<n3; ii++) {
               C.Data[i*C.Ncol+n3+ii] = M1.Data[i*n3+ii]+M2.Data[i*n3+ii];
               C.Data[(i+n1)*C.Ncol+ii] = M4.Data[i*n3+ii]-M3.Data[i*n3+ii];
          }
     }
     return C;
}

Matrix& Matrix::operator /= (const Matrix& A) {
     if (this == &A) {cout<<"Matrix operator /= : wrong divider\n"; return *this;}
     if (A.Ncol != A.Nrow) {cout<<"Matrix operator /= : divider is not square\n"; return *this;}
     if (A.Nrow != Ncol) {cout<<"Matrix operator /= : wrong divider dimension\n"; return *this;}
     unsigned int i, ii;
     Matrix R(Ncol,Nrow), AA(Ncol);
     for (i=0; i<Ncol; i++) for (ii=0; ii<Nrow; ii++) R.Data[i*Nrow+ii] = Data[ii*Ncol+i];
     for (i=0; i<Ncol; i++) for (ii=0; ii<Ncol; ii++) AA.Data[i*Ncol+ii] = A.Data[ii*Ncol+i]; 
     R %= AA;
     for (i=0; i<Nrow; i++) for (ii=0; ii<Ncol; ii++) Data[i*Ncol+ii] = R.Data[ii*Nrow+i];
     return *this;
}

Matrix& Matrix::operator %= (const Matrix& A) {
  if (this == &A) {cout<<"Matrix operator /= : wrong argument\n"; return *this;}
  if (A.Ncol != A.Nrow) {cout<<"Matrix operator /= : argument is not square\n"; return *this;}
  if (A.Nrow != Nrow) {cout<<"Matrix operator /= : wrong argument dimension\n"; return *this;}
     //if (Nrow < 0x200)
     {
          double MaxMod, NewMod;
          const double NearZero = 1.e-15;
          Matrix R(A);
          unsigned int i, ii, j, jj, t = 2*Ncol;
          double *vva, *vvb, *vvc, *vv, va, vb, vc, vd, tv;
          Complex *a = R.Data, *b = Data, *aa, *bb, *aaa = R.Data+Nrow*Nrow-1, ra; 

          for(i=0; i<Nrow-1; i++,a+=Nrow+1,b+=Ncol) {
               MaxMod = abs(*a); j = i; vva = reinterpret_cast <double*>(a+Nrow);
               for (ii=i+1; ii<Nrow; ii++,vva+=2*Nrow) {
                    va = *vva; vb = *(vva+1); NewMod = sqrt(va*va+vb*vb);
                    if (MaxMod < NewMod) {MaxMod = NewMod; j = ii;}
               }
               if (MaxMod <= NearZero) return *this; // exception ???
               jj = 2*(Nrow-i);
               if (j > i) {
                    vva = reinterpret_cast <double*>(a); vvc = reinterpret_cast <double*>(a+(j-i)*Nrow); vv = vvc+jj;
                    while(vvc < vv) {tv = *vva; *vva++ = *vvc; *vvc++ = tv;}
                    vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(b+(j-i)*Ncol); vv = vvc+t;
                    while(vvc < vv) {tv = *vva; *vva++ = *vvc; *vvc++ = tv;}
               }
               jj -= 2; 
               //ra = -1./(*a);
               va = a->real(); vb = a->imag(); tv = va*va+vb*vb; ra = Complex(-va/tv,vb/tv);
               for (aa=a+Nrow,bb=b+Ncol; aa<aaa; aa+=Nrow,bb+=Ncol) {
                    va = ra.real(); vb = ra.imag(); 
                    tv = aa->real(); vd = aa->imag(); vc = va*tv-vb*vd; vd = va*vd+vb*tv;
                    vva = reinterpret_cast <double*>(a+1); vvc = reinterpret_cast <double*>(aa+1); vv = vvc+jj;        
                    while(vvc < vv) {va = *vva++; vb = *vva++; *vvc++ += va*vc-vb*vd; *vvc++ += va*vd+vb*vc;}
                    vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(bb); vv = vvc+t;        
                    while(vvc < vv) {va = *vva++; vb = *vva++; *vvc++ += va*vc-vb*vd; *vvc++ += va*vd+vb*vc;}
               }
          }
          for (i=Nrow-1; i>=0; i--,a-=Nrow+1,b-=Ncol) {
               //ra = 1./(*a); vc = ra.real(); vd = ra.imag(); 
               va = a->real(); vb = a->imag(); tv = va*va+vb*vb; vc = va/tv; vd = -vb/tv;
               vvc = reinterpret_cast <double*>(b); vv = vvc+t;        
               while(vvc < vv) {va = *vvc; vb = *(vvc+1); *vvc++ = va*vc-vb*vd; *vvc++ = va*vd+vb*vc;}
               vvb = reinterpret_cast <double*>(R.Data+i);    
               for (bb=Data; bb<b; vvb+=2*Nrow,bb+=Ncol) {
                    vc = *vvb; vd = *(vvb+1); 
                    vva = reinterpret_cast <double*>(b); vvc = reinterpret_cast <double*>(bb); vv = vvc+t;        
                    while(vvc < vv) {va = *vva++; vb = *vva++; *vvc++ -= va*vc-vb*vd; *vvc++ -= va*vd+vb*vc;}
               }
          }
          return *this;
     }
}
