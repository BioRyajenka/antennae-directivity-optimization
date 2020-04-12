
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

#include <memory.h>
#include "sphereml.h"
#include "spfunc.h"

Vector SphereML::calc_pw(double as, double ap, double th, double ph) {
    int n, m; double /*tv,*/ tp, tt; Complex tc = j_, tcc, te; Vector VA(2*N*N);
     memset(VA.Data,0,2*N*N*sizeof(Complex));
     for (n=1; n<N; ++n) {
          tcc = tc*4.*sqrt(M_PI)/sqrt(2.*n*(n+1.));
          for (m=-n; m<n+1; ++m) {
               te = exp(-j_*double(m)*ph);
               tp = LPin(th,n,m); tt = LTaun(th,n,m);
               VA.Data[n*(n+1)+m] = -tcc*(j_*tp*ap + as*tt)*te;
               VA.Data[N*N+n*(n+1)+m] = -tcc*(j_*tt*ap + as*tp)*te;
          }
          tc *= j_;
     }
     return VA;
}

Vector SphereML::calc_edz(double px, double py, double pz, Complex krz, int in) {
    int n/*, m*/;
     double tv = -0.25/sqrt(M_PI), tvn/*, tv1, tv2*/;
     Complex pp = Complex(px,py), pm = conj(pp), zf, zfd;
     Vector VA(2*N*N);
     memset(VA.Data,0,2*N*N*sizeof(Complex));

     if (in == 1) { // field inside dipole radius
          for (n=1; n<N; ++n) {
               tvn = tv*sqrt(2*n+1.);
               zf = besh1(krz,n); zfd = besh1d(krz,n);
               VA.Data[n*(n+1)-1] = pm*( VA.Data[n*(n+1)+1] = tvn*zf );
               VA.Data[n*(n+1)+1] *= pp;
               VA.Data[N*N+n*(n+1)+0] = -2.*j_*tvn*pz*sqrt(n*(n+1.))*zf/krz;
               VA.Data[N*N+n*(n+1)+1] = -pp*( VA.Data[N*N+n*(n+1)-1] = j_*tvn*(zfd + zf/krz) );
               VA.Data[N*N+n*(n+1)-1] *= pm;
               tv = -tv;
          }
     }
     else { // field outside dipole radius
          for (n=1; n<N; ++n) {
               tvn = tv*sqrt(2*n+1.);
               zf = besj(krz,n); zfd = besjd(krz,n);
               VA.Data[n*(n+1)-1] = pm*( VA.Data[n*(n+1)+1] = tvn*zf );
               VA.Data[n*(n+1)+1] *= pp;
               VA.Data[N*N+n*(n+1)+0] = -2.*j_*tvn*pz*sqrt(n*(n+1.))*zf/krz;
               VA.Data[N*N+n*(n+1)+1] = -pp*( VA.Data[N*N+n*(n+1)-1] = j_*tvn*(zfd + zf/krz) );
               VA.Data[N*N+n*(n+1)-1] *= pm;
               tv = -tv;
          }
     }

     return VA;
}

Vector SphereML::calc_far(const Vector &V, double th, double ph) {
     int m, n, NN = N*N; double tv;
     Complex tc, tc1, tc2, tc3, tc4, *te;
     Vector VE(2);
     te = new Complex [N]; for (m=0; m<N; ++m) te[m] = exp(j_*double(m)*ph);
     tc = -j_; VE.Data[0] = VE.Data[1] = 0.;
     for (n=1; n<N; ++n) {
          tc1 = V(n*(n+1))*LPin(th,n,0); tc2 = V(NN+n*(n+1))*LTaun(th,n,0);
          tc3 = V(n*(n+1))*LTaun(th,n,0); tc4 = V(NN+n*(n+1))*LPin(th,n,0);
          for (m=1; m<n+1; ++m) {
               tc1 += V(n*(n+1)-m)*LPin(th,n,-m)*conj(te[m]) + V(n*(n+1)+m)*LPin(th,n,m)*te[m]; // ae*pi
               tc2 += V(NN+n*(n+1)-m)*LTaun(th,n,-m)*conj(te[m]) + V(NN+n*(n+1)+m)*LTaun(th,n,m)*te[m]; // ah*tau
               tc3 += V(n*(n+1)-m)*LTaun(th,n,-m)*conj(te[m]) + V(n*(n+1)+m)*LTaun(th,n,m)*te[m]; // ae*tau
               tc4 += V(NN+n*(n+1)-m)*LPin(th,n,-m)*conj(te[m]) + V(NN+n*(n+1)+m)*LPin(th,n,m)*te[m]; // ah*pi
          }
          tv = 1./sqrt(n*(n+1.));
          VE.Data[0] -= tc*tv*(tc1 + tc2);
          tc *= -j_; VE.Data[1] += tc*tv*(tc3 + tc4);
     }
     VE.Data[0] *= M_SQRT1_2PI; VE.Data[1] *= M_SQRT1_2PI;
     delete [] te;
     return VE;
}

double SphereML::calc_Psca(const Vector &VS, double tC) {
     int n, NN = N*N; double tv = 0.;
     for (n=1; n<NN; ++n) tv += abs(VS(n)*VS(n)) + abs(VS(NN+n)*VS(NN+n));
     return 0.5*tv/tC;
}

double SphereML::calc_Pext(const Vector &VI, const Vector &VS, double tC) {
     int n, NN = N*N; double tv = 0.;
     for (n=1; n<NN; ++n) tv += (VS(n)*conj(VI(n))).real() + (VS(NN+n)*conj(VI(NN+n))).real();
     return -tv/tC;
}

double SphereML::directivity(const Vector &VS, double th, double ph, double tC) {
     int n, m, nm, NN = N*N;
     double tp, tt;
     Complex tc1, tc2, tc3, tc4, tc = -j_, tcc, tce;
     tc1 = tc2 = tc3 = tc4 = 0.;
     for (n=nm=1; n<N; ++n) {
          tcc = tc/sqrt(double(n*(n+1)));
          for (m=-n; m<n+1; m++,nm++) {
               tp = LPin(th,n,m); tt = LTaun(th,n,m);
               tce = exp(j_*double(m)*ph);
               tc1 += tcc*tce*(VS(nm)*tp + VS(nm+NN)*tt);
               tc2 += tcc*tce*(VS(nm)*tt + VS(nm+NN)*tp);
          }
          tc *= -j_;
     }
     return (tc1*conj(tc1) + tc2*conj(tc2)).real()/calc_Psca(VS,tC)/tC;
}

Matrix SphereML::calc_RT(double kr, Complex e1, Complex e2, Complex m1, Complex m2) {
     Complex kR1, kR2, te, tm, tc;
     kR1 = kr*sqrt(e1*m1); if (arg(kR1) < -1.e-8) kR1 = -kR1;
     kR2 = kr*sqrt(e2*m2); if (arg(kR2) < -1.e-8) kR2 = -kR2;
     te = e1/e2; tm = 1.;//m1/m2;
     Matrix M(4,2*N);//0 1 2 3 -> 00 01 10 11
     for (int n=0; n<N; ++n) {
          tc = 1./(besj(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzj(kR1,n));
          M.Data[0*N+n] = tc*(besh1(kR2,n)*bes_dzh1(kR1,n) - besh1(kR1,n)*bes_dzh1(kR2,n)); // 00e
          M.Data[2*N+n] = tc*j_/kR1; // 01e
          M.Data[6*N+n] = tc*(besj(kR2,n)*bes_dzj(kR1,n) - besj(kR1,n)*bes_dzj(kR2,n)); // 11e
          M.Data[4*N+n] = tc*j_/kR2; // 10e
          tc = 1./(te*besj(kR1,n)*bes_dzh1(kR2,n) - besh1(kR2,n)*bes_dzj(kR1,n));
          M.Data[1*N+n] = tc*(besh1(kR2,n)*bes_dzh1(kR1,n) - te*besh1(kR1,n)*bes_dzh1(kR2,n)); // 00h
          M.Data[3*N+n] = tc*j_/kR2; // 01h
          M.Data[7*N+n] = tc*(besj(kR2,n)*bes_dzj(kR1,n) - te*besj(kR1,n)*bes_dzj(kR2,n)); // 11h
          M.Data[5*N+n] = tc*j_*te/kR1; // 10h
     }
     return M;
}

Matrix SphereML::calc_SML(Matrix **SM, int ns) {
     int n, k; Complex tc;
     Matrix SML1(4,2*N), SML2(4,2*N);
     memcpy(SML2.Data,SM[0]->Data,8*N*sizeof(Complex));
     for (k=1; k<ns; ++k) {
          SML1 = SML2;
          for (n=0; n<N; ++n) {
               tc = 1./(1. - SML1(3,n)*(*SM[k])(0,n));
               SML2.Data[n+0*N] = SML1(0,n) + tc*SML1(1,n)*SML1(2,n)*(*SM[k])(0,n); // 00e
               SML2.Data[n+2*N] = tc*SML1(1,n)*(*SM[k])(1,n); // 01e
               SML2.Data[n+4*N] = tc*SML1(2,n)*(*SM[k])(2,n); // 10e
               SML2.Data[n+6*N] = (*SM[k])(3,n) + tc*(*SM[k])(1,n)*(*SM[k])(2,n)*SML1(3,n); // 11e
               tc = 1./(1. - SML1(3,n+N)*(*SM[k])(0,n+N));
               SML2.Data[n+1*N] = SML1(0,n+N) + tc*SML1(1,n+N)*SML1(2,n+N)*(*SM[k])(0,n+N); // 00h
               SML2.Data[n+3*N] = tc*SML1(1,n+N)*(*SM[k])(1,n+N); // 01h
               SML2.Data[n+5*N] = tc*SML1(2,n+N)*(*SM[k])(2,n+N); // 10h
               SML2.Data[n+7*N] = (*SM[k])(3,n+N) + tc*(*SM[k])(1,n+N)*(*SM[k])(2,n+N)*SML1(3,n+N); // 11h
          }
     }
     return SML2;
}
