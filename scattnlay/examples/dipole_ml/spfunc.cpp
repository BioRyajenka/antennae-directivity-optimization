
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

#include "./spfunc.h"

#include <complex>
#include <algorithm>
#include <math.h>
#define _USE_MATH_DEFINES

#define DG 12

void xyz2rtp(double x, double y, double z, double &r, double &th, double &ph) {
     double tv = x*x + y*y;
     r = sqrt(tv+z*z); tv = sqrt(tv);
     if (r < 1.e-14) {th = ph = 0.; return;}
     th = acos(z/r);
     if (tv < 1.e-14) {ph = 0.; return;}
     ph = (y > 0.) ? acos(x/tv) : (M_PI + acos(x/tv));
}

     // spherical Bessel functions

Complex besj(Complex z, int n) {
     if (n == 0) return besj0(z);
     else if (n == 1) return besj1(z);
//     else if (n < 0) return (n%2) ? besy(z,-n-1) : -besy(z,-n-1);

     if (abs(z) < 0.6) {
          int i, ip;
          double tvv;
          Complex tc, tcc;

          tc = 1.; for (i=1; i<n+1; i++) tc *= z/double(2*i+1); i = 1; tcc = tc;
          ip = int(0.4343*log(abs(tc))) - DG; tvv = (ip < -100) ? 1.e-100 : exp(2.3*ip);
          do {tcc += ( tc *= -0.5*z*z/double(i*(2*(n+i)+1)) ); i++;} while (abs(tc) > tvv);

          return tcc;
     }
     else if (abs(z)/double(n) < 1.) {
          int nn, i, pw;
          double tv1, tv2, tv3, tx, tx1, tx2, ty1, ty2;
          Complex tc, tc1, tc2, tcc;

          if (n > int(abs(z))) pw = abs(int( (n*log(0.5*M_E*abs(z)/n) - 0.5*log(n*abs(z)) - M_LN2)/M_LN10 )) + DG;
          else pw = DG;

          tv1 = abs(z); tv2 = 2./M_E/tv1; tv3 = pw*M_LN10 - M_LN2;
          tx2 = (tx1 = n) + 2;
          ty1 = tx1*log(tv2*tx1) + 0.5*log(tv1*tx1) - tv3;
          ty2 = tx2*log(tv2*tx2) + 0.5*log(tv1*tx2) - tv3;
          do {
               tx = tx2 - ty2*(tx2-tx1)/(ty2-ty1); tx1 = tx2; tx2 = tx;
               ty1 = ty2; ty2 = tx2*log(tv2*tx2) + 0.5*log(tv1*tx2) - tv3;
          } while (fabs(tx2-tx1) > 0.5);
          nn = int(tx2);

          tc1 = 0.; tc2 = exp(-pw*M_LN10);
          for (i=nn; i>0; i--) {
               tc = (2.*i+1.)/z*tc2 - tc1; tc1 = tc2; tc2 = tc;
               if (i == n+1) tcc = tc2;
          }
          return (abs(besj1(z)) > abs(besj0(z))) ? tcc*besj1(z)/tc1 : tcc*besj0(z)/tc2;
     }
     else {
          int i, tm;
          double tvv;
          Complex tc, tp, tq, pp, qq;
          i = -int(0.4343*log(abs(z)/n)) - DG; tvv = (i < -100) ? 1.e-100 : exp(2.3*i);
          tp = pp = 1.; tm = (2*n+1)*(2*n+1); tq = qq = double(tm-1)*(tc = 0.125/z); tc *= tc; i = 1;
          do {
               pp += ( tp *= -tc*double((tm - (4*i-1)*(4*i-1))*(tm - (4*i-3)*(4*i-3)))/double(2*i*(2*i-1)) );
               qq += ( tq *= -tc*double((tm - (4*i+1)*(4*i+1))*(tm - (4*i-1)*(4*i-1)))/double(2*i*(2*i+1)) );
               i++;
          } while ((abs(tp) > tvv) && (abs(tq) > tvv));
          return (pp*cos(z-M_PI_2*(n+1)) - qq*sin(z-M_PI_2*(n+1)))/z;
     }
}

Complex besjd(Complex z, int n) {
     if (n == 0) return besj0d(z);
     else if (n == 1) return besj1d(z);
//     else if (n < 0) return (n%2) ? besy(z,-n-1) : -besy(z,-n-1);

     if (abs(z) < 0.6) {
          int i, ip;
          double tv, tvv;
          Complex tc, tcc;

          tc = 1/3.; for (i=2; i<n+1; i++) tc *= z/double(2*i+1); i = 1; tcc = double(n)*tc;
          ip = int(0.4343*log(abs(tc))) - DG; tvv = (ip < -100) ? 1.e-100 : exp(2.3*ip);
          do {tv = abs( tc *= -0.5*z*z/double(i*(2*(n+i)+1)) ); tcc += double(n+2*i)*tc; i++;} while (tv > tvv);

          return tcc;
     }
     else if (abs(z)/n < 1.) {
          int nn, i, pw;
          double tv1, tv2, tv3, tx, tx1, tx2, ty1, ty2;
          Complex tc, tc1, tc2, tcc;

          if (n > int(abs(z))) pw = abs(int( (n*log(0.5*M_E*abs(z)/n) - 0.5*log(n*abs(z)) - M_LN2)/M_LN10 )) + DG;
          else pw = DG;

          tv1 = abs(z); tv2 = 2./M_E/tv1; tv3 = pw*M_LN10 - M_LN2;
          tx2 = (tx1 = n) + 2;
          ty1 = tx1*log(tv2*tx1) + 0.5*log(tv1*tx1) - tv3;
          ty2 = tx2*log(tv2*tx2) + 0.5*log(tv1*tx2) - tv3;
          do {
               tx = tx2 - ty2*(tx2-tx1)/(ty2-ty1); tx1 = tx2; tx2 = tx;
               ty1 = ty2; ty2 = tx2*log(tv2*tx2) + 0.5*log(tv1*tx2) - tv3;
          } while (fabs(tx2-tx1) > 0.5);
          nn = int(tx2);

          tc1 = 0.; tc2 = exp(-pw*M_LN10);
          for (i=nn; i>0; i--) {
               tc = (2.*i+1.)/z*tc2 - tc1; tc1 = tc2; tc2 = tc;
               if (i == n+1) tcc = double(n)/z*tc2 - tc1;
          }
          return (abs(besj1(z)) > abs(besj0(z))) ? tcc*besj1(z)/tc1 : tcc*besj0(z)/tc2;
     }
     else {
          return (double(n)*besj(z,n-1) - double(n+1)*besj(z,n+1))/double(2*n+1);
     }
}

Complex besy(Complex z, int n) {
     if (n == 0) return besy0(z);
     else if (n == 1) return besy1(z);
//     else if (n < 0) return (abs(n)%2) ? besj(z,n) : -besj(z,n);

     if (abs(z) < 0.6) {
          int i, ip;
          double tv, tvv;
          Complex tc, tcc;

          tc = -1./z; for (i=1; i<n+1; i++) tc *= double(2*i-1)/z; i = 1; tcc = tc;
          ip = int(0.4343*log(abs(tc))) - DG; tvv = (ip < -100) ? 1.e-100 : exp(2.3*ip);
          do {tv = abs( tc *= 0.5*z*z/double(2*(n-i+1)-1)/double(i) ); tcc += tc; i++;} while (tv > tvv);

          return tcc;
     }
     else if (abs(z)/n < 1.) {
          Complex tc, tc1 = besy0(z), tc2 = besy1(z);
          for (int i=2; i<n+1; i++) {
               tc = double(2*i-1)*tc2/z - tc1; tc1 = tc2; tc2 = tc;
          }
          return tc2;
     }
     else {
          int i, tm;
          double tvv;
          Complex tc, tp, tq, pp, qq;
          i = -int(0.4343*log(abs(z))) - DG; tvv = (i < -100) ? 1.e-100 : exp(2.3*i);
          tp = pp = 1.; tm = (2*n+1)*(2*n+1); tq = qq = double(tm-1)*(tc = 0.125/z); tc *= tc; i = 1;
          do {
               pp += ( tp *= -tc*double((tm - (4*i-1)*(4*i-1))*(tm - (4*i-3)*(4*i-3)))/double(2*i*(2*i-1)) );
               qq += ( tq *= -tc*double((tm - (4*i+1)*(4*i+1))*(tm - (4*i-1)*(4*i-1)))/double(2*i*(2*i+1)) );
               i++;
          } while ((abs(tp) > tvv) && (abs(tq) > tvv));
          return (pp*sin(z-M_PI_2*(n+1)) + qq*cos(z-M_PI_2*(n+1)))/z;
     }
}

Complex besyd(Complex z, int n) {
     if (n == 0) return besy0d(z);
     else if (n == 1) return besy1d(z);
//     else if (n < 0) return (abs(n)%2) ? besj(z,n) : -besj(z,n);

     if (abs(z) < 0.6) {
          int i, ip;
          double tvv;
          Complex tc, tcc;
          tc = 1./z/z; for (i=1; i<n+1; i++) tc *= double(2*i-1)/z; i = 1; tcc = double(n+1)*tc;
          ip = int(0.4343*log(abs(tc))) - DG; tvv = (ip < -100) ? 1.e-100 : exp(2.3*ip);
          do {tcc -= double(2*i-n-1)*( tc *= 0.5*z*z/double(2*(n-i+1)-1)/double(i) ); i++;} while (abs(tc) > tvv);
          return tcc;
     }
     else if (abs(z)/n < 1.) {
          Complex tc, tc1 = besy0(z), tc2 = besy1(z);
          for (int i=2; i<n+1; i++) {
               tc = double(2*i-1)*tc2/z - tc1; tc1 = tc2; tc2 = tc;
          }
          return (tc1 - double(n+1)/z*tc2);
     }
     else return (double(n)*besy(z,n-1) - double(n+1)*besy(z,n+1))/double(2*n+1);
}

Complex besh1(Complex z, int n) {
     if (n == 0) return besh10(z);
     else if (n == 1) return besh11(z);

//     if (15.*abs(z) < n)
          return (besj(z,n) + j_*besy(z,n));
//     else {
//          Complex tc, tc1, tc2;
//          int i = 1; tc1 = besh10(z); tc2 = besh11(z);
//          do {tc = double(2*i+1)*tc2/z-tc1; tc1 = tc2; tc2 = tc;} while (++i < n);
//          return tc;
//     }
}

Complex besh2(Complex z, int n) {
     if (n == 0) return besh20(z);
     else if (n == 1) return besh21(z);

//     if (15.*abs(z) < n)
          return (besj(z,n) - j_*besy(z,n));
//     else {
//          Complex tc, tc1, tc2;
//          int i = 1; tc1 = besh20(z); tc2 = besh21(z);
//          do {tc = double(2*i+1)*tc2/z-tc1; tc1 = tc2; tc2 = tc;} while (++i < n);
//          return tc;
//     }
}

Complex besh1d(Complex z, int n) {
     if (n == 0) return besh10d(z);
     else if (n == 1) return besh11d(z);

     return (besjd(z,n) + j_*besyd(z,n));
//     if (15.*abs(z) < n)
//     else {
//          Complex tc, tc1, tc2;
//          int i = 1; tc1 = besh10(z); tc2 = besh11(z);
//          do {tc = double(2*i+1)*tc2/z-tc1; tc1 = tc2; tc2 = tc;} while (++i < n);
//          return (tc1 - double(n+1)/z*tc2);
//     }
}

Complex besh2d(Complex z, int n) {
     if (n == 0) return besh20d(z);
     else if (n == 1) return besh21d(z);

//     if (15.*abs(z) < n)
          return (besjd(z,n) - j_*besyd(z,n));
//     else {
//          Complex tc, tc1, tc2;
//          int i = 1; tc1 = besh20(z); tc2 = besh21(z);
//          do {tc = double(2*i+1)*tc2/z-tc1; tc1 = tc2; tc2 = tc;} while (++i < n);
//          return (tc1 - double(n+1)/z*tc2);
//     }
}

     // Legendre polynomials

double pLegn(double t, int nn) {
     if (nn < 0) return 0;
     else if (fabs(t) < 1.e-14) return sqrt(0.5*(2*nn+1));
     else if (fabs(t-M_PI) < 1.e-14) return (nn%2) ? -sqrt(0.5*(2*nn+1)) : sqrt(0.5*(2*nn+1));
     else switch (nn) {
          case 0: return pLegn0(t);
          case 1: return pLegn1(t);
          default: {
               int n = 1; double tv1 = pLegn0(t), tv2 = pLegn1(t), tv, t2 = sqrt(3.), t1 = 1.;
               do {tv = (t2*tv2*cos(t) - n/t1*tv1)/(n+1.); t1 = t2; t2 = sqrt(2*n+3.); tv1 = tv2; tv2 = tv*t2;} while (++n < nn);
               return tv2;
          }
     }
}

double pLegnd(double t, int nn) {
     switch (nn) {
          case 0: return pLegnd0(t);
          case 1: return pLegnd1(t);
          default: {
               int n = 1; double tv1 = pLegn0(t), tv2 = pLegn1(t), tv, t2 = sqrt(3.), t1 = 1., t0, td0, td1 = 0., td2 = sqrt(1.5);
               do {
                    t0 = t1; t1 = t2; t2 = sqrt(2*n+3.);
                    tv = t2*(t1*tv2*cos(t) - n/t0*tv1)/(n+1.); tv1 = tv2; tv2 = tv*t2;
                    td0 = td1; td1 = td2; td2 = t2*(t1*tv1 + td0/t0);
               } while (++n < nn);
               return td2;
          }
     }
}

double paLegn(double t, int n, int m) {
     if (abs(m) > n) return 0.;
     else if (m == 0) return pLegn(t,n);
     else {
          if (n == 1) return (m > 0) ? paLegn11(t) : paLegn1m1(t);
          else if ((fabs(t) < 1.e-14) || (fabs(t-M_PI) < 1.e-14)) return 0.;
          else if (m < 0) return (abs(m)%2) ? -paLegn(t,n,abs(m)) : paLegn(t,n,abs(m));
          else {
               int i; double tv1, tv2, tv, t1, t2, tc = cos(t);
               tv = log(m+0.5); for (i=2; i<2*m+1; ++i) tv -= log(double(i)); tv2 = exp(0.5*tv);
               tv = 0.; for (i=2; 2*i-1<2*m; ++i) tv += log(double(2*i-1));
               tv2 *= exp(tv)*exp(m*log(sin(t))); //tv2 = exp(tv + m*log(sin(t)));
               if (m == n) return tv2;
               tv1 = t2 = 0.; i = m;
               do {
                    t1 = t2; t2 = sqrt((i+m+1.)*(i-m+1.)/(2*i+1.)/(2*i+3));
                    tv = (tc*tv2 - t1*tv1)/t2; tv1 = tv2; tv2 = tv;
               } while (++i < n);
               return tv2;
          }
     }
}

double paLegnd(double t, int n, int m) {
     if (m < 0) return (abs(m)%2) ? -paLegnd(t,n,abs(m)) : paLegnd(t,n,abs(m)); 
     if (fabs(t) < 1.e-14) {
          switch (m) {
               case 0: return 0.5*n*(n+1)*sqrt(n+0.5);
               case 1: return -std::numeric_limits<double>::infinity();
               case 2: return -0.25*sqrt((n+0.5)*(n-1)*n*(n+1)*(n+2));
               default: return 0.;
          }
     }
     else if (fabs(t-M_PI) < 1.e-14) {
          switch (m) {
               case 0: return (n%2) ? 0.5*n*(n+1)*sqrt(n+0.5) : -0.5*n*(n+1)*sqrt(0.5*(2*n+1.));
               case 1: return (n%2) ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
               case 2: return (n%2) ? -0.25*sqrt((n+0.5)*(n-1)*n*(n+1)*(n+2)) : 0.25*sqrt((n+0.5)*(n-1)*n*(n+1)*(n+2));
               default: return 0.;
          }
     }
     return ((n+1)*cos(t)*paLegn(t,n,m) - sqrt((n-m+1.)*(n+m+1.)*(2*n+1.)/(2*n+3.))*paLegn(t,n+1,m))/sin(t)/sin(t);
     //return (-double(m)*cos(t)*paLegn(t,n,m) + sin(t)*sqrt((n-m)*(n+m+1.))*paLegn(t,n,m+1))/sin(t)/sin(t);
/**
     else {
          int i; double tv1, tv2, tv, t1, t2, tc = cos(t);
          tv = log(m+0.5); for (i=2; i<2*m+1; ++i) tv -= log(double(i));
          tv *= 0.5; for (i=2; 2*i-1<2*m; ++i) tv += log(double(2*i-1));
          tv2 = exp(tv + m*log(sin(t))); if (m == n) return tv2;
          tv1 = t2 = 0.; i = m;
          do {
               t1 = t2; t2 = sqrt((i+m+1.)*(i-m+1.)/(2*i+1.)/(2*i+3));
               tv = (tc*tv2 - t1*tv1)/t2; tv1 = tv2; tv2 = tv;
          } while (++i < (n+1));
          return ((n+1)*cos(t)*tv1 - sqrt((n-m+1.)*(n+m+1.)*(2*n+1.)/(2*n+3.))*tv2)/sin(t)/sin(t);
     }
**/
}

double LPin(double t, int n, int m) {
     if (m < 0) return (abs(m)%2) ? LPin(t,n,abs(m)) : -LPin(t,n,abs(m)); 
     if (fabs(t) < 1.e-14) {
          switch (m) {
               case 1: return 0.5*sqrt(n*(n+1)*(n+0.5));
               case -1: return 0.5*sqrt(n*(n+1)*(n+0.5));
               default: return 0.;
          }
     }
     else if (fabs(t-M_PI) < 1.e-14) {
          switch (m) {
               case 1: return (n%2) ? 0.5*sqrt(n*(n+1)*(n+0.5)) : -0.5*sqrt(n*(n+1)*(n+0.5));
               case -1: return (n%2) ? 0.5*sqrt(n*(n+1)*(n+0.5)) : -0.5*sqrt(n*(n+1)*(n+0.5));
               default: return 0.;
          }
     }
     else {return m*paLegn(t,n,m)/sin(t);}
}

double LTaun(double t, int n, int m) {
     if (m < 0) return (abs(m)%2) ? -LTaun(t,n,abs(m)) : LTaun(t,n,abs(m)); 
     if (fabs(t) < 1.e-14) {
          switch (m) {
               case 1: return 0.5*sqrt(n*(n+1)*(n+0.5));
               case -1: return -0.5*sqrt(n*(n+1)*(n+0.5));
               default: return 0.;
          }
     }
     else if (fabs(t-M_PI) < 1.e-14) {
          switch (m) {
               case 1: return (n%2) ? -0.5*sqrt(n*(n+1)*(n+0.5)) : 0.5*sqrt(n*(n+1)*(n+0.5));
               case -1: return (n%2) ? 0.5*sqrt(n*(n+1)*(n+0.5)) : -0.5*sqrt(n*(n+1)*(n+0.5));
               default: return 0.;
          }
     }
     else {return -sin(t)*paLegnd(t,n,m);}
}

     // spherical vector functions

Vector svfRgM(Complex z, double th, double ph, int n, int m) {
     Vector M(3);
     M.Data[0] = 0.;
     M.Data[1] = j_*LPin(th,n,m)*( M.Data[2] = besj(z,n)*exp(j_*double(m)*ph)*M_SQRT1_2PI/sqrt(n*(n+1.)) );//double(m)*
     M.Data[2] *= -LTaun(th,n,m);
     return M;
}

Vector svfM1(Complex z, double th, double ph, int n, int m) {
     Vector M(3);
     M.Data[0] = 0.;
     M.Data[1] = j_*LPin(th,n,m)*( M.Data[2] = besh1(z,n)*exp(j_*double(m)*ph)*M_SQRT1_2PI/sqrt(n*(n+1.)) );//double(m)*
     M.Data[2] *= -LTaun(th,n,m);
     return M;
}

Vector svfM2(Complex z, double th, double ph, int n, int m) {
     Vector M(3);
     M.Data[0] = 0.;
     M.Data[1] = j_*LPin(th,n,m)*( M.Data[2] = besh2(z,n)*exp(j_*double(m)*ph)*M_SQRT1_2PI/sqrt(n*(n+1.)) );//double(m)*
     M.Data[2] *= -LTaun(th,n,m);
     return M;
}

Vector svfRgMx(Complex z, double th, double ph, int n, int m) {
     double tv = 0.5/sqrt(n*(n+1.)); Complex tc; Vector M(3);
     M.Data[0] = j_*( M.Data[1] = tv*sqrt((n-m)*(n+m+1))*sfRg(z,th,ph,n,m+1) );
     tc = tv*sqrt((n+m)*(n-m+1))*sfRg(z,th,ph,n,m-1);
     M.Data[0] += j_*tc; M.Data[1] -= tc;
     M.Data[2] = -2.*j_*tv*double(m)*sfRg(z,th,ph,n,m);
     return M;
}

Vector svfM1x(Complex z, double th, double ph, int n, int m) {
     double tv = 1./sqrt(n*(n+1.)); Complex tc; Vector M(3);
     M.Data[0] = j_*( M.Data[1] = 0.5*tv*sqrt((n-m)*(n+m+1))*sf1(z,th,ph,n,m+1) );
     tc = 0.5*tv*sqrt((n+m)*(n-m+1))*sf1(z,th,ph,n,m-1);
     M.Data[0] += j_*tc; M.Data[1] -= tc;
     M.Data[2] = -j_*tv*double(m)*sf1(z,th,ph,n,m);
     return M;
}

Vector svfM2x(Complex z, double th, double ph, int n, int m) {
     double tv = 1./sqrt(n*(n+1.)); Complex tc; Vector M(3);
     M.Data[0] = j_*( M.Data[1] = 0.5*tv*sqrt((n-m)*(n+m+1))*sf2(z,th,ph,n,m+1) );
     tc = 0.5*tv*sqrt((n+m)*(n-m+1))*sf2(z,th,ph,n,m-1);
     M.Data[0] += j_*tc; M.Data[1] -= tc;
     M.Data[2] = -j_*tv*double(m)*sf2(z,th,ph,n,m);
     return M;
}

Vector svfRgN(Complex z, double th, double ph, int n, int m) {
     double tv = sqrt(n*(n+1.)); Complex tc, tz; Vector N(3);
     tc = exp(j_*double(m)*ph); tz = besj(z,n)/z;
     N.Data[0] = tv*M_SQRT1_2PI*tz*paLegn(th,n,m)*tc;
     N.Data[2] = j_*LPin(th,n,m)*( N.Data[1] = M_SQRT1_2PI/tv*tc*(tz + besjd(z,n)) );//double(m)*
     N.Data[1] *= LTaun(th,n,m);
     return N;
}

Vector svfN1(Complex z, double th, double ph, int n, int m) {
     double tv = sqrt(n*(n+1.)); Complex tc, tz; Vector N(3);
     tc = exp(j_*double(m)*ph); tz = besh1(z,n)/z;
     N.Data[0] = tv*M_SQRT1_2PI*tz*paLegn(th,n,m)*tc;
     N.Data[2] = j_*LPin(th,n,m)*( N.Data[1] = M_SQRT1_2PI/tv*tc*(tz + besh1d(z,n)) );//double(m)*
     N.Data[1] *= LTaun(th,n,m);
     return N;
}

Vector svfN2(Complex z, double th, double ph, int n, int m) {
     double tv = sqrt(n*(n+1.)); Complex tc, tz; Vector N(3);
     tc = exp(j_*double(m)*ph); tz = besh2(z,n)/z;
     N.Data[0] = tv*M_SQRT1_2PI*tz*paLegn(th,n,m)*tc;
     N.Data[2] = j_*LPin(th,n,m)*( N.Data[1] = M_SQRT1_2PI/tv*tc*(tz + besh2d(z,n)) );//double(m)*
     N.Data[1] *= LTaun(th,n,m);
     return N;
}

Vector svfRgNx(Complex z, double th, double ph, int n, int m) {
     double tv; Complex tc; Vector N(3);
     tv = 1./sqrt(n*(n+1.)*(2*n+1.));
     N.Data[0] = N.Data[1] = n*sqrt((n+m+1.)*(n+m+2.)/(2*n+3.))*sfRg(z,th,ph,n+1,m+1)
          - (n+1)*sqrt((n-m)*(n-m-1.)/(2*n-1.))*sfRg(z,th,ph,n-1,m+1);
     tc = -n*sqrt((n-m+1.)*(n-m+2.)/(2*n+3.))*sfRg(z,th,ph,n+1,m-1)
          + (n+1)*sqrt((n+m)*(n+m-1.)/(2*n-1.))*sfRg(z,th,ph,n-1,m-1);
     N.Data[0] += tc; N.Data[1] -= tc; N.Data[0] *= 0.5*tv; N.Data[1] *= -0.5*j_*tv;
     N.Data[2] = tv*( n*sqrt((n-m+1.)*(n+m+1.)/(2*n+3.))*sfRg(z,th,ph,n+1,m)
          + (n+1)*sqrt((n-m)*(n+m)/(2*n-1.))*sfRg(z,th,ph,n-1,m) );
     return N;
}

Vector svfN1x(Complex z, double th, double ph, int n, int m) {
     double tv; Complex tc; Vector N(3);
     tv = 1./sqrt(n*(n+1.)*(2*n+1.));
     N.Data[0] = N.Data[1] = n*sqrt((n+m+1.)*(n+m+2.)/(2*n+3.))*sf1(z,th,ph,n+1,m+1)
          - (n+1)*sqrt((n-m)*(n-m-1.)/(2*n-1.))*sf1(z,th,ph,n-1,m+1);
     tc = -n*sqrt((n-m+1.)*(n-m+2.)/(2*n+3.))*sf1(z,th,ph,n+1,m-1)
          + (n+1)*sqrt((n+m)*(n+m-1.)/(2*n-1.))*sf1(z,th,ph,n-1,m-1);
     N.Data[0] += tc; N.Data[1] -= tc; N.Data[0] *= 0.5*tv; N.Data[1] *= -0.5*j_*tv;
     N.Data[2] = tv*( n*sqrt((n-m+1.)*(n+m+1.)/(2*n+3.))*sf1(z,th,ph,n+1,m)
          + (n+1)*sqrt((n-m)*(n+m)/(2*n-1.))*sf1(z,th,ph,n-1,m) );
     return N;
}

Vector svfN2x(Complex z, double th, double ph, int n, int m) {
     double tv; Complex tc; Vector N(3);
     tv = 1./sqrt(n*(n+1.)*(2*n+1.));
     N.Data[0] = N.Data[1] = n*sqrt((n+m+1.)*(n+m+2.)/(2*n+3.))*sf2(z,th,ph,n+1,m+1)
          - (n+1)*sqrt((n-m)*(n-m-1.)/(2*n-1.))*sf2(z,th,ph,n-1,m+1);
     tc = -n*sqrt((n-m+1.)*(n-m+2.)/(2*n+3.))*sf2(z,th,ph,n+1,m-1)
          + (n+1)*sqrt((n+m)*(n+m-1.)/(2*n-1.))*sf2(z,th,ph,n-1,m-1);
     N.Data[0] += tc; N.Data[1] -= tc; N.Data[0] *= 0.5*tv; N.Data[1] *= -0.5*j_*tv;
     N.Data[2] = tv*( n*sqrt((n-m+1.)*(n+m+1.)/(2*n+3.))*sf2(z,th,ph,n+1,m)
          + (n+1)*sqrt((n-m)*(n+m)/(2*n-1.))*sf2(z,th,ph,n-1,m) );
     return N;
}
