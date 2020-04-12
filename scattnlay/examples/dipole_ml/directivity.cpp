
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

#include "./matrix.h"
#include "./sphereml.h"

#include <math.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <memory.h>
#include <vector>

Vector evaluate_harmonics(const std::vector<double> &RL,
                          const std::vector< std::complex<double> > &eL_in,
                          const double &Rd, const double &wl,
                          const double &px, const double &py, const double &pz,
                          const int N = 41) {
    int  NL, il, n, m, nm;
    double /*dx, dz,*/ wv/*, kRs, tm*/;
    Complex es, ris, kRd;
    Vector VD1(2*N*N), VD2(2*N*N), VS2(2*N*N), E(3), ED(3), SM(3), SN(3), RC(2*N), VF(2);
    SphereML MS(N);

    wv = 2.*M_PI/wl;

    memset(VS2.Data,0,2*N*N*sizeof(Complex));

        // initialize spherical multilayer:
    // we can place some variables on the stack in order to simplify memory management
    NL = RL.size(); 
    double kRL[NL];
    Complex eL[NL+1];
    Matrix M1(4,2*N), M2(4,2*N), **M;

    for (int i=0; i<NL+1; ++i) eL[i] = eL_in[i];
    for (int i=0; i<NL; ++i) {kRL[i] = wv*RL[i]; eL[i] *= eL[i];}
        // scattering matrices of all spherical interfaces:
    M = new Matrix* [NL];
    for (int i=0; i<NL; ++i) {
        M[i] = new Matrix(4,2*N);
        *M[i] = MS.calc_RT(kRL[i],eL[i],eL[i+1],1.,1.);
    }
    M2 = MS.calc_SML(M,NL); // initial scattering matrix
    il = 0; // initial dipole position (inside the smallest sphere)
    memset(M1.Data,0,8*N*sizeof(Complex));
    
    if (Rd < RL[0]) { // loop for dipole positions inside the smallest sphere
        kRd = wv*Rd*sqrt(eL[0]); VD2 = MS.calc_edz(px,py,pz,kRd,0);
        for (n=1; n<N; n++) for (m=-1; m<2; m+=2) {
            nm = n*(n+1)+m; VS2.Data[nm] = VD2(nm)*M2(1,n); VS2.Data[nm+N*N] = VD2(nm+N*N)*M2(1,n+N);
        }
    } else if (Rd < RL[NL-1]) { // dipole inside multilayer
        while (Rd > RL[il]) il++;
        M1 = MS.calc_SML(M,il); M2 = MS.calc_SML(M+il,NL-il);
        kRd = wv*Rd*sqrt(eL[il]); VD1 = MS.calc_edz(px,py,pz,kRd,1); VD2 = MS.calc_edz(px,py,pz,kRd,0);
        for (n=1; n<N; n++) for (m=-1; m<2; m+=2) {
            nm = n*(n+1)+m;
            VS2.Data[nm] = (VD1(nm)*M1(3,n) + VD2(nm))*M2(1,n)/(1.-M1(3,n)*M2(0,n));
            VS2.Data[nm+N*N] = (VD1(nm+N*N)*M1(3,n+N) + VD2(nm+N*N))*M2(1,n+N)/(1.-M1(3,n+N)*M2(0,n+N));
        }
    } else {         // dipole outside the mutilayer
        M1 = MS.calc_SML(M,NL);
        kRd = wv*Rd*sqrt(eL[NL]);
        VD1 = MS.calc_edz(px,py,pz,kRd,1);
        VD2 = MS.calc_edz(px,py,pz,kRd,0);
        for (n=1; n<N; n++)
            for (m=-1; m<2; m+=2) {
                nm = n*(n+1)+m;
                VS2.Data[nm] = VD1(nm)*M1(3,n) + VD2(nm);
                VS2.Data[nm+N*N] = VD1(nm+N*N)*M1(3,n+N) + VD2(nm+N*N);
            }
    }

    for (int i=0; i<NL; ++i) delete M[i];
    delete[] M; 
    
    return VS2;
}

double evaluate_directivity(const std::vector<double> &RL,
                            const std::vector< std::complex<double> > &eL,
                            const double &Rd, const double &wl,
                            const double &px, const double &py, const double &pz,
                            const double th=M_PI*0., // angle for directivity evaluation
                            const double ph=0.,
                            const int N = 41) {
    const Vector& VS2 = evaluate_harmonics(RL, eL, Rd, wl, px, py, pz, N);
    SphereML MS(N);
    return MS.directivity(VS2,th,ph,1.);
}
