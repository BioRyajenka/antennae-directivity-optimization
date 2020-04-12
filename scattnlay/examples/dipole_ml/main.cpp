
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
#include "./directivity.h"

#include <math.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <memory.h>
#include <vector>
//#include <Windows.h>

#define USE_MATH_DEFINES



int main(int argc, char* argv[]) {
    std::cout.precision(6);
    double px = 1., py = 0., pz = 0.; // dipole amplitudes
    double wl = 0.455;
    int NL = 3;
    
    // std::vector<double> RL = {0.37469538, 0.86627232,  0.91};
    // std::vector< std::complex<double> > eL = {3.07589484, 3.95118756, 2.59045949, 1.};
    // double dRd = 0.4591807*0.1;
    std::vector<double> RL(NL);
    std::vector<Complex> eL(NL+1);
    
    for (int i=0; i<NL; ++i) {eL[i] = 2.*(2.+i) + 0.*j_;} eL[NL] = 1.;
    RL[0] = 0.09; for (int i=1; i<NL; ++i) {RL[i] = RL[i-1] + 0.02;}
    double dRd = 0.001*M_SQRT2;
    for (int i = 0; i < 200; ++i) {
        double Rd = dRd*i;     // dipole position
        std::cout<<Rd<<" "<<evaluate_directivity(RL, eL, Rd, wl, px, py, pz)<<std::endl;
    }
    return 0;
}

