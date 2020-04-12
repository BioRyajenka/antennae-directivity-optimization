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

#ifndef _DIRECTIVITY_H
#define _DIRECTIVITY_H

#include "./matrix.h"
#include "./sphereml.h"

#include <cmath>
#include <complex>
#include <fstream>
#include <math.h>
#include <memory.h>
#include <vector>


#define _USE_MATH_DEFINES
Vector evaluate_harmonics(const std::vector<double> &RL,
                          const std::vector< std::complex<double> > &eL_in,
                          const double &Rd, const double &wl,
                          const double &px, const double &py, const double &pz,
                          const int N = 41);

double evaluate_directivity(const std::vector<double> &RL_in,
                            const std::vector< std::complex<double> > &eL_in,
                            const double &Rd, const double &wl,
                            const double &px, const double &py, const double &pz,
                            const double th=M_PI*0., // angle for directivity evaluation
                            const double ph=0.,
                            const int N = 41);
#endif