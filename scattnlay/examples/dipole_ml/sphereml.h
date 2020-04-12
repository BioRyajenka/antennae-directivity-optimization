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

#ifndef _SPHEREML_H
#define _SPHEREML_H

#include "matrix.h"

#define _USE_MATH_DEFINES
#include <math.h>

class SphereML {
public:
     int N;

     SphereML(int N_) {N = N_;}

     Vector calc_pw(double as, double ap, double th, double ph);
     Vector calc_edz(double px, double py, double pz, Complex krz, int in);

     Vector calc_far(const Vector &V, double th, double ph);
     double calc_Psca(const Vector &VS, double tC);
     double calc_Pext(const Vector &VI, const Vector &VS, double tC);
     double directivity(const Vector &VS, double th, double ph, double tC);

     Matrix calc_RT(double kr, Complex e1, Complex e2, Complex m1, Complex m2);
     Matrix calc_SML(Matrix **SM, int ns);
};

#endif
