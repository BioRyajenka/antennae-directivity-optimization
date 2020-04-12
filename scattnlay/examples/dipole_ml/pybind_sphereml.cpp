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

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>


namespace py = pybind11;


std::vector<double> Py2VectorDouble(const py::array_t<double> &py_x) {
  std::vector<double> c_x(py_x.size());
  std::memcpy(c_x.data(), py_x.data(), py_x.size()*sizeof(double));
  return c_x;
}


std::vector< std::complex<double> > Py2VectorComplex(const py::array_t< std::complex<double> > &py_x){
  std::vector< std::complex<double> > c_x(py_x.size());
  std::memcpy(c_x.data(), py_x.data(), py_x.size()*sizeof( std::complex<double>));
  return c_x;
}

py::array_t< std::complex<double>, py::array::c_style | py::array::forcecast > VectorComplex2Py(const Vector& v) {
  return py::array(v.Nrow, v.Data);
}

py::array_t< std::complex<double>, py::array::c_style | py::array::forcecast > py_evaluate_harmonics(const py::array_t<double, py::array::c_style | py::array::forcecast> &RL,
                             const py::array_t< std::complex<double>, py::array::c_style | py::array::forcecast> &eL,
                             const double Rd, const double wl,
                             const double px, const double py, const double pz,
                             const int N) {
    const auto& c_RL = Py2VectorDouble(RL);
    const auto& c_eL = Py2VectorComplex(eL);
    const Vector& res = evaluate_harmonics(c_RL, c_eL, Rd, wl, px, py, pz, N);
    return VectorComplex2Py(res);
}

double py_evaluate_directivity(const py::array_t<double, py::array::c_style | py::array::forcecast> &RL,
                               const py::array_t< std::complex<double>, py::array::c_style | py::array::forcecast> &eL,
                               const double Rd, const double wl,
                               const double px, const double py, const double pz,
                               const double th, const double ph,
                               const int N) {
    const auto& c_RL = Py2VectorDouble(RL);
    const auto& c_eL = Py2VectorComplex(eL);
    return evaluate_directivity(c_RL, c_eL, Rd, wl, px, py, pz, th, ph, N);
}


PYBIND11_MODULE(sphereml, m) {
    m.doc() = "sphereml evaluates excitation of a multilayerd sphere by a dipole source"; // optional module docstring

    m.def("evaluate_directivity", &py_evaluate_directivity, "evaluate directivity",
          py::arg("RL"), py::arg("eL"),
          py::arg("Rd"), py::arg("wl"),
          py::arg("px")=1., py::arg("py")=0., py::arg("pz")=0.,
          py::arg("th")=0., py::arg("ph")=0.,
          py::arg("N")=41);

    m.def("evaluate_harmonics", &py_evaluate_harmonics, "evaluate harmonics",
          py::arg("RL"), py::arg("eL"),
          py::arg("Rd"), py::arg("wl"),
          py::arg("px")=1., py::arg("py")=0., py::arg("pz")=0.,
          py::arg("N")=41);
}

