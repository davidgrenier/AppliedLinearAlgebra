#include <iostream>
#include <numbers>
#include <vector>
#include <cmath>
#define EIGEN_USE_BLAS
#include <Eigen/Dense>
#include <complex>
// #include "gnuplot-iostream.h"
#include <algorithm>

using namespace Eigen;
using namespace std;
typedef std::vector<double> vec;

constexpr double pi = std::numbers::pi;
constexpr double tau = 2*pi;

int main() {
    auto u = Vector4d::Random();
    auto v = Vector4d::Random();
    auto E = u*v.adjoint();

    for (int i = 0; i < 100; i++) {
        Matrix4d F = E;
        cout << abs(pow((F*F.adjoint()).trace(),0.5) - F.operatorNorm()) << endl;
        cout << abs(F.norm() - F.operatorNorm()) << endl;
        cout << abs(pow((F*F.adjoint()).trace(),0.5) - F.norm()) << endl << endl;
    }
}
