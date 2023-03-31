#include <iostream>
#include <numbers>
#include <cmath>
#define EIGEN_USE_BLAS
#include <Eigen/Dense>
#include <complex>
#include "gnuplot-iostream.h"
#include <algorithm>

using namespace Eigen;
using namespace std;

constexpr double pi = std::numbers::pi;
constexpr double tau = 2*pi;

Vector2d grad(Matrix2d m, Vector2d v, double factor) {
    return factor*m.transpose()*m*v;
}

int main() {
    srand(time(0));
    Matrix2d m = Matrix2d::Random();
    cout << m << endl << endl;
    cout << "Determinant: " << m.determinant() << endl;
    auto hi = make_pair(pi,(m*Vector2d(-1,0)).norm());
    auto low = make_pair(0.0,(m*Vector2d(1,0)).norm());
    
    int n = 100;
    ArrayXd angles = ArrayXd::LinSpaced(n, 0, tau);
    MatrixXd points(2,n);
    points.row(0) = cos(angles);
    points.row(1) = sin(angles);
    MatrixXd target = m*points;

    auto diff = numeric_limits<double>::infinity();
    auto eps = 2*numeric_limits<double>::epsilon();
    Vector2d v(1,0);
    int i = 0;
    int factor = 2;
    while (diff > eps) {
        i++;
        Vector2d next = (v+grad(m,v,factor)).normalized();
        diff = (next-v).norm();
        v = next;
        // factor = (1+factor)/2;
    }

    cout << i << " iterations\n";

    Vector2d u = m*v;
    Matrix2d rot;
    rot << 0, -1, 1, 0;
    Vector2d v2 = rot*v;
    Vector2d u2 = m*v2;

    Gnuplot g;
    g << "set style line 1 lc rgb \"blue\"\n"
        << "set style line 2 lc rgb \"red\"\n"
        << "set style arrow 2 ls 1\n"
        << "set style arrow 3 ls 2\n"
        << "set arrow from 0,0 to " << v(0) << "," << v(1) << "\n"
        << "set arrow from 0,0 to " << u(0) << "," << u(1) << "as 2\n"
        << "set arrow from 0,0 to " << v2(0) << "," << v2(1) << "\n"
        << "set arrow from 0,0 to " << u2(0) << "," << u2(1) << "as 3\n"
        << "set param\n"
        << "plot" << g.file1d(points.transpose()) << "with lines,"
        << g.file1d(target.transpose()) << "with lines\n";
}
