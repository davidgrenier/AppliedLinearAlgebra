#include <iostream>
#include <numbers>
#include <vector>
#include <cmath>
#define EIGEN_USE_BLAS
#include <Eigen/Dense>
#include <complex>
#include "gnuplot-iostream.h"
#include <algorithm>

using namespace Eigen;
using namespace std;
typedef std::vector<double> vec;

constexpr double pi = std::numbers::pi;
constexpr double tau = 2*pi;

ArrayXd norm(ArrayXd xs, ArrayXd ys, double p) {
    return (xs.abs().pow(p)+ys.abs().pow(p)).pow(1/p);
}

int main() {
    Gnuplot g;

    ArrayXd t = ArrayXd::LinSpaced(500,0,tau);
    ArrayXd xs = cos(t);
    ArrayXd ys = sin(t);
    auto l1 = norm(xs,ys,1);
    ArrayXd xs1 = xs/l1;
    ArrayXd ys1 = ys/l1;
    auto l3 = norm(xs,ys,3);
    ArrayXd xs3 = xs/l3;
    ArrayXd ys3 = ys/l3;
    auto l10 = norm(xs,ys,10);
    ArrayXd xs10 = xs/l10;
    ArrayXd ys10 = ys/l10;

    g << "plot [x=-1.1:1.1] [y=-1.1:1.1]" << g.file1d(make_pair(xs,ys)) << "with lines,"
        << g.file1d(make_pair(xs1, ys1)) << "with lines,"
        << g.file1d(make_pair(xs3, ys3)) << "with lines,"
        << g.file1d(make_pair(xs10, ys10)) << "with lines\n";

    Matrix2d A;
    A << 1,2
        ,0,2;

    MatrixXd coords(xs.size(),2);
    coords << xs, ys;
    g << "plot [x=-3:3] [y=-3:3]" << g.file1d(coords*A.transpose()) << "with lines,";
    coords << xs1,ys1;
    g << g.file1d(coords*A.transpose()) << "with lines,";
    coords << xs3,ys3;
    g << g.file1d(coords*A.transpose()) << "with lines,";
    coords << xs10,ys10;
    g << g.file1d(coords*A.transpose()) << "with lines\n";
}
