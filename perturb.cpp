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

int main() {
    Vector4d u(1,1,2,3);
    Vector4d v(1,2,3,1);
    Matrix4d A = Matrix4d::Identity() + u*v.transpose();
    cout << A.eigenvalues() << endl;
    cout << A.jacobiSvd().singularValues() << endl;

    u = Vector4d(1,0,0,0);
    v = Vector4d(-1,1,0,0);
    cout << "Should be -1: " << u.transpose()*v << endl;
    A = Matrix4d::Identity() + u*v.transpose();

    cout << "\nInverse\n" << A.inverse() << endl;
    cout << "\nPartialPivLU\n" << A.partialPivLu().solve(Matrix4d::Identity()) << endl;
    cout << "\nHouseholderQr\n" << A.householderQr().solve(Matrix4d::Identity()) << endl;
    cout << "\ncolvPivHouseholderQr\n" << A.colPivHouseholderQr().solve(Matrix4d::Identity()) << endl;
    cout << "\nfullHouseHolderQr\n" << A.fullPivHouseholderQr().solve(Matrix4d::Identity()) << endl;
    cout << "\nOrtho\n" << A.completeOrthogonalDecomposition().solve(Matrix4d::Identity()) << endl;
    cout << "\nPseudoinverse\n" << (A.transpose()*A).fullPivLu().inverse()*A.transpose() << endl;
    cout << endl << A.eigenvalues() << endl;
}
