#include <iostream>
#include <numbers>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
using std::numbers::pi;

typedef Matrix2d M;
typedef Vector2d V;

int main() {
    M m;
    m << -2, 11, -10, 5;
    cout << "M =" << endl << m << endl << endl;

    BDCSVD<M> svd(m, ComputeFullU | ComputeFullV);
    Matrix2d u = svd.matrixU();
    Matrix2d v = svd.matrixV();
    cout << svd.matrixU() << endl << endl;
    cout << svd.matrixV() << endl << endl;

    Matrix4d c;
    c.col(0) << v.col(0), u.col(0);
    c.col(1) << v.col(1), u.col(1);
    c.col(2) << v.col(0), -u.col(0);
    c.col(3) << v.col(1), -u.col(1);
    c /= sqrt(2);

    cout << "Orthonormal eigenbasis from svd:" << endl << c << endl << endl;

    Matrix4d b;
    b << M::Zero(), m.adjoint(), m, M::Zero();

    cout << "B =" << endl << b << endl << endl;

    SelfAdjointEigenSolver<Matrix4d> solv(b);
    cout << solv.eigenvectors() << endl << endl;
    cout << solv.eigenvalues();
}
