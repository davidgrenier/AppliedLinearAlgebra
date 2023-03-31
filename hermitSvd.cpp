#include <iostream>
// #define EIGEN_USE_BLAS
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main() {
    typedef MatrixXd M;
    M h = Matrix4d::Random(4,4).selfadjointView<Lower>();

    cout << h << endl << endl;

    BDCSVD<M> d(h, ComputeThinU | ComputeThinV);
    M u = d.matrixU();
    M v = d.matrixV();
    M s = d.singularValues().asDiagonal();

    cout << u << endl << endl;
    cout << v << endl << endl;
    cout << h - u*s*v.adjoint() << endl << endl;

    SelfAdjointEigenSolver<M> solver(h);
    cout << solver.eigenvalues() << endl;
    cout << solver.eigenvectors() << endl;
}
