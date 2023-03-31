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

    auto t = atan(72.0/21)/2;
    cout << "θ = " << t << endl;

    M v;
    v.col(1) << cos(t), sin(t);
    v.col(0) << -sin(t), cos(t);

    cout << "V =" << endl << v << endl << endl;

    M u = m*v;

    M s;
    s << u.col(0).norm(), 0, 0, u.col(1).norm();
    cout << "Σ =" << endl << s << endl << endl;

    u.col(0) = u.col(0).normalized();
    u.col(1) = u.col(1).normalized();
    cout << "U =" << endl << u << endl << endl;

    cout << "M-U*Σ*V' =" << endl << m-u*s*v.adjoint() << endl << endl;

    MatrixXd m2 = m;
    BDCSVD<M> svd(m2, ComputeFullU | ComputeFullV);

    cout << "bdcU" << endl << svd.matrixU() << endl << endl;
    MatrixXd bdcS = svd.singularValues().asDiagonal();
    cout << "bdcS" << endl << bdcS << endl;
    cout << "bdcV" << endl << svd.matrixV() << endl << endl;
    cout << "M-bdcU*s*bdcV'" << endl << m-svd.matrixU()*bdcS*svd.matrixV().adjoint() << endl << endl;

    cout << m.eigenvalues() << endl << endl;
}
