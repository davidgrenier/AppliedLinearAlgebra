#include <iostream>
#include <numbers>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
using std::numbers::pi;

typedef Matrix<double, 3, 2> M;
typedef Vector3d V;

M getQ(const M& m) {
    HouseholderQR<M> qr(m);
    return qr.householderQ()*M::Identity();
}

int main() {
    M u;
    u << 1, 0, 0, 1, 1, 0;
    M qu = getQ(u);

    M v;
    v << 1, 2, 0, 1, 1, 0;
    M qv = getQ(v);
    cout << "U =\n" << u << "\n\n";
    cout << "V =\n" << v << "\n\n";
    cout << qu << "\n\n";
    cout << qv << "\n\n";

    Matrix3d m = qu*qu.adjoint()*qv*qv.adjoint();
    HouseholderQR<Matrix3d> qrm(m);
    Matrix3d qm = qrm.householderQ();
    cout << qm << "\n\n";

    Vector3d mu = m.col(0);
    cout << "mu =\n" << mu << "\n\n";
    cout << "m*mu =\n" << m*mu << "\n\n";
}
