#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

using Vx = VectorXd;
using Mx = MatrixXd;
template<int m, int n>
using M = Matrix<double, m, n>;
template<int m>
using V = M<m,1>;
template<int m>
using A = Array<double, m, 1>;

template<typename T>
int sign(T val) {
    return (T(0) < val) - (T(0) > val);
}

template<int m, int n>
pair<M<m,m>, M<m,n>> QRh(M<m,n> r) {
    M<m,m> q = M<m,m>::Identity();

    for (int i = 0; i < n; i++) {
        VectorXd ri = r(lastN(m-i), i);
        ri(0) += sign(ri(0))*ri.norm();
        ri.normalize();
        r.bottomRightCorner(m-i,n-i) -= 2*ri*(ri.adjoint()*r.bottomRightCorner(m-i,n-i));
        q.bottomRightCorner(m-i,m) -= 2*ri*(ri.adjoint()*q.bottomRightCorner(m-i,m));
    }

    return make_pair(q.adjoint(),r);
}

template<int m, int n>
pair<M<m,n>, M<n,n>> thinQRh(const M<m,n>& a) {
    auto [q,r] = QRh(a);
    return make_pair(q*M<m,n>::Identity(), M<n,m>::Identity()*r);
}

pair<Mx, Mx> QRh(Mx r) {
    int m = r.rows(), n = r.cols();
    Mx q = Mx::Identity(m,m);
    for (int i = 0; i < n; i++) {
        VectorXd ri = r(lastN(m-i), i);
        ri(0) += sign(ri(0))*ri.norm();
        ri.normalize();
        r.bottomRightCorner(m-i,n-i) -= 2*ri*(ri.adjoint()*r.bottomRightCorner(m-i,n-i));
        q.bottomRightCorner(m-i,m) -= 2*ri*(ri.adjoint()*q.bottomRightCorner(m-i,m));
    }

    return make_pair(q.adjoint(), r);
}

pair<Mx, Mx> thinQRh(const Mx& a) {
    int m = a.rows(), n = a.cols();
    auto [q,r] = QRh(a);
    return make_pair(q*Mx::Identity(m,n), Mx::Identity(n,m)*r);
}

