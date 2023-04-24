#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

template<int m, int n>
using M = Matrix<double, m, n>;
template<int m>
using V = Vector<double, m>;
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
        M<m,m> qi = M<m,m>::Identity();
        r.bottomRightCorner(m-i,n-i) -= 2*ri*(ri.adjoint()*r.bottomRightCorner(m-i,n-i));
        q.bottomRightCorner(m-i,m) -= 2*ri*(ri.adjoint()*q.bottomRightCorner(m-i,m));
    }

    return make_pair(q.adjoint(),r);
}
