#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

constexpr auto eps = numeric_limits<double>::epsilon();

template<int m, int n>
using M = Matrix<double, m, n>;
template<int m>
using V = Vector<double, m>;
template<int m>
using A = Array<double, m, 1>;

template<int m, int n>
pair<V<m>, V<n>> classic(const V<m>& aj, const M<m,n>& q, int j) {
    V<m> qj = aj;
    V<n> rj = V<n>::Zero();
    for (int i = 0; i < j; i++) {
        rj(i) = q.col(i).dot(aj);
        qj -= rj(i)*q.col(i);
    }
    rj(j) = qj.norm();

    return make_pair(qj.normalized(),rj);
}

template<int m, int n>
pair<V<m>, V<n>> modified(V<m> qj, const M<m,n>& q, int j) {
    V<n> rj = V<n>::Zero();
    for (int i = 0; i < j; i++) {
        rj(i) = q.col(i).dot(qj);
        qj -= rj(i)*q.col(i);
    }
    rj(j) = qj.norm();

    return make_pair(qj.normalized(),rj);
}

template<int m, int n, typename Ortho>
pair<M<m,n>, M<n,n>> QR(const M<m,n>& a, const Ortho& ortho) {
    M<m,n> q = M<m,n>::Zero();
    M<n,n> r = M<n,n>::Zero();

    for (int j = 0; j < n; j++) {
        auto qrj = ortho(a.col(j), q, j);
        r.col(j) = qrj.second;
        // for (int i = 0; qrj.second(j) <= 2*eps && i < n; i++) {
        //     r(j,j) = 0;
        //     qrj = ortho(Vector<double, m>::Unit(i), q, j);
        // }
        q.col(j) = qrj.first;
    }

    return make_pair(q,r);
}

template<int m, int n>
pair<M<m,n>, M<n,n>> QRm(const M<m,n>& a) {
    return QR(a, modified<m,n>);
}

template<int m, int n>
pair<M<m,n>, M<n,n>> QRc(const M<m,n>& a) {
    return QR(a, classic<m,n>);
}

// This is an equivalent implementation to QRm, but the order of projection is different.
// The result is numerically identical.
template<int m, int n>
pair<M<m,n>, M<n,n>> QRmod(const M<m,n>& a) {
    M<m,n> q = a;
    M<n,n> r = M<n,n>::Zero();
    for (int i = 0; i < n; i++) {
        r(i,i) = q.col(i).norm();
        q.col(i).normalize();
        auto qi = q.col(i);
        for (int j = i+1; j < n; j++) {
            auto rij = qi.dot(q.col(j));
            q.col(j) -= rij*qi;
            r(i,j) = rij;
        }
    }
    return make_pair(q,r);
}
