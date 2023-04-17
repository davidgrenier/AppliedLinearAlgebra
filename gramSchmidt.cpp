#include <iostream>
#include <numbers>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
using std::numbers::pi;

constexpr auto eps = numeric_limits<double>::epsilon();

template<int m, int n>
using M = Matrix<double, m, n>;
template<int m>
using V = Vector<double, m>;

template<int m, int n>
pair<V<m>, V<n>> classic(const V<m>& aj, const M<m,n>& q, int j) {
    V<m> qj = aj;
    V<n> rj = V<n>::Zero();
    for (int i = 0; i < j; i++) {
        rj(i) = q.col(i).dot(aj);
        qj -= rj(i)*q.col(i);
    }
    rj(j) = qj.norm();
    qj.normalize();

    return make_pair(qj,rj);
}

template<int m, int n>
pair<V<m>, V<n>> modified(V<m> qj, const M<m,n>& q, int j) {
    V<n> rj = V<n>::Zero();
    for (int i = 0; i < j; i++) {
        rj(i) = q.col(i).dot(qj);
        qj -= rj(i)*q.col(i);
    }
    rj(j) = qj.norm();
    qj.normalize();

    return make_pair(qj,rj);
}

template<int m, int n, typename Ortho>
pair<M<m,n>, M<n,n>> QR(const M<m,n>& a, const Ortho& o) {
    M<m,n> q = M<m,n>::Zero();
    M<n,n> r = M<n,n>::Zero();

    for (int j = 0; j < n; j++) {
        auto qrj = o(a.col(j), q, j);
        r.col(j) = qrj.second;
        for (int i = 0; qrj.second(j) <= 2*eps && i < n; i++) {
            r(j,j) = 0;
            qrj = o(Vector<double, m>::Unit(i), q, j);
        }
        q.col(j) = qrj.first;
    }

    return make_pair(q,r);
}

template<int m, int n>
void showQR(const M<m,n> a, const pair<M<m,n>, M<n,n>>& qr) {
    auto q = qr.first;
    auto r = qr.second;

    cout << "Q =\n" << q << "\n\n";
    cout << "R =\n" << r << "\n\n";
    cout << "Q*Q =\n" << q.adjoint()*q << "\n\n";
    cout << "A-QR =\n" << a-q*r << "\n\n";
    cout << "FrobN of diff = " << (a-q*r).norm() << "\n\n";
}

int main() {
    srand(time(0));
    constexpr int m = 5, n = 3;

    M<m,n> a = M<m,n>::Random();
    cout << "A =\n" << a << "\n\n";

    showQR(a, QR(a, classic<m,n>));
    showQR(a, QR(a, modified<m,n>));
}
