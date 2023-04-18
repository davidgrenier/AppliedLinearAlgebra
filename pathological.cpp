#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include "gnuplot-iostream.h"

using namespace Eigen;
using namespace std;
using std::numbers::pi;

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
pair<M<m,n>, M<n,n>> QRm(const M<m,n>& a) {
    return QR(a, modified<m,n>);
}

template<int m, int n>
pair<M<m,n>, M<n,n>> QRc(const M<m,n>& a) {
    return QR(a, classic<m,n>);
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

template<int m>
M<m,m> randomOrtho() {
    M<m,m> rand1 = M<m,m>::Random();
    HouseholderQR<M<m,m>> qr1(rand1);
    return qr1.householderQ();
}

int main() {
    srand(time(0));
    constexpr int m = 80;

    auto u = randomOrtho<m>();
    auto v = randomOrtho<m>();
    A<80> ix = A<80>::LinSpaced(1,80);
    A<80> d = (A<80>::Ones()*2).pow(-ix);
    M<m,m> s = d.matrix().asDiagonal();
    
    M<m,m> a = u*s*v.adjoint();
    auto qrc = QRc(a);
    auto qrm = QRm(a);

    Gnuplot g;

    g << "plot"
        << g.file1d(make_pair(ix, log(A<80>::Ones()*eps)/log(2))) << "with lines,"
        << g.file1d(make_pair(ix, log(A<80>::Ones()*sqrt(eps))/log(2))) << "with lines,"
        << g.file1d(make_pair(ix, log(qrc.second.diagonal().array())/log(2))) << "with points lc \"red\" pt 7 ps 0.5,"
        << g.file1d(make_pair(ix, log(qrm.second.diagonal().array())/log(2))) << "with points lc \"black\" pt 7 ps 0.5\n";
}
