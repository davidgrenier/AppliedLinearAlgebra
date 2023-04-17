#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include <gnuplot-iostream.h>

using namespace Eigen;
using namespace std;
using std::numbers::pi;

template<int m, int n>
using M = Matrix<double, m, n>;
template<int n>
using V = Vector<double, n>;
template<int n>
using A = Array<double, n, 1>;

template<int m, int n>
pair<M<m,n>, M<n,n>> qr(const M<m,n>& a) {
    HouseholderQR<M<m,n>> aqr(a);
    M<m,n> q = aqr.householderQ()*M<m,n>::Identity();
    M<m,n> r = aqr.matrixQR().template triangularView<Upper>();
    return make_pair(q,r(seqN(0,n),all));
}

int main() {
    constexpr int m = 257;
    constexpr int n = 5;
    A<m> l = A<m>::LinSpaced(-1, 1);
    M<m, n> a;
    for (int i = 0; i < n; i++)
        a.col(i) = l.pow(i);

    auto aqr = qr(a);

    M<m,n> leg = aqr.first;
    M<n,n> norms;
    norms.diagonal() = leg.row(m-1).array().pow(-1);
    leg *= norms;

    Gnuplot g;
    g << "set yrange [-1.25:1.25]\nplot";
    for (int i = 0; i < n-1; i++)
        g << g.file1d(make_pair(l, leg.col(i))) << "with lines,";
    g << g.file1d(make_pair(l, leg.col(n-1))) << "with lines\n";
}
