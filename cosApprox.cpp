#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include <cmath>
#include "myGram.h"
#include "householder.h"
#include "gnuplot-iostream.h"

using std::numbers::pi;

template<int m, int n>
M<m,n> vander(double low, double high) {
    A<m> xs = A<m>::LinSpaced(low,high);
    M<m,n> v;
    for (int i = 0; i < n; i++)
        v.col(i) = xs.pow(i);
    return v;
}

template<int m, int n>
V<n> normal(const M<m,n>& a, const V<m>& b) {
    return (a.adjoint()*a).ldlt().solve(a.adjoint()*b);
}

template<int m, int n>
V<n> eigenSolve(const M<m,n>& a, const V<m>& b) {
    return HouseholderQR<M<m,n>>{a}.solve(b);
}

template<int m, int n>
V<n> hhQRSolve(const M<m,n>& a, const V<m>& b) {
    auto [q,r] = thinQRh(a);
    return r.template triangularView<Upper>().solve(q.adjoint()*b);
}

template<int m, int n>
V<n> normalRRSolve(const M<m,n>& a, const V<m>& b) {
    auto [q,r] = thinQRh(a);
    return (r.adjoint()*r).ldlt().solve(a.adjoint()*b);
}

template<int m, int n>
V<n> luOrthProjSolve(const M<m,n>& a, const V<m>& b) {
    auto [q,r] = thinQRh(a);
    return FullPivLU<M<m,n>>{a}.solve(q*(q.adjoint()*b));
}

template<int m, int n>
V<n> gramQRSolve(const M<m,n>& a, const V<m>& b) {
    auto [q,r] = QRm(a);
    return r.template triangularView<Upper>().solve(q.adjoint()*b);
}

template<int m, int n>
V<n> svdSolve(const M<m,n>& a, const V<m>& b) {
    BDCSVD<Mx> svd(Mx{a}, ComputeThinU | ComputeThinV);
    return svd.solve(b);
}

int main() {
    constexpr int m = 50;
    constexpr int n = 12;

    auto v = vander<m,n>(0,1);

    A<m> xs = v.col(1);
    V<m> b = cos(xs*4);
    
    constexpr int methods = 7;
    M<n,methods> cs = M<n,methods>::Zero();
    cs.col(0) = normal(v, b);
    cs.col(1) = eigenSolve(v, b);
    cs.col(2) = hhQRSolve(v, b);
    cs.col(3) = normalRRSolve(v, b);
    cs.col(4) = luOrthProjSolve(v, b);
    cs.col(5) = gramQRSolve(v, b);
    cs.col(6) = svdSolve(v, b);

    cout.precision(11);
    Array<string,n+1,methods> show;
    show << "normal", "eigen hh", "my hh", "rr(hh)", "qq(lu)", "my gram", "svd"
        , cs.array().unaryExpr([] (double x) { return to_string(x); });
    cout << show << "\n\n";

    Gnuplot g;
    g << "plot";
    g << g.file1d(make_pair(xs, b)) << "with lines title 'cos',";
    for (int i = 0; i < methods; i++) {
        V<m> ys = v*cs.col(i);
        g << g.file1d(make_pair(xs,ys)) << "with lines lw 0.25 title '" << show(0,i) << "',";
        cout << "residual(" << show(0,i) << "):" << (b.array()-ys.array()).abs().pow(2).sum() << "\n";
    }
    
    g << "\n";
}
