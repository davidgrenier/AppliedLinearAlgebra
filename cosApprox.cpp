#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include <cmath>
#include "myGram.h"
#include "householder.h"
#include "gnuplot-iostream.h"

using std::numbers::pi;

int main() {

    constexpr int m = 50;
    constexpr int n = 12;
    A<m> xs = A<m>::LinSpaced(0,1);

    M<m,n> v;
    for (int i = 0; i < n; i++)
        v.col(i) = xs.pow(i);

    V<m> b = cos(xs*4);
    
    constexpr int methods = 7;
    M<n,methods> cs = M<n,methods>::Zero();
    cs.col(0) = (v.adjoint()*v).ldlt().solve(v.adjoint()*b);

    HouseholderQR<M<m,n>> qre{v};
    cs.col(1) = qre.solve(b);

    auto [q,r] = thinQRh(v);
    cs.col(2) = r.triangularView<Upper>().solve(q.adjoint()*b);

    cs.col(3) = (r.adjoint()*r).ldlt().solve(v.adjoint()*b);
    FullPivLU<M<m,n>> lu{v};
    cs.col(4) = lu.solve(q*q.adjoint()*b);


    auto [qg, rg] = QRm(v);
    cs.col(5) = rg.triangularView<Upper>().solve(qg.adjoint()*b);

    Mx forSvd{v};
    BDCSVD<Mx> svd(forSvd, ComputeThinU | ComputeThinV);
    cs.col(6) = svd.solve(b);

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
