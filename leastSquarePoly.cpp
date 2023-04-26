#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include "householder.h"
#include "gnuplot-iostream.h"

void plotPoly(Gnuplot& g, const Vx& p, const string& title) {
    int n = p.rows();
    constexpr int m = 150;
    A<m> xs = A<m>::LinSpaced(-1,1);
    Mx v = Mx::Zero(m,n);
    for (int i = 0; i < n; i++)
        v.col(i) = xs.pow(i);
    A<m> ys = v*p;
    g << g.file1d(make_pair(xs,ys)) << "with lines title '" << title << "',";
}

int main() {
    auto t = time(0);
    cout << t << "\n";
    srand(t);

    constexpr int m = 11;
    V<m> x = V<m>::LinSpaced(-1,1);
    V<m> y = V<m>::Random();
    M<m,m> v;
    for (int i = 0; i < m; i++)
        v.col(i) = x.array().pow(i);

    Gnuplot g;
    g << "plot" << g.file1d(make_pair(x, y)) << "with points title 'data',";

    array<int,5> degs{11,10,9,7,5};
    for (int i = 0; i < 5; i++) {
        auto [q,r] = thinQRh(v(all,seqN(0,degs[i])));
        Vx c = r.triangularView<Upper>().solve(q.adjoint()*y);

        plotPoly(g, c, "degree = " + to_string(degs[i]));
    }

    g << "\n";
}
