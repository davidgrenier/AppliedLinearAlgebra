#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include "myGram.h"
#include "gnuplot-iostream.h"

// Slower
// #include "../cpp/bench.h"
template<int m, int n>
V<m> evalPoly(const A<m>& x, const A<n>& coef) {
    A<m> tot = A<m>::Zero();
    for (int i = n-1; i >= 0; i--) {
        tot *= x;
        tot += coef(i,0);
    }
    return tot;
}

int main() {
    constexpr int m = 256;
    constexpr int n = 4;
    A<m> ix = A<m>::LinSpaced(-1,1.0);
    M<m,n> a;
    for (int i = 0; i < n; i++)
        a.col(i) = ix.pow(i);
    M<m,n> qm = QRm(a).first;
    qm *= qm.row(m-1).array().pow(-1).matrix().asDiagonal();

    M<n,n> legPoly;
    legPoly << 1, 0, -0.5,  0.0
             , 0, 1,  0.0, -1.5
             , 0, 0,  1.5,  0.0
             , 0, 0,  0.0,  2.5;

    M<m,n> realPoly = a*legPoly;
    M<m,n> diff = realPoly-qm;

    Gnuplot g;
    g << "plot";

    for (int i = 0; i < n; i++)
        g << g.file1d(make_pair(ix, realPoly.col(i)-qm.col(i))) << "with lines lc \"blue\" notitle,";

    g << "\n";
}
