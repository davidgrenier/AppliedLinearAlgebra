#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include "gnuplot-iostream.h"
#include "myGram.h"
#include "householder.h"

template<int m, int n>
M<m,n> vander(double low, double high) {
    A<m> xs = A<m>::LinSpaced(low,high);
    M<m,n> v;
    for (int i = 0; i < n; i++)
        v.col(i) = xs.pow(i);
    return v;
}

template<int m, int n>
A<m> evalPoly(const A<m>& x, const V<n>& coef) {
    A<m> tot = A<m>::Zero();
    for (int i = n-1; i >= 0; i--) {
        tot *= x;
        tot += coef(i,0);
    }
    return tot;
}

int main() {
    V<3> p{-10,1,3};
    M<100,3> v = vander<100,3>(-4,4);
    A<2> roots = A<2>{-2,5.0/3};
    cout << "roots are: " << roots.transpose() << endl;
    A<100> xs = v.col(1);

    M<2,2> a;
    double b = p(1)/p(2);
    double c = p(0)/p(2);
    a << -b-1, 1, -b-1-c, 1;

    A<2> eigv = a.eigenvalues().real();
    cout << "estimated roots are: " << eigv.transpose() << endl;
    cout << "2-norm of residual: " << (eigv-roots).matrix().norm() << endl;
    V<3> p2{eigv(0)*eigv(1), -eigv.sum(), 1};

    Gnuplot g;
    g << "set arrow from " << eigv(0) << ",-4 to " << eigv(0) << ",6 nohead linecolor 'gray'\n";
    g << "set arrow from " << eigv(1) << ",-4 to " << eigv(1) << ",6 nohead linecolor 'gray'\n";
    g << "set arrow from -3,0 to 3,0 nohead linecolor 'gray'\n";
    g << "plot"
        << g.file1d(make_pair(xs, v*p)) << "with lines lc 'red' lw 2,";

    g << g.file1d(make_pair(xs, v*p2)) << "with lines lc 'black' lw 1,";

    g << "\n";
}
