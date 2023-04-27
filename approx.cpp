#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include <cmath>
#include "householder.h"
#include "gnuplot-iostream.h"

Mx buildA(const Ax& xs) {
    int m = xs.rows();
    Mx a = Mx::Zero(m,3);
    a << sin(xs), exp(xs), xs.unaryExpr([](double x) { return tgamma(x); });
    return a;
}

template<typename F>
void plotF(Gnuplot& g, string name, const F& f) {
    A<150> xs = A<150>::LinSpaced(0.2, 3);
    Array<double,150,2> pts;
    pts << xs, f(xs);
    g << g.file1d(pts) << "with lines title '" << name << "',";
}

void showF(Gnuplot& g) {
    g << "set arrow from 1,0 to 1,10 nohead linecolor 'gray'\n";
    g << "set arrow from 2,0 to 2,10 nohead linecolor 'gray'\n";
    g << "plot";
    plotF(g, "1/x", [](auto xs) { return xs.pow(-1); });
}

void showCurve(Gnuplot& g, const V<3>& cs) {
    stringstream buffer;
    buffer.precision(3);
    buffer << cs(0) << "sin(x) + " << cs(1) << "e^x + " << cs(2) << "⌈(x)";
    plotF(g, buffer.str(), [&](auto xs) { return buildA(xs)*cs; });
}

template<int n>
double diffNorm(const V<3>& cs, double low, double high) {
    A<n> xs = A<n>::LinSpaced(low, high);
    A<n> ys = xs.pow(-1);
    A<n> est = buildA(xs)*cs;
    V<n> diff = (ys-est).abs();
    A<n-1> traps = (diff(seqN(0,n-1))+diff(lastN(n-1)))/2;
    traps = traps.pow(2)*(high-low)/n;
    return traps.sum();
}

int main() {
    double low = 1, high = 2;
    Gnuplot g;

    showF(g);

    for (int m = 2; m < 10; m++) {
        Ax xs = Ax::LinSpaced(pow(2,m),low,high);
        Vx b = xs.pow(-1);
        auto a = buildA(xs);

        auto [q,r] = thinQRh(a);
        V<3> cs = r.triangularView<Upper>().solve(q.adjoint()*b);
        
        showCurve(g, cs);
        auto norm = diffNorm<200>(cs, 1, 2);
        cout << "Approximate L2 norm of difference for m = " << pow(2,m) << ": " << norm << "\n";
    }

    g << "\n";
}
