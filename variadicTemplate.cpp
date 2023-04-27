#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include <cmath>
#include "householder.h"
#include "gnuplot-iostream.h"

using std::numbers::pi;

template<int m, int n = 1, typename F>
M<m,2> points(double low, double high, F f) {
    A<m> xs = A<m>::LinSpaced(low, high);
    M<m,2> pts;
    pts << xs, f(xs);
    return pts;
}

template<int m, int n, typename G, typename... F>
M<m,n+1> points(double low, double high, G g, F... rest) {
    M<m,n> xs = points<m,n-1>(low, high, rest...);
    M<m,n+1> ys;
    ys << xs, g(xs.col(0));
    return ys;
}

int main() {
    Gnuplot g;
    g << "plot";
    auto pts = points<100,2>(-pi, pi
            , [] (A<100> xs) { return cos(xs); }
            , [] (auto xs) { return xs.sin(); });
    for (int i = 1; i < 3; i++)
        g << g.file1d(pts(all,{0,i})) << "with lines,";
    g << "\n";
}
