#include <iostream>
#include <Eigen/Dense>
#include "gnuplot-iostream.h"
#include "myGram.h"

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
