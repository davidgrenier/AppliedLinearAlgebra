#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "myGram.h"
// #include "gnuplot-iostream.h"
#include "CImg.h"

using namespace cimg_library;

int main() {
    constexpr int m = 15;
    constexpr int n = 40;
    M<8,6> h = M<8,6>::Ones();
    h(all,seqN(2,2)) = M<8,2>::Zero();
    h(seqN(3,2), seqN(2,2)) = M<2,2>::Ones();

    M<8,6> e = M<8,6>::Ones();
    e(2, seqN(2,4)) = M<1,4>::Zero();
    e(5, seqN(2,4)) = M<1,4>::Zero();

    M<8,6> l = M<8,6>::Ones();
    l(seqN(0,6), seqN(2,4)) = M<6,4>::Zero();

    M<m,n> hello = M<m,n>::Zero();
    hello(seqN(2,8), seqN(2,6)) = h;
    hello(seqN(3,8), seqN(10,6)) = e;
    hello(seqN(4,8), seqN(18,6)) = l;
    hello(seqN(5,8), seqN(26,6)) = l;
    hello(seqN(6,8), seqN(34,6)) = M<8,6>::Ones();
    hello(seqN(8,4), seqN(36,2)) = M<4,2>::Zero();

    auto qr = QRm(hello);
    M<m,n-1> highqr = qr.first*M<n,n-1>::Identity();
    hello = highqr*highqr.adjoint()*hello;

    CImg<float> img(n,m);
    cimg_forXY(img, x, y) {
        img(x,y) = hello(y,x);
    }
    img.display();

    // Gnuplot g;
    // g << "plot" << g.file1d(m, "test.tmp") << "\n";
}
