#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "myGram.h"
#include "gnuplot-iostream.h"

#include "CImg.h"
using namespace cimg_library;

using std::numbers::pi;

template<int m, int n>
void show(const M<m,n>& a) {
    CImg<float> img(n,m);
    cimg_forXY(img, x, y) {
        img(x,y) = a(y,x);
    }
    img.display();
}

template<int m>
Matrix<complex<double>, m, m> dft() {
    Matrix<complex<double>, m, m> a = Matrix<complex<double>, m, m>::Zero();
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            a(i,j) = exp(-2*pi*(i*j%m)/m);
    return a;
}

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

    BDCSVD<M<m,n>> svd(hello, ComputeFullU | ComputeFullV);

    Gnuplot g;
    g << "set logscale y\nplot" << g.file1d(make_pair(A<m>::LinSpaced(0,m-1), svd.singularValues())) << "with lines\n";

    M<m,m> U = svd.matrixU();
    M<n,n> V = svd.matrixV();
    M<m,n> S;
    S.diagonal() = svd.singularValues();

    for (int i = 1; i <= 10; i++) {
        auto ns = seqN(0,i);
        // ns = seqN(i-1,11-i);
        M<m,n> happrox = U(all, ns)*S(ns, ns)*V(all, ns).adjoint();
        cout << i << ": " << (hello-happrox).norm() << "\n\n";

        show(happrox);
    }
}
