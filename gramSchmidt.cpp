#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include "myGram.h"

template<int m, int n>
using M = Matrix<double, m, n>;
template<int m>
using V = Vector<double, m>;

template<int m, int n>
void showQR(const M<m,n> a, const pair<M<m,n>, M<n,n>>& qr) {
    auto q = qr.first;
    auto r = qr.second;

    cout << "Q =\n" << q << "\n\n";
    cout << "R =\n" << r << "\n\n";
    cout << "Q*Q =\n" << q.adjoint()*q << "\n\n";
    cout << "A-QR =\n" << a-q*r << "\n\n";
    cout << "FrobN of diff = " << (a-q*r).norm() << "\n\n";
}

int main() {
    srand(time(0));
    constexpr int m = 5, n = 3;

    M<m,n> a = M<m,n>::Random();
    cout << "A =\n" << a << "\n\n";

    showQR(a, QRc(a));
    showQR(a, QRm(a));
}
