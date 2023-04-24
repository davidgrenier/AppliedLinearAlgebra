#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include "myGram.h"
#include "householder.h"

int main() {
    srand(time(0));

    constexpr int m = 6, n = 4;
    M<m,n> a = M<m,n>::Random();
    cout << "A =\n" << a << "\n\n";
    auto qr = QRh(a);
    auto q = qr.first;
    auto r = qr.second;
    cout << "Q =\n" << q << "\n\n";
    cout << "R =\n" << r << "\n\n";
    M<m,n> zer = a-q*r;
    cout << "A-QR =\n" << zer << "\n\n";
    cout << "Q*Q =\n" << q.adjoint()*q << "\n\n";
    cout << "norm(A-QR) = " << zer.norm() << "\n";
    cout << "norm(I-Q*Q) = " << (M<m,m>::Identity()-q.adjoint()*q).norm() << "\n";
}
