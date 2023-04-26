#include <iostream>
#include <Eigen/Dense>
#include "myGram.h"
#include "householder.h"

int main() {
    M<2,2> a {
        {0.700000, 1},
        {0.700001, 1}
    };
    auto [q,r] = QRm(a);
    cout << "Q =\n" << q << "\n\n";
    cout << "R =\n" << r << "\n\n";
    cout << "q1*q2 = " << q.col(0).dot(q.col(1)) << "\n";
    cout << "norm(Q*Q-I) = " << (q.adjoint()*q-M<2,2>::Identity()).norm() << "\n\n";

    auto [qh, rh] = QRh(a);
    cout << "Qh =\n" << qh << "\n\n";
    cout << "Rh =\n" << rh << "\n\n";
    cout << "q1*q2 = " << qh.col(0).dot(qh.col(1)) << "\n";
    cout << "norm(Q*Q-I) = " << (qh.adjoint()*qh-M<2,2>::Identity()).norm() << "\n\n";
}
