#include <iostream>
#include <Eigen/Dense>
#include "myGram.h"
#include "householder.h"

int main() {
    M<2,2> a {
        {0.70000, 1},
        {0.70001, 1}
    };
    auto qra = QRm(a);
    auto q = qra.first;
    cout << q << "\n\n";
    cout << q.col(0).dot(q.col(1)) << "\n";
    cout << (q.adjoint()*q-M<2,2>::Identity()).norm() << "\n\n";

    auto qrh = QRh(a);
    auto qh = qrh.first;
    cout << qh << "\n\n";
    cout << qh.col(0).dot(qh.col(1)) << "\n";
    cout << (qh.adjoint()*qh-M<2,2>::Identity()).norm() << "\n\n";
}
