#include <iostream>
#include <numbers>
#define EIGEN_USE_BLAS
#include <Eigen/Dense>
#include "gnuplot-iostream.h"

template<int m, int n>
using Mat = Eigen::Matrix<double,m,n>;

using M4 = Mat<4,4>;

void e1pa(const M4& B) {
    auto id = M4::Identity();

    M4 dc1 = id;
    dc1(0,0) = 2;

    M4 hr3 = id;
    hr3(2,2) = .5;

    M4 ar3r1 = id;
    ar3r1(0,2) = 1;

    M4 ic14 = id;
    ic14.col(0) << 0,0,0,1;
    ic14.col(3) << 1,0,0,0;

    M4 repc4byc3 = id;
    repc4byc3.col(3) << 0,0,1,0;

    M4 sr2 = id;
    sr2(0,1) = -1;
    sr2(2,1) = -1;
    sr2(3,1) = -1;

    Mat<4,3> dropc1 = id.rightCols(3);

    std::cout << sr2*ar3r1*hr3 << "\n* B *\n";
    std::cout << dc1*ic14*repc4byc3*dropc1 << "\n=\n";
    std::cout << sr2*ar3r1*hr3*B*dc1*ic14*repc4byc3*dropc1 << "\n\n";
}

void e1pb(const M4& B) {
    auto id = M4::Identity();

    Mat<4,3> C = id.rightCols(3);
    C.col(2) << 0,0,1,0;

    Mat<4,4> A;
    A << 1,-1,.5,0,
        0,1,0,0,
        0,-1,.5,0,
        0,-1,0,1;

    std::cout << A << "\n* B *\n";
    std::cout << A*B*C << std::endl;
}

int main() {
    M4 B = M4::Random();
    e1pa(B);
    e1pb(B);
}
