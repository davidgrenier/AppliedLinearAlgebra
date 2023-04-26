#include <iostream>
// #define EIGEN_USE_BLAS
#include <Eigen/Dense>
#include "householder.h"
#include "bench.h"

int main() {
    constexpr int m = 1024, n = 16;
    Mx a = Mx::Random(m,n);

    bench("Eigen's", [&]() {
        HouseholderQR<Mx> h(a);
        Mx q = h.householderQ();
        Mx r = h.matrixQR().template triangularView<Upper>();
    });

    bench("Naive", [&]() {
        auto [q,r] = QRh(a);
    });
}
