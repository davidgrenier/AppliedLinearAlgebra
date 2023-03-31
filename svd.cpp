#include <iostream>
// #define EIGEN_USE_BLAS
#include <Eigen/Eigen>
// #include "gnuplot-iostream.h"

using namespace Eigen;
using namespace std;

int main() {
    typedef Matrix<double,2,1> M;
    M m;
    m << 9,16;

    cout << m << endl << endl;

    BDCSVD<M> svdM(m);

    cout << svdM.singularValues() << endl;

    typedef Matrix<double,1,2> N;
    N n;
    n << 16,9;

    cout << n << endl << endl;

    BDCSVD<N> svdN(n);

    cout << svdN.singularValues() << endl;

}
