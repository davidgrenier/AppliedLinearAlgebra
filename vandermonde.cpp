#include <iostream>
#include <numbers>
#define EIGEN_USE_BLAS
#include <Eigen/Dense>
#include "gnuplot-iostream.h"

typedef std::vector<double> vec;

constexpr double pi = std::numbers::pi;
constexpr double tau = 2*pi;

Eigen::ArrayXd evalP(const Eigen::ArrayXd& poly, const Eigen::ArrayXd& points) {
    auto n = poly.size()-1;
    Eigen::ArrayXd result = points*0+poly[n--];
    while (n >= 0)
        result = result*points+poly[n--];
    return result;
}

Eigen::ArrayXd range(double low, double step, double high) {
    return Eigen::ArrayXd::LinSpaced((high-low)/step, low, high);
}

int main() {
    Gnuplot gp;

    int pdegree = 8;
    int npts = 2*pdegree+1;

    auto xs = Eigen::ArrayXd::LinSpaced(100,-4,4);
    auto ys = (-xs.square()/2).exp();

    auto pts = Eigen::ArrayXd::LinSpaced(npts,-pi,pi);
    auto ypts = (-pts.square()/2).exp().matrix();

    Eigen::MatrixXd A(npts,pdegree);
    for (int i = 0; i < pdegree; i++)
        A.col(i) = pts.pow(i).matrix();

    Eigen::MatrixXd At = A.transpose();
    Eigen::VectorXd projYs = At*ypts;
    Eigen::VectorXd coefs = (At*A).llt().solve(projYs);

    Eigen::MatrixXd AtA = At*A;
    Eigen::VectorXd coefs2 = AtA.inverse()*projYs;
    Eigen::VectorXd coefs3 = AtA.householderQr().solve(At*ypts);
    Eigen::VectorXd coefs4 = A.householderQr().solve(ypts);

    Eigen::MatrixXd S(coefs.size(),4);
    S << coefs, coefs2, coefs3, coefs4;;
    std::cout << std::setprecision(4) << std::scientific;
    std::cout << S << std::endl;

    Eigen::MatrixXd E = Eigen::MatrixXd::Zero(coefs.size(),4);
    E.col(1) = coefs-coefs2;
    E.col(2) = coefs-coefs3;
    E.col(3) = coefs-coefs4;
    std::cout << E << std::endl;

    auto est = evalP(coefs,xs);

    gp << "plot [x=-4:4] [y=-3:2]" << gp.file1d(std::make_pair(xs,ys)) << "with lines,"
        << gp.file1d(std::make_pair(pts,ypts)) << "with points,"
        << gp.file1d(std::make_pair(xs, est)) << "with lines\n";
}
