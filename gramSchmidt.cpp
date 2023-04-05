#include <iostream>
#include <numbers>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
using std::numbers::pi;

typedef Matrix4d M;
typedef Vector4d V;

constexpr auto eps = numeric_limits<double>::epsilon();

pair<V, V> orthogonalize(const V& v, const M& q, int j) {
    V qj = v;
    V rj = V::Zero();
    for (int i = 0; i < j; i++) {
        rj(i) = q.col(i).dot(qj);
        qj -= rj(i)*q.col(i);
    }
    rj(j) = qj.norm();
    qj.normalize();

    return make_pair(qj,rj);
}

int main() {
    srand(time(0));

    M m = M::Random();
    // m.col(2) = m.col(1);
    // m.col(3) = m.col(2);
    cout << "M =\n" << m << endl << endl;

    M q = M::Zero();
    M r = M::Zero();

    for (int j = 0; j < m.cols(); j++) {
        auto qrj = orthogonalize(m.col(j), q, j);
        r.col(j) = qrj.second;
        if (qrj.second(j) > 2*eps) {
            q.col(j) = qrj.first;
        } else {
            r(j,j) = 0;
            for (int i = 0; i < m.cols(); i++) {
                qrj = orthogonalize(V::Unit(i), q, j);
                if (qrj.second(j) > 2*eps) {
                    q.col(j) = qrj.first;
                    break;
                }
            }
        }
    }

    cout << "Q =\n" << q << endl << endl;
    cout << "Q'Q =\n" << q.adjoint()*q << endl << endl;
    cout << "R =\n" << r << endl << endl;
    cout << "M-Q*R =\n" << m-q*r << endl << endl;
    cout << "eps = " << eps << ", Error norm = " << (m-q*r).norm() << endl << endl;
}
