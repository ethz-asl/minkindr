#ifndef MINKINDR_COMMON_GTSAM_H
#define MINKINDR_COMMON_GTSAM_H

#include <Eigen/Core>

typedef Eigen::Matrix<double, 3, 6> Jacobian3x6;
typedef Eigen::Matrix<double, 3, 3> Jacobian3x3;
typedef Eigen::Matrix<double, 6, 3> Jacobian6x3;
typedef Eigen::Matrix<double, 6, 6> Jacobian6x6;

namespace kindr {
namespace minimal {

template <int N>
Eigen::Matrix<double, N, 1> vectorScalingImplementation(const Eigen::Matrix<double, N, 1> & v, double alpha, gtsam::OptionalJacobian<N, N> H) {
  if (H) {
    *H = gtsam::OptionalJacobian<N,N>::Jacobian::Identity()*alpha;
  }
  return v*alpha;
}

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > vectorScaling(const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v, double alpha) {
  return gtsam::Expression<Eigen::Matrix<double, N, 1> >(boost::bind(&vectorScalingImplementation<N>, _1, alpha, _2), v);
}

template <int N>
Eigen::Matrix<double, N, 1> vectorSumImplementation(const Eigen::Matrix<double, N, 1> & v1, const Eigen::Matrix<double, N, 1> & v2,
                                                    gtsam::OptionalJacobian<N, N> H1, gtsam::OptionalJacobian<N, N> H2) {
  if (H1) {
    H1->setIdentity();
  }
  if (H2) {
    H2->setIdentity();
  }
  return v1+v2;
}

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > vectorSum(const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v1,
                                                          const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v2) {
  return gtsam::Expression<Eigen::Matrix<double, N, 1> >(vectorSumImplementation<N>, v1, v2);
}


template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > operator*(const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v, double alpha) {
  return vectorScaling(v, alpha);
}

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > operator/(const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v, double alpha) {
  return vectorScaling(v, 1.0/alpha);
}

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > operator-(const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v) {
  return vectorScaling(v, -1.0);
}

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > operator+(const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v1,
                                                          const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v2) {
  return vectorSum(v1, v2);
}


}
}
#endif // MINKINDR_COMMON_GTSAM_H
