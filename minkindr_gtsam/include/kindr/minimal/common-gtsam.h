#ifndef MINKINDR_COMMON_GTSAM_H
#define MINKINDR_COMMON_GTSAM_H

#include <Eigen/Core>

#include <gtsam/base/Manifold.h>
#include <gtsam/nonlinear/Expression.h>

typedef Eigen::Matrix<double, 3, 6> Jacobian3x6;
typedef Eigen::Matrix<double, 3, 3> Jacobian3x3;
typedef Eigen::Matrix<double, 6, 3> Jacobian6x3;
typedef Eigen::Matrix<double, 6, 6> Jacobian6x6;

namespace kindr {
namespace minimal {

template <int N>
Eigen::Matrix<double, N, 1> vectorScalingImplementation(const Eigen::Matrix<double, N, 1> & v, double alpha,
                                                        gtsam::OptionalJacobian<N, N> H1,
                                                        gtsam::OptionalJacobian<N, 1> H2);

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > vectorScaling(const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v,
                                                              double alpha);

template <int N>
Eigen::Matrix<double, N, 1> vectorSumImplementation(const Eigen::Matrix<double, N, 1> & v1,
                                                    const Eigen::Matrix<double, N, 1> & v2,
                                                    gtsam::OptionalJacobian<N, N> H1,
                                                    gtsam::OptionalJacobian<N, N> H2);

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > vectorSum(const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v1,
                                                          const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v2);

template <int N>
Eigen::Matrix<double, N, 1> vectorDifferenceImplementation(const Eigen::Matrix<double, N, 1> & v1,
                                                           const Eigen::Matrix<double, N, 1> & v2,
                                                           gtsam::OptionalJacobian<N, N> H1,
                                                           gtsam::OptionalJacobian<N, N> H2);

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > vectorDifference(const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v1,
                                                                 const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v2);

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > operator*(const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v,
                                                          double alpha);

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > operator/(const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v, double alpha);

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > operator-(const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v);

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > operator+(const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v1,
                                                          const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v2);

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > operator-(const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v1,
                                                          const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v2);

}
}

#include <kindr/minimal/implementation/common-gtsam-inl.h>

#endif // MINKINDR_COMMON_GTSAM_H
