// Copyright (c) 2015, Autonomous Systems Lab, ETH Zurich
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the <organization> nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#ifndef KINDR_MIN_COMMON_GTSAM_INL_H_
#define KINDR_MIN_COMMON_GTSAM_INL_H_

#include <kindr/minimal/common-gtsam.h>

namespace kindr {
namespace minimal {

template <int N>
Eigen::Matrix<double, N, 1> vectorScalingImplementation(const Eigen::Matrix<double, N, 1> & v, double alpha,
                                                        gtsam::OptionalJacobian<N, N> H1,
                                                        gtsam::OptionalJacobian<N, 1> H2) {
  if (H1) {
    *H1 = gtsam::OptionalJacobian<N,N>::Jacobian::Identity()*alpha;
  }
  if (H2) {
    *H2 = v;
  }
  return v*alpha;
}

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > vectorScaling(const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v, double alpha) {
  return gtsam::Expression<Eigen::Matrix<double, N, 1> >(boost::bind(&vectorScalingImplementation<N>, _1, alpha, _2, boost::none), v);
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
Eigen::Matrix<double, N, 1> vectorDifferenceImplementation(const Eigen::Matrix<double, N, 1> & v1,
                                                           const Eigen::Matrix<double, N, 1> & v2,
                                                           gtsam::OptionalJacobian<N, N> H1, gtsam::OptionalJacobian<N, N> H2) {
  if (H1) {
    H1->setIdentity();
  }
  if (H2) {
    H2->setIdentity();
    *H2 = -*H2;
  }
  return v1-v2;
}

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > vectorDifference(const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v1,
                                                                 const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v2) {
  return gtsam::Expression<Eigen::Matrix<double, N, 1> >(vectorDifferenceImplementation<N>, v1, v2);
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

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > operator-(const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v1,
                                                          const gtsam::Expression<Eigen::Matrix<double, N, 1> >&v2) {
  return vectorDifference(v1, v2);
}

}
}
#endif  // KINDR_MIN_ROTATION_ANGLE_AXIS_INL_H_
