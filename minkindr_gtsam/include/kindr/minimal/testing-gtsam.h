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
#ifndef KINDR_MINIMAL_TESTING_GTSAM_H
#define KINDR_MINIMAL_TESTING_GTSAM_H

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include <gtsam/linear/VectorValues.h>
#include <gtsam/nonlinear/Expression.h>
#include <gtsam/nonlinear/ExpressionFactor.h>
#include <eigen-checks/gtest.h>

namespace gtsam {

// Compute finite difference Jacobians for an expression factor.
// This is extremely useful for unit testing during expression development.
template<typename T>
gtsam::JacobianFactor computeFiniteDifferenceJacobians(gtsam::ExpressionFactor<T> expression_factor,
                                                       const gtsam::Values& values,
                                                       double fd_step) {
  Eigen::VectorXd e = expression_factor.unwhitenedError(values);
  const size_t rows = e.size();

  std::map<gtsam::Key, Eigen::MatrixXd> jacobians;
  typename gtsam::ExpressionFactor<T>::const_iterator key_it = expression_factor.begin();
  gtsam::VectorValues dX = values.zeroVectors();
  for(; key_it != expression_factor.end(); ++key_it) {
    size_t key = *key_it;
    // Compute central differences using the values struct.
    const size_t cols = dX.dim(key);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(rows, cols);
    for(size_t col = 0; col < cols; ++col) {
      Eigen::VectorXd dx = Eigen::VectorXd::Zero(cols);
      dx[col] = fd_step;
      dX[key] = dx;
      gtsam::Values eval_values = values.retract(dX);
      Eigen::VectorXd left = expression_factor.unwhitenedError(eval_values);
      dx[col] = -fd_step;
      dX[key] = dx;
      eval_values = values.retract(dX);
      Eigen::VectorXd right = expression_factor.unwhitenedError(eval_values);
      J.col(col) = (left - right) * (1.0/(2.0 * fd_step));
    }
    jacobians[key] = J;
  }

  // Next step...return JacobianFactor
  return gtsam::JacobianFactor(jacobians, -e);
}

template<typename T>
void testExpressionJacobians(gtsam::Expression<T> expression,
                             const gtsam::Values& values,
                             double fd_step,
                             double tolerance) {
  using namespace gtsam;
  // Create factor
  size_t size = traits<T>::dimension;
  ExpressionFactor<T> f(noiseModel::Unit::Create(size), expression.value(values), expression);

  // Check linearization
  JacobianFactor expected = computeFiniteDifferenceJacobians(f, values, fd_step);
  boost::shared_ptr<GaussianFactor> gf = f.linearize(values);
  boost::shared_ptr<JacobianFactor> jf = //
      boost::dynamic_pointer_cast<JacobianFactor>(gf);

  typedef std::pair<Eigen::MatrixXd, Eigen::VectorXd> Jacobian;
  Jacobian evalJ = jf->jacobianUnweighted();
  Jacobian estJ = expected.jacobianUnweighted();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(evalJ.first, estJ.first, tolerance))
    << "Mismatch in the Jacobian matrix.";
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(evalJ.second, Eigen::VectorXd::Zero(size), tolerance))
    << "Mismatch in the error vector.";
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(estJ.second, Eigen::VectorXd::Zero(size), tolerance))
    << "Mismatch in the error vector.";

}

// Do a full concept check and test the invertibility of
// local() vs. retract().
template<typename T>
void testDefaultChart(const T& value) {
  T other = value;
  gtsam::traits<T>::Print(value, "value");
  gtsam::traits<T>::Print(other, "other");
  EXPECT_TRUE(gtsam::traits<T>::Equals(value, other, 1e-12));

  typedef typename gtsam::traits<T> Traits;
  typedef typename Traits::vector Vector;

  Vector dx = Traits::Local(value, other);
  EXPECT_EQ(Traits::dimension, dx.size());
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Eigen::VectorXd::Zero(dx.size()), dx, 1e-9));

  dx.setRandom();
  T updated = Traits::Retract(value, dx);
  Vector invdx = Traits::Local(value, updated);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(dx, invdx, 1e-9));

  dx = -dx;
  updated = Traits::Retract(value, dx);
  invdx = Traits::Local(value, updated);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(dx, invdx, 1e-9));
}
}  // namespace gtsam

#endif // KINDR_MINIMAL_TESTING_GTSAM_H
