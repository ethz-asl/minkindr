#ifndef KINDR_MINIMAL_TESTING_GTSAM_H
#define KINDR_MINIMAL_TESTING_GTSAM_H

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include <gtsam/linear/VectorValues.h>
#include <gtsam_unstable/nonlinear/Expression.h>
#include <gtsam_unstable/nonlinear/ExpressionFactor.h>
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
  size_t size = traits::dimension<T>::value;
  ExpressionFactor<T> f(noiseModel::Unit::Create(size), expression.value(values), expression);

  // Check linearization
  //JacobianFactor expected(1, eye(3,6), 2, Tval.getRotationMatrix(), -(Tval * pval));
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
}  // namespace gtsam

#endif // KINDR_MINIMAL_TESTING_GTSAM_H
