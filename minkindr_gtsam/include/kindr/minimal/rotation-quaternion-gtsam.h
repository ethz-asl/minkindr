#ifndef MINKINDR_QUATERNION_GTSAM_H
#define MINKINDR_QUATERNION_GTSAM_H

#include <gtsam/base/Manifold.h>
#include <gtsam/nonlinear/Expression.h>

#include <kindr/minimal/rotation-quaternion.h>
#include <kindr/minimal/common.h>

#include "common-gtsam.h"

namespace gtsam {
template<> struct traits<kindr::minimal::RotationQuaternion> {
  // The dimension of the manifold.
  enum {
    dimension = 3
  };

  typedef kindr::minimal::RotationQuaternion type;
  // the "vector" typedef is used by gtsam.
  typedef Eigen::Matrix<double, dimension, 1> vector;

  // Print the type.
  static void Print(const kindr::minimal::RotationQuaternion& q,
                    const std::string& str = "") {
    if (str.size() > 0) {
      std::cout << str << ": ";
    }
    std::cout << "kindr::minimal::RotationQuaternion[" << q.vector().transpose() << "]" << std::endl;
  }

  // Check the equality of two values.
  static bool Equals(const kindr::minimal::RotationQuaternion& q1,
                     const kindr::minimal::RotationQuaternion& q2, double tol) {
    return (q1.getUnique().vector() - q2.getUnique().vector()).array()
        .abs().maxCoeff() < tol;
  }

  static vector Local(const type& origin, const type& other) {
    return (other * origin.inverted()).log();
  }
  static type Retract(const type& origin, const vector& d) {
    return kindr::minimal::RotationQuaternion(d) * origin;
  }
  static int GetDimension(const type& /* origin */) {
    return dimension;
  }
};  // traits

////////////////////////////////////////////////////////////////////////////////
// Convenience functions to make working with expressions easy and fun!

/// \brief Invert a rotation quaternion expression.
Expression<kindr::minimal::RotationQuaternion> invert(
    const Expression<kindr::minimal::RotationQuaternion>& q);

/// \brief Compose two quaternion expressions.
Expression<kindr::minimal::RotationQuaternion>
operator*(const Expression<kindr::minimal::RotationQuaternion>& C1,
          const Expression<kindr::minimal::RotationQuaternion>& C2);

/// \brief Rotate a point.
///
/// This is syntatic sugar to be able to write
/// Expression<Eigen::Vector3d> Cp = C * p;
/// instead of
/// Expression<Eigen::Vector3d> Cp = Expression<Eigen::Vector3d>(&rotate_point, C, p);
gtsam::Expression<Eigen::Vector3d> operator*(
    const gtsam::Expression<kindr::minimal::RotationQuaternion>& C,
    const gtsam::Expression<Eigen::Vector3d>& p);

/// \brief Rotate a point.
gtsam::Expression<Eigen::Vector3d> rotate(
    const gtsam::Expression<kindr::minimal::RotationQuaternion>& C,
    const gtsam::Expression<Eigen::Vector3d>& p);

/// \brief Rotate a point by the inverse of the rotation.
gtsam::Expression<Eigen::Vector3d> inverseRotate(
    const gtsam::Expression<kindr::minimal::RotationQuaternion>& C,
    const gtsam::Expression<Eigen::Vector3d>& p);

/// \brief Expose the rotation log and Jacobian.
Eigen::Vector3d rotationLogImplementation(const kindr::minimal::RotationQuaternion& C,
                                          OptionalJacobian<3, 3> JC);

/// \brief Compute the matrix log of SO3.
gtsam::Expression<Eigen::Vector3d> log(
    const gtsam::Expression<kindr::minimal::RotationQuaternion>& C);

/// \brief Expose the rotation log and Jacobian.
kindr::minimal::RotationQuaternion rotationExpImplementation(const Eigen::Vector3d& p,
                                                             OptionalJacobian<3, 3> Jp);

/// \brief Compute the matrix log of SO3.
gtsam::Expression<kindr::minimal::RotationQuaternion> exp(
    const gtsam::Expression<Eigen::Vector3d>& C);
}  // namespace gtsam

#endif // MINKINDR_QUATERNION_GTSAM_H
