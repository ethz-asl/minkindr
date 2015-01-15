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
} //namespace gtsam

namespace kindr {
namespace minimal {
////////////////////////////////////////////////////////////////////////////////
// Convenience functions to make working with expressions easy and fun!
typedef gtsam::Expression<RotationQuaternion> EQuaternion;
typedef gtsam::Expression<Eigen::Vector3d> EVector3;

/// \brief Invert a rotation quaternion expression.
EQuaternion invert(
    const EQuaternion& q);

/// \brief Compose two quaternion expressions.
EQuaternion
operator*(const EQuaternion& C1,
          const EQuaternion& C2);

/// \brief Rotate a point.
///
/// This is syntatic sugar to be able to write
/// Expression<Eigen::Vector3d> Cp = C * p;
/// instead of
/// Expression<Eigen::Vector3d> Cp = Expression<Eigen::Vector3d>(&rotate_point, C, p);
gtsam::Expression<Eigen::Vector3d> operator*(const EQuaternion& C, const EVector3& p);

/// \brief Rotate a point.
EVector3 rotate(const EQuaternion& C, const EVector3& p);

/// \brief Rotate a point by the inverse of the rotation.
EVector3 inverseRotate(const EQuaternion& C, const EVector3& p);

/// \brief Expose the rotation log and Jacobian.
Eigen::Vector3d rotationLogImplementation(const kindr::minimal::RotationQuaternion& C,
                                          gtsam::OptionalJacobian<3, 3> JC);

/// \brief Compute the matrix log of SO3.
EVector3 quaternionLog(const EQuaternion& C);

/// \brief Expose the rotation log and Jacobian.
kindr::minimal::RotationQuaternion rotationExpImplementation(const Eigen::Vector3d& p,
                                                             gtsam::OptionalJacobian<3, 3> Jp);

/// \brief Compute the matrix log of SO3.
EQuaternion quaternionExp(const EVector3& C);

}  // namespace minimal
}  // namespace kindr

#endif // MINKINDR_QUATERNION_GTSAM_H
