#ifndef MINKINDR_QUATERNION_GTSAM_H
#define MINKINDR_QUATERNION_GTSAM_H

#include <gtsam/base/Manifold.h>
#include <gtsam/nonlinear/Expression.h>

#include <kindr/minimal/rotation-quaternion.h>
#include <kindr/minimal/common.h>

#include "common-gtsam.h"

namespace gtsam {
namespace traits { // Traits define the basic interface required by gtsam.

// The dimension of the manifold.
template<>
struct dimension<kindr::minimal::RotationQuaternion> :
    public boost::integral_constant<int, 3> {
};

// Check the equality of two values.
template<>
struct equals<kindr::minimal::RotationQuaternion> {
  bool operator()(const kindr::minimal::RotationQuaternion& q1,
                  const kindr::minimal::RotationQuaternion& q2, double tol) {
    return (q1.getUnique().vector() - q2.getUnique().vector()).array()
        .abs().maxCoeff() < tol;
  }
};

// Print the type.
template<>
struct print<kindr::minimal::RotationQuaternion> {
  void operator()(const kindr::minimal::RotationQuaternion& q,
                  const std::string& str) {
    if (str.size() > 0) {
      std::cout << str << ": ";
    }
    std::cout << "kindr::minimal::RotationQuaternion[" << q.vector().transpose() << "]" << std::endl;
  }
};

}  // namespace traits

// Chart is a map from T -> vector, retract is its inverse
// This is required by gtsam because our classes don't have the
// correct interface.
template<>
struct DefaultChart<kindr::minimal::RotationQuaternion> {
  // the "type" typedef is used by gtsam.
  typedef kindr::minimal::RotationQuaternion type;
  // the "vector" typedef is used by gtsam.
  typedef Eigen::Matrix<double, traits::dimension<type>::value, 1> vector;

  static inline vector local(const type& origin, const type& other) {
    return (other * origin.inverted()).log();
  }
  static inline type retract(const type& origin, const vector& d) {
    return kindr::minimal::RotationQuaternion(d) * origin;
  }
  static inline int getDimension(const type& /* origin */) {
    return traits::dimension<kindr::minimal::RotationQuaternion>::value;
  }
};

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
