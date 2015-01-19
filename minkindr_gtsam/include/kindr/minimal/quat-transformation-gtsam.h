#ifndef MINKINDR_QUAT_TRANSFORMATION_GTSAM_H
#define MINKINDR_QUAT_TRANSFORMATION_GTSAM_H

#include <gtsam/base/Manifold.h>
#include <gtsam/nonlinear/Expression.h>

#include <kindr/minimal/quat-transformation.h>

#include "common-gtsam.h"
#include "rotation-quaternion-gtsam.h"

namespace gtsam {
template<> struct traits<kindr::minimal::QuatTransformation> {
  // The dimension of the manifold.
  enum {
    dimension = 6
  };

  typedef kindr::minimal::QuatTransformation type;
  // the "vector" typedef is used by gtsam.
  typedef Eigen::Matrix<double, dimension, 1> vector;

  // Print the type.
  static void Print(const kindr::minimal::QuatTransformation& T,
                    const std::string& str) {
    if(str.size() > 0) { std::cout << str << ":\n";}
    std::cout << T.getTransformationMatrix() << std::endl;
  }

  // Check the equality of two values.
  static bool Equals(const kindr::minimal::QuatTransformation& T1,
                     const kindr::minimal::QuatTransformation& T2, double tol) {
    return (T1.getTransformationMatrix() - T2.getTransformationMatrix()).array().abs().maxCoeff() < tol;
  }

  static vector Local(const type& origin, const type& other) {
    return (other * origin.inverted()).log();
  }
  static type Retract(const type& origin, const vector& d) {
    return type(d) * origin;
  }
  static int GetDimension(const type& /* origin */) {
    return dimension;
  }
};  // traits
} // namespace gtsam

namespace kindr {
namespace minimal {
////////////////////////////////////////////////////////////////////////////////
// Convenience functions to make working with expressions easy and fun!

/// \brief Transform a point.
///
/// This is syntactic sugar to be able to write
/// Expression<Eigen::Vector3d> Tp = T * p;
/// instead of
/// Expression<Eigen::Vector3d> Tp = Expression<Eigen::Vector3d>(&transform_point, T, p);
gtsam::Expression<Eigen::Vector3d>
operator*(const gtsam::Expression<kindr::minimal::QuatTransformation>& T,
          const gtsam::Expression<Eigen::Vector3d>& p);

/// \brief Transform a point.
gtsam::Expression<Eigen::Vector3d>
transform(const gtsam::Expression<kindr::minimal::QuatTransformation>& T,
          const gtsam::Expression<Eigen::Vector3d>& p);

/// \brief Build a transformation expression from a rotation expression and a
/// point expression.
gtsam::Expression<kindr::minimal::QuatTransformation> transformationFromComponents(
    const gtsam::Expression<kindr::minimal::RotationQuaternion>& C_A_B,
    const gtsam::Expression<Eigen::Vector3d>& A_t_B);

/// \brief Recover the rotation part of a transformation.
gtsam::Expression<kindr::minimal::RotationQuaternion> rotationFromTransformation(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T);

/// \brief Recover the translation part of a transformation.
gtsam::Expression<Eigen::Vector3d> translationFromTransformation(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T);

/// \brief Transform a point by the inverse of a transformation.
gtsam::Expression<Eigen::Vector3d> inverseTransform(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T,
    const gtsam::Expression<Eigen::Vector3d>& p);

/// \brief Invert a transformation.
gtsam::Expression<kindr::minimal::QuatTransformation> inverse(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T);

/// \brief Compose two transformations.
gtsam::Expression<kindr::minimal::QuatTransformation> compose(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T1,
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T2);

/// \brief Compose two transformations.
inline gtsam::Expression<kindr::minimal::QuatTransformation> operator*(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T1,
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T2) {
  return compose(T1,T2);
}

/// \brief Compose two transformations as inv(T1)*T2.
gtsam::Expression<kindr::minimal::QuatTransformation> invertAndCompose(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T1,
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T2);

/// \brief Recover the matrix log of R^3 x SO3
gtsam::Expression<gtsam::Vector6> log(const gtsam::Expression<kindr::minimal::QuatTransformation>& T);

/// \brief The exponential map of R^3 x SO3
gtsam::Expression<kindr::minimal::QuatTransformation> exp(const gtsam::Expression<gtsam::Vector6>& params);

/// \brief Recover the matrix log of the rotation part of the transformation.
gtsam::Expression<Eigen::Vector3d> rotationLog(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T);

/// \brief Recover the matrix log of the translation part of the transformation.
inline gtsam::Expression<Eigen::Vector3d> translationLog(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T) {
  return translationFromTransformation(T);
}

gtsam::Expression<kindr::minimal::QuatTransformation> slerp(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T0,
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T1,
    double alpha);


}  // namespace minimal
}  // namespace kindr

#endif // MINKINDR_QUAT_TRANSFORMATION_GTSAM_H
