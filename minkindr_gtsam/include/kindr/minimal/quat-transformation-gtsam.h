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
typedef gtsam::Expression<QuatTransformation> ETransformation;
typedef gtsam::Expression<gtsam::Vector6> EVector6;

/// \brief Transform a point.
///
/// This is syntactic sugar to be able to write
/// Expression<Eigen::Vector3d> Tp = T * p;
/// instead of
/// Expression<Eigen::Vector3d> Tp = Expression<Eigen::Vector3d>(&transform_point, T, p);
EVector3
operator*(const ETransformation& T,
          const EVector3& p);

/// \brief Transform a point.
EVector3
transform(const ETransformation& T,
          const EVector3& p);

/// \brief Build a transformation expression from a rotation expression and a
/// point expression.
ETransformation transformationFromComponents(
    const EQuaternion& C_A_B,
    const EVector3& A_t_B);

/// \brief Recover the rotation part of a transformation.
EQuaternion rotationFromTransformation(
    const ETransformation& T);

/// \brief Recover the translation part of a transformation.
EVector3 translationFromTransformation(
    const ETransformation& T);

/// \brief Transform a point by the inverse of a transformation.
EVector3 inverseTransform(
    const ETransformation& T,
    const EVector3& p);

/// \brief Invert a transformation.
ETransformation inverse(
    const ETransformation& T);

/// \brief Compose two transformations.
ETransformation compose(
    const ETransformation& T1,
    const ETransformation& T2);

/// \brief Compose two transformations.
inline ETransformation operator*(
    const ETransformation& T1,
    const ETransformation& T2) {
  return compose(T1,T2);
}

/// \brief Compose two transformations as inv(T1)*T2.
ETransformation invertAndCompose(
    const ETransformation& T1,
    const ETransformation& T2);

/// \brief Recover the matrix log of R^3 x SO3
EVector6 transformationLog(const ETransformation& T);

/// \brief The exponential map of R^3 x SO3
ETransformation transformationExp(const EVector6& params);

/// \brief Recover the matrix log of the rotation part of the transformation.
EVector3 rotationLog(
    const ETransformation& T);

/// \brief Recover the matrix log of the translation part of the transformation.
inline EVector3 translationLog(
    const ETransformation& T) {
  return translationFromTransformation(T);
}

ETransformation slerp(
    const ETransformation& T0,
    const ETransformation& T1,
    double alpha);


}  // namespace minimal
}  // namespace kindr

#endif // MINKINDR_QUAT_TRANSFORMATION_GTSAM_H
