#ifndef MINKINDR_QUATERNION_GTSAM_H
#define MINKINDR_QUATERNION_GTSAM_H

#include <gtsam/base/Manifold.h>
#include <kindr/minimal/rotation-quaternion.h>
#include <kindr/minimal/common.h>
#include <gtsam_unstable/nonlinear/Expression.h>

namespace gtsam {
namespace traits {

// Required by gtsam.
template<>
struct dimension<kindr::minimal::RotationQuaternion> :
    public boost::integral_constant<int, 3> {
};

// Required by gtsam.
template<>
struct equals<kindr::minimal::RotationQuaternion> {
  bool operator()(const kindr::minimal::RotationQuaternion& q1,
                  const kindr::minimal::RotationQuaternion& q2, double tol) {
    return (q1.getUnique().vector() - q2.getUnique().vector()).array()
        .abs().maxCoeff() < tol;
  }
};

// Required by gtsam.
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

  static vector local(const type& origin, const type& other) {
    return (other * origin.inverted()).log();
  }
  static type retract(const type& origin, const vector& d) {
    return type(d) * origin;
  }
  static int getDimension(const type& /* origin */) {
    return traits::dimension<kindr::minimal::RotationQuaternion>::value;
  }
};

typedef Eigen::Matrix<double, 3, 6> Jacobian3x6;
typedef Eigen::Matrix<double, 3, 3> Jacobian3x3;

inline Eigen::Vector3d rotate_point(
    const kindr::minimal::RotationQuaternion& C, const Eigen::Vector3d& p,
    boost::optional<Jacobian3x3&> HC, boost::optional<Jacobian3x3&> Hp) {
  Eigen::Vector3d Cp = C.rotate(p);
  if (HC) {
    kindr::minimal::skewMatrix(-Cp, &(*HC));
  }
  if (Hp) {
    *Hp = C.getRotationMatrix();
  }
  return Cp;
}

inline kindr::minimal::RotationQuaternion invert_rotation_quaternion(
    const kindr::minimal::RotationQuaternion& C,
    boost::optional<Jacobian3x3&> HC) {
  if(HC) {
    *HC = -C.getRotationMatrix().transpose();
  }
  return C.inverted();
}

inline kindr::minimal::RotationQuaternion compose_rotation_quaternion(
    const kindr::minimal::RotationQuaternion& C1,
    const kindr::minimal::RotationQuaternion& C2,
    boost::optional<Jacobian3x3&> HC1,
    boost::optional<Jacobian3x3&> HC2) {
  if(HC1) {
    *HC1 = Eigen::Matrix3d::Identity();
  }
  if(HC2) {
    *HC2 = C1.getRotationMatrix();
  }
  return C1 * C2;
}

inline Expression<kindr::minimal::RotationQuaternion> invert(
    const Expression<kindr::minimal::RotationQuaternion>& q) {
  return Expression<kindr::minimal::RotationQuaternion>(&invert_rotation_quaternion, q);
}

inline Expression<kindr::minimal::RotationQuaternion>
operator*(const Expression<kindr::minimal::RotationQuaternion>& C1,
          const Expression<kindr::minimal::RotationQuaternion>& C2) {
  return Expression<kindr::minimal::RotationQuaternion>(&compose_rotation_quaternion, C1, C2);
}

// This is syntatic sugar to be able to write
// Expression<Eigen::Vector3d> Tp = T * p;
// instead of
// Expression<Eigen::Vector3d> Tp = Expression<Eigen::Vector3d>(&transform_point, T, p);
inline gtsam::Expression<Eigen::Vector3d> operator*(
    const gtsam::Expression<kindr::minimal::RotationQuaternion>& C,
    const gtsam::Expression<Eigen::Vector3d>& p) {
  return Expression<Eigen::Vector3d>(&rotate_point, C, p);
}

}  // namespace gtsam

#endif // MINKINDR_QUATERNION_GTSAM_H
