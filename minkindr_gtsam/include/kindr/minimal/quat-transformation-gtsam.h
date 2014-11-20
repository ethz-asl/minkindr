#ifndef MINKINDR_QUAT_TRANSFORMATION_GTSAM_H
#define MINKINDR_QUAT_TRANSFORMATION_GTSAM_H

#include <gtsam/base/Manifold.h>
#include <kindr/minimal/quat-transformation.h>
#include <gtsam_unstable/nonlinear/Expression.h>

namespace gtsam {
namespace traits {

//template<>
//struct is_group<kindr::minimal::QuatTransformation> : public boost::false_type {
//};

//template<>
//struct is_manifold<kindr::minimal::QuatTransformation> : public boost::false_type {
//};

// Required by gtsam.
template<>
struct dimension<kindr::minimal::QuatTransformation> : public boost::integral_constant<int,6> {
};

// Required by gtsam.
template<>
struct equals<kindr::minimal::QuatTransformation>{
  bool operator()(const kindr::minimal::QuatTransformation& T1,
                  const kindr::minimal::QuatTransformation& T2, double tol) {
    return (T1.getTransformationMatrix() - T2.getTransformationMatrix()).array().abs().maxCoeff() > tol;
  }
};

// Required by gtsam.
template<>
struct print<kindr::minimal::QuatTransformation>{
  void operator()(const kindr::minimal::QuatTransformation& T, const std::string& str) {
    if(str.size() > 0) { std::cout << str << ":\n";}
    std::cout << T.getTransformationMatrix() << std::endl;
  }
};

}  // namespace traits

// Chart is a map from T -> vector, retract is its inverse
// This is required by gtsam because our classes don't have the
// correct interface.
template<>
struct DefaultChart<kindr::minimal::QuatTransformation> {
  // the "type" typedef is used by gtsam.
  typedef kindr::minimal::QuatTransformation type;
  // the "vector" typedef is used by gtsam.
  typedef Eigen::Matrix<double, traits::dimension<type>::value, 1> vector;

  static vector local(const type& origin, const type& other) {
    return (other * origin.inverted()).log();
  }
  static type retract(const type& origin, const vector& d) {
    return type(d) * origin;
  }
  static int getDimension(const type& /* origin */) {
    return 6;
  }
};

typedef Eigen::Matrix<double, 3, 6> Jacobian3x6;
typedef Eigen::Matrix<double, 3, 3> Jacobian3x3;

inline Eigen::Vector3d transform_point(const kindr::minimal::QuatTransformation& T, const Eigen::Vector3d& p,
                                boost::optional<Jacobian3x6&> HT,
                                boost::optional<Jacobian3x3&> Hp) {
  if(HT) {
    // TODO(furgalep) fill in.
  }

  if(Hp) {
    // TODO(furgalep) fill in.
  }

  return T * p;
}

// This is syntatic sugar to be able to write
// Expression<Eigen::Vector3d> Tp = T * p;
// instead of
// Expression<Eigen::Vector3d> Tp = Expression<Eigen::Vector3d>(&transform_point, T, p);
inline gtsam::Expression<Eigen::Vector3d>
operator*(const gtsam::Expression<kindr::minimal::QuatTransformation>& T,
          const gtsam::Expression<Eigen::Vector3d>& p) {
  return Expression<Eigen::Vector3d>(&transform_point, T, p);
}

}  // namespace gtsam

#endif // MINKINDR_QUAT_TRANSFORMATION_GTSAM_H
