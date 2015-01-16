#ifndef KINDR_MIN_ROTATION_QUATERNION_INL_H_
#define KINDR_MIN_ROTATION_QUATERNION_INL_H_
#include <boost/math/special_functions/sinc.hpp>
#include <glog/logging.h>
#include <kindr/minimal/rotation-quaternion.h>
#include <kindr/minimal/angle-axis.h>

namespace kindr {
namespace minimal {

/// \brief initialize to identity
template<typename Scalar>
RotationQuaternionTemplate<Scalar>::RotationQuaternionTemplate() :
    q_A_B_(Implementation::Identity()) {
}

/// \brief initialize from real and imaginary components (real first)
template<typename Scalar>
RotationQuaternionTemplate<Scalar>::RotationQuaternionTemplate(
    Scalar w, Scalar x, Scalar y, Scalar z) :
    q_A_B_(w,x,y,z) {
  CHECK_NEAR(squaredNorm(), static_cast<Scalar>(1.0), static_cast<Scalar>(1e-4));
}

  
/// \brief initialize from real and imaginary components
template<typename Scalar>
RotationQuaternionTemplate<Scalar>::RotationQuaternionTemplate(
    Scalar real,
    const typename RotationQuaternionTemplate<Scalar>::Vector3& imaginary) :
    q_A_B_(real, imaginary[0], imaginary[1], imaginary[2]){
  CHECK_NEAR(squaredNorm(), static_cast<Scalar>(1.0),
             static_cast<Scalar>(1e-4));
}


/// \brief initialize from an Eigen quaternion
template<typename Scalar>
RotationQuaternionTemplate<Scalar>::RotationQuaternionTemplate(
    const Implementation& quaternion) :
    q_A_B_(quaternion){
  CHECK_NEAR(squaredNorm(), static_cast<Scalar>(1.0),
             static_cast<Scalar>(1e-4));
}

namespace detail {

template <typename Scalar_>
inline bool isLessThenEpsilons4thRoot(Scalar_ x){
  static const Scalar_ epsilon4thRoot = pow(std::numeric_limits<Scalar_>::epsilon(), 1.0/4.0);
  return x < epsilon4thRoot;
}

template <typename Scalar_>
inline Scalar_ arcSinXOverX(Scalar_ x) {
  if(isLessThenEpsilons4thRoot(fabs(x))){
    return Scalar_(1.0) + x * x * Scalar_(1.0/6.0);
  }
  return asin(x) / x;
}
}  // namespace detail

/// \brief initialize from axis-scaled angle vector
template<typename Scalar>
RotationQuaternionTemplate<Scalar>::RotationQuaternionTemplate(
    const Vector3& axis_scaled_angle) {
  *this = exp(axis_scaled_angle);
}

/// \brief initialize from a rotation matrix
template<typename Scalar>
RotationQuaternionTemplate<Scalar>::RotationQuaternionTemplate(
    const RotationMatrix& matrix) :
    q_A_B_(matrix) {
  // \todo furgalep check that this was a real rotation matrix
}


template<typename Scalar>
RotationQuaternionTemplate<Scalar>::RotationQuaternionTemplate(
    const AngleAxisTemplate<Scalar>& angleAxis) :
    q_A_B_(angleAxis.toImplementation()){

}


template<typename Scalar>
RotationQuaternionTemplate<Scalar>::~RotationQuaternionTemplate() {

}


/// \brief the real component of the quaternion
template<typename Scalar>
Scalar RotationQuaternionTemplate<Scalar>::w() const {
  return q_A_B_.w();
}

/// \brief the first imaginary component of the quaternion
template<typename Scalar>
Scalar RotationQuaternionTemplate<Scalar>::x() const {
  return q_A_B_.x();
}

/// \brief the second imaginary component of the quaternion
template<typename Scalar>
Scalar RotationQuaternionTemplate<Scalar>::y() const {
  return q_A_B_.y();
}

/// \brief the third imaginary component of the quaternion
template<typename Scalar>
Scalar RotationQuaternionTemplate<Scalar>::z() const {
  return q_A_B_.z();
}

/// \brief assignment operator
template<typename Scalar>
RotationQuaternionTemplate<Scalar>&
RotationQuaternionTemplate<Scalar>::operator=(
    const RotationQuaternionTemplate<Scalar>& rhs) {
  if(this != &rhs) {
    q_A_B_ = rhs.q_A_B_;
  }
  return *this;
}

/// \brief the imaginary components of the quaterion.
template<typename Scalar>
typename RotationQuaternionTemplate<Scalar>::Imaginary
RotationQuaternionTemplate<Scalar>::imaginary() const {
  return Imaginary(q_A_B_.x(),q_A_B_.y(),q_A_B_.z());
}

/// \brief get the components of the quaternion as a vector (real first)
template<typename Scalar>
typename RotationQuaternionTemplate<Scalar>::Vector4
RotationQuaternionTemplate<Scalar>::vector() const {
  return Vector4(q_A_B_.w(), q_A_B_.x(),q_A_B_.y(),q_A_B_.z());
}

/// \brief set the quaternion by its values (real, imaginary)
template<typename Scalar>
void RotationQuaternionTemplate<Scalar>::setValues(Scalar w, Scalar x,
                                                   Scalar y, Scalar z) {
  q_A_B_ = Implementation(w,x,y,z);
  CHECK_NEAR(squaredNorm(), static_cast<Scalar>(1.0),
             static_cast<Scalar>(1e-4));
}

/// \brief set the quaternion by its real and imaginary parts
template<typename Scalar>
void RotationQuaternionTemplate<Scalar>::setParts(Scalar real,
                                                  const Imaginary& imag) {
  q_A_B_ = Implementation(real, imag[0], imag[1], imag[2]);
  CHECK_NEAR(squaredNorm(), static_cast<Scalar>(1.0),
             static_cast<Scalar>(1e-4));
}


/// \brief get a copy of the representation that is unique
template<typename Scalar>
RotationQuaternionTemplate<Scalar>
RotationQuaternionTemplate<Scalar>::getUnique() const {
  if(this->w() > 0) {
    return *this;
  } else if (this->w() < 0){
    return RotationQuaternionTemplate<Scalar>(
        -this->w(),-this->x(),-this->y(),-this->z());
  } 
  // w == 0
  if(this->x() > 0) {
    return *this;
  } else if (this->x() < 0){
    return RotationQuaternionTemplate<Scalar>(
        -this->w(),-this->x(),-this->y(),-this->z());
  }
  // x == 0
  if(this->y() > 0) {
    return *this;
  } else if (this->y() < 0){
    return RotationQuaternionTemplate<Scalar>(
        -this->w(),-this->x(),-this->y(),-this->z());
  } 
  // y == 0
  if(this->z() > 0) { // z must be either -1 or 1 in this case
    return *this;
  } else {
    return RotationQuaternionTemplate<Scalar>(
        -this->w(),-this->x(),-this->y(),-this->z());
  }
}

/// \brief set the quaternion to its unique representation
template<typename Scalar>
RotationQuaternionTemplate<Scalar>&
RotationQuaternionTemplate<Scalar>::setUnique() {
  *this = getUnique();
  return *this;
}

/// \brief set the quaternion to identity
template<typename Scalar>
RotationQuaternionTemplate<Scalar>&
RotationQuaternionTemplate<Scalar>::setIdentity() {
  q_A_B_.setIdentity();
  return *this;
}

/// \brief set to random rotation
template<typename Scalar>
RotationQuaternionTemplate<Scalar>&
RotationQuaternionTemplate<Scalar>::setRandom() {
  Vector4 coeffs;
  coeffs.setRandom().normalize();
  q_A_B_ = Implementation(coeffs(0), coeffs(1), coeffs(2), coeffs(3));
  this->setUnique();
  return *this;
}

/// \brief set to random rotation  with a given angle
template<typename Scalar>
RotationQuaternionTemplate<Scalar>&
RotationQuaternionTemplate<Scalar>::setRandom(Scalar angle_rad) {
  Vector3 rotation_axis;
  rotation_axis.setRandom().normalize();
  q_A_B_ = Implementation(Eigen::AngleAxis<Scalar>(angle_rad, rotation_axis));
  return *this;
}

/// \brief get a copy of the quaternion inverted.
template<typename Scalar>
RotationQuaternionTemplate<Scalar>
RotationQuaternionTemplate<Scalar>::inverted() const {
  return conjugated();
}

/// \brief get a copy of the conjugate of the quaternion.
template<typename Scalar>
RotationQuaternionTemplate<Scalar>
RotationQuaternionTemplate<Scalar>::conjugated() const {
  // Own implementation since Eigen::conjugate does not use the correct
  // scalar type for the greater than zero comparison.
  return RotationQuaternionTemplate(
      Implementation(q_A_B_.w(),-q_A_B_.x(),-q_A_B_.y(),-q_A_B_.z()));
}


/// \brief rotate a vector, v
template<typename Scalar>
typename RotationQuaternionTemplate<Scalar>::Vector3
RotationQuaternionTemplate<Scalar>::rotate(
    const typename RotationQuaternionTemplate<Scalar>::Vector3& v) const {
  return q_A_B_*v;
}


/// \brief rotate a vector, v
template<typename Scalar>
typename RotationQuaternionTemplate<Scalar>::Vector4
RotationQuaternionTemplate<Scalar>::rotate4(
    const typename RotationQuaternionTemplate<Scalar>::Vector4& v) const {
  typename RotationQuaternionTemplate<Scalar>::Vector4 vprime;
  vprime[3] = v[3];
  vprime.template head<3>() = q_A_B_*v.template head<3>();
  return vprime;
}


/// \brief rotate a vector, v
template<typename Scalar>
typename RotationQuaternionTemplate<Scalar>::Vector3
RotationQuaternionTemplate<Scalar>::inverseRotate(
    const typename RotationQuaternionTemplate<Scalar>::Vector3& v) const {
  return q_A_B_.inverse()*v;
}


/// \brief rotate a vector, v
template<typename Scalar>
typename RotationQuaternionTemplate<Scalar>::Vector4
RotationQuaternionTemplate<Scalar>::inverseRotate4(
    const typename RotationQuaternionTemplate<Scalar>::Vector4& v) const {
  typename RotationQuaternionTemplate<Scalar>::Vector4 vprime;
  vprime[3] = v[3];
  vprime.template head<3>() = q_A_B_.inverse()*v.template head<3>();
  return vprime;
}


/// \brief cast to the implementation type
template<typename Scalar>
typename  RotationQuaternionTemplate<Scalar>::Implementation&
RotationQuaternionTemplate<Scalar>::toImplementation() {
  return q_A_B_;
}

/// \brief cast to the implementation type
template<typename Scalar>
const typename RotationQuaternionTemplate<Scalar>::Implementation&
RotationQuaternionTemplate<Scalar>::toImplementation() const {
  return q_A_B_;
}

/// \brief get the norm of the quaternion
template<typename Scalar>
Scalar RotationQuaternionTemplate<Scalar>::norm() const {
  return q_A_B_.norm();
}

/// \brief get the squared norm of the quaternion
template<typename Scalar>
Scalar RotationQuaternionTemplate<Scalar>::squaredNorm() const {
  return q_A_B_.squaredNorm();
}

/// \brief enforce the unit length constraint
template<typename Scalar>
RotationQuaternionTemplate<Scalar>&
RotationQuaternionTemplate<Scalar>::normalize() {
  q_A_B_.normalize();
  return *this;
}

template<typename Scalar>
RotationQuaternionTemplate<Scalar>
RotationQuaternionTemplate<Scalar>::operator*(
    const RotationQuaternionTemplate<Scalar>& rhs) const {
  return RotationQuaternionTemplate<Scalar>(q_A_B_ * rhs.q_A_B_);
}

template<typename Scalar>
RotationQuaternionTemplate<Scalar>
RotationQuaternionTemplate<Scalar>::operator*(
    const AngleAxisTemplate<Scalar>& rhs) const {
  return RotationQuaternionTemplate<Scalar>(q_A_B_ *
                                            rhs.toImplementation());
}

template<typename Scalar>
std::ostream& operator<<(std::ostream& out,
                         const RotationQuaternionTemplate<Scalar>& rhs) {
  out << rhs.vector();
  return out;
}

template<typename Scalar>
typename RotationQuaternionTemplate<Scalar>::RotationMatrix
RotationQuaternionTemplate<Scalar>::getRotationMatrix() const {
  return q_A_B_.matrix();
}

template<typename Scalar>
Scalar RotationQuaternionTemplate<Scalar>::getDisparityAngle(
    const RotationQuaternionTemplate<Scalar>& rhs) const{
  return AngleAxis(rhs * this->inverted()).getUnique().angle();
}

template<typename Scalar>
Scalar RotationQuaternionTemplate<Scalar>::getDisparityAngle(
    const AngleAxisTemplate<Scalar>& rhs) const{
  return AngleAxis(rhs * this->inverted()).getUnique().angle();
}

template<typename Scalar>
typename RotationQuaternionTemplate<Scalar>::Vector3
RotationQuaternionTemplate<Scalar>::log(const RotationQuaternionTemplate<Scalar>& q) {
  const Eigen::Matrix<Scalar, 3, 1> a = q.imaginary();
  const Scalar na = a.norm();
  const Scalar eta = q.w();
  Scalar scale;
  if(fabs(eta) < na){ // use eta because it is more precise than na to calculate the scale. No singularities here.
    if (eta >= 0) {
      scale = acos(eta) / na;
    } else {
      scale = -acos(-eta) / na;
    }
  } else {
    /*
     * In this case more precision is in na than in eta so lets use na only to calculate the scale:
     *
     * assume first eta > 0 and 1 > na > 0.
     *               u = asin (na) / na  (this implies u in [1, pi/2], because na i in [0, 1]
     *    sin (u * na) = na
     *  sin^2 (u * na) = na^2
     *  cos^2 (u * na) = 1 - na^2
     *                              (1 = ||q|| = eta^2 + na^2)
     *    cos^2 (u * na) = eta^2
     *                              (eta > 0,  u * na = asin(na) in [0, pi/2] => cos(u * na) >= 0 )
     *      cos (u * na) = eta
     *                              (u * na in [ 0, pi/2] )
     *                 u = acos (eta) / na
     *
     * So the for eta > 0 it is acos(eta) / na == asin(na) / na.
     * From some geometric considerations (mirror the setting at the hyper plane q==0) it follows for eta < 0 that (pi - asin(na)) / na = acos(eta) / na.
     */
    if(eta > 0) {
      // For asin(na)/ na the singularity na == 0 can be removed. We can ask (e.g. Wolfram alpha) for its series expansion at na = 0. And that is done in the following function.
      scale = detail::arcSinXOverX(na);
    } else {
      // (pi - asin(na))/ na has a pole at na == 0. So we cannot remove this singularity.
      // It is just the cut locus of the unit quaternion manifold at identity and thus the axis angle description becomes necessarily unstable there.
      //scale = (M_PI - asin(na)) / na;
      scale = -detail::arcSinXOverX(na);
    }
  }
  return a * (Scalar(2.0) * scale);
}

template<typename Scalar>
RotationQuaternionTemplate<Scalar>
RotationQuaternionTemplate<Scalar>::exp(const Vector3& dx) {
  // Method of implementing this function that is accurate to numerical precision from
  // Grassia, F. S. (1998). Practical parameterization of rotations using the exponential map. journal of graphics, gpu, and game tools, 3(3):29–48.
  double theta = dx.norm();
  // na is 1/theta sin(theta/2)
  double na;
  if(detail::isLessThenEpsilons4thRoot(theta)){
    static const double one_over_48 = 1.0/48.0;
    na = 0.5 + (theta * theta) * one_over_48;
  } else {
    na = sin(theta*0.5) / theta;
  }
  double ct = cos(theta*0.5);
  return RotationQuaternionTemplate<Scalar>(ct,
                                            dx[0]*na,
                                            dx[1]*na,
                                            dx[2]*na);
}

template<typename Scalar>
typename RotationQuaternionTemplate<Scalar>::Vector3
RotationQuaternionTemplate<Scalar>::log() const {
  return log(*this);
}

} // namespace minimal
} // namespace kindr
#endif  // KINDR_MIN_ROTATION_QUATERNION_INL_H_
