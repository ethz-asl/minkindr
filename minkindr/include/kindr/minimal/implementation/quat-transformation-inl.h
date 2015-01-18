#ifndef KINDR_MINIMAL_QUAT_TRANSFORMATION_H_INL_
#define KINDR_MINIMAL_QUAT_TRANSFORMATION_H_INL_
#include <kindr/minimal/quat-transformation.h>

namespace kindr {
namespace minimal {

template<typename Scalar>
QuatTransformationTemplate<Scalar>::QuatTransformationTemplate() {
  setIdentity();
}

template<typename Scalar>
QuatTransformationTemplate<Scalar>::QuatTransformationTemplate(
    const RotationQuaternionTemplate<Scalar>& q_A_B, const Position& A_t_A_B) :
    q_A_B_(q_A_B),
    A_t_A_B_(A_t_A_B) {

}

template<typename Scalar>
QuatTransformationTemplate<Scalar>::QuatTransformationTemplate(
    const typename Rotation::Implementation& q_A_B,
    const Position& A_t_A_B) :
        q_A_B_(q_A_B),
        A_t_A_B_(A_t_A_B) {

}

template<typename Scalar>
QuatTransformationTemplate<Scalar>::QuatTransformationTemplate(
    const Position& A_t_A_B, const RotationQuaternionTemplate<Scalar>& q_A_B) :
    q_A_B_(q_A_B),
    A_t_A_B_(A_t_A_B) {

}

template<typename Scalar>
QuatTransformationTemplate<Scalar>::QuatTransformationTemplate(
    const Position& A_t_A_B, const typename Rotation::Implementation& q_A_B) :
        q_A_B_(q_A_B),
        A_t_A_B_(A_t_A_B) {

}

template<typename Scalar>
QuatTransformationTemplate<Scalar>::QuatTransformationTemplate(
    const TransformationMatrix& T) :
    q_A_B_(T.template topLeftCorner<3,3>().eval()),
    A_t_A_B_(T.template topRightCorner<3,1>().eval()) {
}

template<typename Scalar>
QuatTransformationTemplate<Scalar>::QuatTransformationTemplate(
     const QuatTransformationTemplate<Scalar>::Vector6& x_t_r) : 
  q_A_B_(x_t_r.template tail<3>().eval()),
  A_t_A_B_(x_t_r.template head<3>().eval()) {
}
                                                              
template<typename Scalar>
QuatTransformationTemplate<Scalar>::~QuatTransformationTemplate() {

}

template<typename Scalar>
void QuatTransformationTemplate<Scalar>::setIdentity() {
  q_A_B_.setIdentity();
  A_t_A_B_.setZero();
}

template<typename Scalar>
typename QuatTransformationTemplate<Scalar>::Position&
QuatTransformationTemplate<Scalar>::getPosition() {
  return A_t_A_B_;
}

template<typename Scalar>
const typename QuatTransformationTemplate<Scalar>::Position&
QuatTransformationTemplate<Scalar>::getPosition() const {
  return A_t_A_B_;
}

template<typename Scalar>
typename QuatTransformationTemplate<Scalar>::Rotation&
QuatTransformationTemplate<Scalar>::getRotation() {
  return q_A_B_;
}

template<typename Scalar>
const typename QuatTransformationTemplate<Scalar>::Rotation&
QuatTransformationTemplate<Scalar>::getRotation() const {
  return q_A_B_;
}

template<typename Scalar>
typename QuatTransformationTemplate<Scalar>::TransformationMatrix
QuatTransformationTemplate<Scalar>::getTransformationMatrix() const {
  TransformationMatrix transformation_matrix;
  transformation_matrix.setIdentity();
  transformation_matrix.template topLeftCorner<3,3>() =
      q_A_B_.getRotationMatrix();
  transformation_matrix.template topRightCorner<3,1>() = A_t_A_B_;
  return transformation_matrix;
}

template<typename Scalar>
Eigen::Matrix<Scalar, 7, 1>
QuatTransformationTemplate<Scalar>::asVector() const {
  return (Eigen::Matrix<Scalar, 7, 1>() <<
      q_A_B_.w(), q_A_B_.x(), q_A_B_.y(), q_A_B_.z(), A_t_A_B_).finished();
}


template<typename Scalar>
typename QuatTransformationTemplate<Scalar>::RotationMatrix
QuatTransformationTemplate<Scalar>::getRotationMatrix() const {
  return q_A_B_.getRotationMatrix();
}

template<typename Scalar>
QuatTransformationTemplate<Scalar>
QuatTransformationTemplate<Scalar>::operator*(
    const QuatTransformationTemplate<Scalar>& rhs) const {
  return QuatTransformationTemplate<Scalar>(q_A_B_ * rhs.q_A_B_, A_t_A_B_ +
                                            q_A_B_.rotate(rhs.A_t_A_B_));
}

template<typename Scalar>
typename QuatTransformationTemplate<Scalar>::Vector3
QuatTransformationTemplate<Scalar>::transform(
    const typename QuatTransformationTemplate<Scalar>::Vector3& rhs) const {
  return q_A_B_.rotate(rhs) + A_t_A_B_;
}

template<typename Scalar>
typename QuatTransformationTemplate<Scalar>::Vector3
QuatTransformationTemplate<Scalar>::operator*(
    const typename QuatTransformationTemplate<Scalar>::Vector3& rhs) const {
  return transform(rhs);
}

template<typename Scalar>
typename QuatTransformationTemplate<Scalar>::Vector4
QuatTransformationTemplate<Scalar>::transform4(
    const typename QuatTransformationTemplate<Scalar>::Vector4& rhs) const {
  QuatTransformationTemplate<Scalar>::Vector4 rval;
  rval[3] = rhs[3];
  rval.template head<3>() =
      q_A_B_.rotate(rhs.template head<3>()) + rhs[3]*A_t_A_B_;
  return rval;
}

template<typename Scalar>
typename QuatTransformationTemplate<Scalar>::Vector3
QuatTransformationTemplate<Scalar>::inverseTransform(
    const Vector3& rhs) const {
  return q_A_B_.inverseRotate(rhs - A_t_A_B_);
}

template<typename Scalar>
typename QuatTransformationTemplate<Scalar>::Vector4
QuatTransformationTemplate<Scalar>::inverseTransform4(
    const typename QuatTransformationTemplate<Scalar>::Vector4& rhs) const {
  typename QuatTransformationTemplate<Scalar>::Vector4 rval;
  rval.template head<3>() = q_A_B_.inverseRotate(rhs.template head<3>() -
                                                 A_t_A_B_*rhs[3]);
  rval[3] = rhs[3];
  return rval;
}

template<typename Scalar>
QuatTransformationTemplate<Scalar>
QuatTransformationTemplate<Scalar>::inverted() const {
  return QuatTransformation(q_A_B_.inverted(), -q_A_B_.inverseRotate(A_t_A_B_));
}

template<typename Scalar>
typename QuatTransformationTemplate<Scalar>::Vector6 
QuatTransformationTemplate<Scalar>::log() const {
  AngleAxisTemplate<Scalar> angleaxis(q_A_B_);
  return (Vector6() << A_t_A_B_, (angleaxis.axis()*angleaxis.angle())).finished();
}

template<typename Scalar>
QuatTransformationTemplate<Scalar>&
QuatTransformationTemplate<Scalar>::setRandom() {
  q_A_B_.setRandom();
  A_t_A_B_.setRandom();
  return *this;
}

template<typename Scalar>
QuatTransformationTemplate<Scalar>&
QuatTransformationTemplate<Scalar>::setRandom(Scalar norm_translation) {
  setRandom();
  A_t_A_B_.normalize();
  A_t_A_B_ *= norm_translation;
  return *this;
}

template<typename Scalar>
QuatTransformationTemplate<Scalar>&
QuatTransformationTemplate<Scalar>::setRandom(Scalar norm_translation, Scalar angle_rad) {
  q_A_B_.setRandom(angle_rad);
  A_t_A_B_.setRandom().normalize();
  A_t_A_B_ *= norm_translation;
  return *this;
}

template<typename Scalar>
std::ostream & operator<<(std::ostream & out,
                          const QuatTransformationTemplate<Scalar>& pose) {
  out << pose.getTransformationMatrix();
  return out;
}

template<typename Scalar>
bool QuatTransformationTemplate<Scalar>::operator==(
    const QuatTransformationTemplate<Scalar>& rhs) const {
  return q_A_B_ == rhs.q_A_B_ && A_t_A_B_ == rhs.A_t_A_B_;
}

} // namespace minimal
} // namespace kindr
#endif  // KINDR_MINIMAL_QUAT_TRANSFORMATION_H_INL_
