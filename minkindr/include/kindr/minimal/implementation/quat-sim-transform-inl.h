#ifndef KINDR_MINIMAL_IMPLEMENTATION_QUAT_SIM_TRANSFORM_INL_H_
#define KINDR_MINIMAL_IMPLEMENTATION_QUAT_SIM_TRANSFORM_INL_H_

namespace kindr {
namespace minimal {

template <typename Scalar>
QuatSimTransformTemplate<Scalar>::QuatSimTransformTemplate() : scale_A_B_(1.) {}

template <typename Scalar>
QuatSimTransformTemplate<Scalar>::QuatSimTransformTemplate(
    const Transform& T_A_B, const Scalar scale_A_B_)
    : T_A_B_(T_A_B), scale_A_B_(scale_A_B_) {}

template <typename Scalar>
typename QuatSimTransformTemplate<Scalar>::Vectors
QuatSimTransformTemplate<Scalar>::operator*(const Vectors& rhs) const {
  return T_A_B_.transformVectorized(scale_A_B_ * rhs);
}

template <typename Scalar>
QuatSimTransformTemplate<Scalar>
QuatSimTransformTemplate<Scalar>::inverse() const {
  return QuatSimTransformTemplate<Scalar>(
      QuatTransformationTemplate<Scalar>(
          T_A_B_.getRotation().inverse(), -T_A_B_.getRotation().inverseRotate(
              T_A_B_.getPosition() / scale_A_B_)), 1. / scale_A_B_);
}

template<typename Scalar>
std::ostream & operator<<(std::ostream & out,
                          const QuatSimTransformTemplate<Scalar>& sim_3) {
  out << "Transform:" << std::endl << sim_3.T_A_B_.getTransformationMatrix() <<
      std::endl;
  out << "Scale: " << sim_3.scale_A_B_;
  return out;
}

}  // namespace minimal
}  // namespace kindr

#endif  // KINDR_MINIMAL_IMPLEMENTATION_QUAT_SIM_TRANSFORM_INL_H_
