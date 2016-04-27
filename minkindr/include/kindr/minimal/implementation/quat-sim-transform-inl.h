#ifndef KINDR_MINIMAL_IMPLEMENTATION_QUAT_SIM_TRANSFORM_INL_H_
#define KINDR_MINIMAL_IMPLEMENTATION_QUAT_SIM_TRANSFORM_INL_H_

namespace kindr {
namespace minimal {

template <typename Scalar>
QuatSimTransformTemplate<Scalar>::QuatSimTransformTemplate() : scale_A_B_(1) {}

template <typename Scalar>
QuatSimTransformTemplate<Scalar>::QuatSimTransformTemplate(
    const Transform& T_A_B, const Scalar scale_A_B_)
    : T_A_B_(T_A_B), scale_A_B_(scale_A_B_) {}

template <typename Scalar>
QuatSimTransformTemplate<Scalar>::Vectors
QuatSimTransformTemplate<Scalar>::operator*(const Vectors& rhs) const {
  return T_A_B_.transformVectorized(rhs) * scale_A_B_;
}

}  // namespace minimal
}  // namespace kindr

#endif  // KINDR_MINIMAL_IMPLEMENTATION_QUAT_SIM_TRANSFORM_INL_H_
