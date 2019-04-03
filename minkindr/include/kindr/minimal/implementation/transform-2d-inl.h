#include <limits>

#include <glog/logging.h>

namespace kindr {
namespace minimal {

template <typename Scalar>
Transformation2DTemplate<Scalar>::Transformation2DTemplate() {
  setIdentity();
}

template <typename Scalar>
Transformation2DTemplate<Scalar>::Transformation2DTemplate(
    const Rotation& r_A_B, const Position& A_t_A_B)
    : r_A_B_(r_A_B), A_t_A_B_(A_t_A_B) {}

template <typename Scalar>
Transformation2DTemplate<Scalar>::Transformation2DTemplate(
    const Position& A_t_A_B, const Rotation& r_A_B)
    : Transformation2DTemplate<Scalar>(r_A_B, A_t_A_B) {}

template <typename Scalar>
Transformation2DTemplate<Scalar>::Transformation2DTemplate(
    const TransformationMatrix& T)
    : Transformation2DTemplate<Scalar>(
          Rotation2D(T.template topLeftCorner<2, 2>().eval()),
          T.template topRightCorner<2, 1>().eval()) {
  constexpr Scalar kScalarOne = static_cast<Scalar>(1.0);
  constexpr Scalar kEpsilon = std::numeric_limits<Scalar>::epsilon();
  CHECK_LE((T(2, 2) - kScalarOne), kEpsilon);
  const Eigen::Matrix<Scalar, 2, 2> rotation_matrix =
      T.template topLeftCorner<2, 2>().eval();
  CHECK_LE(rotation_matrix.determinant() - kScalarOne, kEpsilon);
}

template <typename Scalar>
void Transformation2DTemplate<Scalar>::setIdentity() {
  r_A_B_ = Rotation2D(static_cast<Scalar>(0.0));
  A_t_A_B_.setZero();
}

template <typename Scalar>
typename Transformation2DTemplate<Scalar>::Position&
Transformation2DTemplate<Scalar>::getPosition() {
  return A_t_A_B_;
}

template <typename Scalar>
const typename Transformation2DTemplate<Scalar>::Position&
Transformation2DTemplate<Scalar>::getPosition() const {
  return A_t_A_B_;
}

template <typename Scalar>
typename Transformation2DTemplate<Scalar>::Rotation&
Transformation2DTemplate<Scalar>::getRotation() {
  return r_A_B_;
}

template <typename Scalar>
const typename Transformation2DTemplate<Scalar>::Rotation&
Transformation2DTemplate<Scalar>::getRotation() const {
  return r_A_B_;
}

template <typename Scalar>
typename Transformation2DTemplate<Scalar>::RotationMatrix
Transformation2DTemplate<Scalar>::getRotationMatrix() const {
  return r_A_B_.toRotationMatrix();
}

template <typename Scalar>
typename Transformation2DTemplate<Scalar>::TransformationMatrix
Transformation2DTemplate<Scalar>::getTransformationMatrix() const {
  TransformationMatrix transformation_matrix;
  transformation_matrix.setIdentity();
  transformation_matrix.template topLeftCorner<2, 2>() =
      r_A_B_.toRotationMatrix();
  transformation_matrix.template topRightCorner<2, 1>() = A_t_A_B_;
  return transformation_matrix;
}

template <typename Scalar>
Eigen::Matrix<Scalar, 3, 1> Transformation2DTemplate<Scalar>::asVector() const {
  return (Eigen::Matrix<Scalar, 3, 1>() << r_A_B_.angle(), A_t_A_B_).finished();
}

template <typename Scalar>
Transformation2DTemplate<Scalar> Transformation2DTemplate<Scalar>::operator*(
    const Transformation2DTemplate<Scalar>& rhs) const {
  return Transformation2DTemplate<Scalar>(
      r_A_B_ * rhs.r_A_B_, A_t_A_B_ + r_A_B_ * rhs.A_t_A_B_);
}

template <typename Scalar>
typename Transformation2DTemplate<Scalar>::Vector2
    Transformation2DTemplate<Scalar>::operator*(const Vector2& rhs) const {
  return transform(rhs);
}

template <typename Scalar>
typename Transformation2DTemplate<Scalar>::Vector2
Transformation2DTemplate<Scalar>::transform(const Vector2& rhs) const {
  return r_A_B_ * rhs + A_t_A_B_;
}

template <typename Scalar>
typename Transformation2DTemplate<Scalar>::Matrix2X
Transformation2DTemplate<Scalar>::transformVectorized(
    const Matrix2X& rhs) const {
  CHECK_GT(rhs.cols(), 0);
  return (r_A_B_.toRotationMatrix() * rhs).colwise() + A_t_A_B_;
}

template <typename Scalar>
Transformation2DTemplate<Scalar> Transformation2DTemplate<Scalar>::inverse()
    const {
  return Transformation2DTemplate<Scalar>(
      r_A_B_.inverse(), -(r_A_B_.inverse() * A_t_A_B_));
}

template <typename Scalar>
bool Transformation2DTemplate<Scalar>::operator==(
    const Transformation2DTemplate<Scalar>& rhs) const {
  return r_A_B_.angle() == rhs.r_A_B_.angle() && A_t_A_B_ == rhs.A_t_A_B_;
}

template <typename Scalar>
bool Transformation2DTemplate<Scalar>::operator!=(
    const Transformation2DTemplate<Scalar>& rhs) const {
  return !(*this == rhs);
}

template <typename Scalar>
template <typename ScalarAfterCast>
Transformation2DTemplate<ScalarAfterCast>
Transformation2DTemplate<Scalar>::cast() const {
  return Transformation2DTemplate<ScalarAfterCast>(
      getRotation().template cast<ScalarAfterCast>(),
      getPosition().template cast<ScalarAfterCast>());
}

}  // namespace minimal
}  // namespace kindr
