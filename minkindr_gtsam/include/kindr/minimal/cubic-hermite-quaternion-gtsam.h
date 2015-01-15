#ifndef MINKINDR_QUATERNION_CUBIC_HERMITE_GTSAM_H
#define MINKINDR_QUATERNION_CUBIC_HERMITE_GTSAM_H

#include <gtsam/base/Manifold.h>
#include <gtsam/nonlinear/Expression.h>

#include <kindr/minimal/rotation-quaternion.h>
#include <kindr/minimal/common.h>

#include "common-gtsam.h"

namespace kindr {
namespace minimal {
////////////////////////////////////////////////////////////////////////////////
// Convenience functions to make working with expressions easy and fun!
typedef Expression<RotationQuaternion> EQuaternion;
typedef gtsam::Expression<Eigen::Vector3d> EVector3;

/// cubic hermite interpolation on the unit interval for quaternions
/// param[in] qa, expression for the quaternion at the begining of the interval
/// param[in] qb, expression for the quaternion at the end of the interval
/// param[in] wa, expression for the rotational velocity at the begining of the interval (in world coordinate frame)
/// param[in] wb, expression for the rotational velocity at the end of the interval (in world coordinate frame)
EQuaternion hermite_interpolation(const EQuaternion& qa, const EQuaternion& qb, const EVector3& wa, const EVector3& wb, double alpha) {

  /// Equations for the unit interval:
    // Let qa, qb denote the control point values (=unit quaternions) and va, vb
    // the control derivatives (=angular speeds in R³) at the interval's boundaries.
    // Then (w == omega, b == beta):
    // w_1 = va / 3
    // w_2 = log[ exp(w_1)^{-1} * qa^{-1} * qb * exp(w_3)^{-1} ]
    // w_3 = vb / 3
    // b_1 = t³-3t²+3t, b_2 = -2t³+3t², b_3 = t³
    // Spline equation:
    // q(t) = p_1 * exp(w_1*b_1) * exp(w_2*b_2) * exp(w_3*b_3)

  EVector3 w1 = inverseRotate(qa, wa * (1 / 3.0));
  EVector3 w3 = inverseRotate(qb, wb * (1 / 3.0));
  EVector3 w2 = log( exp(w1 * -1.0) * invert(qa) * qb * exp(w3 * -1.0));
  double alpha2 = alpha * alpha;
  double alpha3 = alpha2 * alpha;
  double beta1 = alpha3 - 3 * alpha2 + 3 * alpha;
  double beta2 = -2.0 * alpha3 + 3 * alpha2;
  double beta3 = alpha3;

  return qa * exp(w1 * beta1) * exp(w2 * beta2) * exp(w3 * beta3);
}

/// derivative of cubic hermite interpolation on the unit interval for quaternions
/// param[in] qa, expression for the quaternion at the begining of the interval
/// param[in] qb, expression for the quaternion at the end of the interval
/// param[in] wa, expression for the rotational velocity at the begining of the interval (in world coordinate frame)
/// param[in] wb, expression for the rotational velocity at the end of the interval (in world coordinate frame)
/// output is rotational velocity in world coordinate frame.
EVector3 hermite_interpolation_derivative(const EQuaternion& qa, const EQuaternion& qb, const EVector3& wa, const EVector3& wb, double alpha)  {
  // In order to obtain the spline's derivative apply the chain rule. Pro memoria the spline equation:
  // q(t) = p_1 * exp(w_1*b_1) * exp(w_2*b_2) * exp(w_3*b_3)
  // Thus, the derivatives of b_1, b_2, b_3 have to be calculated.

  EVector3 w1 = inverseRotate(qa, wa * (1 / 3.0));
  EVector3 w3 = inverseRotate(qb, wb * (1 / 3.0));
  EVector3 w2 = log( exp(w1 * -1.0) * invert(qa) * qb * exp(w3 * -1.0));

  double alpha2 = alpha * alpha;
  double alpha3 = alpha2 * alpha;
  double beta1 = alpha3 - 3.0 * alpha2 + 3.0 * alpha;
  double beta2 = -2.0 * alpha3 + 3.0 * alpha2;
  double beta3 = alpha3;

  // derivatives of b_1, b_2, b_3
  double dotBeta1 = 3.0 * alpha2 - 6.0 * alpha + 3.0;
  double dotBeta2 = -6.0 * alpha2 + 6.0 * alpha;
  double dotBeta3 = 3.0 * alpha2;

  EQuaternion qa_exp_w1beta1 = qa * exp(w1 * beta1);

  return
      qa * w1 * dotBeta1 +
      qa_exp_w1beta1 * w2 * dotBeta2 +
      qa_exp_w1beta1 * exp(w2 * beta2) * w3 * dotBeta3;
}

}  // namespace minimal
}  // namespace kindr

#endif // MINKINDR_QUATERNION_CUBIC_HERMITE_GTSAM_H
