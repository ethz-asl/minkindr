#ifndef MINKINDR_CUBIC_HERMITE_TRANSFORMATION_GTSAM_H
#define MINKINDR_CUBIC_HERMITE_TRANSFORMATION_GTSAM_H

#include <gtsam/base/Manifold.h>
#include <gtsam/nonlinear/Expression.h>

#include <kindr/minimal/common.h>

#include "common-gtsam.h"


namespace kindr {
namespace minimal {
typedef Eigen::Matrix<double, 6, 1> Vector6d;
////////////////////////////////////////////////////////////////////////////////
// Convenience functions to make working with expressions easy and fun!

/// cubic hermite interpolation on the unit interval for TRANSLATIONS
/// param[in] t_W_A, expression for the translation at the beginning of the interval
/// param[in] t_W_B, expression for the translation at the end of the interval
/// param[in] W_v_W_A, expression for the velocity at the beginning of the interval (in world coordinate frame)
/// param[in] W_v_W_B, expression for the velocity at the end of the interval (in world coordinate frame)
/// param[in] alpha, the interpolation coefficient [0 .. 1]
/// param[in] dt, the time difference between A and B
EVector3 hermiteTranslationInterpolation(const EVector3& t_W_A, const EVector3& W_v_W_A,
                                         const EVector3& t_W_B, const EVector3& W_v_W_B,
                                         double alpha, double dt) {

  /// Equations for the unit interval:
  // Let t_W_A, t_W_B denote the control point values (=translations) and W_v_W_A, W_v_W_B
  // the control derivatives (=velocities in R³) at the interval's boundaries.
  // Then (b == beta):
  //
  // p0 = t_W_A, p1 = t_W_B, p2 = W_v_W_A, p3 = W_v_W_B
  // b0 = 2t³-3t²+1, b_1 = -2t³+3t², b_2 = t³-2t²+t, b_3 = t³-t²
  // Spline equation:
  // p(t) = p0 * b0 + p1 * b1 + p2 * b2 + p3 + b3
  double alpha2 = alpha * alpha;
  double alpha3 = alpha2 * alpha;

  double beta0 = 2.0 * alpha3 - 3.0 * alpha2 + 1.0;
  double beta1 = -2.0 * alpha3 + 3.0 * alpha2;
  double beta2 = alpha3 - 2.0 * alpha2 + alpha;
  double beta3 = alpha3 - alpha2;

  return t_W_A * beta0 + t_W_B * beta1 + W_v_W_A * beta2 * dt + W_v_W_B * beta3 * dt;
}

/// derivative of cubic hermite interpolation on the unit interval for translations
/// param[in] t_W_A, expression for the translation at the beginning of the interval
/// param[in] t_W_B, expression for the translation at the end of the interval
/// param[in] W_v_W_A, expression for the velocity at the beginning of the interval (in world coordinate frame)
/// param[in] W_v_W_B, expression for the velocity at the end of the interval (in world coordinate frame)
/// param[in] alpha, the interpolation coefficient [0 .. 1]
/// param[in] dt, the time difference between A and B
/// output is velocity in world coordinate frame.
EVector3 hermiteTranslationInterpolationDerivative(const EVector3& t_W_A, const EVector3& W_v_W_A,
                                                   const EVector3& t_W_B, const EVector3& W_v_W_B,
                                                   double alpha, double dt)  {
  // In order to obtain the spline's derivative apply the chain rule. Pro memoria the spline equation:
  // p(t) = p_0 * b_0 + p_1 * b_1 + p_2 * b_2 + p_3 * b_3
  // Thus, the derivatives of b_0, b_1, b_2, b_3 have to be calculated.

  double alpha2 = alpha * alpha;

  // derivatives of b_0, b_1, b_2, b_3
  double dotBeta0 = 6.0 * alpha2 - 6.0 * alpha;
  double dotBeta1 = -6.0 * alpha2 + 6.0 * alpha;
  double dotBeta2 = 3.0 * alpha2 - 4.0 * alpha + 1.0;
  double dotBeta3 = 3.0 * alpha2 - 2.0 * alpha;

  return t_W_A * dotBeta0 + t_W_B * dotBeta1 + W_v_W_A * dotBeta2 * dt + W_v_W_B * dotBeta3 * dt;
}

}  // namespace minimal
}  // namespace kindr

#endif // MINKINDR_CUBIC_HERMITE_TRANSFORMATION_GTSAM_H
