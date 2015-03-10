// Copyright (c) 2015, Autonomous Systems Lab, ETH Zurich
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the <organization> nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#ifndef MINKINDR_QUATERNION_CUBIC_HERMITE_GTSAM_H
#define MINKINDR_QUATERNION_CUBIC_HERMITE_GTSAM_H

#include <gtsam/base/Manifold.h>
#include <gtsam/nonlinear/Expression.h>

#include <kindr/minimal/rotation-quaternion.h>
#include <kindr/minimal/common.h>

#include "common-gtsam.h"
#include "kindr/minimal/rotation-quaternion-gtsam.h"


namespace kindr {
namespace minimal {
////////////////////////////////////////////////////////////////////////////////
// Convenience functions to make working with expressions easy and fun!

/// cubic hermite interpolation on the unit interval for quaternions
/// param[in] quat_W_A, expression for the quaternion at the begining of the interval
/// param[in] quat_W_B, expression for the quaternion at the end of the interval
/// param[in] W_omega_W_A, expression for the rotational velocity at the begining of the interval (in world coordinate frame)
/// param[in] W_omega_W_B, expression for the rotational velocity at the end of the interval (in world coordinate frame)
/// param[in] alpha, the interpolation coefficient [0 .. 1]
EQuaternion hermiteInterpolation(const EQuaternion& quat_W_A, const EVector3& W_omega_W_A,
                                  const EQuaternion& quat_W_B, const EVector3& W_omega_W_B,
                                  double alpha) {

  /// Equations for the unit interval:
    // Let quat_W_A, quat_W_B denote the control point values (=unit quaternions) and va, vb
    // the control derivatives (=angular speeds in R³) at the interval's boundaries.
    // Then (w == omega, b == beta):
    // w_1 = va / 3
    // w_2 = log[ exp(w_1)^{-1} * quat_W_A^{-1} * quat_W_B * exp(w_3)^{-1} ]
    // w_3 = vb / 3
    // b_1 = t³-3t²+3t, b_2 = -2t³+3t², b_3 = t³
    // Spline equation:
    // q(t) = p_1 * exp(w_1*b_1) * exp(w_2*b_2) * exp(w_3*b_3)
  const double one_third = 1.0 / 3.0;
  EVector3 w1 = inverseRotate(quat_W_A, W_omega_W_A * one_third);
  EVector3 w3 = inverseRotate(quat_W_B, W_omega_W_B * one_third);
  EVector3 w2 = quaternionLog( quaternionExp(-w1) * invert(quat_W_A) * quat_W_B * quaternionExp(-w3));
  double alpha2 = alpha * alpha;
  double alpha3 = alpha2 * alpha;
  double beta1 = alpha3 - 3.0 * alpha2 + 3.0 * alpha;
  double beta2 = -2.0 * alpha3 + 3.0 * alpha2;
  double beta3 = alpha3;

  return quat_W_A * quaternionExp(w1 * beta1) * quaternionExp(w2 * beta2) * quaternionExp(w3 * beta3);
}

/// derivative of cubic hermite interpolation on the unit interval for quaternions
/// param[in] quat_W_A, expression for the quaternion at the begining of the interval
/// param[in] quat_W_B, expression for the quaternion at the end of the interval
/// param[in] W_omega_W_A, expression for the rotational velocity at the begining of the interval (in world coordinate frame)
/// param[in] W_omega_W_B, expression for the rotational velocity at the end of the interval (in world coordinate frame)
/// param[in] alpha, the interpolation coefficient [0 .. 1]
/// output is rotational velocity in world coordinate frame.
EVector3 hermiteInterpolationDerivative(const EQuaternion& quat_W_A, const EVector3& W_omega_W_A,
                                          const EQuaternion& quat_W_B, const EVector3& W_omega_W_B,
                                          double alpha)  {
  // In order to obtain the spline's derivative apply the chain rule. Pro memoria the spline equation:
  // q(t) = p_1 * exp(w_1*b_1) * exp(w_2*b_2) * exp(w_3*b_3)
  // Thus, the derivatives of b_1, b_2, b_3 have to be calculated.
  const double one_third = 1.0 / 3.0;
  EVector3 w1 = inverseRotate(quat_W_A, W_omega_W_A * one_third);
  EVector3 w3 = inverseRotate(quat_W_B, W_omega_W_B * one_third);
  EVector3 w2 = quaternionLog( quaternionExp(-w1) * invert(quat_W_A) * quat_W_B * quaternionExp(-w3));

  double alpha2 = alpha * alpha;
  double alpha3 = alpha2 * alpha;
  double beta1 = alpha3 - 3.0 * alpha2 + 3.0 * alpha;
  double beta2 = -2.0 * alpha3 + 3.0 * alpha2;
  //double beta3 = alpha3;

  // derivatives of b_1, b_2, b_3
  double dotBeta1 = 3.0 * alpha2 - 6.0 * alpha + 3.0;
  double dotBeta2 = -6.0 * alpha2 + 6.0 * alpha;
  double dotBeta3 = 3.0 * alpha2;

  EQuaternion qa_exp_w1beta1 = quat_W_A * quaternionExp(w1 * beta1);

  return
      quat_W_A * w1 * dotBeta1 +
      qa_exp_w1beta1 * w2 * dotBeta2 +
      qa_exp_w1beta1 * quaternionExp(w2 * beta2) * w3 * dotBeta3;
}

}  // namespace minimal
}  // namespace kindr

#endif // MINKINDR_QUATERNION_CUBIC_HERMITE_GTSAM_H
