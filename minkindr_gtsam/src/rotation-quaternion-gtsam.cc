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
#include <kindr/minimal/rotation-quaternion-gtsam.h>

using namespace gtsam;

namespace kindr {
namespace minimal {
Eigen::Vector3d rotate_point(
    const RotationQuaternion& C, const Eigen::Vector3d& p,
    OptionalJacobian<3, 3> HC, OptionalJacobian<3, 3> Hp) {
  Eigen::Vector3d Cp = C.rotate(p);
  if (HC) {
    kindr::minimal::skewMatrix(-Cp, &(*HC));
  }
  if (Hp) {
    *Hp = C.getRotationMatrix();
  }
  return Cp;
}

Eigen::Vector3d inverse_rotate_point(
    const RotationQuaternion& C, const Eigen::Vector3d& p,
    OptionalJacobian<3, 3> HC, OptionalJacobian<3, 3> Hp) {
  Eigen::Vector3d Ctp = C.inverseRotate(p);
  if(HC || Hp) {
    Eigen::Matrix3d Ct = C.getRotationMatrix().transpose();
    if (HC) {
      *HC = Ct * kindr::minimal::skewMatrix(p);
    }
    if (Hp) {
      *Hp = Ct;
    }
  }
  return Ctp;
}

RotationQuaternion invert_rotation_quaternion(
    const RotationQuaternion& C,
    OptionalJacobian<3, 3> HC) {
  if(HC) {
    *HC = -C.getRotationMatrix().transpose();
  }
  return C.inverse();
}

RotationQuaternion compose_rotation_quaternion(
    const RotationQuaternion& C1,
    const RotationQuaternion& C2,
    OptionalJacobian<3, 3> HC1,
    OptionalJacobian<3, 3> HC2) {
  if(HC1) {
    *HC1 = Eigen::Matrix3d::Identity();
  }
  if(HC2) {
    *HC2 = C1.getRotationMatrix();
  }
  return C1 * C2;
}

EQuaternion invert(const EQuaternion& q) {
  return EQuaternion(&invert_rotation_quaternion, q);
}

EQuaternion
operator*(const EQuaternion& C1, const EQuaternion& C2) {
  return EQuaternion(&compose_rotation_quaternion, C1, C2);
}

EVector3 operator*(const EQuaternion& C, const EVector3& p) {
  return EVector3(&rotate_point, C, p);
}
EVector3 rotate(const EQuaternion& C, const EVector3& p) {
  return EVector3(&rotate_point, C, p);
}
EVector3 inverseRotate(const EQuaternion& C, const EVector3& p) {
  return EVector3(&inverse_rotate_point, C, p);
}

template <typename Scalar_ = double>
inline bool isLessThanEpsilons4thRoot(Scalar_ x){
  static const Scalar_ epsilon4thRoot = pow(std::numeric_limits<Scalar_>::epsilon(), 1.0/4.0);
  return x < epsilon4thRoot;
}

Eigen::Vector3d rotationLogImplementation(const RotationQuaternion& C,
                                          OptionalJacobian<3, 3> JC) {
  Eigen::Vector3d aa = C.log();
  if(JC) {
    double phi = aa.norm();
    if(phi == 0) {
      JC->setIdentity();
    } else {

      double phiAbs = fabs(phi);
      Eigen::Matrix3d px = kindr::minimal::skewMatrix(aa);

      double a;
      if(!isLessThanEpsilons4thRoot(phiAbs)) {
        double phiHalf = 0.5 * phi;
        a = ((1 - phiHalf / tan(phiHalf))/phi/phi);
      } else {
        a = 1.0 / 12 * (1 + 1.0 / 60 * phi * phi);
      }

      (*JC) = Eigen::Matrix3d::Identity() - 0.5 * px + a * px * px;
    }
  }
  return aa;
}

EVector3 quaternionLog(const EQuaternion& C) {
  return EVector3(&rotationLogImplementation, C);
}

/// \brief Expose the rotation log and Jacobian.
RotationQuaternion rotationExpImplementation(const Eigen::Vector3d& p,
                                             OptionalJacobian<3, 3> Jp) {
  if(Jp) {
    double angle = p.norm();
    if(angle < 1e-14) {
      Jp->setIdentity();
    } else {
      double recip_angle = 1.0/angle;
      Eigen::Vector3d axis = p * recip_angle;
      double st2 = sin(angle * 0.5);
      double st  = sin(angle);

      double c1 = -2.0 * st2 * st2 * recip_angle;
      double c2 = (angle - st) * recip_angle;
      Eigen::Matrix3d crossA = kindr::minimal::skewMatrix(axis);

      (*Jp) = Eigen::Matrix3d::Identity() - (c1 * crossA) + (c2 * crossA * crossA);
    }
  }
  return RotationQuaternion::exp(p);
}

/// \brief Compute the matrix log of SO3.
EQuaternion quaternionExp(const EVector3& C) {
  return EQuaternion(&rotationExpImplementation, C);
}

}  // namespace minimal
}  // namespace kindr
