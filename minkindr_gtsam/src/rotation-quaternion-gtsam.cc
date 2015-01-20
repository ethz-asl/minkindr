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
  return C.inverted();
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

Eigen::Vector3d rotationLogImplementation(const RotationQuaternion& C,
                                          OptionalJacobian<3, 3> JC) {
  Eigen::Vector3d aa = C.log();
  if(JC) {
    double phi = aa.norm();
    if(phi < 1e-10) {
      JC->setIdentity();
    } else {
      double cot = - std::sin(phi)/(std::cos(phi)-1);
      double a1 = 1/(phi*phi) * (1.0 - 0.5*phi*cot);
      Eigen::Matrix3d px = kindr::minimal::skewMatrix(aa);
      (*JC) = Eigen::Matrix3d::Identity() - 0.5 * px + a1 * px * px;
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
