#include <kindr/minimal/quat-transformation-gtsam.h>

using namespace gtsam;

namespace kindr {
namespace minimal {
Eigen::Vector3d transform_point(const QuatTransformation& T,
                                const Eigen::Vector3d& p,
                                OptionalJacobian<3, 6> HT,
                                OptionalJacobian<3, 3> Hp) {
  Eigen::Vector3d Tp = T * p;
  if(HT) {
    // TODO(furgalep) fill in.
    HT->leftCols<3>().setIdentity();
    HT->rightCols<3>() = -kindr::minimal::skewMatrix(Tp);
  }

  if(Hp) {
    *Hp = T.getRotationMatrix();
  }

  return Tp;
}

// This is syntatic sugar to be able to write
// EVector3 Tp = T * p;
// instead of
// EVector3 Tp = EVector3(&transform_point, T, p);
EVector3
operator*(const ETransformation& T,
          const EVector3& p) {
  return EVector3(&transform_point, T, p);
}

EVector3
transform(const ETransformation& T,
          const EVector3& p) {
  return EVector3(&transform_point, T, p);
}

QuatTransformation combine_components(
    const RotationQuaternion& C_A_B, const Eigen::Vector3d& A_t_B,
    OptionalJacobian<6,3> HC, OptionalJacobian<6,3> Hp) {
  if(HC) {
    HC->topRows<3>() = kindr::minimal::skewMatrix(A_t_B);
    HC->bottomRows<3>() = Eigen::Matrix3d::Identity();
  }

  if(Hp) {
    Hp->topRows<3>() = Eigen::Matrix3d::Identity();
    Hp->bottomRows<3>().setZero();
  }

  return QuatTransformation(C_A_B, A_t_B);
}

// Build a transformation expression from a rotation expression and a point expression.
ETransformation transformationFromComponents(
    const EQuaternion& C_A_B,
    const EVector3& A_t_B) {
  return ETransformation(&combine_components, C_A_B, A_t_B);
}

RotationQuaternion rotationFromTransformationImplementation(
    const QuatTransformation& T, OptionalJacobian<3, 6> HT) {
  if(HT) {
    HT->leftCols<3>().setZero();
    HT->rightCols<3>().setIdentity();
  }
  return T.getRotation();
}

EQuaternion rotationFromTransformation(
    const ETransformation& T) {
  return EQuaternion(
      &rotationFromTransformationImplementation, T);
}

Eigen::Vector3d translationFromTransformationImplementation(
    const QuatTransformation& T, OptionalJacobian<3, 6> HT) {
  if(HT) {
    HT->leftCols<3>().setIdentity();
    HT->rightCols<3>() = -kindr::minimal::skewMatrix(T.getPosition());
  }
  return T.getPosition();
}

EVector3 translationFromTransformation(
    const ETransformation& T) {
  return EVector3(
      &translationFromTransformationImplementation, T);
}

Eigen::Vector3d inverseTransformImplementation(
    const QuatTransformation& T, const Eigen::Vector3d& p,
    OptionalJacobian<3, 6> HT, OptionalJacobian<3, 3> Hp) {
  Eigen::Vector3d Tp = T.inverseTransform(p);
  if(HT || Hp) {
    Eigen::Matrix3d invC = T.getRotation().inverted().getRotationMatrix();
    if(HT) {
      HT->leftCols<3>() = -invC;
      HT->rightCols<3>() = invC * kindr::minimal::skewMatrix(p);
    }
    if(Hp) {
      (*Hp) = invC;
    }
  }
  return Tp;
}

EVector3 inverseTransform(
    const ETransformation& T,
    const EVector3& p) {
  return EVector3(&inverseTransformImplementation, T, p);
}

QuatTransformation inverseImplementation(
    const QuatTransformation& T, OptionalJacobian<6, 6> HT) {
  QuatTransformation invT = T.inverted();
  if(HT) {
    Eigen::Matrix3d ninvC = -invT.getRotationMatrix();
    HT->topLeftCorner<3,3>() = ninvC;
    HT->bottomRightCorner<3,3>() = ninvC;
    HT->bottomLeftCorner<3,3>().setZero();
    HT->topRightCorner<3,3>() = kindr::minimal::skewMatrix(invT.getPosition())*ninvC;
  }
  return invT;
}

ETransformation inverse(
    const ETransformation& T) {
  return ETransformation(&inverseImplementation, T);
}

QuatTransformation composeImplementation(
    const QuatTransformation& T1,
    const QuatTransformation& T2,
    OptionalJacobian<6, 6> HT1,
    OptionalJacobian<6, 6> HT2) {
  QuatTransformation T1T2 = T1 * T2;
  if(HT1) {
    HT1->setIdentity();
  }

  if(HT2) {
    Eigen::Matrix3d C1 = T1.getRotationMatrix();
    HT2->topLeftCorner<3,3>() = C1;
    HT2->bottomRightCorner<3,3>() = C1;
    HT2->bottomLeftCorner<3,3>().setZero();
    HT2->topRightCorner<3,3>() = kindr::minimal::skewMatrix(T1.getPosition()) * C1;
  }

  return T1T2;
}

ETransformation compose(
    const ETransformation& T1,
    const ETransformation& T2) {
  return ETransformation(&composeImplementation, T1, T2);
}

Vector6 transformationLogImplementation(const QuatTransformation& T,
                                        OptionalJacobian<6, 6> HT) {
  if(HT) {
    Vector6 logT;
    Eigen::Vector3d logC;
    Jacobian3x3 HC;
    logT.tail<3>() = rotationLogImplementation(T.getRotation(), HC);
    logT.head<3>() = T.getPosition();
    HT->topRightCorner<3,3>() = kindr::minimal::skewMatrix(-T.getPosition());
    HT->bottomLeftCorner<3,3>().setZero();
    HT->topLeftCorner<3,3>().setIdentity();
    HT->bottomRightCorner<3,3>() = HC;
    return logT;
  } else {
    return T.log();
  }
}

EVector6 transformationLog(const ETransformation& T) {
  return EVector6(&transformationLogImplementation, T);
}

Eigen::Vector3d rotationFromTransformationLogImplementation(const QuatTransformation& T,
                                                            OptionalJacobian<3, 6> HT) {

  if(HT) {
    Eigen::Vector3d logC;
    Jacobian3x3 HC;
    logC = rotationLogImplementation(T.getRotation(), HC);
    HT->leftCols<3>().setZero();
    HT->rightCols<3>() = HC;
    return logC;
  } else {
    return T.getRotation().log();
  }
}

EVector3 rotationLog(
    const ETransformation& T) {
  return EVector3(&rotationFromTransformationLogImplementation, T);
}

QuatTransformation invertAndComposeImplementation(
    const QuatTransformation& T1,
    const QuatTransformation& T2,
    OptionalJacobian<6, 6> HT1,
    OptionalJacobian<6, 6> HT2) {
  QuatTransformation invT1 = inverseImplementation(T1, HT1);

  return composeImplementation(invT1, T2, boost::none, HT2);
}

/// \brief Compose two transformations as inv(T1)*T2.
ETransformation invertAndCompose(
    const ETransformation& T1,
    const ETransformation& T2) {
  return ETransformation(&invertAndComposeImplementation, T1, T2);
}

ETransformation slerp(
    const ETransformation& T0,
    const ETransformation& T1,
    double alpha) {
  return compose(T0, transformationExp(vectorScaling(transformationLog(invertAndCompose(T0, T1)), alpha)));
}

QuatTransformation transformationExpImplementation(
    const Vector6& params,
    OptionalJacobian<6, 6> Hp) {
  if(Hp) {
    Eigen::Matrix3d S;
    RotationQuaternion q = rotationExpImplementation(params.tail<3>(), S);
    Hp->topLeftCorner<3,3>().setIdentity();
    Hp->bottomLeftCorner<3,3>().setZero();
    Hp->topRightCorner<3,3>() = kindr::minimal::skewMatrix(params.head<3>()) * S;
    Hp->bottomRightCorner<3,3>() = S;
    return QuatTransformation(q, params.head<3>());
  } else {
    return QuatTransformation::exp(params);
  }
}

ETransformation transformationExp(const EVector6& params) {
  return ETransformation(
      &transformationExpImplementation, params);
}

}  // namespace minimal
}  // namespace kindr
