#include <kindr/minimal/quat-transformation-gtsam.h>
namespace kindr {
namespace minimal {
Eigen::Vector3d transform_point(const kindr::minimal::QuatTransformation& T,
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
// Expression<Eigen::Vector3d> Tp = T * p;
// instead of
// Expression<Eigen::Vector3d> Tp = Expression<Eigen::Vector3d>(&transform_point, T, p);
gtsam::Expression<Eigen::Vector3d>
operator*(const gtsam::Expression<kindr::minimal::QuatTransformation>& T,
          const gtsam::Expression<Eigen::Vector3d>& p) {
  return Expression<Eigen::Vector3d>(&transform_point, T, p);
}

gtsam::Expression<Eigen::Vector3d>
transform(const gtsam::Expression<kindr::minimal::QuatTransformation>& T,
          const gtsam::Expression<Eigen::Vector3d>& p) {
  return Expression<Eigen::Vector3d>(&transform_point, T, p);
}

kindr::minimal::QuatTransformation combine_components(
    const kindr::minimal::RotationQuaternion& C_A_B, const Eigen::Vector3d& A_t_B,
    OptionalJacobian<6,3> HC, OptionalJacobian<6,3> Hp) {
  if(HC) {
    HC->topRows<3>() = kindr::minimal::skewMatrix(A_t_B);
    HC->bottomRows<3>() = Eigen::Matrix3d::Identity();
  }

  if(Hp) {
    Hp->topRows<3>() = Eigen::Matrix3d::Identity();
    Hp->bottomRows<3>().setZero();
  }

  return kindr::minimal::QuatTransformation(C_A_B, A_t_B);
}

// Build a transformation expression from a rotation expression and a point expression.
gtsam::Expression<kindr::minimal::QuatTransformation> transformationFromComponents(
    const gtsam::Expression<kindr::minimal::RotationQuaternion>& C_A_B,
    const gtsam::Expression<Eigen::Vector3d>& A_t_B) {
  return Expression<kindr::minimal::QuatTransformation>(&combine_components, C_A_B, A_t_B);
}

kindr::minimal::RotationQuaternion rotationFromTransformationImplementation(
    const kindr::minimal::QuatTransformation& T, OptionalJacobian<3, 6> HT) {
  if(HT) {
    HT->leftCols<3>().setZero();
    HT->rightCols<3>().setIdentity();
  }
  return T.getRotation();
}

Expression<kindr::minimal::RotationQuaternion> rotationFromTransformation(
    const Expression<kindr::minimal::QuatTransformation>& T) {
  return Expression<kindr::minimal::RotationQuaternion>(
      &rotationFromTransformationImplementation, T);
}

Eigen::Vector3d translationFromTransformationImplementation(
    const kindr::minimal::QuatTransformation& T, OptionalJacobian<3, 6> HT) {
  if(HT) {
    HT->leftCols<3>().setIdentity();
    HT->rightCols<3>() = -kindr::minimal::skewMatrix(T.getPosition());
  }
  return T.getPosition();
}

gtsam::Expression<Eigen::Vector3d> translationFromTransformation(
    const gtsam::Expression<kindr::minimal::QuatTransformation>& T) {
  return Expression<Eigen::Vector3d>(
      &translationFromTransformationImplementation, T);
}

Eigen::Vector3d inverseTransformImplementation(
    const kindr::minimal::QuatTransformation& T, const Eigen::Vector3d& p,
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

Expression<Eigen::Vector3d> inverseTransform(
    const Expression<kindr::minimal::QuatTransformation>& T,
    const Expression<Eigen::Vector3d>& p) {
  return Expression<Eigen::Vector3d>(&inverseTransformImplementation, T, p);
}

kindr::minimal::QuatTransformation inverseImplementation(
    const kindr::minimal::QuatTransformation& T, OptionalJacobian<6, 6> HT) {
  kindr::minimal::QuatTransformation invT = T.inverted();
  if(HT) {
    Eigen::Matrix3d ninvC = -invT.getRotationMatrix();
    HT->topLeftCorner<3,3>() = ninvC;
    HT->bottomRightCorner<3,3>() = ninvC;
    HT->bottomLeftCorner<3,3>().setZero();
    HT->topRightCorner<3,3>() = kindr::minimal::skewMatrix(invT.getPosition())*ninvC;
  }
  return invT;
}

Expression<kindr::minimal::QuatTransformation> inverse(
    const Expression<kindr::minimal::QuatTransformation>& T) {
  return Expression<kindr::minimal::QuatTransformation>(&inverseImplementation, T);
}

kindr::minimal::QuatTransformation composeImplementation(
    const kindr::minimal::QuatTransformation& T1,
    const kindr::minimal::QuatTransformation& T2,
    OptionalJacobian<6, 6> HT1,
    OptionalJacobian<6, 6> HT2) {
  kindr::minimal::QuatTransformation T1T2 = T1 * T2;
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

Expression<kindr::minimal::QuatTransformation> compose(
    const Expression<kindr::minimal::QuatTransformation>& T1,
    const Expression<kindr::minimal::QuatTransformation>& T2) {
  return Expression<kindr::minimal::QuatTransformation>(&composeImplementation, T1, T2);
}

Vector6 transformationLogImplementation(const kindr::minimal::QuatTransformation& T,
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

Expression<Vector6> log(const Expression<kindr::minimal::QuatTransformation>& T) {
  return Expression<Vector6>(&transformationLogImplementation, T);
}

Eigen::Vector3d rotationFromTransformationLogImplementation(const kindr::minimal::QuatTransformation& T,
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

Expression<Eigen::Vector3d> rotationLog(
    const Expression<kindr::minimal::QuatTransformation>& T) {
  return Expression<Eigen::Vector3d>(&rotationFromTransformationLogImplementation, T);
}

kindr::minimal::QuatTransformation invertAndComposeImplementation(
    const kindr::minimal::QuatTransformation& T1,
    const kindr::minimal::QuatTransformation& T2,
    OptionalJacobian<6, 6> HT1,
    OptionalJacobian<6, 6> HT2) {
  kindr::minimal::QuatTransformation invT1 = inverseImplementation(T1, HT1);

  return composeImplementation(invT1, T2, boost::none, HT2);
}

/// \brief Compose two transformations as inv(T1)*T2.
Expression<kindr::minimal::QuatTransformation> invertAndCompose(
    const Expression<kindr::minimal::QuatTransformation>& T1,
    const Expression<kindr::minimal::QuatTransformation>& T2) {
  return Expression<kindr::minimal::QuatTransformation>(&invertAndComposeImplementation, T1, T2);
}

Vector6 vectorScalingImplementation(const Vector6& v, double alpha, OptionalJacobian<6, 6> H) {
  if (H) {
    *H = OptionalJacobian<6,6>::Jacobian::Identity()*alpha;
  }
  return v*alpha;
}

Expression<Vector6> vectorScaling(const Expression<Vector6>& v, double alpha) {
  return Expression<Vector6>(boost::bind(&vectorScalingImplementation, _1, alpha, _2), v);
}

Expression<kindr::minimal::QuatTransformation> slerp(
    const Expression<kindr::minimal::QuatTransformation>& T0,
    const Expression<kindr::minimal::QuatTransformation>& T1,
    double alpha) {
  return compose(T0, exp(vectorScaling(log(invertAndCompose(T0, T1)), alpha)));
}

kindr::minimal::QuatTransformation transformationExpImplementation(
    const Vector6& params,
    OptionalJacobian<6, 6> Hp) {
  if(Hp) {
    Eigen::Matrix3d S;
    kindr::minimal::RotationQuaternion q = rotationExpImplementation(params.tail<3>(), S);
    Hp->topLeftCorner<3,3>().setIdentity();
    Hp->bottomLeftCorner<3,3>().setZero();
    Hp->topRightCorner<3,3>() = kindr::minimal::skewMatrix(params.head<3>()) * S;
    Hp->bottomRightCorner<3,3>() = S;
    return kindr::minimal::QuatTransformation(q, params.head<3>());
  } else {
    return kindr::minimal::QuatTransformation::exp(params);
  }
}

Expression<kindr::minimal::QuatTransformation> exp(const Expression<Vector6>& params) {
  return Expression<kindr::minimal::QuatTransformation>(
      &transformationExpImplementation, params);
}

}  // namespace minimal
}  // namespace kindr
