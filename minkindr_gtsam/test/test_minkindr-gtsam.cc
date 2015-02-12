
#include <kindr/minimal/quat-transformation-gtsam.h>
#include <kindr/minimal/rotation-quaternion-gtsam.h>
#include <kindr/minimal/cubic-hermite-quaternion-gtsam.h>
#include <kindr/minimal/common-gtsam.h>
#include <kindr/minimal/testing-gtsam.h>
#include "../include/kindr/minimal/cubic-hermite-translation-gtsam.h"

typedef kindr::minimal::QuatTransformation Transformation;
typedef kindr::minimal::RotationQuaternion Quaternion;

using kindr::minimal::ETransformation;
using kindr::minimal::EQuaternion;
using kindr::minimal::EVector3;
using kindr::minimal::EVector6;

TEST(MinkindrGtsamTests, testRotateVector1) {
  using gtsam::Expression;
  Quaternion Cval;
  Eigen::Vector3d vval;
  Cval.setRandom();
  vval.setRandom();

  // Create some values
  gtsam::Values values;
  values.insert(1, Cval);
  values.insert(2, vval);

  EQuaternion C(1);
  EVector3 v(2);
  EVector3 Cv = C * v;

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Cv.value(values), Cval.rotate(vval), 1e-9));

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(Cv, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testRotateVector2) {
  using gtsam::Expression;
  Quaternion Cval;
  Eigen::Vector3d vval;
  Cval.setRandom();
  vval.setRandom();

  // Create some values
  gtsam::Values values;
  values.insert(1, Cval);
  values.insert(2, vval);

  EQuaternion C(1);
  EVector3 v(2);
  EVector3 Cv = rotate(C,v);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Cv.value(values), Cval.rotate(vval), 1e-9));

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(Cv, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testInverseRotateVector) {
  using gtsam::Expression;
  Quaternion Cval;
  Eigen::Vector3d vval;
  Cval.setRandom();
  vval.setRandom();

  // Create some values
  gtsam::Values values;
  values.insert(1, Cval);
  values.insert(2, vval);

  EQuaternion C(1);
  EVector3 v(2);
  EVector3 Ctv = inverseRotate(C,v);

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Ctv.value(values), Cval.inverseRotate(vval), 1e-9));

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(Ctv, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testRotationQuaternionChart) {
  Quaternion Cval;
  Cval.setRandom();
  SCOPED_TRACE("Testing Default Chart.");
  gtsam::testDefaultChart(Cval);
}

TEST(MinkindrGtsamTests, testTransformPoint) {
  using namespace gtsam;

  Transformation Tval;
  Eigen::Vector3d pval;
  Tval.setRandom();
  pval.setRandom();

  // Create some values
  Values values;
  values.insert(1, Tval);
  values.insert(2, pval);

  ETransformation T(1);
  EVector3 p(2);

  EVector3 Tp = T * p;

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(Tp, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testInverseTransformPoint) {
  using namespace gtsam;

  Transformation Tval;
  Eigen::Vector3d pval;
  Tval.setRandom();
  pval.setRandom();

  // Create some values
  Values values;
  values.insert(1, Tval);
  values.insert(2, pval);

  ETransformation T(1);
  EVector3 p(2);

  EVector3 invTp = inverseTransform(T,p);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(invTp, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testRotationExpression) {
  Quaternion Cval;
  Cval.setRandom();
  // Create some values
  gtsam::Values values;
  values.insert(1, Cval);
  EQuaternion C(1);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(C, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testRotationInverseExpression) {
  Quaternion Cval;
  Cval.setRandom();
  // Create some values
  gtsam::Values values;
  values.insert(1, Cval);
  EQuaternion C(1);
  EQuaternion invC = invert(C);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(invC, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testComposeRotationExpression) {
  Quaternion C1val;
  C1val.setRandom();
  Quaternion C2val;
  C2val.setRandom();
  // Create some values
  gtsam::Values values;
  values.insert(1, C1val);
  values.insert(2, C2val);
  EQuaternion C1(1);
  EQuaternion C2(1);
  EQuaternion C1C2 = C1 * C2;

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(C1C2, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testCombineRotationTranslation) {
  Quaternion Cval;
  Cval.setRandom();
  Eigen::Vector3d tval;
  tval.setRandom();
  // Create some values
  gtsam::Values values;
  values.insert(1, Cval);
  values.insert(2, tval);
  EQuaternion C(1);
  EVector3 t(2);
  ETransformation T = transformationFromComponents(C,t);

  Transformation Tval(Cval, tval);
  Transformation Teval = T.value(values);
  Transformation eye = Tval * Teval.inverted();

  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Tval.getTransformationMatrix(), Teval.getTransformationMatrix(), 1e-9));

  const double fd_step = 1e-5;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(T, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testSO3Log) {
  using gtsam::Expression;
  Quaternion Cval;
  Cval.setRandom();

  // Create some values
  gtsam::Values values;
  values.insert(1, Cval);

  EQuaternion C(1);
  EVector3 logC = kindr::minimal::quaternionLog(C);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(logC, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testSO3Exp) {
  using gtsam::Expression;
  Eigen::Vector3d pval;
  pval.setRandom();

  // Create some values
  gtsam::Values values;
  values.insert(1, pval);

  EVector3 p(1);
  EQuaternion C = kindr::minimal::quaternionExp(p);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(C, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testSE3Exp) {
  using gtsam::Expression;
  gtsam::Vector6 pval;
  pval.setRandom();

  // Create some values
  gtsam::Values values;
  values.insert(1, pval);

  Expression<gtsam::Vector6> p(1);
  ETransformation T = kindr::minimal::transformationExp(p);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(T, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testSO3LogIdentity) {
  using gtsam::Expression;
  Quaternion Cval;

  // Create some values
  gtsam::Values values;
  values.insert(1, Cval);

  EQuaternion C(1);
  EVector3 logC = kindr::minimal::quaternionLog(C);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(logC, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testRotationFromTransformation) {
  using namespace gtsam;
  using ::Quaternion;
  Transformation Tval;
  Tval.setRandom();

  // Create some values
  Values values;
  values.insert(1, Tval);

  ETransformation T(1);

  EQuaternion q = rotationFromTransformation(T);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(q, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testTranslationFromTransformation) {
  using namespace gtsam;

  Transformation Tval;
  Tval.setRandom();

  // Create some values
  Values values;
  values.insert(1, Tval);

  ETransformation T(1);

  EVector3 t = translationFromTransformation(T);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(t, values, fd_step, tolerance);
}


TEST(MinkindrGtsamTests, testTransform) {
  using namespace gtsam;

  Transformation Tval;
  Tval.setRandom();

  // Create some values
  Values values;
  values.insert(1, Tval);

  ETransformation T(1);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(T, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testInverseTransform) {
  using namespace gtsam;

  Transformation Tval;
  Tval.setRandom();

  // Create some values
  Values values;
  values.insert(1, Tval);

  ETransformation T(1);
  ETransformation invT = inverse(T);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(invT, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testComposeTransform1) {
  using namespace gtsam;

  Transformation T1val;
  T1val.setRandom();
  Transformation T2val;
  T2val.setRandom();

  // Create some values
  Values values;
  values.insert(1, T1val);
  values.insert(2, T2val);

  ETransformation T1(1);
  ETransformation T2(2);
  ETransformation T1T2 = compose(T1,T2);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(T1T2, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testComposeTransform2) {
  using namespace gtsam;

  Transformation T1val;
  T1val.setRandom();
  Transformation T2val;
  T2val.setRandom();

  // Create some values
  Values values;
  values.insert(1, T1val);
  values.insert(2, T2val);

  ETransformation T1(1);
  ETransformation T2(2);
  ETransformation T1T2 = T1 * T2;

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(T1T2, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testTransformationLog) {
  using namespace gtsam;

  Transformation Tval;
  Tval.setRandom();

  // Create some values
  Values values;
  values.insert(1, Tval);

  ETransformation T(1);
  Expression<Vector6> logT = kindr::minimal::transformationLog(T);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(logT, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testTransformationTranslationLog) {
  using namespace gtsam;

  Transformation Tval;
  Tval.setRandom();

  // Create some values
  Values values;
  values.insert(1, Tval);

  ETransformation T(1);
  EVector3 logt = translationLog(T);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(logt, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testTransformationRotationLog) {
  using namespace gtsam;

  Transformation Tval;
  Tval.setRandom();

  // Create some values
  Values values;
  values.insert(1, Tval);

  ETransformation T(1);
  EVector3 logC = rotationLog(T);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(logC, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testInvertAndComposeTransform) {
  using namespace gtsam;

  Transformation T1val;
  T1val.setRandom();
  Transformation T2val;
  T2val.setRandom();

  // Create some values
  Values values;
  values.insert(1, T1val);
  values.insert(2, T2val);

  ETransformation T1(1);
  ETransformation T2(2);
  ETransformation invT1T2 = invertAndCompose(T1,T2);

  const double fd_step = 1e-6;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(invT1T2, values, fd_step, tolerance);
}

TEST(MinkindrGtsamTests, testSlerp) {
  using namespace gtsam;

  Transformation T1val;
  T1val.setRandom();
  Transformation T2val;
  T2val.setRandom();

  // Create some values
  Values values;
  values.insert(1, T1val);
  values.insert(2, T2val);

  ETransformation T1(1);
  ETransformation T2(2);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;

  {
    ETransformation slerpT1 = slerp(T1, T2, 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(slerpT1, values, fd_step, tolerance);
  }
  {
    ETransformation slerpT2 = slerp(T1, T2, 1.0 - 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(slerpT2, values, fd_step, tolerance);
  }
  {
    ETransformation slerpTa = slerp(T1, T2, 0.25);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(slerpTa, values, fd_step, tolerance);
  }

  {
    ETransformation slerpTb = slerp(T1, T2, 0.75);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(slerpTb, values, fd_step, tolerance);
  }
}

TEST(MinkindrGtsamTests, testCubicHermiteQuaternion) {
  using namespace gtsam;

  Quaternion qaVal;  qaVal.setRandom();
  Quaternion qbVal;  qbVal.setRandom();
  Eigen::Vector3d waVal; waVal.setRandom();
  Eigen::Vector3d wbVal; wbVal.setRandom();

  // Create some values
  Values values;
  values.insert(1, qaVal);
  values.insert(2, qbVal);
  values.insert(3, waVal);
  values.insert(4, wbVal);

  EQuaternion qA(1), qB(2);
  EVector3 wA(3), wB(4);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;

  {
    EQuaternion interpT1 = hermiteInterpolation(qA, wA, qB, wB, 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpT1, values, fd_step, tolerance);
  }
  {
    EQuaternion interpT2 = hermiteInterpolation(qA, wA, qB, wB, 1.0 - 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpT2, values, fd_step, tolerance);
  }
  {
    EQuaternion interpTa = hermiteInterpolation(qA, wA, qB, wB, 0.25);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpTa, values, fd_step, tolerance);
  }
  {
    EQuaternion interpTb = hermiteInterpolation(qA, wA, qB, wB, 0.75);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpTb, values, fd_step, tolerance);
  }

  {
    EVector3 interpT1 = hermiteInterpolationDerivative(qA, wA, qB, wB, 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpT1, values, fd_step, tolerance);
  }
  {
    EVector3 interpT2 = hermiteInterpolationDerivative(qA, wA, qB, wB, 1.0 - 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpT2, values, fd_step, tolerance);
  }
  {
    EVector3 interpTa = hermiteInterpolationDerivative(qA, wA, qB, wB, 0.25);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpTa, values, fd_step, tolerance);
  }
  {
    EVector3 interpTb = hermiteInterpolationDerivative(qA, wA, qB, wB, 0.75);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpTb, values, fd_step, tolerance);
  }
}

TEST(MinkindrGtsamTests, testCubicHermiteQuaternionDerivative) {
  using namespace gtsam;

  Quaternion qaVal;  qaVal.setRandom();
  Quaternion qbVal;  qbVal.setRandom();
  Eigen::Vector3d waVal; waVal.setRandom() *= 10.0;
  Eigen::Vector3d wbVal; wbVal.setRandom() *= 10.0;

  // Create some values
  Values values;
  values.insert(1, qaVal);
  values.insert(2, qbVal);
  values.insert(3, waVal);
  values.insert(4, wbVal);

  EQuaternion qA(1), qB(2);
  EVector3 wA(3), wB(4);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;

  {
    EQuaternion interpQ = hermiteInterpolation(qA, wA, qB, wB, 1e-5);
    EVector3 interpV = hermiteInterpolationDerivative(qA, wA, qB, wB, 1e-5);

    Quaternion qI = interpQ.value(values);
    Eigen::Vector3d vI = interpV.value(values);

    EXPECT_TRUE(EIGEN_MATRIX_NEAR(vI,waVal,1e-3));
    Eigen::Vector3d dq = (qI * qaVal.inverted()).log()/1e-5;
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(vI,dq,1e-3));

  }

  {
    const int N = 100;
    for (int i = 0; i < N; ++i) {
      double alpha = double(i)/(N-1);
      EQuaternion interpQ = hermiteInterpolation(qA, wA, qB, wB, alpha);
      EQuaternion interpQp = hermiteInterpolation(qA, wA, qB, wB, alpha+fd_step);
      EVector3 interpV = hermiteInterpolationDerivative(qA, wA, qB, wB, alpha);

      Eigen::Vector3d dq = (interpQp.value(values) * interpQ.value(values).inverted()).log()/fd_step;
      Eigen::Vector3d v = interpV.value(values);

      EXPECT_TRUE(EIGEN_MATRIX_NEAR(v,dq,tolerance*5));
    }

  }

}

TEST(MinkindrGtsamTests, testCubicHermiteTranslation) {
  using namespace gtsam;

  Eigen::Vector3d taVal; taVal.setRandom();
  Eigen::Vector3d tbVal; tbVal.setRandom();
  Eigen::Vector3d vaVal; vaVal.setRandom();
  Eigen::Vector3d vbVal; vbVal.setRandom();

  // Create some values
  Values values;
  values.insert(1, taVal);
  values.insert(2, tbVal);
  values.insert(3, vaVal);
  values.insert(4, vbVal);

  EVector3 tA(1), tB(2);
  EVector3 vA(3), vB(4);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;

  {
    EVector3 interpT1 = kindr::minimal::hermiteTranslationInterpolation(tA, vA, tB, vB, 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpT1, values, fd_step, tolerance);
  }
  {
    EVector3 interpT2 = kindr::minimal::hermiteTranslationInterpolation(tA, vA, tB, vB, 1.0 - 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpT2, values, fd_step, tolerance);
  }
  {
    EVector3 interpTa = kindr::minimal::hermiteTranslationInterpolation(tA, vA, tB, vB, 0.25);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpTa, values, fd_step, tolerance);
  }
  {
    EVector3 interpTb = kindr::minimal::hermiteTranslationInterpolation(tA, vA, tB, vB, 0.75);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpTb, values, fd_step, tolerance);
  }

  {
    EVector3 interpT1 = kindr::minimal::hermiteTranslationInterpolationDerivative(tA, vA, tB, vB, 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpT1, values, fd_step, tolerance);
  }
  {
    EVector3 interpT2 = kindr::minimal::hermiteTranslationInterpolationDerivative(tA, vA, tB, vB, 1.0 - 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpT2, values, fd_step, tolerance);
  }
  {
    EVector3 interpTa = kindr::minimal::hermiteTranslationInterpolationDerivative(tA, vA, tB, vB, 0.25);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpTa, values, fd_step, tolerance);
  }
  {
    EVector3 interpTb = kindr::minimal::hermiteTranslationInterpolationDerivative(tA, vA, tB, vB, 0.75);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(interpTb, values, fd_step, tolerance);
  }
}

TEST(MinkindrGtsamTests, testCubicHermiteTranslationDerivative) {
  using namespace gtsam;


  Eigen::Vector3d  taVal;  taVal.setRandom();
  Eigen::Vector3d  tbVal;  tbVal.setRandom();
  Eigen::Vector3d vaVal; vaVal.setRandom() *= 10.0;
  Eigen::Vector3d vbVal; vbVal.setRandom() *= 10.0;

  // Create some values
  Values values;
  values.insert(1, taVal);
  values.insert(2, tbVal);
  values.insert(3, vaVal);
  values.insert(4, vbVal);

  EVector3 tA(1), tB(2);
  EVector3 vA(3), vB(4);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;

  {
    EVector3 interpTrans = kindr::minimal::hermiteTranslationInterpolation(tA, vA, tB, vB, 1e-5);
    EVector3 interpV = kindr::minimal::hermiteTranslationInterpolationDerivative(tA, vA, tB, vB, 1e-5);

    Eigen::Vector3d tI = interpTrans.value(values);
    Eigen::Vector3d vI = interpV.value(values);

    EXPECT_TRUE(EIGEN_MATRIX_NEAR(vI,vaVal,1e-3));
    Eigen::Vector3d dtrans = (tI - taVal)/1e-5;
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(vI,dtrans,1e-3));

  }

  {
    const int N = 100;
    for (int i = 0; i < N; ++i) {
      double alpha = double(i)/(N-1);
      EVector3 interpTrans = kindr::minimal::hermiteTranslationInterpolation(tA, vA, tB, vB, alpha);
      EVector3 interpTransp = kindr::minimal::hermiteTranslationInterpolation(tA, vA, tB, vB, alpha+fd_step);
      EVector3 interpV = kindr::minimal::hermiteTranslationInterpolationDerivative(tA, vA, tB, vB, alpha);

      Eigen::Vector3d dtrans = (interpTransp.value(values) - interpTrans.value(values))/fd_step;
      Eigen::Vector3d v = interpV.value(values);

      EXPECT_TRUE(EIGEN_MATRIX_NEAR(v,dtrans,tolerance*5));
    }

  }

}

TEST(MinkindrGtsamTests, testVectorScaling) {
  using namespace gtsam;
  Eigen::Vector3d v; v.setRandom();
  double a = (double)rand()/RAND_MAX*10;

  // Create some values
  Values values;
  values.insert(1, v);

  EVector3 w(1);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;

  {
    EVector3 wScaled = kindr::minimal::vectorScaling(w, a);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(wScaled.value(values),v*a,tolerance));
    testExpressionJacobians(wScaled, values, fd_step, tolerance);
  }
}

TEST(MinkindrGtsamTests, testVectorSumAndDifference) {
  using namespace gtsam;
  Eigen::Vector3d vA; vA.setRandom();
  Eigen::Vector3d vB; vB.setRandom();

  // Create some values
  Values values;
  values.insert(1, vA);
  values.insert(2, vB);

  EVector3 wA(1), wB(2);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;

  {
    EVector3 sum = kindr::minimal::vectorSum(wA, wB);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(sum.value(values),vA+vB,tolerance));
    testExpressionJacobians(sum, values, fd_step, tolerance);
  }

  {
    EVector3 difference = kindr::minimal::vectorDifference(wA, wB);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(difference.value(values),vA-vB,tolerance));
    testExpressionJacobians(difference, values, fd_step, tolerance);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, false);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  return RUN_ALL_TESTS();
}
