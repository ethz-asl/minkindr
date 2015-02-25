
#include <kindr/minimal/quat-transformation-gtsam.h>
#include <kindr/minimal/rotation-quaternion-gtsam.h>
#include <kindr/minimal/cubic-hermite-quaternion-gtsam.h>
#include <kindr/minimal/common-gtsam.h>
#include <kindr/minimal/testing-gtsam.h>

#define N_TEST_ITERATIONS 10000

typedef kindr::minimal::QuatTransformation Transformation;
typedef kindr::minimal::RotationQuaternion Quaternion;

using kindr::minimal::ETransformation;
using kindr::minimal::EQuaternion;
using kindr::minimal::EVector3;
using kindr::minimal::EVector6;

const double tolerance = 1e-5;
const double fd_step = 1e-9;

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

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    Tval.setRandom();
    pval.setRandom();
    values.update(1, Tval);
    values.update(2, pval);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(Tp, values, fd_step, tolerance);
  }

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

  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(logC, values, fd_step, tolerance);

  for(int i = 0; i < N_TEST_ITERATIONS; ++i) {
    Cval.setRandom();
    values.update(1, Cval);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(logC, values, fd_step, tolerance);
  }

  // Testing at values near identity quaternion
  Eigen::Vector3d aa(0,0,0);
  for(int i = 0; i < N_TEST_ITERATIONS; ++i) {
    Cval = Quaternion(aa);
    values.update(1, Cval);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(logC, values, fd_step, tolerance);
    aa(0) += 0.5e-11;
    aa(1) -= 0.4e-12;
    aa(2) += 0.7e-12;
  }
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

  for(int i = 0; i < N_TEST_ITERATIONS; ++i) {
    pval.setRandom();
    values.update(1, pval);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(T, values, fd_step, tolerance);
  }

}

TEST(MinkindrGtsamTests, testSO3LogIdentity) {
  using gtsam::Expression;
  Quaternion Cval;

  // Create some values
  gtsam::Values values;
  values.insert(1, Cval);

  EQuaternion C(1);
  EVector3 logC = kindr::minimal::quaternionLog(C);

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
  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    Tval.setRandom();
    values.update(1, Tval);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(invT, values, fd_step, tolerance);
  }

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

  for(int i=0; i < N_TEST_ITERATIONS; ++i) {
    T1val.setRandom();
    T2val.setRandom();
    values.update(1, T1val);
    values.update(2, T2val);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(T1T2, values, fd_step, tolerance);
  }

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

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    Tval.setRandom();
    values.update(1, Tval);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(logT, values, fd_step, tolerance);
  }

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

  for(int i = 0; i < N_TEST_ITERATIONS; ++i) {
    T1val.setRandom();
    T2val.setRandom();
    values.update(1, T1val);
    values.update(2, T2val);

    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(invT1T2, values, fd_step, tolerance);
  }

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

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    T1val.setRandom();
    T2val.setRandom();
    values.update(1, T1val);
    values.update(2, T2val);
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

TEST(MinkindrGtsamTests, testVectorScaling) {
  using namespace gtsam;
  Eigen::Vector3d v; v.setRandom();
  double a = (double)rand()/RAND_MAX*10;

  // Create some values
  Values values;
  values.insert(1, v);

  EVector3 w(1);

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    v.setRandom();
    values.update(1, v);
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

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    vA.setRandom();
    vB.setRandom();
    values.update(1, vA);
    values.update(2, vB);
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

}

//This test simulates having a match constraint between two pairs of
//transformations. The result is a complex expression tree
TEST(MinkindrGtsamTests, expressionCombination) {
  using namespace gtsam;

  Transformation T1aVal, T1bVal, T2aVal, T2bVal, Tconst1Val, Tconst2Val;
  Eigen::Vector3d v1Val, v2Val;
  Tconst1Val.setRandom();
  Tconst2Val.setRandom();
  v1Val.setRandom();
  v2Val.setRandom();

  // Create some values
  Values values;
  values.insert(1, T1aVal);
  values.insert(2, T1bVal);
  values.insert(3, T2aVal);
  values.insert(4, T2bVal);


  ETransformation T1a(1);
  ETransformation T1b(2);
  ETransformation T2a(3);
  ETransformation T2b(4);

  ETransformation Tconst1(Tconst1Val);
  ETransformation Tconst2(Tconst2Val);

  EVector3 v1(v1Val);
  EVector3 v2(v2Val);

  ETransformation slerpT1 = slerp(T1a, T1b, 0.2154);
  ETransformation slerpT2 = slerp(T2a, T2b, 0.3658);

  ETransformation composedT1 = compose(slerpT1, Tconst1);
  ETransformation composedT2 = compose(slerpT2, Tconst2);

  EVector3 transformedV1 = transform(composedT1, v1);
  EVector3 transformedV2 = transform(composedT2, v2);

  EVector3 diff = kindr::minimal::vectorDifference(transformedV1, transformedV2);

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    T1aVal.setRandom();
    T1bVal.setRandom();
    T2aVal.setRandom();
    T2bVal.setRandom();

    values.update(1, T1aVal);
    values.update(2, T1bVal);
    values.update(3, T2aVal);
    values.update(4, T2bVal);

    testExpressionJacobians(diff, values, fd_step, tolerance);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, false);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  return RUN_ALL_TESTS();
}
