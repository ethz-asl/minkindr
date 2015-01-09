
#include <kindr/minimal/quat-transformation-gtsam.h>
#include <kindr/minimal/rotation-quaternion-gtsam.h>
#include <kindr/minimal/testing-gtsam.h>

typedef kindr::minimal::QuatTransformation Transformation;
typedef kindr::minimal::RotationQuaternion Quaternion;

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

  Expression<Quaternion> C(1);
  Expression<Eigen::Vector3d> v(2);
  Expression<Eigen::Vector3d> Cv = C * v;

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

  Expression<Quaternion> C(1);
  Expression<Eigen::Vector3d> v(2);
  Expression<Eigen::Vector3d> Cv = rotate(C,v);

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

  Expression<Quaternion> C(1);
  Expression<Eigen::Vector3d> v(2);
  Expression<Eigen::Vector3d> Ctv = inverseRotate(C,v);

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

  Expression<Transformation> T(1);
  Expression<Eigen::Vector3d> p(2);

  Expression<Eigen::Vector3d> Tp = T * p;

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

  Expression<Transformation> T(1);
  Expression<Eigen::Vector3d> p(2);

  Expression<Eigen::Vector3d> invTp = inverseTransform(T,p);

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
  gtsam::Expression<Quaternion> C(1);

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
  gtsam::Expression<Quaternion> C(1);
  gtsam::Expression<Quaternion> invC = invert(C);

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
  gtsam::Expression<Quaternion> C1(1);
  gtsam::Expression<Quaternion> C2(1);
  gtsam::Expression<Quaternion> C1C2 = C1 * C2;

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
  gtsam::Expression<Quaternion> C(1);
  gtsam::Expression<Eigen::Vector3d> t(2);
  gtsam::Expression<Transformation> T = transformationFromComponents(C,t);

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

  Expression<Quaternion> C(1);
  Expression<Eigen::Vector3d> logC = log(C);

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

  Expression<Eigen::Vector3d> p(1);
  Expression<Quaternion> C = exp(p);

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
  Expression<Transformation> T = exp(p);

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

  Expression<Quaternion> C(1);
  Expression<Eigen::Vector3d> logC = log(C);

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

  Expression<Transformation> T(1);

  Expression<Quaternion> q = rotationFromTransformation(T);

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

  Expression<Transformation> T(1);

  Expression<Eigen::Vector3d> t = translationFromTransformation(T);

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

  Expression<Transformation> T(1);

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

  Expression<Transformation> T(1);
  Expression<Transformation> invT = inverse(T);

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

  Expression<Transformation> T1(1);
  Expression<Transformation> T2(2);
  Expression<Transformation> T1T2 = compose(T1,T2);

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

  Expression<Transformation> T1(1);
  Expression<Transformation> T2(2);
  Expression<Transformation> T1T2 = T1 * T2;

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

  Expression<Transformation> T(1);
  Expression<Vector6> logT = log(T);

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

  Expression<Transformation> T(1);
  Expression<Eigen::Vector3d> logt = translationLog(T);

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

  Expression<Transformation> T(1);
  Expression<Eigen::Vector3d> logC = rotationLog(T);

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

  Expression<Transformation> T1(1);
  Expression<Transformation> T2(2);
  Expression<Transformation> invT1T2 = invertAndCompose(T1,T2);

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

  Expression<Transformation> T1(1);
  Expression<Transformation> T2(2);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;

  {
    Expression<Transformation> slerpT1 = slerp(T1, T2, 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(slerpT1, values, fd_step, tolerance);
  }
  {
    Expression<Transformation> slerpT2 = slerp(T1, T2, 1.0 - 1e-5);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(slerpT2, values, fd_step, tolerance);
  }
  {
    Expression<Transformation> slerpTa = slerp(T1, T2, 0.25);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(slerpTa, values, fd_step, tolerance);
  }

  {
    Expression<Transformation> slerpTb = slerp(T1, T2, 0.75);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(slerpTb, values, fd_step, tolerance);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, false);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  return RUN_ALL_TESTS();
}
