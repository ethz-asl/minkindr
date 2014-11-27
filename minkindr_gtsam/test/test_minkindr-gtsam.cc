#include <gtest/gtest.h>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtsam_unstable/nonlinear/Expression.h>
#include <gtsam_unstable/nonlinear/ExpressionFactor.h>
#include <kindr/minimal/quat-transformation-gtsam.h>
#include <kindr/minimal/rotation-quaternion-gtsam.h>
#include <gtsam/linear/VectorValues.h>
#include <gtsam/geometry/Rot3.h>
#include <eigen-checks/gtest.h>

// Compute finite difference Jacobians for an expression factor.
// This is extremely useful for unit testing during expression development.
template<typename T>
gtsam::JacobianFactor computeFiniteDifferenceJacobians(gtsam::ExpressionFactor<T> expression_factor,
                                                       const gtsam::Values& values,
                                                       double fd_step) {
  Eigen::VectorXd e = expression_factor.unwhitenedError(values);
  const size_t rows = e.size();

  std::map<gtsam::Key, Eigen::MatrixXd> jacobians;
  typename gtsam::ExpressionFactor<T>::const_iterator key_it = expression_factor.begin();
  gtsam::VectorValues dX = values.zeroVectors();
  for(; key_it != expression_factor.end(); ++key_it) {
    size_t key = *key_it;
    // Compute central differences using the values struct.
    const size_t cols = dX.dim(key);
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(rows, cols);
    for(size_t col = 0; col < cols; ++col) {
      Eigen::VectorXd dx = Eigen::VectorXd::Zero(cols);
      dx[col] = fd_step;
      dX[key] = dx;
      gtsam::Values eval_values = values.retract(dX);
      Eigen::VectorXd left = expression_factor.unwhitenedError(eval_values);
      dx[col] = -fd_step;
      dX[key] = dx;
      eval_values = values.retract(dX);
      Eigen::VectorXd right = expression_factor.unwhitenedError(eval_values);
      J.col(col) = (left - right) * (1.0/(2.0 * fd_step));
    }
    jacobians[key] = J;
  }

  // Next step...return JacobianFactor
  return gtsam::JacobianFactor(jacobians, -e);
}

template<typename T>
void testExpressionJacobians(gtsam::Expression<T> expression,
                             const gtsam::Values& values,
                             double fd_step,
                             double tolerance) {
  using namespace gtsam;
  // Create factor
  size_t size = traits::dimension<T>::value;
  ExpressionFactor<T> f(noiseModel::Unit::Create(size), expression.value(values), expression);

  // Check linearization
  //JacobianFactor expected(1, eye(3,6), 2, Tval.getRotationMatrix(), -(Tval * pval));
  JacobianFactor expected = computeFiniteDifferenceJacobians(f, values, fd_step);
  boost::shared_ptr<GaussianFactor> gf = f.linearize(values);
  boost::shared_ptr<JacobianFactor> jf = //
      boost::dynamic_pointer_cast<JacobianFactor>(gf);

  typedef std::pair<Eigen::MatrixXd, Eigen::VectorXd> Jacobian;
  Jacobian evalJ = jf->jacobianUnweighted();
  Jacobian estJ = expected.jacobianUnweighted();
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(evalJ.first, estJ.first, tolerance))
    << "Mismatch in the Jacobian matrix.";
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(evalJ.second, Eigen::VectorXd::Zero(size), tolerance))
    << "Mismatch in the error vector.";
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(estJ.second, Eigen::VectorXd::Zero(size), tolerance))
    << "Mismatch in the error vector.";

}


// Do a full concept check and test the invertibility of
// local() vs. retract().
template<typename T>
void testDefaultChart(const T& value) {
  T other = value;
  gtsam::traits::print<T>()(value, "value");
  gtsam::traits::print<T>()(other, "other");
  EXPECT_TRUE(gtsam::traits::equals<T>()(value, other, 1e-12));

  typedef typename gtsam::DefaultChart<T> Chart;
  typedef typename Chart::vector Vector;

  Vector dx = Chart::local(value, other);
  EXPECT_EQ(Chart::getDimension(value), dx.size());
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Eigen::VectorXd::Zero(dx.size()), dx, 1e-9));

  dx.setRandom();
  T updated = Chart::retract(value, dx);
  Vector invdx = Chart::local(value, updated);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(dx, invdx, 1e-9));

  dx = -dx;
  updated = Chart::retract(value, dx);
  invdx = Chart::local(value, updated);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(dx, invdx, 1e-9));
}


typedef kindr::minimal::QuatTransformation Transformation;
typedef kindr::minimal::RotationQuaternion Quaternion;

TEST(MinkindrGtsamTests, testRot3Expression) {

  gtsam::Rot3 Cval(gtsam::Rot3::RzRyRx(0.1,0.2,0.3));

  // Create some values
  gtsam::Values values;
  values.insert(1, Cval);

  gtsam::Expression<gtsam::Rot3> C(1);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(C, values, fd_step, tolerance);
}

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
  testDefaultChart(Cval);
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

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, false);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  return RUN_ALL_TESTS();
}
