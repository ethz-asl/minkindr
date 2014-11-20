#include <gtest/gtest.h>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <kindr/minimal/quat-transformation-gtsam.h>
#include <gtsam_unstable/nonlinear/Expression.h>
#include <gtsam_unstable/nonlinear/ExpressionFactor.h>

typedef kindr::minimal::QuatTransformation Transformation;

TEST(MinkindrGtsamTests, testThis) {
  using namespace gtsam;
  Expression<Transformation> T(1);
  Expression<Eigen::Vector3d> p(2);

  Transformation Tval;
  Eigen::Vector3d pval;
  Tval.setRandom();
  pval.setRandom();

  // Create some values
  Values values;
  values.insert(1, Tval);
  values.insert(2, pval);

  Expression<Eigen::Vector3d> Tp = T * p;
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, false);
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  return RUN_ALL_TESTS();
}
