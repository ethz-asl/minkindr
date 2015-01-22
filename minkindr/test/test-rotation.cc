#include <gtest/gtest.h>
#include <kindr/minimal/rotation-quaternion.h>
#include <kindr/minimal/angle-axis.h>
#include <cmath>
#include <eigen-checks/gtest.h>

#ifndef TEST
#define TEST(a, b) int Test_##a##_##b()
#endif

TEST(MinKindrTests,testQuatAxisAngle) {
  using namespace kindr::minimal;
  Eigen::Vector4d v(0.64491714, 0.26382416,  0.51605132,  0.49816637);
  Eigen::Matrix3d C;

  // This is data generated from a trusted python implementation.
  C << -0.02895739, -0.37025845,  0.92847733,
        0.91484567,  0.36445417,  0.17386938,
       -0.40276404,  0.85444827,  0.3281757;

  Eigen::Vector4d aax(1.7397629325686206, 0.34520549,  0.67523668,  0.65183479);

  RotationQuaternion q1(v[0], v[1], v[2], v[3]);
  RotationQuaternion q2(C);
  RotationQuaternion q3(v[0],v.tail<3>());
  RotationQuaternion q4(v[0], v[1], v[2], v[3]);
  RotationQuaternion q5(-v[0], -v[1], -v[2], -v[3]);


  AngleAxis a1(aax[0],aax.tail<3>()); 
  AngleAxis a2(aax[0] + 2*M_PI,aax.tail<3>()); 
  AngleAxis a3(aax[0] + 4*M_PI,aax.tail<3>()); 
  AngleAxis a4(q1);
  RotationQuaternion q6(a1);

  EXPECT_NEAR(q1.getDisparityAngle(q2), 0.0, 1e-3);
  EXPECT_NEAR(q1.getDisparityAngle(q3), 0.0, 1e-3);
  EXPECT_NEAR(q1.getDisparityAngle(q4), 0.0, 1e-3);
  EXPECT_NEAR(q1.getDisparityAngle(q5), 0.0, 1e-3);
  EXPECT_NEAR(q1.getDisparityAngle(q6), 0.0, 1e-3);
  EXPECT_NEAR(q1.getDisparityAngle(a1), 0.0, 1e-3);
  EXPECT_NEAR(q1.getDisparityAngle(a2), 0.0, 1e-3);
  EXPECT_NEAR(q1.getDisparityAngle(a2), 0.0, 1e-3);
  EXPECT_NEAR(q1.getDisparityAngle(a3), 0.0, 1e-3);
  EXPECT_NEAR(q1.getDisparityAngle(a4), 0.0, 1e-3);
  EXPECT_NEAR(a1.getDisparityAngle(q2), 0.0, 1e-3);
  EXPECT_NEAR(a1.getDisparityAngle(q3), 0.0, 1e-3);
  EXPECT_NEAR(a1.getDisparityAngle(q4), 0.0, 1e-3);
  EXPECT_NEAR(a1.getDisparityAngle(q5), 0.0, 1e-3);
  EXPECT_NEAR(a1.getDisparityAngle(q6), 0.0, 1e-3);
  EXPECT_NEAR(a1.getDisparityAngle(q1), 0.0, 1e-3);
  EXPECT_NEAR(a1.getDisparityAngle(a2), 0.0, 1e-3);
  EXPECT_NEAR(a1.getDisparityAngle(a2), 0.0, 1e-3);
  EXPECT_NEAR(a1.getDisparityAngle(a3), 0.0, 1e-3);
  EXPECT_NEAR(a1.getDisparityAngle(a4), 0.0, 1e-3);

  EXPECT_NEAR(q1.inverted().getDisparityAngle(q2.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(q1.inverted().getDisparityAngle(q3.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(q1.inverted().getDisparityAngle(q4.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(q1.inverted().getDisparityAngle(q5.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(q1.inverted().getDisparityAngle(q6.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(q1.inverted().getDisparityAngle(a1.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(q1.inverted().getDisparityAngle(a2.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(q1.inverted().getDisparityAngle(a2.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(q1.inverted().getDisparityAngle(a3.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(q1.inverted().getDisparityAngle(a4.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(a1.inverted().getDisparityAngle(q2.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(a1.inverted().getDisparityAngle(q3.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(a1.inverted().getDisparityAngle(q4.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(a1.inverted().getDisparityAngle(q5.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(a1.inverted().getDisparityAngle(q6.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(a1.inverted().getDisparityAngle(q1.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(a1.inverted().getDisparityAngle(a2.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(a1.inverted().getDisparityAngle(a2.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(a1.inverted().getDisparityAngle(a3.inverted()), 0.0, 1e-3);
  EXPECT_NEAR(a1.inverted().getDisparityAngle(a4.inverted()), 0.0, 1e-3);


}

TEST(MinKindrTests,testComposition) {
  using namespace kindr::minimal;
  Eigen::Vector4d v(0.64491714, 0.26382416,  0.51605132,  0.49816637);
  Eigen::Matrix3d C, Csquared;

  // This is data generated from a trusted python implementation.
  C << -0.02895739, -0.37025845,  0.92847733,
        0.91484567,  0.36445417,  0.17386938,
       -0.40276404,  0.85444827,  0.3281757;

  Csquared << -0.71184809,  0.66911533,  0.21344081,
      0.23689944, -0.05734011,  0.96984059,
      0.66117392,  0.74094318, -0.1176956;

  Eigen::Vector4d aax(1.7397629325686206, 0.34520549,  0.67523668,  0.65183479);
  
  
  RotationQuaternion q1(v[0], v[1], v[2], v[3]);
  AngleAxis a1(aax[0],aax.tail<3>()); 


  RotationQuaternion qsquared(Csquared); 
  
  EXPECT_NEAR((q1*q1).getDisparityAngle(qsquared), 0.0, 1e-3);
  EXPECT_NEAR((q1*a1).getDisparityAngle(qsquared), 0.0, 1e-3);
  EXPECT_NEAR((a1*q1).getDisparityAngle(qsquared), 0.0, 1e-3);
  EXPECT_NEAR((a1*a1).getDisparityAngle(qsquared), 0.0, 1e-3);
  EXPECT_NEAR((a1.inverted()*a1.inverted()).getDisparityAngle(qsquared.inverted()), 0.0, 1e-3);
  EXPECT_NEAR((q1.inverted()*a1.inverted()).getDisparityAngle(qsquared.inverted()), 0.0, 1e-3);
  EXPECT_NEAR((a1.inverted()*q1.inverted()).getDisparityAngle(qsquared.inverted()), 0.0, 1e-3);
  EXPECT_NEAR((q1.inverted()*q1.inverted()).getDisparityAngle(qsquared.inverted()), 0.0, 1e-3);

}

TEST(MinKindrTests, testQuaternionInitialization) {
 Eigen::Vector4d q_coeffs;
 q_coeffs << 1, 0, 0, 0;
 kindr::minimal::RotationQuaternionTemplate<double> q_from_coeffs(
     q_coeffs[0], q_coeffs[1], q_coeffs[2], q_coeffs[3]);
 Eigen::Quaterniond q;
 q.setIdentity();

 for(int i = 0; i < 4; ++i) {
    EXPECT_NEAR(q_from_coeffs.toImplementation().coeffs()[i],
                q.coeffs()[i], 1e-10);
 }
}

TEST(MinKindrTests, testRotate) {
  using namespace kindr::minimal;
  Eigen::Vector4d q(0.64491714, 0.26382416,  0.51605132,  0.49816637);
  Eigen::Matrix3d C, Csquared;

  // This is data generated from a trusted python implementation.
  C << -0.02895739, -0.37025845,  0.92847733,
        0.91484567,  0.36445417,  0.17386938,
       -0.40276404,  0.85444827,  0.3281757;

  Eigen::Vector4d aax(1.7397629325686206, 0.34520549,  0.67523668,  0.65183479);
  
  
  RotationQuaternion q1(q[0], q[1], q[2], q[3]);
  AngleAxis a1(aax[0],aax.tail<3>()); 

  Eigen::Vector3d v(4.67833851,  8.52053031,  6.71796159);
  Eigen::Vector3d Cv(2.94720425,  8.55334831,  7.60075758);
  Eigen::Vector3d Ctv(4.95374446,  7.11329907,  8.02986227);
  Eigen::Vector3d Cv1 = q1.rotate(v);
  Eigen::Vector3d Cv2 = a1.rotate(v);
  Eigen::Vector3d Ctv1 = q1.inverseRotate(v);
  Eigen::Vector3d Ctv2 = a1.inverseRotate(v);
  Eigen::Vector3d Ctv3 = q1.inverted().rotate(v);
  Eigen::Vector3d Ctv4 = a1.inverted().rotate(v);


  for(int i = 0; i < 3; ++i) {
     EXPECT_NEAR(Cv[i], Cv1[i], 1e-4);
     EXPECT_NEAR(Cv[i], Cv2[i], 1e-4);
     EXPECT_NEAR(Ctv[i], Ctv1[i], 1e-4);
     EXPECT_NEAR(Ctv[i], Ctv2[i], 1e-4);
     EXPECT_NEAR(Ctv[i], Ctv3[i], 1e-4);
     EXPECT_NEAR(Ctv[i], Ctv4[i], 1e-4);
  }
}

TEST(MinKindrTests, testRotationExpLog) {
  using namespace kindr::minimal;
  for(int i = 0; i < 10; ++i) {
    RotationQuaternion C1;
    C1.setRandom();
    RotationQuaternion::Vector3 v = C1.log();
    RotationQuaternion C2 = RotationQuaternion::exp(v);
    Eigen::Matrix3d CC1 = C1.getRotationMatrix();
    Eigen::Matrix3d CC2 = C2.getRotationMatrix();
    for(int r = 0; r < 3; ++r) {
      for(int c = 0; c < 3; ++c) {
        EXPECT_NEAR(CC1(r,c), CC2(r,c), 1e-6) << "Failed at (" << r << "," << c << ")";
      }
    }
  }

  RotationQuaternion::Vector3 axis;
  axis << 0.0,0.0,1.0;
  axis /= axis.norm();
  for (double angle = -M_PI; angle <= M_PI; angle+=M_PI/100) {
    RotationQuaternion::Vector3 axisAngle = axis*angle;
    RotationQuaternion C2 = RotationQuaternion::exp(axisAngle);
    RotationQuaternion::Vector3 v = C2.log();
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(v,axisAngle,1e-6));
  }

  Eigen::Matrix<double,4,1> q;
  for(int i=0; i < 1000; ++i) {
    RotationQuaternion C1;
    C1.setRandom();
    RotationQuaternion::Vector3 v1 = C1.log();
    C1.setParts(-C1.w(), -C1.imaginary());
    RotationQuaternion::Vector3 v2 = C1.log();
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(v1,v2,1e-6));
  }
}
