#ifndef MINKINDR_MINIMAL_COMMON_H
#define MINKINDR_MINIMAL_COMMON_H

#include <Eigen/Core>
#include <glog/logging.h>

namespace kindr {
namespace minimal {

void skewMatrix(const Eigen::Vector3d& v, Eigen::Matrix3d * skew) {
  CHECK_NOTNULL(skew);
  skew->setZero();
  (*skew)(0,1) = -v[2];
  (*skew)(1,0) =  v[2];
  (*skew)(0,2) =  v[1];
  (*skew)(2,0) = -v[1];
  (*skew)(1,2) = -v[0];
  (*skew)(2,1) =  v[0];
}

Eigen::Matrix3d skewMatrix(const Eigen::Vector3d& v) {
  Eigen::Matrix3d skew;
  skewMatrix(v, &skew);
  return skew;
}

}  // namespace minimal
}  // namespace kindr

#endif // MINKINDR_MINIMAL_COMMON_H
