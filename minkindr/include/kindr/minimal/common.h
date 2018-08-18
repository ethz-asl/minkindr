// Copyright (c) 2015, Autonomous Systems Lab, ETH Zurich
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the <organization> nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#ifndef MINKINDR_MINIMAL_COMMON_H
#define MINKINDR_MINIMAL_COMMON_H

#include <math.h>
#include <Eigen/Core>
#include <glog/logging.h>

namespace kindr {
namespace minimal {

inline void skewMatrix(const Eigen::Vector3d& v, Eigen::Matrix3d * skew) {
  CHECK_NOTNULL(skew);
  skew->setZero();
  (*skew)(0,1) = -v[2];
  (*skew)(1,0) =  v[2];
  (*skew)(0,2) =  v[1];
  (*skew)(2,0) = -v[1];
  (*skew)(1,2) = -v[0];
  (*skew)(2,1) =  v[0];
}

inline void skewMatrix(const Eigen::Vector3d& v, Eigen::Map<Eigen::Matrix3d> * skew) {
  CHECK_NOTNULL(skew);
  skew->setZero();
  (*skew)(0,1) = -v[2];
  (*skew)(1,0) =  v[2];
  (*skew)(0,2) =  v[1];
  (*skew)(2,0) = -v[1];
  (*skew)(1,2) = -v[0];
  (*skew)(2,1) =  v[0];
}

inline Eigen::Matrix3d skewMatrix(const Eigen::Vector3d& v) {
  Eigen::Matrix3d skew;
  skewMatrix(v, &skew);
  return skew;
}

inline Eigen::Matrix3d eulerAnglesYawPitchRollToRotationMatrix(const double& yaw, const double& pitch, const double& roll) {
    Eigen::Matrix3d C;
    double cx = cos(roll);
    double sx = sin(roll);
    double cy = cos(pitch);
    double sy = sin(pitch);
    double cz = cos(yaw);
    double sz = sin(yaw);
      //[cos(z)*cos(y), -sin(z)*cos(x)+cos(z)*sin(y)*sin(x),  sin(z)*sin(x)+cos(z)*sin(y)*cos(x)]
      //[sin(z)*cos(y),  cos(z)*cos(x)+sin(z)*sin(y)*sin(x), -cos(z)*sin(x)+sin(z)*sin(y)*cos(x)]
      //[      -sin(y),                       cos(y)*sin(x),                       cos(y)*cos(x)]
      C <<
	cz*cy,  -sz*cx+cz*sy*sx,   sz*sx+cz*sy*cx,
	sz*cy,   cz*cx+sz*sy*sx,  -cz*sx+sz*sy*cx,
        -sy,         cy*sx,            cy*cx;

      return C;
    }

}  // namespace minimal
}  // namespace kindr

#endif // MINKINDR_MINIMAL_COMMON_H
