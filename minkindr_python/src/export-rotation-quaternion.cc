#include <glog/logging.h>
#include <kindr/minimal/rotation-quaternion.h>
#include <numpy_eigen/boost_python_headers.hpp>

using namespace boost::python;

typedef kindr::minimal::RotationQuaternionTemplate<double> Quaternion;

Eigen::Vector4d getQuaternionXYZW(const Quaternion* quaternion) {
  return CHECK_NOTNULL(quaternion)->toImplementation().coeffs();
}

Quaternion createQuaternionFromXYZW(const Eigen::Vector4d& xyzw) {
  return Quaternion(xyzw(3), xyzw(0), xyzw(1), xyzw(2));
}

Quaternion createQuaternionFromApproximateRotationMatrix(
    const Eigen::Matrix3d& R) {
 return Quaternion::constructAndRenormalize(R);
}

void exportRotationQuaternion() {
  using namespace boost::python;

  class_< Quaternion >( "Quaternion", init<>() )
    .def(init<const Eigen::Matrix3d&>())
    .def(init<const double, const double, const double, const double>())
    .def("w", &Quaternion::w)
    .def("x", &Quaternion::x)
    .def("y", &Quaternion::y)
    .def("z", &Quaternion::z)
    .def("getRotationMatrix", &Quaternion::getRotationMatrix)
    .def("getQuaternionXYZW", getQuaternionXYZW)
    ;

  def("createQuaternionFromXYZW", createQuaternionFromXYZW, "...");
  def("createQuaternionFromApproximateRotationMatrix", createQuaternionFromApproximateRotationMatrix, "...");
}
