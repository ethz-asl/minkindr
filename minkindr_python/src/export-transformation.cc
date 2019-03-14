#include <glog/logging.h>
#include <kindr/minimal/quat-transformation.h>
#include <numpy_eigen/boost_python_headers.hpp>

using namespace boost::python;

using Transformation = kindr::minimal::QuatTransformationTemplate<double>;

Eigen::Vector3d getPosition(const Transformation* transformation) {
  return CHECK_NOTNULL(transformation)->getPosition();
}

kindr::minimal::RotationQuaternionTemplate<double> getRotation(const Transformation* transformation) {
  return CHECK_NOTNULL(transformation)->getRotation();
}

Transformation interpolate(
    const Transformation& T_A, const int64_t time_a_ns,
    const Transformation& T_B, const int64_t time_b_ns,
    const int64_t target_time_ns) {
  CHECK_LE(time_a_ns, target_time_ns);
  CHECK_LE(target_time_ns, time_b_ns);
  CHECK_NE(time_a_ns, time_b_ns);

  if (target_time_ns == time_a_ns) {
    return T_A;
  }

  if (target_time_ns == time_b_ns) {
    return T_B;
  }

  const double f = static_cast<double>(target_time_ns - time_a_ns) /
      static_cast<double>(time_b_ns - time_a_ns);
  CHECK_GT(f, 0.0);
  CHECK_LT(f, 1.0);

  const Eigen::Vector3d p_int =
      T_A.getPosition() + f * (T_B.getPosition() - T_A.getPosition());

  const Eigen::Quaterniond q_int =
      T_A.getEigenQuaternion().slerp(f, T_B.getEigenQuaternion());

  return Transformation(q_int, p_int);
}

void exportTransformation() {
  using namespace boost::python;

  class_< Transformation, boost::shared_ptr<Transformation> >( "Transformation", init<>() )
    .def(init<const Eigen::Matrix4d&>())
    .def(init<const Transformation::Rotation&, const Transformation::Position&>())
    .def("getTransformationMatrix", &Transformation::getTransformationMatrix)
    .def("getRotation", getRotation)
    .def("getRotationMatrix", &Transformation::getRotationMatrix)
    .def("getPosition", getPosition)
    .def("inverse", &Transformation::inverse)
    .def(self * self)
    ;

  def("interpolate", &interpolate);
}
