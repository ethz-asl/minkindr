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

void exportTransformation() {
  using namespace boost::python;

  class_< Transformation, boost::shared_ptr<Transformation>>("Transformation", init<>())
    .def(init<const Eigen::Matrix4d&>())
    .def(init<const Transformation::Rotation&, const Transformation::Position&>())
    .def("getTransformationMatrix", &Transformation::getTransformationMatrix)
    .def("getRotation", getRotation)
    .def("getRotationMatrix", &Transformation::getRotationMatrix)
    .def("getPosition", getPosition)
    .def("inverse", &Transformation::inverse)
    .def(self * self)
    ;

  def("interpolateLinearly", &kindr::minimal::interpolateComponentwise<double>);
}
