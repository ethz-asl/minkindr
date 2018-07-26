#include <glog/logging.h>
#include <kindr/minimal/quat-transformation.h>
#include <numpy_eigen/boost_python_headers.hpp>

using namespace boost::python;

typedef kindr::minimal::QuatTransformationTemplate<double> Transformation;

Eigen::Vector3d getPosition(const Transformation* transformation) {
  return CHECK_NOTNULL(transformation)->getPosition();
}

void exportTransformation() {
  using namespace boost::python;

  class_< Transformation >( "Transformation", init<>() )
    .def(init<const Eigen::Matrix4d&>())
    .def("getTransformationMatrix", &Transformation::getTransformationMatrix)
    .def("getRotationMatrix", &Transformation::getRotationMatrix)
    .def("getPosition", getPosition)
    ;
}
