cmake_minimum_required(VERSION 2.8.3)
project(minkindr_python)

find_package(catkin_simple REQUIRED)
catkin_simple()

add_definitions(-std=c++0x -D__STRICT_ANSI__)

# Set up the python exports.
SET(PY_PROJECT_NAME minkindr_python)
SET(PY_PACKAGE_DIR python/minkindr)

add_python_export_library(${PY_PROJECT_NAME} ${PY_PACKAGE_DIR}
  src/module.cpp
  src/export-rotation-quaternion.cc
  src/export-transformation.cc
)

cs_install()
cs_export()