cmake_minimum_required(VERSION 2.8.3)
project(minkindr)

find_package(catkin_simple REQUIRED)
catkin_simple()

find_package (Threads)

# add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=0)
catkin_add_gtest(${PROJECT_NAME}_tests 
  test/test-main.cc
  test/test-rotation.cc
  test/test-transformation.cc
)

# Add all files to show up in QtCreator.
add_custom_target(${PROJECT_NAME} SOURCES
  include/kindr/minimal/angle-axis.h
  include/kindr/minimal/common.h
  include/kindr/minimal/position.h
  include/kindr/minimal/quat-transformation.h
  include/kindr/minimal/rotation-quaternion.h
  include/kindr/minimal/implementation/angle-axis-inl.h
  include/kindr/minimal/implementation/quat-transformation-inl.h
  include/kindr/minimal/implementation/rotation-quaternion-inl.h
)

target_link_libraries(${PROJECT_NAME}_tests ${catkin_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

cs_install()
cs_export()

