cmake_minimum_required(VERSION 2.8.3)
project(minkindr_gtsam)

find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

# cs_add_library(${PROJECT_NAME} src/dummy.cc)
# cs_add_executable(${PROJECT_NAME} src/main.cc)
# find_package(Boost REQUIRED COMPONENTS system serialization filesystem)

# add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=0)
catkin_add_gtest(${PROJECT_NAME}_tests test/test_minkindr-gtsam.cc)
target_link_libraries(${PROJECT_NAME}_tests ${catkin_LIBRARIES})

cs_install()
cs_export()

