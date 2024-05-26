cmake_minimum_required(VERSION 3.28)
project(hmsh)

set(CMAKE_CXX_STANDARD 17)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(eigen)

add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
target_include_directories(${PROJECT_NAME} INTERFACE ./include/hmsh)
