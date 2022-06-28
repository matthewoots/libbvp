cmake_minimum_required(VERSION 3.0.2)
project(obvp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-declarations")

file(GLOB INCLUDE_FILES "include/*.h")
file(GLOB SRC_FILES "src/*.cpp")

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    PX4-Matrix/matrix
    mathlib
    eigen
)

# add_executable(${PROJECT_NAME}_bvp
#     src/bvp_test.cpp
#     src/geo.cpp
# )

# add_executable(${PROJECT_NAME}_precision_landing
#     src/landing_test.cpp
#     src/geo.cpp
# )

add_executable(${PROJECT_NAME}_precision_landing_eigen
    src/landing_eigen_test.cpp
    src/geo.cpp
)

add_library(${PROJECT_NAME} STATIC ${SRC_FILES} ${INCLUDE_FILES})
target_include_directories (${PROJECT_NAME} PUBLIC 
    ${PROJECT_SOURCE_DIR} 
    ${PROJECT_SOURCE_DIR}/include)