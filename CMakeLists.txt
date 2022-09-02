cmake_minimum_required(VERSION 3.0.2)
project(obvp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-declarations")

file(GLOB INCLUDE_FILES "include/*.h")
file(GLOB SRC_FILES "src/*.cpp")
# file(COPY "src/parameters.yaml" DESTINATION ${CMAKE_BINARY_DIR})

find_package(PythonLibs REQUIRED)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${NLOPT_INCLUDE_DIR}
    ${PYTHON_INCLUDE_DIRS}
    include/PX4-Matrix/matrix
    include/mathlib
    include/eigen
    include/matplotlib-cpp
)

add_subdirectory(yaml-cpp)

# add_executable(${PROJECT_NAME}_bvp
#     src/bvp_test.cpp
#     src/geo.cpp
# )

add_executable(${PROJECT_NAME}_precision_landing
    src/landing_test.cpp
    src/geo.cpp
)

# add_executable(${PROJECT_NAME}_precision_landing_eigen
#     src/landing_eigen_test.cpp
#     src/geo.cpp
# )

add_executable(${PROJECT_NAME}_opt_landing
    src/opt_landing.cpp
    src/geo.cpp
)
target_link_libraries(${PROJECT_NAME}_opt_landing 
    ${PYTHON_LIBRARIES}
    yaml-cpp
    nlopt
)

add_library(${PROJECT_NAME} STATIC ${SRC_FILES} ${INCLUDE_FILES})
target_include_directories (${PROJECT_NAME} 
    PUBLIC 
        ${PROJECT_SOURCE_DIR} 
        ${PROJECT_SOURCE_DIR}/include
        yaml-cpp
)
# https://stackoverflow.com/a/34800230
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_SOURCE_DIR}/src/parameters.yaml
            ${CMAKE_CURRENT_BINARY_DIR}/parameters.yaml)