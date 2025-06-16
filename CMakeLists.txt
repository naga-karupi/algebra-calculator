project(algebra CXX)
cmake_minimum_required(VERSION 3.10)

add_executable(determinant src/det.cpp)
add_executable(eigenvalue src/eigen-value.cpp)
add_executable(multiplication src/multi-mat.cpp)

target_compile_features(determinant PRIVATE cxx_std_20)
target_compile_features(eigenvalue PRIVATE cxx_std_20)
target_compile_features(multiplication PRIVATE cxx_std_20)

target_compile_options(determinant PRIVATE -Wall -Wextra -pedantic -O2)
target_compile_options(eigenvalue PRIVATE -Wall -Wextra -pedantic -O2)
target_compile_options(multiplication PRIVATE -Wall -Wextra -pedantic -O2)

target_include_directories(determinant PRIVATE /usr/include/eigen3 inc)
target_include_directories(eigenvalue PRIVATE /usr/include/eigen3 inc)
target_include_directories(multiplication PRIVATE /usr/include/eigen3 inc)

enable_testing()
find_package(GTest REQUIRED)
add_executable(test_matrix test/test.cpp)

target_link_libraries(test_matrix GTest::GTest GTest::Main)
target_compile_features(test_matrix PRIVATE cxx_std_20)
include(GoogleTest)
