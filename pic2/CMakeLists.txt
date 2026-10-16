cmake_minimum_required(VERSION 3.16)
project(boris_pusher CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-O2 -Wall -Wextra)

# C++: pusher + trajectory drivers + phase-1 PPC (moments); tests in Python.
add_library(boris_lib src/boris.cpp src/ppc.cpp)
target_include_directories(boris_lib PUBLIC src)

add_executable(boris_traj tools/boris_traj.cpp)
target_link_libraries(boris_traj PRIVATE boris_lib)

add_executable(boris_traj_n tools/boris_traj_n.cpp)
target_link_libraries(boris_traj_n PRIVATE boris_lib)

add_executable(boris_traj_field tools/boris_traj_field.cpp)
target_link_libraries(boris_traj_field PRIVATE boris_lib)

add_executable(boris_traj_grid tools/boris_traj_grid.cpp)
target_link_libraries(boris_traj_grid PRIVATE boris_lib)

# Python bindings (phase-1 moments + PPC vocabulary; skipped if no pybind11)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
if(NOT pybind11_DIR)
  execute_process(COMMAND ${Python_EXECUTABLE} -m pybind11 --cmakedir
                  OUTPUT_VARIABLE pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET)
endif()
find_package(pybind11 CONFIG QUIET)
if(pybind11_FOUND)
  pybind11_add_module(pic2_ext bindings/pic2_ext.cpp)
  target_link_libraries(pic2_ext PRIVATE boris_lib)
  target_link_options(pic2_ext PRIVATE -static-libstdc++ -static-libgcc)
else()
  message(STATUS "pybind11 not found: skipping pic2_ext (pip install pybind11)")
endif()
