cmake_minimum_required(VERSION 3.16)
project(boris_pusher CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-O2 -Wall -Wextra)

# C++: pusher + trajectory drivers; tests in trajectories.ipynb.
add_library(boris_lib src/boris.cpp)
target_include_directories(boris_lib PUBLIC src)

add_executable(boris_traj tools/boris_traj.cpp)
target_link_libraries(boris_traj PRIVATE boris_lib)

add_executable(boris_traj_n tools/boris_traj_n.cpp)
target_link_libraries(boris_traj_n PRIVATE boris_lib)

add_executable(boris_traj_field tools/boris_traj_field.cpp)
target_link_libraries(boris_traj_field PRIVATE boris_lib)

add_executable(boris_traj_grid tools/boris_traj_grid.cpp)
target_link_libraries(boris_traj_grid PRIVATE boris_lib)
