include(${PROJECT_SOURCE_DIR}/cmake_helpers.cmake)

# Libfruit provides the machinery for unit testing in Fortran
add_Fortran_library(fruit
    fruit.f90)
