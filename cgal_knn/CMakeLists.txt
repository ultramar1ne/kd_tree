cmake_minimum_required(VERSION 3.1...3.22)
project(Spatial_searching_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)


# include for local directory

# include for local package

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program("knn_3d.cpp")


