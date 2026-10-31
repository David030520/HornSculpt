find_package(CGAL 6.2.1 REQUIRED)
find_package(Threads REQUIRED)
add_library(hs_meshing src/mesh_quality.cpp src/mesh_intersections.cpp src/axisymmetric_mesh.cpp
                      src/surface_sizing_field.cpp src/surface_charts.cpp
                      src/frontal_delaunay.cpp src/frontal_mesher.cpp)
target_include_directories(hs_meshing PUBLIC include)
target_link_libraries(hs_meshing PUBLIC hs_core PRIVATE CGAL::CGAL Threads::Threads)
target_compile_features(hs_meshing PUBLIC cxx_std_20)
