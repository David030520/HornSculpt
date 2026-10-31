find_package(Qt6 6.8 REQUIRED COMPONENTS Core Network Concurrent)
add_library(hs_freecad STATIC
    freecad_link.cpp
    cad_snapshot.cpp
    cad_validation.cpp
    include/hornsculpt/freecad_link.hpp)
set_target_properties(hs_freecad PROPERTIES AUTOMOC ON)
target_include_directories(hs_freecad PUBLIC include)
target_link_libraries(hs_freecad PUBLIC hs_core hs_cad Qt6::Core PRIVATE hs_meshing Qt6::Network Qt6::Concurrent)
target_compile_features(hs_freecad PUBLIC cxx_std_20)
install(FILES bridge.py workbench/Init.py workbench/InitGui.py workbench/hornsculpt_commands.py
    DESTINATION bin/integrations/freecad/workbench)
