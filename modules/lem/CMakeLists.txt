add_library(hs_lem STATIC src/lem.cpp)
target_include_directories(hs_lem PUBLIC include)
target_link_libraries(hs_lem PUBLIC hs_core)
target_compile_features(hs_lem PUBLIC cxx_std_20)
