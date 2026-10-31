add_library(hs_geometry STATIC
    draft_geometry.cpp
    draft_instances.cpp
    axisymmetric_draft.cpp
    longitudinal_profile.cpp
    shape_transition.cpp
    profiles/osse_profile.cpp
    curves/superellipse.cpp
    curves/rounded_rectangle.cpp
    curves/superformula.cpp
    curves/physical_contour.cpp
    transforms/mouth_morph.cpp
    transforms/diffraction_slot.cpp)
target_include_directories(hs_geometry PUBLIC include)
target_link_libraries(hs_geometry PUBLIC hs_core)
target_compile_features(hs_geometry PUBLIC cxx_std_20)
