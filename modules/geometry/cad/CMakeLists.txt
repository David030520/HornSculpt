# The SDK and its runtime directory are explicit inputs; no dependency downloads.
set(HORNSCULPT_OCCT_ROOT "" CACHE PATH "Open CASCADE 7.9.3 SDK")
set(HORNSCULPT_OCCT_RUNTIME_DIR "" CACHE PATH "Open CASCADE third-party runtime DLLs")
find_path(HS_OCCT_INCLUDE Standard_Version.hxx HINTS "${HORNSCULPT_OCCT_ROOT}/inc" REQUIRED)
file(STRINGS "${HS_OCCT_INCLUDE}/Standard_Version.hxx" hs_occt_version
  REGEX "^#define OCC_VERSION_COMPLETE ")
if(NOT hs_occt_version MATCHES "\"7\\.9\\.3\"")
  message(FATAL_ERROR "HornSculpt requires the pinned Open CASCADE SDK 7.9.3")
endif()
set(hs_occt_libraries TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo
  TKTopAlgo TKShHealing TKMesh TKPrim TKBO TKDESTEP TKXSBase TKDE)
set(HS_CAD_RUNTIME_FILES "")
foreach(name IN LISTS hs_occt_libraries)
  find_library(HS_OCCT_${name}_LIBRARY NAMES ${name}
    HINTS "${HORNSCULPT_OCCT_ROOT}/win64/vc14/lib" REQUIRED)
  add_library(hs_occt_${name} SHARED IMPORTED GLOBAL)
  set_target_properties(hs_occt_${name} PROPERTIES
    IMPORTED_IMPLIB "${HS_OCCT_${name}_LIBRARY}"
    IMPORTED_LOCATION "${HORNSCULPT_OCCT_ROOT}/win64/vc14/bin/${name}.dll")
endforeach()
# Runtime closure of the pinned Windows STEP/meshing toolkits. STEP's XDE
# toolkit imports visualization support even for headless shape transfer.
set(hs_occt_runtime TKernel TKMath TKG2d TKG3d TKGeomBase TKBRep TKGeomAlgo TKTopAlgo
  TKShHealing TKMesh TKPrim TKBO TKDESTEP TKXSBase TKDE TKCAF TKCDF TKLCAF TKHLR
  TKService TKV3d TKVCAF TKXCAF)
foreach(name IN LISTS hs_occt_runtime)
  list(APPEND HS_CAD_RUNTIME_FILES "${HORNSCULPT_OCCT_ROOT}/win64/vc14/bin/${name}.dll")
endforeach()
foreach(name tbb12 jemalloc avcodec-57 avformat-57 avutil-55 swscale-4 FreeImage freetype openvr_api)
  find_file(HS_CAD_${name}_RUNTIME NAMES "${name}.dll"
    HINTS "${HORNSCULPT_OCCT_RUNTIME_DIR}" REQUIRED)
  list(APPEND HS_CAD_RUNTIME_FILES "${HS_CAD_${name}_RUNTIME}")
endforeach()
set(HS_CAD_RUNTIME_FILES "${HS_CAD_RUNTIME_FILES}" PARENT_SCOPE)
function(hs_deploy_cad_runtime target)
  add_custom_command(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E
    copy_if_different ${HS_CAD_RUNTIME_FILES} "$<TARGET_FILE_DIR:${target}>"
    COMMAND_EXPAND_LISTS VERBATIM)
endfunction()
add_library(hs_cad STATIC cad_geometry.cpp cad_axisymmetric.cpp cad_boundary_corners.cpp cad_parameterization.cpp
    cad_revolved_projection.cpp draft_step_export.cpp)
target_include_directories(hs_cad PUBLIC ../include)
target_include_directories(hs_cad SYSTEM PRIVATE "${HS_OCCT_INCLUDE}")
target_link_libraries(hs_cad PUBLIC hs_core hs_geometry)
foreach(name IN LISTS hs_occt_libraries)
  target_link_libraries(hs_cad PRIVATE hs_occt_${name})
endforeach()

if(BUILD_TESTING)
  add_executable(draft_step_export_tests ../../../tests/draft_step_export_tests.cpp)
  target_include_directories(draft_step_export_tests SYSTEM PRIVATE "${HS_OCCT_INCLUDE}")
  target_link_libraries(draft_step_export_tests PRIVATE hs_cad hs_occt_TKDESTEP hs_occt_TKBRep
    hs_occt_TKTopAlgo hs_occt_TKG3d hs_occt_TKMath hs_occt_TKernel)
  hs_deploy_cad_runtime(draft_step_export_tests)
  add_test(NAME draft_step_export COMMAND draft_step_export_tests)
  set_tests_properties(draft_step_export PROPERTIES TIMEOUT 180)
endif()
