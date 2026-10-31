cmake_minimum_required(VERSION 3.25)
project(HornSculpt VERSION 0.1.60 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CTest)
if(MSVC)
  include(InstallRequiredSystemLibraries)
endif()
option(HORNSCULPT_BUILD_DESKTOP "Build the native Qt desktop application" ON)
option(HORNSCULPT_BUILD_CUDA_BEM "Build the native CUDA BEM solver" ON)

add_library(hs_options INTERFACE)
if(MSVC)
  target_compile_options(hs_options INTERFACE $<$<COMPILE_LANGUAGE:CXX>:/W4;/permissive-;/utf-8;/Zc:__cplusplus>)
else()
  target_compile_options(hs_options INTERFACE -Wall -Wextra -Wpedantic)
endif()

add_library(hs_core INTERFACE)
target_include_directories(hs_core INTERFACE modules/core/include)
target_link_libraries(hs_core INTERFACE hs_options)

add_subdirectory(modules/geometry)
add_subdirectory(modules/geometry/cad)
add_subdirectory(modules/analysis)
add_subdirectory(modules/lem)
add_subdirectory(modules/meshing)

if(HORNSCULPT_BUILD_CUDA_BEM)
  enable_language(CUDA)
  find_package(CUDAToolkit 13.0 REQUIRED)
endif()
add_subdirectory(modules/bem)

install(DIRECTORY licenses/ DESTINATION licenses)

if(HORNSCULPT_BUILD_DESKTOP)
  if(NOT HORNSCULPT_BUILD_CUDA_BEM)
    message(FATAL_ERROR "The desktop application includes the CUDA solver; enable HORNSCULPT_BUILD_CUDA_BEM.")
  endif()
  find_package(Qt6 6.8 REQUIRED COMPONENTS Core Gui Widgets Network Concurrent)
  add_subdirectory(modules/visualization)
  add_subdirectory(integrations/freecad)
  add_subdirectory(apps/desktop)
endif()

if(BUILD_TESTING)
  add_executable(acoustic_partition_tests tests/acoustic_partition_tests.cpp)
  target_link_libraries(acoustic_partition_tests PRIVATE hs_meshing)
  add_test(NAME acoustic_partition COMMAND acoustic_partition_tests)
  add_executable(symmetry_sector_tests tests/symmetry_sector_tests.cpp)
  target_link_libraries(symmetry_sector_tests PRIVATE hs_meshing)
  add_test(NAME symmetry_sector COMMAND symmetry_sector_tests)
  set_tests_properties(acoustic_partition PROPERTIES TIMEOUT 120)
  add_executable(gmsh_mesh_tests tests/gmsh_mesh_tests.cpp)
  target_link_libraries(gmsh_mesh_tests PRIVATE hs_meshing)
  add_test(NAME gmsh_mesh COMMAND gmsh_mesh_tests)
  set_tests_properties(gmsh_mesh PROPERTIES TIMEOUT 30)
  add_executable(lem_tests tests/lem_tests.cpp)
  target_link_libraries(lem_tests PRIVATE hs_lem)
  add_test(NAME lem COMMAND lem_tests)
  add_executable(surface_geometry_tests tests/surface_geometry_tests.cpp)
  target_include_directories(surface_geometry_tests PRIVATE modules/meshing/src)
  target_link_libraries(surface_geometry_tests PRIVATE hs_meshing)
  add_test(NAME surface_geometry COMMAND surface_geometry_tests)
  add_executable(mesh_orientation_tests tests/mesh_orientation_tests.cpp)
  target_link_libraries(mesh_orientation_tests PRIVATE hs_meshing hs_bem_galerkin)
  add_test(NAME mesh_orientation COMMAND mesh_orientation_tests)
  add_executable(cad_geometry_tests tests/cad_geometry_tests.cpp)
  target_include_directories(cad_geometry_tests SYSTEM PRIVATE "${HS_OCCT_INCLUDE}")
  target_link_libraries(cad_geometry_tests PRIVATE hs_cad hs_meshing hs_occt_TKPrim hs_occt_TKBO
    hs_occt_TKDESTEP hs_occt_TKBRep hs_occt_TKMath hs_occt_TKernel)
  add_test(NAME cad_geometry COMMAND cad_geometry_tests)
  add_executable(cad_parameter_tests tests/cad_parameter_tests.cpp)
  target_include_directories(cad_parameter_tests SYSTEM PRIVATE "${HS_OCCT_INCLUDE}")
  target_link_libraries(cad_parameter_tests PRIVATE hs_cad hs_occt_TKPrim hs_occt_TKBO
    hs_occt_TKBRep hs_occt_TKMath hs_occt_TKernel)
  hs_deploy_cad_runtime(cad_parameter_tests)
  add_test(NAME cad_parameters COMMAND cad_parameter_tests)
  set_tests_properties(cad_parameters PROPERTIES TIMEOUT 180)
  add_executable(parameterized_atlas_tests tests/parameterized_atlas_tests.cpp)
  target_include_directories(parameterized_atlas_tests PRIVATE modules/meshing/src)
  target_include_directories(parameterized_atlas_tests SYSTEM PRIVATE "${HS_OCCT_INCLUDE}")
  target_link_libraries(parameterized_atlas_tests PRIVATE hs_cad hs_meshing
    hs_occt_TKPrim hs_occt_TKBO hs_occt_TKBRep hs_occt_TKMath hs_occt_TKernel)
  target_compile_definitions(parameterized_atlas_tests PRIVATE
    HS_SPLINE_PARAMETER_FIXTURE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/spline_boundary_projection.brep"
    HS_CUSP_PARAMETER_FIXTURE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/spline_boundary_cusp.brep")
  hs_deploy_cad_runtime(parameterized_atlas_tests)
  add_test(NAME parameterized_atlas COMMAND parameterized_atlas_tests)
  set_tests_properties(parameterized_atlas PROPERTIES TIMEOUT 60)
  add_custom_command(TARGET cad_geometry_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
    copy_if_different ${HS_CAD_RUNTIME_FILES} "$<TARGET_FILE_DIR:cad_geometry_tests>"
    COMMAND_EXPAND_LISTS VERBATIM)
  add_executable(meshing_corner_tests tests/meshing_corner_tests.cpp)
  target_include_directories(meshing_corner_tests SYSTEM PRIVATE "${HS_OCCT_INCLUDE}")
  target_link_libraries(meshing_corner_tests PRIVATE hs_cad hs_meshing
    hs_occt_TKPrim hs_occt_TKBO hs_occt_TKBRep hs_occt_TKMath hs_occt_TKernel)
  target_compile_definitions(meshing_corner_tests PRIVATE
    HS_CUSP_CORNER_FIXTURE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/spline_boundary_cusp.brep")
  hs_deploy_cad_runtime(meshing_corner_tests)
  add_test(NAME meshing_corners COMMAND meshing_corner_tests)
  set_tests_properties(meshing_corners PROPERTIES TIMEOUT 30)
  add_executable(geometry_tests tests/geometry_tests.cpp)
  add_executable(axisymmetric_geometry_tests tests/axisymmetric_geometry_tests.cpp)
  target_include_directories(axisymmetric_geometry_tests SYSTEM PRIVATE "${HS_OCCT_INCLUDE}")
  target_link_libraries(axisymmetric_geometry_tests PRIVATE hs_cad hs_geometry hs_meshing hs_occt_TKPrim)
  hs_deploy_cad_runtime(axisymmetric_geometry_tests)
  add_test(NAME axisymmetric_geometry COMMAND axisymmetric_geometry_tests)
  set_tests_properties(axisymmetric_geometry PROPERTIES TIMEOUT 90)
  target_link_libraries(geometry_tests PRIVATE hs_geometry hs_analysis)
  add_test(NAME geometry COMMAND geometry_tests)
  add_executable(diffraction_slot_tests tests/diffraction_slot_tests.cpp)
  target_link_libraries(diffraction_slot_tests PRIVATE hs_geometry hs_analysis hs_meshing)
  add_test(NAME diffraction_slot_geometry COMMAND diffraction_slot_tests)
  add_test(NAME diffraction_slot_mesh COMMAND diffraction_slot_tests --mesh)
  add_test(NAME diffraction_slot_default_mesh COMMAND diffraction_slot_tests --default-mesh)
  add_test(NAME diffraction_slot_radii_default_mesh COMMAND diffraction_slot_tests --radii-default-mesh)
  set_tests_properties(diffraction_slot_radii_default_mesh PROPERTIES TIMEOUT 180 LABELS "meshing")
  set_tests_properties(diffraction_slot_default_mesh PROPERTIES TIMEOUT 180 LABELS "meshing")
  set_tests_properties(diffraction_slot_mesh PROPERTIES TIMEOUT 180 LABELS "meshing")
  if(HORNSCULPT_BUILD_CUDA_BEM)
    target_link_libraries(diffraction_slot_tests PRIVATE hs_bem)
    target_compile_definitions(diffraction_slot_tests PRIVATE HORNSCULPT_SLOT_CUDA_TEST)
    add_test(NAME diffraction_slot_bem COMMAND diffraction_slot_tests --bem)
    set_tests_properties(diffraction_slot_bem PROPERTIES TIMEOUT 600 RUN_SERIAL TRUE LABELS "cuda")
    if(WIN32)
      add_custom_command(TARGET diffraction_slot_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:diffraction_slot_tests>"
        COMMAND_EXPAND_LISTS VERBATIM)
    endif()
  endif()
  add_executable(acoustic_analysis_tests tests/acoustic_analysis_tests.cpp)
  target_link_libraries(acoustic_analysis_tests PRIVATE hs_analysis)
  add_test(NAME acoustic_analysis COMMAND acoustic_analysis_tests)
  add_executable(meshing_tests tests/meshing_tests.cpp)
  target_link_libraries(meshing_tests PRIVATE hs_meshing hs_geometry)
  add_test(NAME meshing COMMAND meshing_tests)
  add_executable(draft_meshing_tests tests/draft_meshing_tests.cpp)
  target_link_libraries(draft_meshing_tests PRIVATE hs_geometry hs_meshing)
  add_test(NAME draft_meshing COMMAND draft_meshing_tests)
  add_executable(planar_meshing_tests tests/planar_meshing_tests.cpp)
  target_link_libraries(planar_meshing_tests PRIVATE hs_meshing)
  add_test(NAME planar_meshing COMMAND planar_meshing_tests)
  add_executable(surface_charts_tests tests/surface_charts_tests.cpp)
  target_link_libraries(surface_charts_tests PRIVATE hs_meshing hs_geometry)
  add_test(NAME surface_charts COMMAND surface_charts_tests)
  add_executable(frontal_delaunay_tests tests/frontal_delaunay_tests.cpp)
  target_link_libraries(frontal_delaunay_tests PRIVATE hs_meshing)
  add_test(NAME frontal_delaunay COMMAND frontal_delaunay_tests)
  set_tests_properties(meshing draft_meshing surface_charts frontal_delaunay planar_meshing
                      PROPERTIES TIMEOUT 90)
  add_executable(mesh_embedding_tests tests/mesh_embedding_tests.cpp)
  target_link_libraries(mesh_embedding_tests PRIVATE hs_meshing)
  add_test(NAME mesh_embedding COMMAND mesh_embedding_tests)
  add_executable(galerkin_space_tests tests/galerkin_space_tests.cpp)
  target_link_libraries(galerkin_space_tests PRIVATE hs_bem_galerkin)
  add_test(NAME galerkin_space COMMAND galerkin_space_tests)
  add_executable(axisymmetric_space_tests tests/axisymmetric_space_tests.cpp)
  target_include_directories(axisymmetric_space_tests PRIVATE modules/bem)
  target_link_libraries(axisymmetric_space_tests PRIVATE hs_bem_galerkin hs_meshing)
  add_test(NAME axisymmetric_space COMMAND axisymmetric_space_tests)
  set_tests_properties(axisymmetric_space PROPERTIES LABELS "numerical" TIMEOUT 60)
  add_executable(axisymmetric_kernel_tests tests/axisymmetric_kernel_tests.cpp)
  target_link_libraries(axisymmetric_kernel_tests PRIVATE hs_bem_galerkin)
  add_test(NAME axisymmetric_kernel COMMAND axisymmetric_kernel_tests)
  set_tests_properties(axisymmetric_kernel PROPERTIES LABELS "numerical" TIMEOUT 240)
  add_executable(field_evaluation_tests tests/field_evaluation_tests.cpp)
  target_link_libraries(field_evaluation_tests PRIVATE hs_bem_field hs_bem_galerkin)
  add_test(NAME field_evaluation COMMAND field_evaluation_tests)
  set_tests_properties(field_evaluation PROPERTIES LABELS "numerical" TIMEOUT 60)
  add_executable(bem_support_tree_tests tests/bem_support_tree_tests.cpp)
  target_link_libraries(bem_support_tree_tests PRIVATE hs_bem_galerkin)
  add_test(NAME bem_support_tree COMMAND bem_support_tree_tests)
  set_tests_properties(bem_support_tree PROPERTIES LABELS "numerical" TIMEOUT 120)
  add_executable(bem_far_quadrature_tests tests/bem_far_quadrature_tests.cpp)
  target_link_libraries(bem_far_quadrature_tests PRIVATE hs_bem_galerkin)
  add_test(NAME bem_far_quadrature COMMAND bem_far_quadrature_tests)
  set_tests_properties(bem_far_quadrature PROPERTIES LABELS "numerical" TIMEOUT 120)
  add_executable(bem_hf_reference_tests tests/bem_hf_reference_tests.cpp)
  add_executable(bem_support_near_tests tests/bem_support_near_tests.cpp)
  target_link_libraries(bem_support_near_tests PRIVATE hs_bem_galerkin)
  add_test(NAME bem_support_near COMMAND bem_support_near_tests)
  set_tests_properties(bem_support_near PROPERTIES LABELS "numerical" TIMEOUT 120)
  target_link_libraries(bem_hf_reference_tests PRIVATE hs_bem_galerkin)
  add_test(NAME bem_hf_reference COMMAND bem_hf_reference_tests)
  set_tests_properties(bem_hf_reference PROPERTIES LABELS "numerical" TIMEOUT 120)
  add_executable(scaled_fixture_tests tests/scaled_fixture_tests.cpp)
  target_link_libraries(scaled_fixture_tests PRIVATE hs_meshing)
  add_test(NAME scaled_fixtures COMMAND scaled_fixture_tests)
  add_executable(pair_quadrature_tests tests/pair_quadrature_tests.cpp)
  target_link_libraries(pair_quadrature_tests PRIVATE hs_bem_galerkin)
  add_test(NAME pair_quadrature COMMAND pair_quadrature_tests)
  add_executable(exterior_domain_tests tests/exterior_domain_tests.cpp)
  target_link_libraries(exterior_domain_tests PRIVATE hs_bem_galerkin)
  add_test(NAME exterior_domain COMMAND exterior_domain_tests)
  add_executable(interaction_plan_tests tests/interaction_plan_tests.cpp)
  target_link_libraries(interaction_plan_tests PRIVATE hs_bem_galerkin)
  add_test(NAME interaction_plan COMMAND interaction_plan_tests)
  add_executable(baffled_plan_cache_tests tests/baffled_plan_cache_tests.cpp)
  target_include_directories(baffled_plan_cache_tests PRIVATE modules/bem)
  target_link_libraries(baffled_plan_cache_tests PRIVATE hs_bem_galerkin hs_meshing)
  add_test(NAME baffled_plan_cache COMMAND baffled_plan_cache_tests)
  set_tests_properties(baffled_plan_cache PROPERTIES LABELS "numerical" TIMEOUT 60)
  add_executable(workspace_size_tests tests/workspace_size_tests.cpp)
  target_include_directories(workspace_size_tests PRIVATE modules/bem)
  target_link_libraries(workspace_size_tests PRIVATE hs_core hs_options)
  add_test(NAME workspace_size COMMAND workspace_size_tests)
  if(HORNSCULPT_BUILD_CUDA_BEM)
    add_executable(bem_regular_quadrature_cuda_tests tests/bem_regular_quadrature_cuda_tests.cu)
    target_include_directories(bem_regular_quadrature_cuda_tests PRIVATE modules/bem tests)
    target_link_libraries(bem_regular_quadrature_cuda_tests PRIVATE hs_bem_galerkin CUDA::cudart)
    set_target_properties(bem_regular_quadrature_cuda_tests PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON)
    target_compile_options(bem_regular_quadrature_cuda_tests PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/utf-8,/Zc:__cplusplus,/wd4211>)
    add_test(NAME bem_regular_quadrature_cuda COMMAND bem_regular_quadrature_cuda_tests)
    set_tests_properties(bem_regular_quadrature_cuda PROPERTIES LABELS "cuda;numerical" RUN_SERIAL TRUE SKIP_RETURN_CODE 77 TIMEOUT 120)
    if(WIN32)
      add_custom_command(TARGET bem_regular_quadrature_cuda_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_regular_quadrature_cuda_tests>" VERBATIM)
    endif()
    add_executable(bem_adaptive_patch_cuda_tests tests/bem_adaptive_patch_cuda_tests.cpp)
    target_link_libraries(bem_adaptive_patch_cuda_tests PRIVATE hs_bem hs_bem_galerkin CUDA::cudart)
    add_test(NAME bem_adaptive_patch_cuda COMMAND bem_adaptive_patch_cuda_tests)
    set_tests_properties(bem_adaptive_patch_cuda PROPERTIES LABELS "cuda;numerical" RUN_SERIAL TRUE
      SKIP_RETURN_CODE 77 TIMEOUT 120)
    if(WIN32)
      add_custom_command(TARGET bem_adaptive_patch_cuda_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_adaptive_patch_cuda_tests>" VERBATIM)
    endif()
    add_executable(bem_support_direct_audit_tests tests/bem_support_direct_audit_tests.cpp)
    target_link_libraries(bem_support_direct_audit_tests PRIVATE hs_bem hs_bem_galerkin CUDA::cudart)
    add_test(NAME bem_support_direct_audit COMMAND bem_support_direct_audit_tests)
    set_tests_properties(bem_support_direct_audit PROPERTIES LABELS "cuda;numerical" RUN_SERIAL TRUE
      SKIP_RETURN_CODE 77 TIMEOUT 120)
    if(WIN32)
      add_custom_command(TARGET bem_support_direct_audit_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_support_direct_audit_tests>" VERBATIM)
    endif()
    add_executable(bem_support_fmm_solve_tests tests/bem_support_fmm_solve_tests.cpp)
    target_link_libraries(bem_support_fmm_solve_tests PRIVATE hs_bem hs_bem_galerkin CUDA::cudart)
    add_test(NAME bem_support_fmm_solve COMMAND bem_support_fmm_solve_tests)
    set_tests_properties(bem_support_fmm_solve PROPERTIES LABELS "cuda;numerical" RUN_SERIAL TRUE
      SKIP_RETURN_CODE 77 TIMEOUT 120)
    if(WIN32)
      add_custom_command(TARGET bem_support_fmm_solve_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_support_fmm_solve_tests>" VERBATIM)
    endif()
    add_executable(bem_device_gmres_tests tests/bem_device_gmres_tests.cpp)
    target_link_libraries(bem_device_gmres_tests PRIVATE hs_bem hs_bem_galerkin CUDA::cudart CUDA::cublas)
    add_test(NAME bem_device_gmres COMMAND bem_device_gmres_tests)
    set_tests_properties(bem_device_gmres PROPERTIES LABELS "cuda;numerical" RUN_SERIAL TRUE
      SKIP_RETURN_CODE 77 TIMEOUT 120)
    if(WIN32)
      add_custom_command(TARGET bem_device_gmres_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_device_gmres_tests>" VERBATIM)
    endif()
    add_executable(bem_support_preconditioner_tests tests/bem_support_preconditioner_tests.cpp)
    target_link_libraries(bem_support_preconditioner_tests PRIVATE hs_bem hs_bem_galerkin CUDA::cudart)
    add_test(NAME bem_support_preconditioner COMMAND bem_support_preconditioner_tests)
    set_tests_properties(bem_support_preconditioner PROPERTIES LABELS "cuda;numerical" RUN_SERIAL TRUE
      SKIP_RETURN_CODE 77 TIMEOUT 120)
    if(WIN32)
      add_custom_command(TARGET bem_support_preconditioner_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_support_preconditioner_tests>" VERBATIM)
    endif()
    add_executable(bem_support_near_cuda_tests tests/bem_support_near_cuda_tests.cpp)
    target_link_libraries(bem_support_near_cuda_tests PRIVATE hs_bem hs_bem_galerkin CUDA::cudart)
    add_test(NAME bem_support_near_cuda COMMAND bem_support_near_cuda_tests)
    set_tests_properties(bem_support_near_cuda PROPERTIES LABELS "cuda;numerical" RUN_SERIAL TRUE
      SKIP_RETURN_CODE 77 TIMEOUT 180)
    if(WIN32)
      add_custom_command(TARGET bem_support_near_cuda_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_support_near_cuda_tests>" VERBATIM)
    endif()
    add_executable(bem_wideband_cuda_tests tests/bem_wideband_cuda_tests.cpp tests/bem_scaled_regular_cuda.cu tests/bem_diagonal_test_operator.cu)
    set_target_properties(bem_wideband_cuda_tests PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON)
    target_compile_options(bem_wideband_cuda_tests PRIVATE
      $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/utf-8,/Zc:__cplusplus,/wd4211>)
    target_link_libraries(bem_wideband_cuda_tests PRIVATE hs_bem hs_bem_galerkin CUDA::cudart)
    add_test(NAME bem_wideband_cuda COMMAND bem_wideband_cuda_tests)
    set_tests_properties(bem_wideband_cuda PROPERTIES LABELS "cuda;numerical" RUN_SERIAL TRUE
      SKIP_RETURN_CODE 77 TIMEOUT 120)
    if(WIN32)
      add_custom_command(TARGET bem_wideband_cuda_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_wideband_cuda_tests>" VERBATIM)
    endif()
    add_executable(bem_rcr_queue_cuda_tests tests/bem_rcr_queue_cuda_tests.cpp)
    target_link_libraries(bem_rcr_queue_cuda_tests PRIVATE hs_bem hs_bem_galerkin CUDA::cudart)
    add_test(NAME bem_rcr_queue_cuda COMMAND bem_rcr_queue_cuda_tests)
    set_tests_properties(bem_rcr_queue_cuda PROPERTIES LABELS "cuda;numerical;fmm"
      RUN_SERIAL TRUE SKIP_RETURN_CODE 77 TIMEOUT 120)
    if(WIN32)
      add_custom_command(TARGET bem_rcr_queue_cuda_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES}
        "$<TARGET_FILE_DIR:bem_rcr_queue_cuda_tests>" VERBATIM)
    endif()
    add_executable(bem_fmm_tests tests/bem_fmm_tests.cpp)
    target_link_libraries(bem_fmm_tests PRIVATE hs_bem hs_bem_galerkin CUDA::cudart)
    target_compile_definitions(bem_fmm_tests PRIVATE
      HS_PARITY_REFERENCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/beat-sphere-512.csv")
    foreach(fmm_case IN ITEMS sphere independent-rhs sparse-rhs convergence lifecycle operator lem)
      add_test(NAME bem_fmm_${fmm_case} COMMAND bem_fmm_tests --${fmm_case})
      set_tests_properties(bem_fmm_${fmm_case} PROPERTIES LABELS "cuda;numerical;fmm"
        RUN_SERIAL TRUE TIMEOUT 600)
    endforeach()
    if(WIN32)
      add_custom_command(TARGET bem_fmm_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_fmm_tests>" VERBATIM)
    endif()
    add_executable(bem_small_mesh_parity_tests tests/bem_small_mesh_parity_tests.cpp)
    target_link_libraries(bem_small_mesh_parity_tests PRIVATE hs_bem hs_bem_galerkin CUDA::cudart)
    foreach(parity_case IN ITEMS graded horn wideband)
      add_test(NAME bem_small_mesh_parity_${parity_case}
        COMMAND bem_small_mesh_parity_tests --${parity_case})
      set_tests_properties(bem_small_mesh_parity_${parity_case} PROPERTIES
        LABELS "cuda;numerical;fmm;parity" RUN_SERIAL TRUE SKIP_RETURN_CODE 77 TIMEOUT 600)
    endforeach()
    if(WIN32)
      add_custom_command(TARGET bem_small_mesh_parity_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES}
        "$<TARGET_FILE_DIR:bem_small_mesh_parity_tests>" VERBATIM)
    endif()
    add_executable(bem_singular_port_tests tests/bem_singular_port_tests.cpp)
    target_link_libraries(bem_singular_port_tests PRIVATE hs_bem CUDA::cudart)
    target_compile_definitions(bem_singular_port_tests PRIVATE
      HS_SINGULAR_PORT_REFERENCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/beat-scalene-tetra-singular12.csv")
    add_test(NAME bem_singular_ports COMMAND bem_singular_port_tests)
    set_tests_properties(bem_singular_ports PROPERTIES LABELS "cuda;numerical;fmm" RUN_SERIAL TRUE
      SKIP_RETURN_CODE 77 TIMEOUT 120)
    if(WIN32)
      add_custom_command(TARGET bem_singular_port_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_singular_port_tests>" VERBATIM)
    endif()
    add_executable(bem_cuda_tests tests/bem_cuda_tests.cpp)
    target_link_libraries(bem_cuda_tests PRIVATE hs_bem hs_geometry hs_meshing CUDA::cudart)
    target_include_directories(bem_cuda_tests PRIVATE modules/bem)
    target_compile_definitions(bem_cuda_tests PRIVATE
      HS_PARITY_REFERENCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/tests/data/beat-sphere-512.csv")
    add_test(NAME bem_cuda COMMAND bem_cuda_tests)
    add_test(NAME bem_cuda_compact_quadrature COMMAND bem_cuda_tests --compact-quadrature)
    set_tests_properties(bem_cuda_compact_quadrature PROPERTIES LABELS "cuda;numerical" TIMEOUT 120)
    set_tests_properties(bem_cuda PROPERTIES LABELS "cuda;numerical" TIMEOUT 180)
    add_test(NAME bem_cuda_low_frequency COMMAND bem_cuda_tests --low-frequency)
    set_tests_properties(bem_cuda_low_frequency PROPERTIES LABELS "cuda;numerical" TIMEOUT 180)
    add_test(NAME bem_cuda_far_phase COMMAND bem_cuda_tests --far-phase)
    set_tests_properties(bem_cuda_far_phase PROPERTIES LABELS "cuda;numerical" TIMEOUT 180)
    add_test(NAME bem_cuda_coarse_draft COMMAND bem_cuda_tests --coarse-draft)
    set_tests_properties(bem_cuda_coarse_draft PROPERTIES LABELS "cuda;numerical" TIMEOUT 180)
    add_test(NAME bem_cuda_memory_preflight COMMAND bem_cuda_tests --memory-preflight)
    set_tests_properties(bem_cuda_memory_preflight PROPERTIES LABELS "cuda" TIMEOUT 180)
    add_test(NAME bem_cuda_mapped_quadrature COMMAND bem_cuda_tests --mapped-quadrature)
    set_tests_properties(bem_cuda_mapped_quadrature PROPERTIES LABELS "cuda;numerical" TIMEOUT 180)
    add_test(NAME bem_cuda_tight_residual COMMAND bem_cuda_tests --tight-residual)
    set_tests_properties(bem_cuda_tight_residual PROPERTIES LABELS "cuda;numerical" TIMEOUT 180)
    add_test(NAME bem_cuda_thin_body COMMAND bem_cuda_tests --thin-body)
    set_tests_properties(bem_cuda_thin_body PROPERTIES LABELS "cuda;numerical" TIMEOUT 180)
    add_executable(bem_resonance_tests tests/bem_resonance_tests.cpp)
    target_link_libraries(bem_resonance_tests PRIVATE hs_bem)
    add_test(NAME bem_resonance COMMAND bem_resonance_tests)
    set_tests_properties(bem_resonance PROPERTIES LABELS "cuda;numerical" TIMEOUT 180)
    if(WIN32)
      add_custom_command(TARGET bem_cuda_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_cuda_tests>" VERBATIM)
      add_custom_command(TARGET bem_resonance_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_resonance_tests>" VERBATIM)
    endif()
  endif()
  if(HORNSCULPT_BUILD_DESKTOP)
    find_package(Qt6 6.8 REQUIRED COMPONENTS Test)
    add_executable(lem_editor_tests tests/lem_editor_tests.cpp)
    target_link_libraries(lem_editor_tests PRIVATE hs_desktop_controls Qt6::Test)
    add_test(NAME lem_editor COMMAND lem_editor_tests)
    set_tests_properties(lem_editor PROPERTIES LABELS "desktop" TIMEOUT 30
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    if(WIN32)
      add_custom_command(TARGET lem_editor_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different $<TARGET_RUNTIME_DLLS:lem_editor_tests> "$<TARGET_FILE_DIR:lem_editor_tests>"
        COMMAND_EXPAND_LISTS VERBATIM)
    endif()
    add_executable(realtime_simulation_tests tests/realtime_simulation_tests.cpp
      apps/desktop/angular_response_widget.cpp
      apps/desktop/phase_response_widget.cpp apps/desktop/acoustic_field_widget.cpp
      apps/desktop/simulation_dock.cpp apps/desktop/acoustic_frame_dialog.cpp apps/desktop/baffle_dialog.cpp apps/desktop/live_trace.cpp apps/desktop/main_window.cpp
      apps/desktop/parameter_panel.cpp apps/desktop/freecad_dock.cpp apps/desktop/mesh_panel.cpp
      apps/desktop/edit_history.cpp apps/desktop/parameter_sweep.cpp)
    add_executable(angular_response_widget_tests tests/angular_response_widget_tests.cpp
      apps/desktop/angular_response_widget.cpp)
    add_executable(field_phase_tests tests/field_phase_tests.cpp apps/desktop/phase_response_widget.cpp
      apps/desktop/acoustic_field_widget.cpp)
    target_include_directories(field_phase_tests PRIVATE apps/desktop)
    target_link_libraries(field_phase_tests PRIVATE hs_bem_field hs_visualization hs_desktop_controls Qt6::Widgets)
    add_test(NAME field_phase COMMAND field_phase_tests)
    set_tests_properties(field_phase PROPERTIES LABELS "desktop;numerical" TIMEOUT 30
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    if(WIN32)
      add_custom_command(TARGET field_phase_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different $<TARGET_RUNTIME_DLLS:field_phase_tests>
        "$<TARGET_FILE_DIR:field_phase_tests>" COMMAND_EXPAND_LISTS VERBATIM)
    endif()
    target_include_directories(angular_response_widget_tests PRIVATE apps/desktop)
    target_link_libraries(angular_response_widget_tests PRIVATE hs_project hs_visualization
      hs_desktop_controls hs_options Qt6::Widgets)
    add_test(NAME angular_response_widget COMMAND angular_response_widget_tests)
    set_tests_properties(angular_response_widget PROPERTIES LABELS "desktop" TIMEOUT 30
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    if(WIN32)
      add_custom_command(TARGET angular_response_widget_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different $<TARGET_RUNTIME_DLLS:angular_response_widget_tests>
        "$<TARGET_FILE_DIR:angular_response_widget_tests>" COMMAND_EXPAND_LISTS VERBATIM)
      hs_deploy_cad_runtime(angular_response_widget_tests)
    endif()
    target_include_directories(realtime_simulation_tests PRIVATE apps/desktop)
    target_include_directories(realtime_simulation_tests SYSTEM PRIVATE "${HS_OCCT_INCLUDE}")
    target_link_libraries(realtime_simulation_tests PRIVATE hs_occt_TKPrim)
    target_compile_definitions(realtime_simulation_tests PRIVATE HORNSCULPT_VERSION="${PROJECT_VERSION}")
    target_link_libraries(realtime_simulation_tests PRIVATE hs_project hs_bem hs_bem_field hs_visualization hs_analysis hs_bem_galerkin
      hs_desktop_controls hs_desktop_export hs_study hs_result_comparison Qt6::Widgets Qt6::Test Qt6::Network)
    add_test(NAME subdomain_lem_gui COMMAND realtime_simulation_tests --subdomain-lem-gui-only)
    set_tests_properties(subdomain_lem_gui PROPERTIES LABELS "cuda;desktop" TIMEOUT 120 RUN_SERIAL TRUE
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_test(NAME subdomain_gui COMMAND realtime_simulation_tests --subdomain-gui-only)
    add_test(NAME symmetry_baffle_gui COMMAND realtime_simulation_tests --symmetry-baffle-gui-only)
    set_tests_properties(symmetry_baffle_gui PROPERTIES LABELS "cuda;desktop" TIMEOUT 180 RUN_SERIAL TRUE
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_test(NAME gmsh_import_gui COMMAND realtime_simulation_tests --gmsh-source-only)
    set_tests_properties(subdomain_gui PROPERTIES LABELS "cuda;desktop" TIMEOUT 120 RUN_SERIAL TRUE
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(gmsh_import_gui PROPERTIES LABELS "cuda;desktop" RUN_SERIAL TRUE TIMEOUT 120
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_test(NAME realtime_simulation COMMAND realtime_simulation_tests)
    add_test(NAME draft_preview_gui COMMAND realtime_simulation_tests --draft-preview-only)
    set_tests_properties(draft_preview_gui PROPERTIES LABELS "desktop" TIMEOUT 60
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_test(NAME draft_assembly_gui COMMAND realtime_simulation_tests --assembly-only)
    add_executable(guided_ui_tests tests/guided_ui_tests.cpp
      apps/desktop/angular_response_widget.cpp apps/desktop/phase_response_widget.cpp
      apps/desktop/acoustic_field_widget.cpp apps/desktop/simulation_dock.cpp
      apps/desktop/acoustic_frame_dialog.cpp apps/desktop/baffle_dialog.cpp apps/desktop/live_trace.cpp
      apps/desktop/main_window.cpp apps/desktop/parameter_panel.cpp apps/desktop/freecad_dock.cpp
      apps/desktop/mesh_panel.cpp apps/desktop/edit_history.cpp apps/desktop/parameter_sweep.cpp)
    target_include_directories(guided_ui_tests PRIVATE apps/desktop)
    target_include_directories(guided_ui_tests SYSTEM PRIVATE "${HS_OCCT_INCLUDE}")
    target_compile_definitions(guided_ui_tests PRIVATE HORNSCULPT_VERSION="${PROJECT_VERSION}")
    target_link_libraries(guided_ui_tests PRIVATE hs_occt_TKPrim hs_project hs_bem hs_bem_field
      hs_visualization hs_analysis hs_bem_galerkin hs_desktop_controls hs_desktop_export hs_study
      hs_result_comparison Qt6::Widgets Qt6::Test Qt6::Network)
    add_test(NAME guided_recovery_gui COMMAND guided_ui_tests --recovery)
    add_test(NAME guided_calculation_gui COMMAND guided_ui_tests --calculation)
    add_test(NAME guided_result_focus COMMAND guided_ui_tests --focus)
    add_test(NAME guided_late_cancel_gui COMMAND guided_ui_tests --late-cancel)
    add_test(NAME guided_recovery_reject_gui COMMAND guided_ui_tests --recovery-reject)
    add_test(NAME guided_comparison_epoch_gui COMMAND guided_ui_tests --comparison-epoch)
    add_test(NAME study_idle_gui COMMAND guided_ui_tests --study-idle)
    add_test(NAME study_handoff_gui COMMAND guided_ui_tests --study-handoff)
    add_test(NAME study_cached_gui COMMAND guided_ui_tests --study-cached)
    add_test(NAME study_epoch_gui COMMAND guided_ui_tests --study-epoch)
    set_tests_properties(study_idle_gui study_handoff_gui study_epoch_gui PROPERTIES LABELS "desktop" TIMEOUT 90
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(study_cached_gui PROPERTIES LABELS "cuda;desktop" RUN_SERIAL TRUE TIMEOUT 240
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(guided_late_cancel_gui guided_recovery_reject_gui guided_comparison_epoch_gui
      PROPERTIES LABELS "desktop" TIMEOUT 60
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(guided_result_focus PROPERTIES LABELS "desktop" TIMEOUT 120
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(guided_recovery_gui PROPERTIES LABELS "desktop" TIMEOUT 120
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(guided_calculation_gui PROPERTIES LABELS "cuda;desktop" RUN_SERIAL TRUE TIMEOUT 240
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_custom_command(TARGET guided_ui_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:guided_ui_tests> ${HS_CUDA_RUNTIME_FILES}
      "$<TARGET_FILE_DIR:guided_ui_tests>" COMMAND_EXPAND_LISTS VERBATIM)
    hs_deploy_cad_runtime(guided_ui_tests)
    add_test(NAME usability_gui COMMAND realtime_simulation_tests --usability-only)
    add_test(NAME diffraction_slot_gui COMMAND realtime_simulation_tests --slot-only)
    add_test(NAME diffraction_slot_radii_gui COMMAND realtime_simulation_tests --slot-radii-only)
    set_tests_properties(diffraction_slot_radii_gui PROPERTIES LABELS "desktop" TIMEOUT 120
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(diffraction_slot_gui PROPERTIES LABELS "desktop" TIMEOUT 90
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(usability_gui PROPERTIES LABELS "desktop" TIMEOUT 60
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(draft_assembly_gui PROPERTIES LABELS "desktop" TIMEOUT 90
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_test(NAME fmm_controls_gui COMMAND realtime_simulation_tests --fmm-controls-only)
    add_test(NAME baffle_controls_gui COMMAND realtime_simulation_tests --baffle-controls-only)
    add_test(NAME angular_response_gui COMMAND realtime_simulation_tests --angular-controls-only)
    add_test(NAME interface_controls_gui COMMAND realtime_simulation_tests --interface-controls-only)
    add_test(NAME axisymmetric_controls_gui COMMAND realtime_simulation_tests --axisymmetric-controls-only)
    set_tests_properties(axisymmetric_controls_gui PROPERTIES LABELS "desktop" TIMEOUT 60
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_test(NAME axisymmetric_integration_gui COMMAND realtime_simulation_tests --axisymmetric-only)
    set_tests_properties(axisymmetric_integration_gui PROPERTIES LABELS "desktop;cuda" RUN_SERIAL TRUE TIMEOUT 180
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(interface_controls_gui PROPERTIES LABELS "desktop" TIMEOUT 45
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(angular_response_gui PROPERTIES LABELS "desktop" TIMEOUT 30
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(baffle_controls_gui PROPERTIES LABELS "desktop" TIMEOUT 30
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(fmm_controls_gui PROPERTIES LABELS "desktop" TIMEOUT 30
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_executable(fmm_gui_tests tests/fmm_gui_tests.cpp apps/desktop/simulation_dock.cpp
      apps/desktop/angular_response_widget.cpp
      apps/desktop/phase_response_widget.cpp apps/desktop/acoustic_field_widget.cpp
      apps/desktop/acoustic_frame_dialog.cpp apps/desktop/baffle_dialog.cpp apps/desktop/live_trace.cpp)
    target_include_directories(fmm_gui_tests PRIVATE apps/desktop)
    target_link_libraries(fmm_gui_tests PRIVATE hs_project hs_visualization hs_bem hs_bem_field hs_bem_galerkin
      hs_desktop_controls Qt6::Widgets)
    add_test(NAME fmm_integration_gui COMMAND fmm_gui_tests)
    set_tests_properties(fmm_integration_gui PROPERTIES LABELS "cuda;desktop;fmm" RUN_SERIAL TRUE TIMEOUT 240
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    if(WIN32)
      add_custom_command(TARGET fmm_gui_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different $<TARGET_RUNTIME_DLLS:fmm_gui_tests> ${HS_CUDA_RUNTIME_FILES}
        "$<TARGET_FILE_DIR:fmm_gui_tests>" COMMAND_EXPAND_LISTS VERBATIM)
      hs_deploy_cad_runtime(fmm_gui_tests)
    endif()
    add_test(NAME lem_integration_gui COMMAND realtime_simulation_tests --lem-only)
    set_tests_properties(lem_integration_gui PROPERTIES LABELS "desktop;cuda" TIMEOUT 180
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_test(NAME face_orientation_gui COMMAND realtime_simulation_tests --orientation-only)
    add_test(NAME job_lifecycle_gui COMMAND realtime_simulation_tests --lifecycle-only)
    set_tests_properties(job_lifecycle_gui PROPERTIES LABELS "desktop" TIMEOUT 60
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(face_orientation_gui PROPERTIES LABELS "desktop" TIMEOUT 60
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_executable(live_trace_tests tests/live_trace_tests.cpp apps/desktop/live_trace.cpp)
    target_include_directories(live_trace_tests PRIVATE apps/desktop)
    target_link_libraries(live_trace_tests PRIVATE hs_options)
    add_test(NAME live_trace COMMAND live_trace_tests)
    set_tests_properties(realtime_simulation PROPERTIES LABELS "cuda;desktop" TIMEOUT 120
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_custom_command(TARGET realtime_simulation_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:realtime_simulation_tests> ${HS_CUDA_RUNTIME_FILES}
      "$<TARGET_FILE_DIR:realtime_simulation_tests>" COMMAND_EXPAND_LISTS VERBATIM)
    add_executable(geometry_mesh_worker_tests tests/geometry_mesh_worker_tests.cpp)
    target_include_directories(geometry_mesh_worker_tests PRIVATE apps/desktop)
    target_link_libraries(geometry_mesh_worker_tests PRIVATE Qt6::Core hs_options)
    add_test(NAME geometry_mesh_worker COMMAND geometry_mesh_worker_tests)
    set_tests_properties(geometry_mesh_worker PROPERTIES LABELS desktop TIMEOUT 30)
    add_executable(numeric_input_tests tests/numeric_input_tests.cpp)
    add_executable(parameter_panel_tests tests/parameter_panel_tests.cpp
      apps/desktop/parameter_panel.cpp)
    add_executable(step_export_dialog_tests tests/step_export_dialog_tests.cpp)
    target_link_libraries(step_export_dialog_tests PRIVATE hs_desktop_export Qt6::Widgets hs_options)
    add_test(NAME step_export_dialog COMMAND step_export_dialog_tests)
    set_tests_properties(step_export_dialog PROPERTIES LABELS "desktop" TIMEOUT 180 ENVIRONMENT
      "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_custom_command(TARGET step_export_dialog_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:step_export_dialog_tests>
      "$<TARGET_FILE_DIR:step_export_dialog_tests>" COMMAND_EXPAND_LISTS VERBATIM)
    hs_deploy_cad_runtime(step_export_dialog_tests)
    target_include_directories(parameter_panel_tests PRIVATE apps/desktop)
    target_link_libraries(parameter_panel_tests PRIVATE hs_desktop_controls hs_geometry Qt6::Test)
    add_test(NAME parameter_panel COMMAND parameter_panel_tests)
    set_tests_properties(parameter_panel PROPERTIES LABELS "desktop" TIMEOUT 30 ENVIRONMENT
      "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_custom_command(TARGET parameter_panel_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:parameter_panel_tests>
      "$<TARGET_FILE_DIR:parameter_panel_tests>" COMMAND_EXPAND_LISTS VERBATIM)
    add_executable(acoustic_frame_dialog_tests tests/acoustic_frame_dialog_tests.cpp
      apps/desktop/acoustic_frame_dialog.cpp)
    target_link_libraries(acoustic_frame_dialog_tests PRIVATE hs_desktop_controls hs_analysis Qt6::Test)
    add_test(NAME acoustic_frame_dialog COMMAND acoustic_frame_dialog_tests)
    set_tests_properties(acoustic_frame_dialog PROPERTIES LABELS "desktop" TIMEOUT 30 ENVIRONMENT
      "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_custom_command(TARGET acoustic_frame_dialog_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:acoustic_frame_dialog_tests>
      "$<TARGET_FILE_DIR:acoustic_frame_dialog_tests>" COMMAND_EXPAND_LISTS VERBATIM)
    target_link_libraries(numeric_input_tests PRIVATE hs_desktop_controls Qt6::Test)
    add_test(NAME numeric_input COMMAND numeric_input_tests)
    add_executable(mesh_panel_tests tests/mesh_panel_tests.cpp apps/desktop/mesh_panel.cpp)
    target_link_libraries(mesh_panel_tests PRIVATE hs_desktop_controls hs_meshing Qt6::Widgets)
    add_test(NAME mesh_panel COMMAND mesh_panel_tests)
    set_tests_properties(mesh_panel PROPERTIES TIMEOUT 30 ENVIRONMENT
      "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_custom_command(TARGET mesh_panel_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:mesh_panel_tests> "$<TARGET_FILE_DIR:mesh_panel_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
    set_tests_properties(numeric_input PROPERTIES ENVIRONMENT
      "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_custom_command(TARGET numeric_input_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:numeric_input_tests>
      "$<TARGET_FILE_DIR:numeric_input_tests>" COMMAND_EXPAND_LISTS VERBATIM)
    add_executable(guided_recovery_tests tests/guided_recovery_tests.cpp)
    target_include_directories(guided_recovery_tests PRIVATE apps/desktop)
    target_link_libraries(guided_recovery_tests PRIVATE hs_project Qt6::Core hs_options)
    add_test(NAME guided_recovery COMMAND guided_recovery_tests)
    set_tests_properties(guided_recovery PROPERTIES TIMEOUT 50)
    add_custom_command(TARGET guided_recovery_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:guided_recovery_tests> "$<TARGET_FILE_DIR:guided_recovery_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
    hs_deploy_cad_runtime(guided_recovery_tests)
    add_executable(result_comparison_tests tests/result_comparison_tests.cpp)
    target_link_libraries(result_comparison_tests PRIVATE hs_result_comparison hs_options)
    add_test(NAME result_comparison COMMAND result_comparison_tests)
    set_tests_properties(result_comparison PROPERTIES TIMEOUT 30
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_custom_command(TARGET result_comparison_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:result_comparison_tests> "$<TARGET_FILE_DIR:result_comparison_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
    hs_deploy_cad_runtime(result_comparison_tests)
    add_executable(project_file_tests tests/project_file_tests.cpp)
    add_executable(draft_assembly_tests tests/draft_assembly_tests.cpp)
    target_link_libraries(draft_assembly_tests PRIVATE hs_project Qt6::Core)
    target_compile_definitions(draft_assembly_tests PRIVATE
      HS_TEST_CAD="${CMAKE_CURRENT_SOURCE_DIR}/tests/fixtures/freecad-box-20x30x40mm.step")
    add_test(NAME draft_assembly COMMAND draft_assembly_tests)
    set_tests_properties(draft_assembly PROPERTIES TIMEOUT 90)
    add_custom_command(TARGET draft_assembly_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:draft_assembly_tests> "$<TARGET_FILE_DIR:draft_assembly_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
    hs_deploy_cad_runtime(draft_assembly_tests)
    target_link_libraries(project_file_tests PRIVATE hs_project Qt6::Core)
    add_test(NAME project_file COMMAND project_file_tests)
    add_test(NAME axisymmetric_archive COMMAND project_file_tests --axisymmetric-only)
    set_tests_properties(axisymmetric_archive PROPERTIES LABELS "numerical" TIMEOUT 90)
    add_executable(parameter_sweep_tests tests/parameter_sweep_tests.cpp)
    target_link_libraries(parameter_sweep_tests PRIVATE hs_study Qt6::Core)
    add_test(NAME parameter_sweep COMMAND parameter_sweep_tests)
    add_executable(edit_history_tests tests/edit_history_tests.cpp apps/desktop/edit_history.cpp)
    target_link_libraries(edit_history_tests PRIVATE hs_project Qt6::Widgets)
    add_test(NAME edit_history COMMAND edit_history_tests)
    add_custom_command(TARGET edit_history_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:edit_history_tests> "$<TARGET_FILE_DIR:edit_history_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
    add_executable(parameter_sweep_gui_tests tests/parameter_sweep_gui_tests.cpp apps/desktop/parameter_sweep.cpp)
    target_link_libraries(parameter_sweep_gui_tests PRIVATE hs_study hs_bem hs_bem_galerkin hs_analysis hs_visualization
      hs_desktop_controls Qt6::Widgets Qt6::Test)
    add_test(NAME parameter_sweep_gui COMMAND parameter_sweep_gui_tests)
    add_test(NAME parameter_sweep_orientation COMMAND parameter_sweep_gui_tests --orientation-only)
    add_test(NAME parameter_sweep_ui COMMAND parameter_sweep_gui_tests --sweep-ui-only)
    set_tests_properties(parameter_sweep_ui PROPERTIES LABELS "desktop" TIMEOUT 90
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(parameter_sweep_orientation PROPERTIES LABELS "desktop" TIMEOUT 90
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    set_tests_properties(parameter_sweep_gui PROPERTIES LABELS "cuda;desktop" TIMEOUT 240
      ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
    add_custom_command(TARGET parameter_sweep_gui_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:parameter_sweep_gui_tests> ${HS_CUDA_RUNTIME_FILES}
      "$<TARGET_FILE_DIR:parameter_sweep_gui_tests>" COMMAND_EXPAND_LISTS VERBATIM)
    add_custom_command(TARGET parameter_sweep_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:parameter_sweep_tests> "$<TARGET_FILE_DIR:parameter_sweep_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
    add_executable(response_io_tests tests/response_io_tests.cpp)
    target_link_libraries(response_io_tests PRIVATE hs_project Qt6::Core)
    add_test(NAME response_io COMMAND response_io_tests)
    add_custom_command(TARGET response_io_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:response_io_tests> "$<TARGET_FILE_DIR:response_io_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
    add_custom_command(TARGET project_file_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:project_file_tests> "$<TARGET_FILE_DIR:project_file_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
    add_executable(export_near_plan tests/export_near_plan.cpp)
    target_link_libraries(export_near_plan PRIVATE hs_bem_galerkin hs_meshing Qt6::Core)
    add_custom_command(TARGET export_near_plan POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:export_near_plan> "$<TARGET_FILE_DIR:export_near_plan>"
      COMMAND_EXPAND_LISTS VERBATIM)
    add_executable(export_near_plan_tests tests/export_near_plan_tests.cpp)
    target_link_libraries(export_near_plan_tests PRIVATE hs_bem_galerkin hs_meshing Qt6::Core)
    add_test(NAME export_near_plan COMMAND export_near_plan_tests $<TARGET_FILE:export_near_plan>)
    set_tests_properties(export_near_plan PROPERTIES TIMEOUT 120)
    add_dependencies(export_near_plan_tests export_near_plan)
    add_executable(solver_benchmark tests/solver_benchmark.cpp)
    target_link_libraries(solver_benchmark PRIVATE hs_bem hs_bem_galerkin hs_project Qt6::Core CUDA::cudart)
    add_custom_command(TARGET solver_benchmark POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:solver_benchmark> ${HS_CUDA_RUNTIME_FILES}
      "$<TARGET_FILE_DIR:solver_benchmark>" COMMAND_EXPAND_LISTS VERBATIM)
    hs_deploy_cad_runtime(solver_benchmark)
    add_executable(export_bem_fixture tests/export_bem_fixture.cpp)
    target_link_libraries(export_bem_fixture PRIVATE hs_geometry hs_meshing Qt6::Core)
    add_custom_command(TARGET export_bem_fixture POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:export_bem_fixture> "$<TARGET_FILE_DIR:export_bem_fixture>"
      COMMAND_EXPAND_LISTS VERBATIM)
    add_executable(freecad_tests tests/freecad_tests.cpp)
    target_link_libraries(freecad_tests PRIVATE hs_freecad hs_project hs_meshing Qt6::Network)
    add_test(NAME freecad_protocol COMMAND freecad_tests)
    add_custom_command(TARGET freecad_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:freecad_tests> "$<TARGET_FILE_DIR:freecad_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
    foreach(target realtime_simulation_tests project_file_tests parameter_sweep_tests
        edit_history_tests parameter_sweep_gui_tests response_io_tests freecad_tests)
      hs_deploy_cad_runtime(${target})
    endforeach()
  endif()
endif()
