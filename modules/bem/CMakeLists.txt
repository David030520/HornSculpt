find_package(Threads REQUIRED)
add_library(hs_bem_api INTERFACE)
target_include_directories(hs_bem_api INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(hs_bem_api INTERFACE hs_core hs_lem)
add_library(hs_bem_galerkin galerkin/space.cpp galerkin/pair_quadrature.cpp
  galerkin/exterior_domain.cpp galerkin/interaction_plan.cpp galerkin/axisymmetric_space.cpp
  galerkin/axisymmetric_quadrature.cpp surface_measures.cpp axisymmetric_reduction.cpp baffle_exterior.cpp
  fmm/support_tree.cpp fmm/support_near_plan.cpp
  fmm/wideband_transfer.cpp fmm/scaled_translation.cpp
  fmm/spherical_rotation.cpp fmm/translation_plan_cache.cpp fmm/support_level_plan.cpp)
target_include_directories(hs_bem_galerkin PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(hs_bem_galerkin PUBLIC hs_core hs_bem_api PRIVATE Threads::Threads hs_meshing)

add_library(hs_bem_field STATIC field_evaluation.cpp axisymmetric_field.cpp)
target_link_libraries(hs_bem_field PUBLIC hs_bem_api PRIVATE hs_bem_galerkin)

if(NOT HORNSCULPT_BUILD_CUDA_BEM)
  return()
endif()

add_library(hs_bem cuda/cuda_solver.cu cuda/subdomain_solver.cu cuda/baffled_solver.cu cuda/axisymmetric_assembly.cu cuda/device_gmres.cu
  cuda/wideband_angular.cu cuda/wideband_axial.cu cuda/wideband_rotation.cu
  cuda/support_near_operator.cu cuda/adaptive_patch_estimator.cu cuda/support_preconditioner.cu cuda/wideband_leaf.cu
  cuda/support_fmm_operator.cu cuda/support_fmm_solve.cu cuda/support_direct_audit.cu lem_coupling.cpp)
target_link_libraries(hs_bem PUBLIC hs_bem_api PRIVATE hs_bem_galerkin hs_bem_field hs_meshing CUDA::cudart CUDA::cublas CUDA::cusolver CUDA::cufft)
set_target_properties(hs_bem PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON)
target_compile_options(hs_bem PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/utf-8,/Zc:__cplusplus,/wd4211>)
# NVCC-generated launch stubs redeclare internal kernel symbols (MSVC C4211).
if(WIN32)
  foreach(runtime_name IN ITEMS cudart64_13 cublas64_13 cublasLt64_13 cusolver64_12 cusparse64_12 cufft64_12)
    find_file(HS_${runtime_name}_RUNTIME "${runtime_name}.dll"
      HINTS "${CUDAToolkit_BIN_DIR}/x64" "${CUDAToolkit_BIN_DIR}" REQUIRED)
    list(APPEND HS_CUDA_RUNTIME_FILES "${HS_${runtime_name}_RUNTIME}")
  endforeach()
  set(HS_CUDA_RUNTIME_FILES "${HS_CUDA_RUNTIME_FILES}" PARENT_SCOPE)
endif()

if(BUILD_TESTING)
  foreach(test_name IN ITEMS subdomain_lem subdomain_accuracy)
    add_executable(${test_name}_tests "${PROJECT_SOURCE_DIR}/tests/${test_name}_tests.cpp")
    target_link_libraries(${test_name}_tests PRIVATE hs_bem hs_bem_field hs_meshing hs_bem_galerkin)
    add_test(NAME ${test_name} COMMAND ${test_name}_tests)
    set_tests_properties(${test_name} PROPERTIES LABELS "cuda;numerical" TIMEOUT 600 RUN_SERIAL TRUE)
    if(WIN32)
      add_custom_command(TARGET ${test_name}_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:${test_name}_tests>" VERBATIM)
    endif()
  endforeach()
  add_executable(subdomain_bem_tests "${PROJECT_SOURCE_DIR}/tests/subdomain_bem_tests.cpp")
  target_link_libraries(subdomain_bem_tests PRIVATE hs_bem hs_bem_field hs_meshing)
  add_test(NAME subdomain_bem COMMAND subdomain_bem_tests)
  set_tests_properties(subdomain_bem PROPERTIES LABELS "cuda;numerical" TIMEOUT 600 RUN_SERIAL TRUE)
  if(WIN32)
    add_custom_command(TARGET subdomain_bem_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:subdomain_bem_tests>" VERBATIM)
  endif()
  add_executable(exterior_baffle_tests "${PROJECT_SOURCE_DIR}/tests/exterior_baffle_tests.cpp")
  target_link_libraries(exterior_baffle_tests PRIVATE hs_bem hs_bem_galerkin hs_bem_field hs_meshing)
  add_test(NAME exterior_baffle COMMAND exterior_baffle_tests)
  set_tests_properties(exterior_baffle PROPERTIES LABELS "cuda;numerical" TIMEOUT 300)
  if(WIN32)
    add_custom_command(TARGET exterior_baffle_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:exterior_baffle_tests>" VERBATIM)
  endif()
  add_executable(baffled_bem_tests "${PROJECT_SOURCE_DIR}/tests/baffled_bem_tests.cpp")
  target_link_libraries(baffled_bem_tests PRIVATE hs_bem hs_bem_field hs_meshing)
  add_test(NAME baffled_bem COMMAND baffled_bem_tests)
  set_tests_properties(baffled_bem PROPERTIES LABELS "cuda;numerical" TIMEOUT 600)
  add_test(NAME baffled_plan_cache_cuda COMMAND baffled_bem_tests --plan-cache-only)
  set_tests_properties(baffled_plan_cache_cuda PROPERTIES LABELS "cuda;numerical"
    RUN_SERIAL TRUE TIMEOUT 120)
  if(WIN32)
    add_custom_command(TARGET baffled_bem_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:baffled_bem_tests>" VERBATIM)
  endif()
  add_executable(bem_lem_tests "${PROJECT_SOURCE_DIR}/tests/bem_lem_tests.cpp")
  target_link_libraries(bem_lem_tests PRIVATE hs_bem)
  add_test(NAME bem_lem COMMAND bem_lem_tests)
  set_tests_properties(bem_lem PROPERTIES LABELS "cuda;numerical" TIMEOUT 180)
  if(WIN32)
    add_custom_command(TARGET bem_lem_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different ${HS_CUDA_RUNTIME_FILES} "$<TARGET_FILE_DIR:bem_lem_tests>" VERBATIM)
  endif()
endif()
