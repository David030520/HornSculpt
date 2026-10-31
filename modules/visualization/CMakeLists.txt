find_package(Qt6 6.8 REQUIRED GLOBAL COMPONENTS OpenGLWidgets)

add_library(hs_visualization STATIC src/mesh_viewport.cpp src/plot_widget.cpp src/plot_axes.cpp
  src/polar_plot_widget.cpp src/heatmap_widget.cpp src/mesh_section_view.cpp src/chart_export.cpp
  src/plot_scale.cpp src/chart_scaling.cpp src/pressure_field_plane.cpp src/radial_heatmap.cpp)
target_include_directories(hs_visualization PUBLIC include)
target_link_libraries(hs_visualization PUBLIC hs_core hs_analysis Qt6::Widgets Qt6::OpenGLWidgets)
target_compile_features(hs_visualization PUBLIC cxx_std_20)

if(BUILD_TESTING)
  add_executable(overlay_layout_tests tests/overlay_layout_tests.cpp)
  target_link_libraries(overlay_layout_tests PRIVATE Qt6::Widgets hs_options)
  add_test(NAME overlay_layout COMMAND overlay_layout_tests)
  set_tests_properties(overlay_layout PROPERTIES TIMEOUT 30 ENVIRONMENT
    "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
  if(WIN32)
    add_custom_command(TARGET overlay_layout_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:overlay_layout_tests>
      "$<TARGET_FILE_DIR:overlay_layout_tests>" COMMAND_EXPAND_LISTS VERBATIM)
  endif()
  add_executable(plot_scaling_tests "${PROJECT_SOURCE_DIR}/tests/plot_scaling_tests.cpp")
  target_link_libraries(plot_scaling_tests PRIVATE hs_visualization)
  add_test(NAME plot_scaling COMMAND plot_scaling_tests)
  set_tests_properties(plot_scaling PROPERTIES TIMEOUT 30 ENVIRONMENT
    "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
  if(WIN32)
    add_custom_command(TARGET plot_scaling_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:plot_scaling_tests>
      "$<TARGET_FILE_DIR:plot_scaling_tests>" COMMAND_EXPAND_LISTS VERBATIM)
  endif()
  add_executable(mesh_section_tests tests/mesh_section_tests.cpp)
  target_link_libraries(mesh_section_tests PRIVATE hs_visualization hs_geometry)
  add_test(NAME mesh_section COMMAND mesh_section_tests)
  if(WIN32)
    add_custom_command(TARGET mesh_section_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:mesh_section_tests>
      "$<TARGET_FILE_DIR:mesh_section_tests>" COMMAND_EXPAND_LISTS VERBATIM)
  endif()
  add_executable(render_benchmark tests/render_benchmark.cpp)
  target_link_libraries(render_benchmark PRIVATE hs_visualization hs_geometry)
  add_test(NAME interface_render COMMAND render_benchmark --interface-only)
  set_tests_properties(interface_render PROPERTIES LABELS "desktop;opengl" TIMEOUT 60
    ENVIRONMENT "QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
  add_test(NAME orientation_render COMMAND render_benchmark --frames 2)
  set_tests_properties(orientation_render PROPERTIES LABELS "desktop;opengl" TIMEOUT 60
    ENVIRONMENT "QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
  if(WIN32)
    add_custom_command(TARGET render_benchmark POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:render_benchmark>
      "$<TARGET_FILE_DIR:render_benchmark>" COMMAND_EXPAND_LISTS VERBATIM)
  endif()
endif()
