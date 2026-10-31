add_library(hs_desktop_controls STATIC numeric_input.cpp lem_editor.cpp draft_assembly_dialog.cpp draft_wizard.cpp)
target_include_directories(hs_desktop_controls PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(hs_desktop_controls PUBLIC Qt6::Widgets hs_lem hs_geometry PRIVATE hs_options)

add_library(hs_desktop_export STATIC step_export_dialog.cpp)
target_include_directories(hs_desktop_export PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(hs_desktop_export PUBLIC hs_geometry Qt6::Widgets
    PRIVATE hs_cad hs_desktop_controls hs_options)

add_library(hs_project STATIC project_file.cpp project_json.cpp project_cad_json.cpp step_import.cpp
    project_simulation_json.cpp project_lem_json.cpp simulation_request.cpp measured_response.cpp result_export.cpp
    angular_response_settings.cpp project_assembly.cpp)
target_include_directories(hs_project PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(hs_project PUBLIC hs_geometry hs_bem_api hs_freecad hs_meshing hs_analysis
    PRIVATE hs_bem_galerkin Qt6::Core)

add_library(hs_result_comparison STATIC result_comparison.cpp)
target_include_directories(hs_result_comparison PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(hs_result_comparison PUBLIC hs_project hs_visualization Qt6::Widgets
    PRIVATE hs_analysis hs_options)

add_library(hs_study STATIC sweep_model.cpp)
target_link_libraries(hs_study PUBLIC hs_project)
target_link_libraries(hs_result_comparison PRIVATE hs_study)
target_link_libraries(hs_result_comparison PRIVATE hs_desktop_export)

qt_add_executable(HornSculpt main.cpp main_window.cpp parameter_panel.cpp
    simulation_dock.cpp angular_response_widget.cpp phase_response_widget.cpp acoustic_field_widget.cpp acoustic_frame_dialog.cpp baffle_dialog.cpp live_trace.cpp freecad_dock.cpp mesh_panel.cpp edit_history.cpp parameter_sweep.cpp
    branding.cpp branding.qrc qt_localization.cpp qt_localization.qrc)
set_target_properties(HornSculpt PROPERTIES AUTORCC ON)
if(WIN32)
  enable_language(RC)
  target_sources(HornSculpt PRIVATE branding.rc)
endif()

target_link_libraries(HornSculpt PRIVATE hs_geometry hs_analysis hs_meshing hs_visualization hs_bem_galerkin
    hs_bem hs_bem_field hs_freecad hs_project hs_study hs_result_comparison hs_desktop_controls hs_desktop_export Qt6::Widgets)
target_compile_features(HornSculpt PRIVATE cxx_std_20)
target_compile_definitions(HornSculpt PRIVATE HORNSCULPT_VERSION="${PROJECT_VERSION}")
set_target_properties(HornSculpt PROPERTIES WIN32_EXECUTABLE TRUE)
add_custom_command(TARGET HornSculpt POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory "$<TARGET_FILE_DIR:HornSculpt>/integrations/freecad/workbench"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${PROJECT_SOURCE_DIR}/integrations/freecad/bridge.py"
    "${PROJECT_SOURCE_DIR}/integrations/freecad/workbench/Init.py"
    "${PROJECT_SOURCE_DIR}/integrations/freecad/workbench/InitGui.py"
    "${PROJECT_SOURCE_DIR}/integrations/freecad/workbench/hornsculpt_commands.py"
    "$<TARGET_FILE_DIR:HornSculpt>/integrations/freecad/workbench" VERBATIM)
install(TARGETS HornSculpt RUNTIME DESTINATION bin)
install(FILES "${PROJECT_SOURCE_DIR}/Start-HornSculpt.cmd" DESTINATION .)
if(WIN32)
  add_custom_command(TARGET HornSculpt POST_BUILD COMMAND ${CMAKE_COMMAND} -E
    copy_if_different ${HS_CUDA_RUNTIME_FILES} ${HS_CAD_RUNTIME_FILES}
    "$<TARGET_FILE_DIR:HornSculpt>" VERBATIM)
  install(FILES ${HS_CUDA_RUNTIME_FILES} ${HS_CAD_RUNTIME_FILES} DESTINATION bin)
endif()
qt_generate_deploy_app_script(TARGET HornSculpt OUTPUT_SCRIPT deployment_script
    NO_TRANSLATIONS NO_UNSUPPORTED_PLATFORM_ERROR)
install(SCRIPT ${deployment_script})
if(WIN32)
  if(NOT TARGET Qt6::QOffscreenIntegrationPlugin)
    message(FATAL_ERROR "Qt offscreen platform plugin is required for installed CLI verification")
  endif()
  install(FILES "$<TARGET_FILE:Qt6::QOffscreenIntegrationPlugin>" DESTINATION plugins/platforms)
endif()

if(BUILD_TESTING)
  add_executable(qt_localization_tests "${PROJECT_SOURCE_DIR}/tests/qt_localization_tests.cpp"
    qt_localization.cpp qt_localization.qrc)
  set_target_properties(qt_localization_tests PROPERTIES AUTORCC ON
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
  target_include_directories(qt_localization_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(qt_localization_tests PRIVATE Qt6::Widgets hs_options)
  add_test(NAME qt_localization COMMAND qt_localization_tests)
  set_tests_properties(qt_localization PROPERTIES LABELS "desktop" TIMEOUT 30
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
  if(WIN32)
    add_custom_command(TARGET qt_localization_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:qt_localization_tests> "$<TARGET_FILE_DIR:qt_localization_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
  endif()
  add_executable(draft_wizard_tests "${PROJECT_SOURCE_DIR}/tests/draft_wizard_tests.cpp")
  target_link_libraries(draft_wizard_tests PRIVATE hs_desktop_controls hs_options)
  set_target_properties(draft_wizard_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
  add_test(NAME draft_wizard COMMAND draft_wizard_tests)
  set_tests_properties(draft_wizard PROPERTIES LABELS "desktop" TIMEOUT 30
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
  if(WIN32)
    add_custom_command(TARGET draft_wizard_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:draft_wizard_tests> "$<TARGET_FILE_DIR:draft_wizard_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
  endif()
  add_executable(branding_tests "${PROJECT_SOURCE_DIR}/tests/branding_tests.cpp"
    branding.cpp branding.qrc)
  set_target_properties(branding_tests PROPERTIES AUTORCC ON
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
  target_include_directories(branding_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(branding_tests PRIVATE Qt6::Widgets hs_options)
  target_compile_definitions(branding_tests PRIVATE HORNSCULPT_VERSION="${PROJECT_VERSION}")
  add_test(NAME branding COMMAND branding_tests)
  set_tests_properties(branding PROPERTIES LABELS "desktop" TIMEOUT 30
    ENVIRONMENT "QT_QPA_PLATFORM=offscreen;QT_PLUGIN_PATH=$<TARGET_FILE_DIR:Qt6::Core>/../plugins")
  if(WIN32)
    add_custom_command(TARGET branding_tests POST_BUILD COMMAND ${CMAKE_COMMAND} -E
      copy_if_different $<TARGET_RUNTIME_DLLS:branding_tests> "$<TARGET_FILE_DIR:branding_tests>"
      COMMAND_EXPAND_LISTS VERBATIM)
  endif()
endif()


add_executable(cad_model_project "${PROJECT_SOURCE_DIR}/apps/tools/create_cad_model_project.cpp")
target_link_libraries(cad_model_project PRIVATE hs_project hs_bem hs_options)
target_include_directories(cad_model_project PRIVATE "${PROJECT_SOURCE_DIR}/integrations/freecad")
set_target_properties(cad_model_project PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
if(WIN32)
  add_custom_command(TARGET cad_model_project POST_BUILD COMMAND ${CMAKE_COMMAND} -E
    copy_if_different $<TARGET_RUNTIME_DLLS:cad_model_project> ${HS_CUDA_RUNTIME_FILES} ${HS_CAD_RUNTIME_FILES}
    "$<TARGET_FILE_DIR:cad_model_project>" COMMAND_EXPAND_LISTS VERBATIM)
endif()
