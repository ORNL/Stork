file(GLOB APP_SOURCES *.cpp)
foreach(_app_file ${APP_SOURCES})
    # Name exe by source file.
    get_filename_component(_app ${_app_file} NAME_WE)
    add_executable(${_app} ${_app_file})
    target_link_libraries(${_app} PRIVATE Stork)
    install(TARGETS ${_app} DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach()
