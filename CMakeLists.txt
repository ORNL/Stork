cmake_minimum_required(VERSION 3.16)
project(Stork LANGUAGES CXX)
set(PROJECT_VERSION "0.1.0-dev")

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
    cmake_policy(SET CMP0135 NEW)
endif()

# Use standard install directories.
include(GNUInstallDirs)

find_package(Kokkos REQUIRED)

# Define the package configuration.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "StorkConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/StorkConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/StorkConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/Stork)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/StorkConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/StorkConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/Stork)

add_subdirectory(include)
add_subdirectory(apps)
