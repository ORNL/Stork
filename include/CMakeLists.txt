add_library(Stork INTERFACE)
# Add alias for downstream use.
add_library(Stork::Stork ALIAS Stork)

set_target_properties(Stork PROPERTIES INTERFACE_COMPILE_FEATURES cxx_std_17)

target_include_directories(Stork INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Link against dependencies.
target_link_libraries(Stork INTERFACE Kokkos::kokkos)

# Install the core library.
install(TARGETS Stork EXPORT StorkTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install public headers
install(DIRECTORY ./
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

# Export targets with namespace.
install(EXPORT StorkTargets
    FILE StorkTargets.cmake
    NAMESPACE Stork::
    DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/Stork
)
