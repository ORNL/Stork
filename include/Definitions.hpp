#pragma once

// Kokkos includes
#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>

// For better coding autocompletion
#ifdef __INTELLISENSE__
#define KOKKOS_INLINE_FUNCTION inline
#define KOKKOS_LAMBDA []
#endif

namespace Stork {

    // Space definitions
    using host_space = Kokkos::DefaultHostExecutionSpace;
    using device_space = Kokkos::DefaultExecutionSpace;

    // Kokkos host memory and execution spaces
    using host_memory = Kokkos::DefaultHostExecutionSpace::memory_space;
    using host_exe = Kokkos::DefaultHostExecutionSpace::execution_space;

    // Kokkos device memory and exectution spaces
    using device_memory = Kokkos::DefaultExecutionSpace::memory_space;
    using device_exe = Kokkos::DefaultExecutionSpace::execution_space;

    // Set memory access layout
    typedef typename device_exe::array_layout layout;

    // Set quick access to views
    typedef Kokkos::View<uint16_t*, layout, device_memory> uint16_deviceView;
    typedef Kokkos::View<uint16_t*, layout, host_memory> uint16_hostView;

    typedef Kokkos::View<uint32_t*, layout, device_memory> uint32_deviceView;
    typedef Kokkos::View<uint32_t*, layout, host_memory> uint32_hostView;

    typedef Kokkos::View<float*, layout, device_memory> float_deviceView;
    typedef Kokkos::View<float*, layout, host_memory> float_hostView;

    typedef Kokkos::View<double*, layout, device_memory> double_deviceView;
    typedef Kokkos::View<double*, layout, host_memory> double_hostView;

}
