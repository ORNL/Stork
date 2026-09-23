---
title: Couple in Memory
parent: Tutorials
nav_order: 2
---

# Couple in Memory

This minimal program creates one coarse cell with a heating interval followed
by a cooling interval. Stork refines the cell by a factor of two, detects the
critical-temperature crossings, and writes 27 fine-grid RDF histories.

```cpp
#include <Stork_Core.hpp>

#include <cstdint>

int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
        using Float = double;

        Stork::Structs::SRDF_Dual<Float> srdf;
        srdf.numSnaps = 2;
        srdf.T_critical = 1500.0;

        // One coarse cell has two nodes along each axis.
        srdf.host_header.global_i0() = 0;
        srdf.host_header.global_j0() = 0;
        srdf.host_header.global_k0() = 0;
        srdf.host_header.local_inum() = 2;
        srdf.host_header.local_jnum() = 2;
        srdf.host_header.local_knum() = 2;
        srdf.host_header.global_x0() = 0.0;
        srdf.host_header.global_y0() = 0.0;
        srdf.host_header.global_z0() = 0.0;
        srdf.host_header.gridResolution() = 20.0e-6;

        srdf.Make_Data_Views<Stork::host_space>(srdf.numSnaps);

        // Both records refer to the cell whose low-index corner is p = 0.
        srdf.host_data.p(0) = 0;
        srdf.host_data.t_prev(0) = 0.0;
        srdf.host_data.t_cur(0) = 1.0e-4;

        srdf.host_data.p(1) = 0;
        srdf.host_data.t_prev(1) = 1.0e-4;
        srdf.host_data.t_cur(1) = 2.0e-4;

        for (std::uint32_t corner = 0; corner < 8; ++corner) {
            // Record 0 heats through the critical temperature.
            srdf.host_data.thermals_view(corner) = 1400.0;
            srdf.host_data.thermals_view(8 + corner) = 1600.0;

            // Record 1 cools back through the critical temperature.
            srdf.host_data.thermals_view(16 + corner) = 1600.0;
            srdf.host_data.thermals_view(24 + corner) = 1400.0;
        }

        Stork::Structs::RDF_Dual<Float> rdf =
            Stork::Run::Interpolate_SRDF_to_RDF<
                Float,
                Stork::host_space,
                Float,
                Stork::host_space>(srdf, 2);

        Stork::IO::Output_RDF_csv(rdf, "minimal-rdf");
    }
    Kokkos::finalize();
    return 0;
}
```

Compile this as a CMake consumer using the target shown in
[Installation]({% link docs/installation.md %}). Running it produces
`minimal-rdf.csv`.

## Move interpolation to the device

When Kokkos is configured with an accelerator backend, change the target space:

```cpp
constexpr std::uint8_t fineFactor = 2;
Stork::Structs::RDF_Dual<float> rdf =
    Stork::Run::Interpolate_SRDF_to_RDF<
        double,
        Stork::host_space,
        float,
        Stork::device_space>(srdf, fineFactor);
```

Before host-side output, allocate mirrors and copy the result:

```cpp
rdf.Make_Data_Mirrors<Stork::device_space, Stork::host_space>();
rdf.Copy_All<Stork::device_space, Stork::host_space>();
Stork::IO::Output_RDF_csv(rdf, "device-rdf");
```

In a fully coupled application, pass `rdf.device_data` directly to the
microstructure stage and avoid that copy.
