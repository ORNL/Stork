---
title: Integrate Condor and Toucan
parent: Tutorials
nav_order: 3
---

# Integrate Condor and Toucan

This page is based on the checked-out `newCondor` and `newToucan` source trees.
It documents two concrete pipelines:

```text
file pipeline:      Condor -> SRDF.stork -> Stork interpolation -> Toucan
in-memory pipeline: Condor -> SRDF_Dual  -> Stork interpolation -> RDF_Dual -> Toucan
```

Condor can also produce RDF directly. The SRDF route is the Stork workflow that
decouples a coarse thermal grid from a finer microstructure grid.

{: .warning }
Keep the three repositories on compatible revisions. The checked-out Toucan
adapter still names the former `Stork::Run::Norm_Times_For_Toucan` function,
while the current Stork tree exposes `Normalize_SRDF_Times`. Update that call in
Toucan before building these exact revisions. Setting `NormalizeTimes` to
`false` changes runtime behavior but does not necessarily avoid the compile-time
name lookup. The replacement normalization implementation also has known
algorithmic issues described below, so keep normalization disabled until those
are corrected and validated.

## Link the projects

Condor and Toucan both consume the exported `Stork::Stork` CMake target. A
standalone producer can link Condor, which publicly carries Stork:

```cmake
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(Kokkos REQUIRED)
find_package(Stork REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(Condor CONFIG REQUIRED)

add_executable(make_thermal make_thermal.cpp)
target_link_libraries(make_thermal PRIVATE Condor::Condor)
```

These are the imported dependencies named by the current Condor
`CMakeLists.txt`; finding them before Condor also ensures its exported target's
references exist. Configure with install prefixes for Kokkos, Stork, Condor,
and nlohmann_json on `CMAKE_PREFIX_PATH`. OpenMP and MPI must be discoverable by
CMake. All projects must resolve the same Kokkos and MPI builds and backend
configuration.

## Condor as an SRDF producer

The public Condor integration is `Condor::Run::Coupled`. It requires Condor's
`Interface` mode and fills a caller-owned Stork container. This example is
adapted from `newCondor/apps/Make_SRDF.cpp`:

```cpp
#include <Condor_Core.hpp>
#include <mpi.h>
#include <string>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);
    {
        using ThermalFloat = double;
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const std::string input =
            argc > 1 ? argv[1] : "ParamInput.json";

        Stork::Structs::SRDF_Dual<ThermalFloat> srdf;
        Condor::Run::Coupled<ThermalFloat>(srdf, input);

        // Condor leaves the authoritative SRDF in device_space.
        srdf.Make_Data_Mirrors<Stork::device_space,
                               Stork::host_space>();
        srdf.Copy_All<Stork::device_space, Stork::host_space>();
        const std::string output =
            "layer-0-rank-" + std::to_string(rank);
        Stork::IO::Output_SRDF_binary(srdf, output);
    }
    Kokkos::finalize();
    MPI_Finalize();
}
```

Every MPI rank writes a distinct stem. Never let multiple ranks call a Stork
writer on the same path; the writer is not collective and does not coordinate
concurrent access.

`Condor::Run::Coupled` disables normal Interface file output and ordinary JSON
hooks for a coupled RDF/SRDF run. It throws when the selected Condor mode is not
`Interface`. The SRDF overload constructs Condor's MPI wrapper with one
requested decomposition dimension; use the behavior of the exact Condor
revision when planning rank-local files.

Condor also provides an RDF overload:

```cpp
Stork::Structs::RDF_Dual<double> rdf;
Condor::Run::Coupled<double>(rdf, "ParamInput.json");
```

That route records thermal solidification events directly and does not perform
Stork's SRDF-to-fine-RDF interpolation.

## Convert Condor SRDF in memory

When the microstructure grid is finer, leave Condor's SRDF on the device and
interpolate there:

```cpp
using ThermalFloat = double;
using MicroFloat = float;

Stork::Structs::SRDF_Dual<ThermalFloat> srdf;
Condor::Run::Coupled<ThermalFloat>(srdf, "ParamInput.json");

constexpr std::uint8_t fineFactor = 4;
auto rdf = Stork::Run::Interpolate_SRDF_to_RDF<
    ThermalFloat,
    Stork::device_space,
    MicroFloat,
    Stork::device_space>(srdf, fineFactor);

// rdf.device_header and rdf.device_data can now be consumed in place.
```

This is the pattern in `newCondor/apps/Test_SRDF.cpp`. It avoids an intermediate
file and a device-to-host round trip.

Keep every user-supplied fine factor in `[1,255]`. Stork accepts `uint8_t`; a
wider integer passed by a caller is narrowed without validation.

## Match an MPI-decomposed Toucan domain

Toucan's adapters avoid overlap between neighboring x-y rank domains by
interpolating and retaining inclusive global fine-grid bounds. For an `I x J`
rank grid at coordinates `(i,j)`, their current calculation is:

```cpp
const uint32_t iMin =
    srdf.host_header.global_i0() * f + (f / 2) * (i != 0);
const uint32_t iMax =
    (srdf.host_header.global_i0() + srdf.host_header.local_inum() - 1) * f
    - ((f - 1) / 2) * (i != I - 1);
const uint32_t jMin =
    srdf.host_header.global_j0() * f + (f / 2) * (j != 0);
const uint32_t jMax =
    (srdf.host_header.global_j0() + srdf.host_header.local_jnum() - 1) * f
    - ((f - 1) / 2) * (j != J - 1);

auto rdf = Stork::Run::Interpolate_And_Trim_SRDF_to_RDF<
    ThermalFloat, Stork::device_space,
    float, Stork::device_space>(
        srdf, f, iMin, iMax, jMin, jMax);

if (rdf.numEvents > 0) {
    Stork::Run::Trim_RDF_in_Z<float, Stork::device_space>(rdf);
    rdf.Copy_Header<Stork::device_space, Stork::host_space>();
}
```

The `iMin`/`iMax` and `jMin`/`jMax` formulas are integration policy from the
current Toucan adapter, not hidden behavior inside Stork. Validate them for a
different decomposition. Guard z trimming when a rank has no events.

{: .warning }
The hand-written example above has the required empty-event guard, but both
checked-out Toucan adapters (`FileSource` and `CondorSource`) currently call
`Trim_RDF_in_Z` unconditionally. A rank/layer with zero retained events is
therefore unsafe in those adapters until they add an equivalent guard.

## Let Toucan run Condor directly

Toucan's `Thermal::CondorSource` wraps the preceding pipeline. The following
configuration is adapted from
`newCondor/Tests/CondorStorkToucan-01/ToucanSettings-Top.json`:

```json
{
  "Thermal": {
    "Condor": {
      "Files": [
        "thermal/ParamInput-A.json",
        "thermal/ParamInput-B.json"
      ],
      "Precision": "double",
      "SRDF": {
        "FineFactor": 4,
        "NormalizeTimes": false,
        "TimeScale": 0.01,
        "Reuse": true,
        "Cache": false
      }
    }
  }
}
```

Toucan cycles through `Files` by layer. Exactly one of `RDF` and `SRDF` must be
present:

- `SRDF` runs Condor into SRDF, optionally normalizes time, interpolates by
  `FineFactor`, trims to the rank's x-y region, then trims unused z depth.
  Here `Precision` selects the Condor/SRDF `float` or `double` type, and Toucan
  converts to its own float type during interpolation.
- `RDF` asks Condor for RDF directly using `ToucanFloat`. Although the adapter
  parses `Precision`, the direct RDF path does not use it. Its settings support
  `Reuse` and `Cache`.
- `Reuse` retains generated Stork objects for repeated layers.
- `Cache` stores rank/decomposition-specific `.stork` files under `cache/`.

The current JSON readers narrow `FineFactor` to `uint8_t` without a complete
upper-bound check (and the live Condor adapter does not reject zero). Supply an
integer in `[1,255]` yourself.

The name mismatch must be fixed before the adapter builds, but the replacement
normalization functions are also known-problematic and should not be enabled
without correcting and validating their cluster assignment and centering. See
[Contracts and limitations]({% link reference/contracts.md %}). `TimeScale`
would be passed as `minGap`; Stork does not infer its units.

## Let Toucan read `.stork` files

Toucan's `Thermal::FileSource` accepts RDF or SRDF. A single-layer SRDF example
is:

```json
{
  "Thermal": {
    "File": {
      "String": "thermal/layer-0-rank-RANK.stork",
      "Layers": 1,
      "Ranks": 4,
      "Precision": "double",
      "SRDF": {
        "FineFactor": 4,
        "NormalizeTimes": false,
        "Reuse": true
      }
    }
  }
}
```

`RANK` and `LAYER` are numeric filename tokens. `Layers` and `Ranks` may be an
integer count starting at zero or an inclusive `[start, stop]` pair. Toucan
indexes matching files during setup, accepts names with or without the
`.stork` suffix, and requires exactly one file for every requested layer on the
current file rank. It rejects ambiguous matches.

For `FileSource`, `Precision` is the floating type stored in either an RDF or
SRDF file. RDF is converted to `ToucanFloat` when necessary; SRDF retains the
stored type until interpolation produces Toucan's output type. For an SRDF,
provide `FineFactor` in `[1,255]`; the parser rejects nonpositive values but
still narrows values above 255.

For multiple layers, for example:

```json
"String": "thermal/rank-RANK/layer-LAYER.srdf.stork",
"Layers": [0, 15],
"Ranks": 4
```

For precomputed RDF, replace the `SRDF` block with:

```json
"RDF": { "Reuse": true }
```

RDF input bypasses interpolation. SRDF input follows the same interpolation and
trimming path as live Condor coupling. Both current SRDF adapter paths
unconditionally z trim and therefore assume at least one retained event per
rank/layer.

## What is Stork's responsibility?

Stork owns the shared grid/data representations, memory-space movement, serialization interpolation, event formation, alongside optional time-normalization and trimming. 
