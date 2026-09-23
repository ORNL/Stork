# Stork

Stork seeks to serve as performance-portable middleware between additive
manufacturing simulation codes. Here, **middleware** means a software layer
that sits between otherwise independent simulation tools and provides shared
data structures, transformations, and data movement so those tools can exchange
results efficiently across CPUs and accelerators.

The current implementation supports one coupling path: thermal simulation data
is represented sparsely in the Super Reduced Data Format (SRDF), then converted
to fine-grid Reduced Data Format (RDF) phase-change events for use by a
microstructure simulation. Stork uses Kokkos-parallel quad-linear interpolation
to perform this conversion without constructing a large intermediate thermal
field.

The accompanying study reports substantial reductions in thermal-data
generation time and storage while retaining comparable microstructure results
for its validation cases. See
[Stump and Coleman (2026)](https://doi.org/10.1016/j.commatsci.2026.114832).

## Documentation

The full user guide, data-format reference, API guide, and tutorials are on the
[Stork documentation site](https://ornl.github.io/Stork/).

## Installation

Stork depends on Kokkos. See the
[installation guide](https://ornl.github.io/Stork/docs/installation/) for
details. With Kokkos installed, build and install Stork with:

```sh
./make.sh
```

Set `Kokkos_DIR` or `CMAKE_PREFIX_PATH` if CMake cannot locate Kokkos.

A downstream CMake project can import the installed target with:

```cmake
find_package(Kokkos CONFIG REQUIRED)
find_package(Stork CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE Stork::Stork)
```

## App usage

Installation includes the `ReadAndInterpolate` utility. It converts a binary
SRDF file to an interpolated RDF CSV file:

```sh
ReadAndInterpolate <filename-stem> <fineFactor> <precisionIn> <precisionOut>
```

- `<filename-stem>`: SRDF path without the `.stork` suffix
- `<fineFactor>`: integer refinement factor
- `<precisionIn>`: input precision, either `float` or `double`
- `<precisionOut>`: output precision, either `float` or `double`

The output is named `<filename-stem>.interp.<fineFactor>.csv`.

## Contributors

- [Benjamin Stump](https://www.ornl.gov/staff-profile/benjamin-c-stump)
- [John Coleman](https://www.ornl.gov/staff-profile/john-s-coleman)

## Contributing

Contributions are welcome. In particular, Stork welcomes new coupling formats,
transformations, backends, validation cases, and integrations that expand its
role as performance-portable middleware for additive manufacturing simulation
codes.
