---
title: Installation
parent: User Guide
nav_order: 1
---

# Installation

Stork is a header-only C++17 library built on Kokkos. Building the repository
also creates the `ReadAndInterpolate` command-line utility.

## Requirements

- A C++17 compiler supported by your Kokkos configuration
- CMake 3.16 or newer
- An installed Kokkos package
- A Kokkos backend appropriate for the target system, such as Serial, OpenMP,
  CUDA, or HIP

The current Stork library does not call MPI directly. A coupled application may
use MPI to distribute local grid regions and construct one Stork object per
rank.

## Build Kokkos

Build and install Kokkos for the hardware you intend to use. This CPU-only
example enables OpenMP:

```sh
git clone https://github.com/kokkos/kokkos.git
cmake -S kokkos -B kokkos/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PWD/kokkos/install" \
  -DKokkos_ENABLE_OPENMP=ON
cmake --build kokkos/build --parallel
cmake --install kokkos/build
```

Kokkos architecture flags are machine-specific. Consult the
[Kokkos build documentation](https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html)
before enabling an accelerator backend.

## Build Stork

From the Stork repository root, the build helper configures in `build/` and
installs Stork under `install/`:

```sh
./make.sh
```

The script uses `Kokkos_DIR` when set and otherwise checks
`$HOME/kokkos/build`. You can also pass additional CMake definitions directly:

```sh
Kokkos_DIR=/path/to/kokkos/install/lib/cmake/Kokkos \
  ./make.sh -D CMAKE_CXX_COMPILER=g++
```

The equivalent CMake commands are:

```sh
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="/path/to/kokkos/install" \
  -DCMAKE_INSTALL_PREFIX="$PWD/install"
cmake --build build --parallel
cmake --install build
```

The library headers, CMake package, and standalone converter are installed
under `install/`:

```sh
./install/bin/ReadAndInterpolate
```

With no arguments, it prints the expected command syntax.

The helper can be redirected without editing it:

```sh
STORK_BUILD_DIR=/tmp/stork-build \
STORK_INSTALL_DIR="$PWD/stork-install" \
CMAKE_BUILD_TYPE=RelWithDebInfo \
CMAKE_BUILD_PARALLEL_LEVEL=8 \
Kokkos_DIR=/path/to/kokkos/lib/cmake/Kokkos \
  ./make.sh
```

Additional arguments after `./make.sh` are passed to the CMake configure step.

## Consume Stork from CMake

After installation, a downstream project can use the exported target:

```cmake
cmake_minimum_required(VERSION 3.16)
project(stork_consumer LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Kokkos CONFIG REQUIRED)
find_package(Stork CONFIG REQUIRED)

add_executable(my_coupler main.cpp)
target_link_libraries(my_coupler PRIVATE Stork::Stork)
```

Configure the consumer with both installations on `CMAKE_PREFIX_PATH`:

```sh
cmake -S . -B build \
  -DCMAKE_PREFIX_PATH="/path/to/kokkos/install;/path/to/stork/install"
cmake --build build --parallel
```

## Verify the selected backend

Kokkos chooses Stork's `device_space` from `Kokkos::DefaultExecutionSpace`.
Confirm that the consumer and Stork use the same Kokkos installation and
configuration. Mixing headers or package files from different Kokkos builds can
produce compile, link, or runtime failures.

Stork does not initialize MPI. Applications may use
MPI around Stork, but the application owns `MPI_Init`/`MPI_Finalize` and its
domain-decomposition policy.
