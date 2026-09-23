#!/usr/bin/env bash

set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${STORK_BUILD_DIR:-${project_dir}/build}"
install_dir="${STORK_INSTALL_DIR:-${project_dir}/install}"
build_type="${CMAKE_BUILD_TYPE:-Release}"
parallel_jobs="${CMAKE_BUILD_PARALLEL_LEVEL:-$(nproc)}"

cmake_args=(
    -S "${project_dir}"
    -B "${build_dir}"
    -D "CMAKE_BUILD_TYPE=${build_type}"
    -D "CMAKE_INSTALL_PREFIX=${install_dir}"
)

if [[ -n "${Kokkos_DIR:-}" ]]; then
    cmake_args+=(-D "Kokkos_DIR=${Kokkos_DIR}")
elif [[ -f "${HOME}/kokkos/build/KokkosConfig.cmake" ]]; then
    cmake_args+=(-D "Kokkos_DIR=${HOME}/kokkos/build")
fi

cmake "${cmake_args[@]}" "$@"
cmake --build "${build_dir}" --parallel "${parallel_jobs}"
cmake --install "${build_dir}"

echo "Stork installed to ${install_dir}"
