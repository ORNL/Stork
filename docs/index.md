---
title: Home
nav_order: 1
description: Stork couples coarse thermal simulations to fine microstructure simulations.
---

![Stork](assets/images/stork-wordmark.png){: .hero-logo }

Stork seeks to serve as performance-portable middleware between additive
manufacturing simulation codes. Middleware is the software layer between
independent simulation tools: it supplies shared data structures,
transformations, and data movement so codes can exchange results efficiently
across CPUs and accelerators.

The current implementation supports thermal-to-microstructure coupling for
microstructure codes that use the Reduced Data Format (RDF). It can either
directly ingest a provided RDF or use a coarse thermal history in the sparse
Super Reduced Data Format (SRDF) and interpolate it into a fine-grid RDF.

[Install Stork](docs/installation/){: .btn .btn-primary }
[Quick start](docs/quick-start/){: .btn }
[Read the paper](https://doi.org/10.1016/j.commatsci.2026.114832){: .btn }

## Why Stork?

Thermal and microstructure solvers operate at different spatial and temporal
scales. Passing a complete fine-grid temperature history between them can make
storage and data movement more expensive than either simulation. Stork bridges
the two grids with:

- A sparse **Super Reduced Data Format (SRDF)** for active coarse cells.
- Kokkos-parallel spatial and temporal interpolation.
- Kokkos-parallel time-normalization.
- A compact **Reduced Data Format (RDF)** containing phase-change events.
- Host/device mirrors for file-based or in-memory workflows.

The accompanying paper reports more than two orders of magnitude reduction in
thermal-data generation time and file size for the studied workflows, while
preserving grain morphology and texture for lower interpolation factors. Those
findings describe the paper's configurations; users should validate the
interpolation factor for their own process conditions.

## Requirements

Stork requires C++17, CMake 3.16 or newer, and Kokkos. The library is
header-only; the selected Kokkos backend determines whether interpolation runs
on a CPU or accelerator.

## At a glance

| Need | Start here |
|:--|:--|
| Build and install the library | [Installation](docs/installation/) |
| Convert an SRDF file | [Quick start](docs/quick-start/) |
| Understand SRDF and RDF | [Data model](reference/data-model/) |
| Integrate Stork in C++ | [Library API](reference/api/) |
| Understand interpolation | [Interpolation](reference/interpolation/) |
| Check preconditions and limitations | [Contracts and limitations](reference/contracts/) |
| Read the methodology paper | [Paper](about/paper/) |

## Contributors

- [Benjamin Stump](https://www.ornl.gov/staff-profile/benjamin-c-stump)
- [John Coleman](https://www.ornl.gov/staff-profile/john-s-coleman)
