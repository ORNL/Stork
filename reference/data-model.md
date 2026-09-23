---
title: Data Model
parent: Reference
nav_order: 1
---

# Data Model

Stork uses a common regular-grid header with two sparse representations:
**SRDF** at the thermal grid scale and **RDF** at the microstructure grid scale.

## Regular-grid header

`Stork::Structs::RegularGrid_Header<FloatType, MemorySpace>` describes the local
region of a logically global Cartesian grid.

| Field | Meaning |
|:--|:--|
| `global_i0`, `global_j0`, `global_k0` | Global grid index of the local region's first node |
| `local_inum`, `local_jnum`, `local_knum` | Number of local nodes along each axis |
| `global_x0`, `global_y0`, `global_z0` | Coordinates of global index `(0,0,0)` |
| `gridResolution` | Uniform node spacing |

The underlying floating header view contains one additional slot with no named
meaning in the current implementation. It is preserved by copies and binary
I/O; leave it unused.

Coordinates follow

\\[
(x,y,z)=(x_0,y_0,z_0)+\Delta x(i,j,k).
\\]

Local point indices use `k` as the fastest-varying dimension:

\\[
p=i(N_jN_k)+jN_k+k.
\\]

The class supplies conversions among local point indices, local/global grid
indices, and global coordinates.

## SRDF

`Stork::Structs::SRDF_Dual<FloatType>` stores sparse thermal snapshots on the
coarse grid.

Object-level values:

| Value | Type | Meaning |
|:--|:--|:--|
| `numSnaps` | `uint32_t` | Number of sparse space-time records |
| `T_critical` | `FloatType` | Threshold used to detect melt and solidification crossings |
| `host_header`, `device_header` | regular-grid headers | Coarse-grid geometry in each memory space |

Each record `n` contains:

| Accessor/storage | Count | Meaning |
|:--|--:|:--|
| `p(n)` | 1 | Local point index of the coarse cell's low-index corner |
| `t_prev(n)`, `t_cur(n)` | 2 | Ordered time endpoints, with `t_cur > t_prev` |
| `T_prev(n)` | 8 | First of eight temperatures at `t_prev` |
| `T_cur(n)` | 8 | First of eight temperatures at `t_cur` |

The eight corner values at each time are contiguous and ordered with `k`
fastest, then `j`, then `i`:

```text
(0,0,0) (0,0,1) (0,1,0) (0,1,1)
(1,0,0) (1,0,1) (1,1,0) (1,1,1)
```

Offsets are relative to the coarse cell identified by `p(n)`.

An SRDF point index identifies the cell's `(i,j,k)` low corner. A producer must
not emit a record rooted on the last node of any axis because eight cell-corner
temperatures would not exist there. Times and temperatures have no units stored
with them; all coupled codes must agree on a consistent unit system.

{: .important }
SRDF sparsity is established by the producer. Stork does not reconstruct cells
that were omitted from the input.

## RDF

`Stork::Structs::RDF_Dual<FloatType>` stores fine-grid phase-change histories.
Its `numEvents` records contain:

| Accessor | Meaning |
|:--|:--|
| `p(n)` | Local fine-grid point index |
| `tm(n)` | Time at which temperature crosses upward through `T_critical` |
| `tl(n)` | Time at which temperature crosses downward through `T_critical` |
| `cr(n)` | Positive cooling-rate magnitude for the downward crossing interval |

The RDF header uses spacing `coarseResolution / fineFactor` and retains the
coarse grid's physical origin.

The positive `cr(n)` convention assumes forward SRDF intervals
(`t_cur > t_prev`) and a downward temperature crossing.

Multiple RDF records may use the same `p(n)` when a location remelts. Stork
sorts cycles at a point by time. Consumers should process only
`0 <= n < numEvents`; view capacity can be larger after trimmed interpolation.

## Dual memory-space objects

Both formats own host and device headers/data. Callers explicitly allocate
mirrors, copy values, and release unneeded data:

```cpp
object.Make_Data_Mirrors<source_space, destination_space>();
object.Copy_All<source_space, destination_space>();
object.Free_Data_Memory<memory_space>();
```

For a dual container, the supported template choices are exactly
`Stork::host_space` and `Stork::device_space`. Its accessors route any type that
is not exactly `host_space` to the single device member; arbitrary third spaces
are therefore unsupported.

When host and device execution spaces are the same, Stork's same-space helpers
can assign the device handles from the host handles instead of deep-copying.
That is not a permanent alias guarantee: binary readers reallocate only host
data views. Explicitly call `Make_Data_Mirrors` after a read before device-side
access.

Header views are allocated when the header objects are constructed. Data views
are empty until a reader or `Make_Data_Views` populates them. Allocation uses
`Kokkos::ViewAllocateWithoutInitializing`, so unfilled elements are undefined.
