---
title: Library API
parent: Reference
nav_order: 5
---

# Library API

Stork is a header-only C++17 library. Include the supported aggregate header:

```cpp
#include <Stork_Core.hpp>
```

Initialize Kokkos before constructing or using Stork objects, and destroy those
objects before `Kokkos::finalize()`:

```cpp
Kokkos::initialize(argc, argv);
{
    // Stork work
}
Kokkos::finalize();
```

The sections below inventory the caller-facing symbols in the current source
tree. Names under a namespace named `impl` are implementation details and are
not a stable API.

## Execution-space aliases

`Definitions.hpp` defines:

| Alias | Definition |
|:--|:--|
| `Stork::host_space` | `Kokkos::DefaultHostExecutionSpace` |
| `Stork::device_space` | `Kokkos::DefaultExecutionSpace` |
| `Stork::host_memory`, `Stork::device_memory` | memory spaces associated with those execution spaces |
| `Stork::host_exe`, `Stork::device_exe` | the corresponding execution spaces |
| `Stork::layout` | array layout selected by the default device execution space |

Eight low-level, one-dimensional Kokkos view aliases are also public:

```text
uint16_hostView   uint16_deviceView
uint32_hostView   uint32_deviceView
float_hostView    float_deviceView
double_hostView   double_deviceView
```

They are implementation conveniences without container invariants or a
separate compatibility guarantee. Application code should normally use the
typed Stork containers.

{: .note }
The template parameter named `memory_space` throughout the current API is used
as a Kokkos execution/device argument. The dual containers dispatch exactly on
`Stork::host_space`; every other type is routed to their one
`Stork::device_space` member. Their supported choices are therefore exactly
`Stork::host_space` and `Stork::device_space`, not an arbitrary compatible
Kokkos space and not the `host_memory`/`device_memory` aliases.

## Regular-grid header

```cpp
Stork::Structs::RegularGrid_Header<FloatType, Space>
```

The header owns a six-element integer view and a five-element floating-point
view. It exposes references through:

```text
global_i0()  global_j0()  global_k0()
local_inum() local_jnum() local_knum()
global_x0()  global_y0()  global_z0()
gridResolution()
```

The fifth floating-point slot is serialized but has no named accessor in the
current implementation. Callers should not assign meaning to it.

The header also provides `KOKKOS_INLINE_FUNCTION` conversions:

| Function | Conversion |
|:--|:--|
| `LOCAL_ijk_to_LOCAL_p` | local `(i,j,k)` to local flat index |
| `LOCAL_p_to_LOCAL_ijk` | local flat index to local `(i,j,k)` |
| `LOCAL_ijk_to_GLOBAL_ijk` | add the global index offset |
| `GLOBAL_ijk_to_LOCAL_ijk` | subtract the global index offset |
| `LOCAL_p_to_GLOBAL_ijk` | local flat index to global `(i,j,k)` |
| `GLOBAL_ijk_to_LOCAL_p` | global `(i,j,k)` to local flat index |
| `GLOBAL_ijk_to_GLOBAL_xyz` | global grid index to coordinates |
| `LOCAL_p_to_GLOBAL_xyz` | local flat index to coordinates |
| `LOCAL_ijk_to_GLOBAL_xyz` | local grid index to coordinates |

These functions do not perform bounds checking. Unsigned subtraction in a
global-to-local conversion can wrap when the requested global index precedes
the local origin.

```cpp
Stork::Structs::IndexBounds<Space>
```

`IndexBounds` exposes `imin()`, `imax()`, `jmin()`, `jmax()`, `kmin()`, and
`kmax()`. Its default-constructed `index_view` is unallocated, so callers would
have to allocate it before using an accessor. The type is currently referenced
only by an internal trimming helper; the supported top-level trimmed
interpolation API accepts four explicit global x-y bounds instead.

## Single-space data views

```cpp
Stork::Structs::SRDF_Data<FloatType, Space>
```

Owns `cellNum_view`, `times_view`, and `thermals_view`. For record `n`, use
`p(n)`, `t_prev(n)`, `t_cur(n)`, `T_prev(n)`, and `T_cur(n)`. `T_prev(n)` and
`T_cur(n)` return references to the first of eight contiguous corner values;
index `thermals_view` directly to access all corners.

```cpp
Stork::Structs::RDF_Data<FloatType, Space>
```

Owns `cellNum_view` and `solInfo_view`. For event `n`, use `p(n)`, `tm(n)`,
`tl(n)`, and `cr(n)`.

Accessors return mutable references even through a `const` data object. Treat
that as an implementation characteristic, not as permission to mutate data
through a logically const interface.

## Dual-space containers

```cpp
Stork::Structs::SRDF_Dual<FloatType>
Stork::Structs::RDF_Dual<FloatType>
```

`SRDF_Dual` contains `numSnaps`, `T_critical`, `host_header`, `device_header`,
`host_data`, and `device_data`. `RDF_Dual` contains `numEvents` and the same
host/device header and data pattern.

Both types provide:

```cpp
object.get_header<Space>();
object.get_data<Space>();
object.Make_Data_Views<Space>(count);
object.Make_Data_Mirrors<SourceSpace, DestinationSpace>();
object.Copy_Header<SourceSpace, DestinationSpace>();
object.Copy_Data<SourceSpace, DestinationSpace>();
object.Copy_All<SourceSpace, DestinationSpace>();
object.Free_Data_Memory<Space>();
```

They also expose `Mirror_Header_When_Same_Space()` and
`Mirror_Data_When_Same_Space()`. These assign device view handles from host
view handles only when `host_space` and `device_space` are the same type. They
are low-level synchronization helpers; they do not copy data between distinct
spaces.

`Make_Data_Views` allocates without initialization. Fill every used element
before reading it. `Make_Data_Mirrors` allocates destination **data** views with
the source extents; headers are allocated by their constructors. It does not
copy values, so follow it with `Copy_Data` or `Copy_All`. Copy operations assume
matching view extents.

`Free_Data_Memory` shrinks data views to zero only when host and device spaces
are distinct; it does not reset `numSnaps` or `numEvents`. When the two spaces
are the same type, the method intentionally does nothing. Do not assume the two
handles always remain aliased: binary readers reallocate only the host data
views. Call `Make_Data_Mirrors` explicitly after a read before accessing the
device handles, even in a host-only configuration.

## Interpolation

```cpp
template<class FloatIn, class SourceSpace,
         class FloatOut, class TargetSpace>
Stork::Structs::RDF_Dual<FloatOut>
Stork::Run::Interpolate_SRDF_to_RDF(
    Stork::Structs::SRDF_Dual<FloatIn>& srdf,
    uint8_t fineFactor);
```

Both `SourceSpace` and `TargetSpace` must be exactly `Stork::host_space` or
`Stork::device_space`. The function first reads `srdf.host_header` to construct
the fine RDF header regardless of `SourceSpace`. It always invokes the SRDF
mirror/copy methods, so the non-const input may be mutated. When `SourceSpace`
and `TargetSpace` differ, those calls allocate/overwrite the target-side SRDF
data handles and deep-copy the source header/data. When they are the same type,
the same-space branch is effectively a no-op on a distinct-host/device build
and may refresh aliases on a host-only build. Interpolation then runs in
`TargetSpace` and returns RDF event data there. `FloatIn` and `FloatOut` may
differ.

```cpp
template<class FloatIn, class SourceSpace,
         class FloatOut, class TargetSpace>
Stork::Structs::RDF_Dual<FloatOut>
Stork::Run::Interpolate_And_Trim_SRDF_to_RDF(
    Stork::Structs::SRDF_Dual<FloatIn>& srdf,
    uint8_t fineFactor,
    uint32_t globalIMin, uint32_t globalIMax,
    uint32_t globalJMin, uint32_t globalJMax);
```

Performs the same interpolation, keeps events in the inclusive global fine-grid
x-y bounds, and remaps their local point indices to the trimmed header. See
[Interpolation]({% link reference/interpolation.md %}) for preconditions and
the capacity caveat. It has the same host-header read, non-const input, and
conditional mirror/copy behavior as the untrimmed function.

## Additional transformations

```cpp
template<class FloatType, class Space>
void Stork::Run::Trim_RDF_in_Z(
    Stork::Structs::RDF_Dual<FloatType>& rdf);
```

Finds the minimum local `k` used by an RDF event, advances `global_z0`, reduces
`local_knum`, and remaps every event's local point index. The RDF header must be
valid in both host and `Space`; event data must be valid in `Space`. Call
`Copy_Header<Space, host_space>()` afterward if downstream host code needs the
updated header. `numEvents` must be greater than zero.

The current source exposes both of these experimental functions:

```cpp
Stork::Run::Normalize_RDF_Times<FloatType, Space>(srdf, minGap);
Stork::Run::Normalize_SRDF_Times<FloatType, Space>(srdf, minGap);
```

Despite the first name, both currently accept `SRDF_Dual<FloatType>&` and have
the same implementation. They mutate only the data in `Space` and require at
least one snapshot.

{: .warning }
Do not rely on the current functions for correct cluster separation or exact
zero-centering. The scan stores `clusterID(n)` before applying the gap between
records `n-1` and `n`; consequently, the first record after a gap stays in the
prior cluster and a trailing gap can create an empty counted cluster. The
`baseShift` formula also includes `numClusters * minGap` even though there are
only `numClusters - 1` inserted gaps, shifting the intended overall center by
`-minGap/2` rather than centering it at zero.

## Input and output

```cpp
Stork::IO::Input_SRDF_binary(srdf, fileStem);
Stork::IO::Output_SRDF_binary(srdf, fileStem);
Stork::IO::Input_RDF_binary(rdf, fileStem);
Stork::IO::Output_RDF_binary(rdf, fileStem);
Stork::IO::Output_RDF_csv(rdf, fileStem);
```

All binary functions append `.stork`; CSV output appends `.csv`. Reads populate
host views. Writes consume host views, so explicitly mirror/copy device results
first. Binary readers reject a mismatched RDF/SRDF tag or `float`/`double` tag.
See [File formats]({% link reference/file-formats.md %}) for the wire layout and
portability limits.

`Stork::IO::FileType` has values `RDF` and `SRDF`;
`Stork::IO::DataType` has values `FLOAT32` and `FLOAT64`. Applications normally
encounter these only through serialization validation.

## Complete interpolation example

```cpp
Stork::Structs::SRDF_Dual<double> srdf;
Stork::IO::Input_SRDF_binary(srdf, "thermal-layer");

auto rdf = Stork::Run::Interpolate_SRDF_to_RDF<
    double,
    Stork::host_space,
    float,
    Stork::device_space>(srdf, 4);

rdf.Make_Data_Mirrors<Stork::device_space, Stork::host_space>();
rdf.Copy_All<Stork::device_space, Stork::host_space>();
Stork::IO::Output_RDF_csv(rdf, "thermal-layer.interp.4");
```
