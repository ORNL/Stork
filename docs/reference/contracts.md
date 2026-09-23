---
title: Contracts and Limitations
parent: Reference
nav_order: 4
---

# Contracts and Limitations

Stork is a low-level coupling library. Its current API favors direct Kokkos
views and explicit movement over defensive validation. Producers and consumers
must uphold these contracts.

## Object invariants

- Initialize Kokkos before constructing Stork headers or dual containers.
- `numSnaps` and `numEvents` are logical counts maintained separately from
  view extents. Keep them synchronized when constructing objects manually.
- A regular-grid header describes nodes. SRDF `p(n)` identifies the low corner
  of a cell, so all eight `(0/1,0/1,0/1)` neighbors must exist.
- Local flat indexing is `p = i*(N_j*N_k) + j*N_k + k`.
- Origins, resolution, times, temperatures, and cooling rates carry no unit
  metadata. Every component must use one consistent unit system.
- The current implementation uses 32-bit unsigned indices and counts. Grid
  products, refined extents, and event counts must fit in `uint32_t`.
- Interpolation temporarily stores refined local flat point IDs in `int`.
  Consequently, the maximum refined local `p` must also fit in `INT_MAX`, even
  though the final RDF field is `uint32_t`.

## Interpolation preconditions

- `numSnaps > 0`.
- `fineFactor` is in `[1,255]`.
- Coarse grid extents and every sparse point index describe valid cells.
- Every crossing interval has `t_cur > t_prev` and distinct endpoint
  temperatures. Forward time is required for chronological crossings and the
  documented positive cooling-rate magnitude.
- The sparse history contains pairable melting and solidification crossings.
- Source views are populated in the template's declared source execution
  space; both source and target are exactly `Stork::host_space` or
  `Stork::device_space`.
- `srdf.host_header` is valid even when `SourceSpace` is `device_space`, because
  interpolation reads the host header before any source-to-target copy.

Stork detects only unequal **total** melt and solidification counts. It does not
separately validate matching per-point cycle counts after sorting. Validate
producer output before relying on event pairing.

## Memory-space contract

`Make_Data_Mirrors` allocates but does not copy. `Copy_Header`, `Copy_Data`, and
`Copy_All` copy but assume compatible destination allocations. A common
host-to-device sequence is:

```cpp
srdf.Make_Data_Mirrors<Stork::host_space, Stork::device_space>();
srdf.Copy_All<Stork::host_space, Stork::device_space>();
```

The interpolation functions always invoke the source-to-target mirror/copy
methods for SRDF, so their input is non-const and **may** be mutated. When
`SourceSpace` and `TargetSpace` differ, Stork allocates/overwrites target-side
SRDF data handles and deep-copies the source header/data. When the types are
the same, it takes the same-space branch: that branch is effectively a no-op
when configured host/device types are distinct, while a host-only
configuration may refresh view aliases. The functions return RDF data in the
target space and do not automatically make a host copy for output.

Even when `host_space` and `device_space` are the same type, do not assume every
operation preserves handle aliasing. Binary readers reallocate and populate
only host data views; call `Make_Data_Mirrors` afterward to refresh the device
handles before using them.

## Current edge cases

| Operation | Current behavior |
|:--|:--|
| Empty interpolation input | Unsupported; the scan reads its last element |
| Unequal total melt/solidification counts | Writes `DEBUG-STORK-NUMEVENTS.csv` in the working directory, then throws |
| Empty RDF z trim | Unsupported; `Trim_RDF_in_Z` assumes an event exists |
| Empty time normalization input | Unsupported |
| Reversed/out-of-domain trim bounds | Not explicitly validated |
| Trimmed interpolation | `numEvents` is reduced, but allocated view extents remain at the untrimmed size |
| `Free_Data_Memory` | Does not reset the logical record count |
| Host-only backend | Same-space helper calls can alias handles and freeing is a no-op, but a binary read refreshes only host data handles |

## I/O behavior

- All binary I/O is host-side and appends `.stork` to the supplied stem.
- Binary reads validate format and floating precision.
- Binary storage is native-endian/native-representation and unversioned.
- `Output_SRDF_binary` size depends on `numSnaps`
- `Output_RDF_binary` size depends on `numEvents`