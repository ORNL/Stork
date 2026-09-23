---
title: Interpolation
parent: Reference
nav_order: 2
---

# Interpolation

The primary algorithm is
`Stork::Run::Interpolate_SRDF_to_RDF`. It converts sparse coarse-grid
temperature intervals into fine-grid phase-change events.

## Spatial interpolation

For normalized position \((x_d,y_d,z_d)\) inside a coarse cell, Stork applies
trilinear interpolation to its eight corner temperatures. The operation is
performed once at `t_prev` and once at `t_cur`.

The fine factor `f` defines normalized positions `0, 1/f, ..., (f-1)/f` inside
each cell. The positive boundary is included only for the last coarse cell on
an axis, preventing duplicate fine-grid points where neighboring cells meet.

## Crossing classification

Let `T_prev` and `T_cur` be the spatially interpolated endpoint temperatures and
`T_critical` the SRDF threshold.

| Condition | Event |
|:--|:--|
| `T_prev < T_critical` and `T_cur >= T_critical` | melting |
| `T_prev >= T_critical` and `T_cur < T_critical` | solidification |
| otherwise | no event |

Stork linearly interpolates the crossing time within the interval. For a
solidification crossing it also calculates the positive cooling-rate magnitude.

## Assembly

The implementation uses this parallel sequence:

1. Count melt and solidification crossings contributed by each SRDF record.
2. Inclusive-scan the counts to calculate non-overlapping output offsets.
3. Fill separate melt and solidification event arrays.
4. Sort both arrays by fine-grid point and then crossing time.
5. Pair sorted events into RDF records.

The total melt and solidification counts must match. Stork raises an error and
writes `DEBUG-STORK-NUMEVENTS.csv` when they differ.

The two sorted arrays are paired by position, not by a second explicit
per-point consistency check. A producer must therefore supply well-formed
cycles: after sorting, every melt event must correspond to the solidification
event at the same point and cycle number. Matching only the global totals is
not sufficient to establish that condition.

## Trimmed interpolation

`Interpolate_And_Trim_SRDF_to_RDF` performs interpolation and retains events
whose global fine-grid `i` and `j` indices fall within inclusive bounds. The
returned RDF header is updated to describe the trimmed x-y region.

The bounds are global **fine-grid** indices, not coarse SRDF indices. They must
be ordered (`min <= max`) and lie inside the interpolated grid. The function
sets `numEvents` to the retained count but currently leaves the event views at
the pre-trim allocation size. Consumers must iterate to `numEvents`, not to
`view.extent(0)`.

{: .warning }
`Output_RDF_binary` currently derives its record count from the host view
extent. Do not serialize an RDF returned by trimmed interpolation unless its
views have first been compacted to `numEvents`. The Toucan integration consumes
the trimmed object in memory and iterates the logical count.

`Trim_RDF_in_Z` is a separate operation that removes unused fine-grid depth
below the lowest event and remaps local point indices. It requires at least one
event. It changes the header in the selected execution space and the host
header used to calculate the new origin; synchronize the header afterward when
those spaces differ.

## Numerical considerations

- `fineFactor` is an unsigned 8-bit integer and must be in `[1, 255]`. The
  public API does not reject zero before dividing by it.
- `numSnaps` must be greater than zero, and each grid extent must describe at
  least one coarse cell (normally at least two nodes per axis).
- Every interval used for crossing-time calculation must have
  `t_cur > t_prev`. Merely nonzero duration is insufficient: reversed time
  breaks chronological crossing interpretation and the positive cooling-rate
  convention.
- A crossing interval must have different endpoint temperatures.
- The sparse producer must retain both melting and corresponding solidification
  intervals.
- `p(n)` must identify the low-index corner of a valid coarse cell, never a
  high-boundary node with no cell above it.
- Refined local flat point IDs must fit in `INT_MAX`; the interpolation event
  structures temporarily store `p` as signed `int` before RDF stores it as
  `uint32_t`.
- Source data must be initialized in the declared `SourceSpace`, which must be
  exactly `Stork::host_space` or `Stork::device_space`; `TargetSpace` has the
  same two supported choices.
- `SRDF.host_header` must be initialized regardless of `SourceSpace`; both
  interpolation entry points read it before copying to `TargetSpace`.
- The non-const SRDF input may be mutated. Both functions invoke mirror/copy
  methods; target-side handles are allocated/overwritten and data are
  deep-copied only when `SourceSpace` and `TargetSpace` differ. The same-space
  branch is a no-op for distinct configured host/device types and may refresh
  aliases in a host-only configuration.
- The result is not automatically copied to the host when `TargetSpace` is a
  device execution space.
- Accuracy depends on coarse-grid spacing, time sampling, thermal gradients,
  and the behavior between samples—not only the fine factor.

The interpolation path allocates count, scan, melt, solidification, and output
arrays proportional to the sparse input and detected events. It performs two
passes over all candidate fine points plus parallel scans and sorts. Increasing
the refinement factor therefore increases both work and possible event count
approximately with refined volume.

The paper's factor-of-16 result is a validation result for its LPBF studies. A
new workflow should perform its own convergence study.
