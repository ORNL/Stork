---
title: How Stork Works
parent: User Guide
nav_order: 3
---

# How Stork Works

Process-scale thermal models and grain-scale microstructure models resolve very
different length and time scales. Storing every temperature at the fine scale
is wasteful because only a small region near a moving melt pool is actively
melting or solidifying at a given time.

Stork moves information between the two models in three stages.

## 1. Select active coarse cells

The thermal producer identifies coarse cells and time intervals that can contain
a crossing of the critical temperature. Each retained SRDF record stores:

- the coarse cell index;
- the beginning and ending time;
- eight corner temperatures at the beginning; and
- eight corner temperatures at the end.

This sparse selection belongs to the producer. Stork assumes its input contains
the records needed to recover every relevant melt and solidification event.

## 2. Interpolate space and time

For each fine point inside a retained coarse cell, Stork performs trilinear
interpolation at both time endpoints. If those endpoint temperatures straddle
the critical temperature, linear interpolation in time locates the crossing.
Together, the three spatial dimensions and time form the paper's quad-linear
interpolation.

For endpoint temperatures \\(T_0\\) and \\(T_1\\), endpoint times \\(t_0\\) and
\\(t_1\\), and critical temperature \\(T_c\\), Stork requires \\(t_1>t_0\\). The
crossing time is

\\[
t_c = t_0 + \frac{T_c-T_0}{T_1-T_0}(t_1-t_0).
\\]

During a forward-time cooling interval, Stork records the positive cooling-rate
magnitude

\\[
\dot{T}_{cool} = \frac{T_0-T_1}{t_1-t_0}.
\\]

## 3. Assemble phase-change histories

Stork detects upward crossings as melting events and downward crossings as
solidification events. It sorts each event list first by fine-grid point and
then by crossing time, pairs corresponding events, and emits one RDF record per
melt/solidification cycle.

The current interpolation implementation requires the total numbers of melting
and solidification events to match. A mismatch usually indicates an incomplete
thermal interval—for example, a retained cell that melts but whose later
solidification was not provided.

## Performance portability

The interpolation loop, scans, sorting, and data arrays use Kokkos. The same
public function can therefore target the configured host or device execution
space without maintaining separate CPU and GPU implementations. `RDF_Dual` and
`SRDF_Dual` provide explicit host/device mirror and copy operations so a caller
controls data movement.

In this documentation, performance portability means that Stork expresses the
same data-parallel algorithm through Kokkos and obtains its concrete backend
from the application's Kokkos build. It does not mean that one compiled binary
automatically contains every backend, nor that every backend has equal
performance. Build and tune Kokkos for the intended CPU or accelerator.

An in-memory device path is:

```text
thermal producer writes SRDF device views
                  |
                  v
Stork counts, scans, interpolates, and sorts in device_space
                  |
                  v
microstructure consumer reads RDF device views
```

No host copy is required in that path. File I/O, CSV inspection, and host-side
consumers require an explicit mirror and copy. On a CPU-only Kokkos build,
`host_space` and `device_space` can be identical, allowing Stork's same-space
helpers to share handles. A binary reader still reallocates only host data
views, so refresh the device handles explicitly afterward.

## Integration boundary

The thermal producer decides which coarse cells and intervals enter SRDF and
must provide complete phase-change cycles. Stork creates fine-grid events but
does not decide an MPI decomposition. The consumer decides which events belong
to each microstructure subdomain and how `tm`, `tl`, and `cr` drive its model.

## Relationship to the paper

The [Stork paper]({% link about/paper.md %}) demonstrates the method by coupling
the semi-analytic thermal solver 3DThesis to the time-parallel cellular automata
solver Toucan. The paper's performance and validation results apply to those
studies.
