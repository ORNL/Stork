---
title: Quick Start
parent: User Guide
nav_order: 2
---

# Quick Start

The repository includes `ReadAndInterpolate`, a small application that reads an
SRDF binary file, interpolates it, and writes the resulting RDF events to CSV.

## Command

Build Stork, then run:

```sh
./install/bin/ReadAndInterpolate \
  /path/to/thermal-layer \
  4 \
  double \
  float
```

The arguments are:

| Position | Meaning |
|:--|:--|
| `thermal-layer` | File stem; Stork reads `thermal-layer.stork` |
| `4` | Fine factor; the output spacing is coarse spacing divided by 4 |
| `double` | Floating-point type stored in the SRDF file |
| `float` | Floating-point type used in the output RDF object |

Both precision strings must be exactly `float` or `double`. Supply a fine factor
in `[1,255]` and all four arguments after the executable name.

The command writes:

```text
/path/to/thermal-layer.interp.4.csv
```

Each output row contains `x,y,z,tm,tl,cr`: spatial coordinates, melting time,
liquidus-crossing time during solidification, and cooling rate.

The utility interpolates the complete SRDF, copies its RDF result to the host,
and writes CSV. It does not write binary RDF, apply x-y/z trimming, normalize
times, or consume the result in a microstructure model. Use the library API for
those workflows.

{: .important }
The file argument is a stem. Do not include the `.stork` suffix; the reader adds
it automatically.

{: .warning }
The current executable's argument guard is less strict than its implementation.
Always supply all four arguments; invoking it with only a filename and fine
factor can access missing precision arguments. It also parses the fine factor
as an unrestricted signed `int` and then narrows it to `uint8_t` without range
validation. Values outside `[1,255]` can wrap; zero can lead to division by
zero. The caller must enforce the range.

## Choose a fine factor

For a coarse grid spacing \\(\Delta x_c\\) and integer fine factor \\(f\\), Stork
creates a fine spacing

\\[
\Delta x_f = \frac{\Delta x_c}{f}.
\\]

For an axis containing \\(N_c\\) coarse nodes, the corresponding fine-grid extent
is \\(f(N_c-1)+1\\) nodes. Larger factors increase interpolation work and the
number of potential phase-change events approximately with the refined volume.

The Stork paper reports preserved morphology and texture through a factor of 16
for its numerical studies. Treat that result as evidence for the studied cases,
not as a universal accuracy limit. Compare the RDF events and downstream
microstructure against a sufficiently resolved thermal reference for each new
material and process regime.

## In-memory use

Production coupling does not require an intermediate `.stork` file. A thermal
solver can fill `Stork::Structs::SRDF_Dual<T>` and call
`Stork::Run::Interpolate_SRDF_to_RDF` directly. See the
[library API]({% link reference/api.md %}) and
[in-memory tutorial]({% link tutorials/in-memory.md %}).
