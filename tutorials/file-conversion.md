---
title: Convert an SRDF File
parent: Tutorials
nav_order: 1
---

# Convert an SRDF File

Use the repository application to inspect an SRDF file or create CSV input for
analysis tools.

## 1. Build the converter

Follow [Installation]({% link docs/installation.md %}). The executable is
installed at `install/bin/ReadAndInterpolate`.

## 2. Identify the stored precision

The producer chooses either 32-bit `float` or 64-bit `double`. The binary file
contains a type tag and Stork rejects a mismatched command-line type.

Assume the producer wrote `layer-0007.stork` using `double`.

## 3. Interpolate

To refine every coarse-grid interval by a factor of 8 and store RDF values as
single precision:

```sh
./install/bin/ReadAndInterpolate layer-0007 8 double float
```

Supply an integer fine factor in `[1,255]`. The current utility narrows its
parsed `int` to `uint8_t` without validating that range.

Expected terminal output includes the number of RDF events:

```text
NumEvents: <event-count>
```

The generated `layer-0007.interp.8.csv` contains:

```text
x,y,z,tm,tl,cr
```

## 4. Check the result

Before passing a new dataset to a microstructure model, verify:

- event coordinates fall inside the expected physical domain;
- `tm <= tl` for every event;
- crossing times fall inside the thermal simulation window;
- cooling rates use the expected time and temperature units; and
- event count and spatial distribution converge as the coarse thermal data is
  refined.

{: .note }
Stork does not attach unit metadata. Coordinates, times, temperatures, and
cooling rates retain the producer's consistent unit system.

## Common failures

`Invalid data type`
: The `precisionIn` argument does not match the binary type tag.

`Number of melting and solidification events are not equal`
: The sparse history contains an unmatched transition. Check that every
  retained melt cycle includes the subsequent cooling interval.

`Failed to open file`
: Pass the file stem and run from a directory where the relative path resolves.
  Do not include `.stork` in the argument.
