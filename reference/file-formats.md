---
title: File Formats
parent: Reference
nav_order: 3
---

# File Formats

Stork supports binary serialization for SRDF and RDF objects and CSV output for
RDF events.

## Binary `.stork`

All binary readers and writers accept a file stem and append `.stork`.

Every file begins with two one-byte tags:

| Byte | Values |
|:--|:--|
| file type | `0` = RDF, `1` = SRDF |
| data type | `0` = 32-bit float, `1` = 64-bit float |

### SRDF sequence

```text
FileType
DataType
numSnaps                  uint32
T_critical                FloatType
header indices            6 × uint32
header floating values    5 × FloatType
cell indices              numSnaps × uint32
time endpoints            2 × numSnaps × FloatType
corner temperatures       16 × numSnaps × FloatType
```

### RDF sequence

```text
FileType
DataType
numEvents                 uint32
header indices            6 × uint32
header floating values    5 × FloatType
cell indices              numEvents × uint32
event values              3 × numEvents × FloatType
```

The three RDF event values are `tm`, `tl`, and `cr`.

{: .warning }
The current binary representation has no magic number, schema version, byte
order conversion, or integrity checksum. It writes native primitive values.
Treat `.stork` as an efficient interchange format between compatible builds,
not as a permanent archival format. Retain the producing Stork revision and
floating-point type with archived datasets.

The reader verifies both the file-type and floating-point tags. Reading a
double-precision file into `SRDF_Dual<float>`, for example, fails rather than
silently converting values.

Readers reallocate and populate only host data views; they do not create or
refresh device data handles. This remains true when the configured host and
device execution spaces are the same type, so do not rely on an earlier alias
surviving a read. Call `Make_Data_Mirrors<host_space, device_space>()` and
`Copy_All<host_space, device_space>()` before device access. Writers likewise
read host views; mirror/copy device results to the host first.

`Output_SRDF_binary` writes exactly `numSnaps` records even if its host views
have spare capacity. `Output_RDF_binary`, by contrast, currently uses the host
cell-index view extent as its stored count. For ordinary RDF construction those
values match; after trimmed interpolation they may not. See
[Interpolation]({% link reference/interpolation.md %}).

Readers validate stream success only after allocating according to the stored
record count. They do not validate semantic grid bounds, event ordering, or
reasonable allocation sizes. Read only trusted `.stork` files.

## RDF CSV

`Output_RDF_csv` writes:

```text
x,y,z,tm,tl,cr
0.0001,0.0002,-0.00003,0.0011,0.0014,120000
```

Coordinates are calculated from the RDF header and local point index. CSV is
useful for inspection and plotting, but it does not retain the complete grid
header or numeric type and is less efficient than in-memory or binary coupling.

`Output_RDF_csv` catches failures thrown while opening or writing and prints a
message. Its exception-enabled `close()` occurs after that catch and can itself
throw, so callers must not assume every failure merely prints and returns.
Ensure the parent directory exists and verify the resulting file in automated
workflows.
