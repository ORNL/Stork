---
title: Reference
nav_order: 3
has_children: true
---

# Reference

This section describes the interfaces implemented by the current Stork source.

| Reference | Contents |
|:--|:--|
| [Data model](data-model/) | Regular-grid header, SRDF records, and RDF events |
| [Interpolation](interpolation/) | Refinement, crossing detection, ordering, and validation |
| [File formats](file-formats/) | Binary `.stork` and RDF CSV layouts |
| [Contracts and limitations](contracts/) | Preconditions, memory ownership, errors, and current edge cases |
| [Library API](api/) | Public types and functions in `Stork_Core.hpp` |

Stork is template-based and header-only. Include `Stork_Core.hpp` for the
supported aggregate interface or individual headers when minimizing includes is
important.
