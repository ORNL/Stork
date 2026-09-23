---
title: Tutorials
nav_order: 4
has_children: true
---

# Tutorials

The tutorials show the two supported coupling styles.

| Tutorial | Use it when |
|:--|:--|
| [Convert an SRDF file](file-conversion/) | A thermal producer and microstructure consumer run as separate programs |
| [Couple in memory](in-memory/) | Both stages can share Stork objects in one application or library stack |

The examples deliberately use small synthetic data or a single layer. A
production thermal solver is responsible for selecting a complete sparse set of
active intervals.
