---
title: User Guide
nav_order: 2
has_children: true
---

# User Guide

Use this guide to build Stork, convert sparse thermal data, and understand where
the library fits in a coupled simulation.

1. [Install Stork](installation/) against a configured Kokkos installation.
2. Follow the [quick start](quick-start/) to inspect the converter workflow.
3. Read [How Stork works](how-stork-works/) before choosing an interpolation
   factor or implementing a new producer.

Stork's public interface is provided by `Stork_Core.hpp`. The
[reference](../reference/) documents the data structures, algorithms, file
layout, and callable functions exposed by the current source tree.
