---
title: Paper
nav_order: 5
---

# Paper

## Methodology paper

The mathematical framework, performance study, and microstructure validation
are presented in:

> Benjamin Stump and John Coleman. “A scalable framework for efficient coupling
> of thermal and microstructural simulations in additive manufacturing.”
> *Computational Materials Science* 272 (2026), article 114832.
> [doi:10.1016/j.commatsci.2026.114832](https://doi.org/10.1016/j.commatsci.2026.114832)

```bibtex
@article{stump2026stork,
  title   = {A scalable framework for efficient coupling of thermal and
             microstructural simulations in additive manufacturing},
  author  = {Stump, Benjamin and Coleman, John},
  journal = {Computational Materials Science},
  volume  = {272},
  pages   = {114832},
  year    = {2026},
  doi     = {10.1016/j.commatsci.2026.114832}
}
```

## Reported findings

The paper demonstrates Stork by coupling the semi-analytic thermal solver
3DThesis with the time-parallel cellular automata solver Toucan. For the studied
workflows it reports:

- more than two orders of magnitude reduction in thermal-data generation time
  and file size relative to the prior workflow;
- file-based CPU and GPU-resident in-memory implementations; and
- preservation of grain morphology and crystallographic texture for
  interpolation factors through 16 in the reported LPBF studies.

These are results for the published configurations rather than guarantees for
every thermal field or microstructure model. Document the convergence study
used for your application.

## License

Stork source and these original documentation assets are distributed under the
repository's BSD 3-Clause License.
