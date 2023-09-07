# Learning Dictionaries from Physical-Based Interpolation for Water Network Leak Localization

Implementation and experimental data for the paper

> P. Irofti, L. Romero-Ben, F. Stoican, and V. Puig,
“Learning Dictionaries from Physical-Based Interpolation
for Water Network Leak Localization,”
pp. 1--12, 2023.

If you use our work in your research, please cite as:
```
@article{IRSP23_awgsi,
  author = {Irofti, P. and Romero-Ben, L. and Stoican, F. and Puig, V.},
  title = {Learning Dictionaries from Physical-Based Interpolation
           for Water Network Leak Localization},
  year = {2023},
  pages = {1-12},
  eprint = {2304.10932},
  archiveprefix = {arXiv},
}
```

## Prerequisite
[Dictionary Learning Toolbox](https://github.com/pirofti/dl-box) implementation
by Paul Irofti.

## Usage
1. **Initialization:** Run `setupDL.m` from the Dictionary Learning Toolbox
  (required only once per session).

2. **Interpolation:** obtain the interpolated hydraulic heads and residuals
   with [launcher_interpolation.m](code/launcher_interpolation.m) 

2. **Learning:** perform dictionary learning for FDI
   with [launcher_DL.m](code/launcher_DL.m)
