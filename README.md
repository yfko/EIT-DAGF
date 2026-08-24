# Measurement-domain decomposition for heart/lung separation in EIT

MATLAB code accompanying:

> Lin Y-D, Ko Y-F. *A quantitative framework for comparing measurement-domain
> decompositions in electrical impedance tomography: principal component
> analysis, empirical mode decomposition and data-adaptive Gaussian average
> filtering.* Current Medical Imaging (submitted).

Three decompositions — PCA, EMD and DAGAF — are reduced to one contract and
compared under identical conditions on a finite element phantom and on a public
human recording. Everything reported in the paper is produced by the scripts
here.

---

## What is and is not included

Included: the phantom, the three decompositions, the component-selection rule,
the evaluation, and the scripts that generate each figure and table.

Not included, and why:

* **The RCADA EMD implementation.** `emd_rcada.m` is a thin wrapper; the
  algorithm it calls is `eemd.m` by Zhaohua Wu, distributed from
  rcada.ncu.edu.tw. That site is offline and the package carries no licence we
  can redistribute under. See *Obtaining the EMD implementation* below.
* **The prior study's scripts.** The 2022 study's code has authors who are not
  authors of this paper.
* **EIDORS and Netgen.** Installed separately; see *Requirements*.

## Requirements

| | Version used | Notes |
|---|---|---|
| MATLAB | R2026a | No toolboxes required. `gausswin` is inlined in `dagaf2.m` so the Signal Processing Toolbox is not needed. |
| EIDORS | v3.8 | v3.8 specifically: the prior study's scripts were written against it, and the return order of `eval_GREIT_fig_merit` (AR, PE, RES, SD, RNG) was verified there. |
| Netgen | 6.1 | Called by EIDORS as an external binary. On Apple Silicon the 6.1 build is x86_64 and needs Rosetta 2. |

On Apple Silicon, EIDORS v3.8 ships no `mexmaca64` build of `eidors_var_id`, and
MATLAB requires a full Xcode installation to compile one. A pure-MATLAB
replacement is described in `docs/eidors_apple_silicon.md`.

## Getting started

```matlab
setup_eit_env('/path/to/eidors')   % sets NETGENDIR, starts EIDORS, smoke-tests both
eit_compare_table1                  % Tables 3 and 4, and Figure 2
```

`setup_eit_env` verifies that Netgen can actually build a mesh and that the
DAGAF decomposition closes to numerical precision before reporting success, so a
broken dependency surfaces immediately rather than three hundred lines later.

## Obtaining the EMD implementation

`emd_rcada.m` expects `eemd.m` and its helpers (`extrema.m`, `ifndq.m`,
`dist_value.m`, `significance.m`) on the MATLAB path, in a subfolder named
`rcada_emd`. They were published at `rcada.ncu.edu.tw`, which is no longer
reachable; archived copies of the 2009-06-21 snapshot can be retrieved from the
Internet Archive. Place them in `src/rcada_emd/`.

The wrapper calls `eemd(Y, 0, 1)`, which reduces ensemble EMD to plain EMD, and
returns the modes as rows with the residual separated, so that the rows sum
exactly to the input. For a 208-point vector the package returns six modes,
since it fixes the count at `fix(log2(N)) - 1`.

## Reproducing the paper

| Script | Produces |
|---|---|
| `eit_show_phantom` | Figure 1 — mesh, electrode plane, conductivity maps, centre cuts |
| `eit_compare_table1` | Table 2 (component screening), Table 3 (heart), Table 4 (lung), Figure 2 |
| `eit_dagaf_images` | Figure 3 — DAGAF decomposition of one frame |
| `eit_dagaf_imf1_weight` | Figure 4 — the three weight states, and the sweep that fixes the boost at 1.06 |
| `eit_partial_weight` | Figure 5 — weight sweep, both axes logarithmic |
| `eit_dagaf_ringing` | Table 5 — window cap against the figures of merit |
| `eit_axis_test` | Section 3.6, hypothesis 2 — the injection-axis test |
| `eit_lung_table` | Section 3.4 and the mask diagnosis of Section 2.6 |
| `eit_sigma_check` | Section 4.4 — heart conductivity 1.5 against 2.0 |
| `eit_fig11_15_repro` | Figures 6 and 7 — human data |

All phantom scripts use σ_heart = 2.0, the value stated in the prior study's
Methods. `eit_sigma_check` shows that the table published there corresponds to
1.5; see Section 4.4 of the paper.

## The contract

Each decomposition returns a `K × 208` matrix whose rows sum exactly to the
measured vector. Separation is then a weight vector: zero on the lung-dominant
row for heart imaging, above one for lung imaging. Because the three methods
share this interface, the comparison isolates the decomposition and nothing
else.

`dagaf2.m` is a corrected re-implementation of the published DAGAF algorithm.
Three differences from the reference release are deliberate and are documented
in the file header: α defaults to 3 as in the paper rather than 4.0728; the
energy-ratio stop criterion follows the published equation; and the boundary
extension argument is honoured, which it is not in the reference release. A
`maskCap` argument was added, and it is required here — the published adaptive
window rule degenerates on a 208-point measurement vector, as Section 2.2 of the
paper explains.

## Human data

`montreal_data_1995` ships with EIDORS. `eit_fig11_15_repro` loads it directly;
no download is needed.

## Licence

[TODO — choose. MIT or BSD-3 is conventional for MATLAB research code and is
compatible with EIDORS, which is GPL-2/3. Note that the code does not
incorporate EIDORS source, it calls it, so a permissive licence is available.]

## Citation

[TODO — add once the paper has a DOI.]
