# Measurement-domain component removal in thoracic EIT

MATLAB code accompanying:

> Lin Y-D, Ko Y-F. *Measurement-domain component removal in thoracic electrical
> impedance tomography reconstructs the removed component, not the heart: a
> phantom study with a human recording.* Submitted to *Physiological
> Measurement*, 2026.

Three decompositions of a single EIT frame — PCA, EMD and DAGAF — are reduced
to one contract (components that sum exactly to the measured voltages, one of
which is removed) and the output of a linear one-step reconstructor is split
into the unprocessed difference image and the image of the removed component.
The paper shows that the central object reported as a heart image is the
second term. This repository holds the phantom, the three decompositions, the
component-selection rule, the figures of merit and the human split; the split
measurements on the phantom were pre-registered and run in a companion
repository (below).

An earlier version of this code accompanied a submission of the same study to
*Current Medical Imaging* (August 2026, desk-rejected). Commits before
October 2026 and the file names `eit_fig11_15_repro`, `eit_revision_r1`,
`eit_revision_r4_lung` refer to that draft's numbering; the table below gives
the mapping to the *Physiological Measurement* manuscript.

---

## Two repositories

| | Holds | Where |
|---|---|---|
| **This repository** (`EIT-DAGF`) | phantom, decompositions (`dagaf2.m`, `emd_rcada.m`, PCA inline), component selection, GREIT figures of merit with the negative-contrast adjustments, weight sweeps, the human reproduction figures, the human split, the signed Table 2 correlations, the assembly of Figures 2 and 9, and the literature-search record of Section 1.3 | here, MIT |
| **Companion repository** `EIT_noise_geometry_pre-registered` | the pre-registration, deviations D1–D9 and their registered results; the phantom split scripts `p2_separating_component*.m`, the displaced heart, the back-projection and Gauss–Newton reconstructors, the thorax, cardiac-reference and total-variation checks, the window sweep, `p2_probeA_blockscaling.m`; a `vendor/` copy of this repository's decomposition code | https://github.com/yfko/EIT_noise_geometry_pre-registered, **tag v1.0, commit 1e068214** |

Two scripts here read the companion repository's `results/*.mat` files
(`eit_table2_sign.m`, `make_fig2_fig9.m`). Clone it next to this repository
(`../EIT_noise_geometry_pre-registered`) or set the environment variable
`EIT_NOISE_GEOMETRY` to its root; both scripts resolve the path through the
local function `companion_root` and stop with a message if it is not found.

## What is and is not included

Included: the phantom, the three decompositions, the component-selection rule,
the evaluation, the human split, and the scripts that generate each figure and
table of this repository's part of the paper (mapping below).

Not included, and why:

* **The RCADA EMD implementation.** `emd_rcada.m` is a thin wrapper; the
  algorithm it calls is `eemd.m` by Zhaohua Wu, distributed from
  rcada.ncu.edu.tw. That site is offline and the package carries no licence we
  can redistribute under. See *Obtaining the EMD implementation* below.
* **The 2022 study's scripts.** The code of Cheng, Su and Ko (*Curr Med
  Imaging* 2022) has authors who are not authors of this paper. The phantom
  geometry and the human settings are reproduced from its Methods in
  `eit_show_phantom.m` and `eit_fig11_15_repro.m`.
* **EIDORS and Netgen.** Installed separately; see *Requirements*.
* **Results files of the companion repository.** Read from that repository,
  not copied here.

## Requirements

| | Version used | Notes |
|---|---|---|
| MATLAB | R2026a | No toolboxes required. `gausswin` is inlined in `dagaf2.m` so the Signal Processing Toolbox is not needed. |
| EIDORS | v3.8 | v3.8 specifically: the human settings follow the 2022 study, which was written against it, and the return order of `eval_GREIT_fig_merit` (AR, PE, RES, SD, RNG) was verified there. `montreal_data_1995` ships with it. |
| Netgen | 6.1 | Called by EIDORS as an external binary. On Apple Silicon the 6.1 build is x86_64 and needs Rosetta 2. |

On Apple Silicon, EIDORS v3.8 ships no `mexmaca64` build of `eidors_var_id`, and
MATLAB requires a full Xcode installation to compile one. A pure-MATLAB
replacement is described in `docs/eidors_apple_silicon.md`.

## Getting started

```matlab
setup_eit_env('/path/to/eidors')   % sets NETGENDIR, starts EIDORS, smoke-tests both
eit_compare_table1                  % Table 1, Supplementary Table S6 (reference mesh), Figure 3
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

Numbering is that of the *Physiological Measurement* manuscript. "Companion"
means the script lives in `EIT_noise_geometry_pre-registered`.

| Item | Script | Notes |
|---|---|---|
| Figure 1 (phantom) | `eit_show_phantom` | `figures/fig_phantom.png`; panel (a) is drawn from the mesh and also written as `figures/fig1a_phantom.png` (see the note below the table) |
| Figure 2 (the split, ten rows) | `make_fig2_fig9` | assembles companion `results/p2_separating_component*.mat` and `human_separating_component.mat` from `eit_human_separating_component`; writes `fig2_split_composite.png` and a consistency log |
| Figure 3 (contract output and centre cuts) | `eit_compare_table1` | `figures/fig_compare_profiles.png` |
| Figure 4 (DAGAF decomposition of one frame) | `eit_dagaf_images` | `figures/fig_dagaf_components.png` |
| Figure 5 (three weight states) | `eit_dagaf_imf1_weight` | `figures/fig_dagaf_imf1_states.png`; the sweep that fixes the boost at 1.06 |
| Figure 6 (weight swept 1 → 0) | `eit_partial_weight` | `figures/fig_partial_weight.png` |
| Figures 7 and 8 (human heart and lung directions) | `eit_fig11_15_repro` | `figures/fig11_repro_abcd.png`, `figures/fig15_repro_abcd.png`; reproduces the 2022 study's Figures 11 and 15 and adds DAGAF |
| Figure 9 (window sweep) | `make_fig2_fig9` | from companion `results/p2_separating_component_window.mat`; writes `fig9_window_sweep.png` |
| Table 1 (shape deformation after removing each component) | `eit_compare_table1` | |
| Table 2 (signed correlations with the lung-only and heart-only signatures) | `eit_revision_r1` (R1a, magnitudes), `eit_table2_sign` (signs) | reads companion `results/p2_env.mat` |
| Table 3, phantom and control rows | companion `p2_separating_component.m`, `p2_separating_component_ref.m` | |
| Table 3, human rows | `eit_human_separating_component` | 33 frames of `montreal_data_1995` breath-hold; writes `human_separating_component.mat`, `.txt`, `fig_human_separating_component.png` |
| Table 3, thorax rows; Table S9; Figure S6 | companion `p2_separating_component_thorax.m`, `_thorax_ref.m` | |
| Table 4 (displaced heart); Figure S1 | companion `p2_separating_component_offctr.m` | |
| Table S2 (lung figures of merit, with mesh-density ranges) | `eit_lung_table`, `eit_revision_r4_lung` | the negative-contrast adjustments of Section 2.6 |
| Table S3 (window cap) | `eit_dagaf_ringing`, `eit_revision_r1` (R2) | |
| Table S4 (fixed linear high-pass controls) | `eit_revision_r1` (R3) | |
| Table S5 (collinearity thresholds) | companion | |
| Table S6 (heart figures of merit, with mesh-density ranges) | `eit_compare_table1`, `eit_revision_r1` (R4) | |
| Table S7, Figure S2 (four reconstructors) | companion `p2_separating_component_recon.m` | |
| Table S8 (transfer) | companion | |
| Figure S4 (cardiac-reference check), Figure S5 (all sixteen window images) | companion `p2_separating_component_ecg.m`, `_window.m` | |
| Section 3.4, hypothesis 2 (injection axis) | `eit_axis_test`, `dagaf_axis` | |
| Section 4.2 (window-rule probe) | companion `p2_probeA_blockscaling.m` | |
| Heart conductivity 1.5 against 2.0 | `eit_sigma_check` | all phantom scripts use σ_heart = 2.0, the value stated in the 2022 study's Methods; this check shows the table published there corresponds to 1.5 |
| Section 1.3 (literature search) | `docs/literature_search_2026-08-24.md` | the two PubMed queries, the OpenAlex forward-citation screening lists, and a re-run of 31 October 2026 (`docs/literature_search_rerun_2026-10-31.json`) |

`eit_lung_table` also contains the mask diagnosis of Section 2.6, and
`eit_graphical_abstract` draws a graphical abstract that is not part of the
manuscript.

**Figure 1, panel (a).** In the submitted manuscript this panel is a 3-D
render first published as Fig. (2) of Cheng, Su and Ko 2022 (*Curr Med
Imaging* 18:1396). Since 31 October 2026 this repository no longer carries that
image: `eit_show_phantom` draws panel (a) from the phantom mesh it builds, and
`eit_graphical_abstract` uses the same file. The model is the same (cylinder
r = 2, h = 2, lungs r = 0.6 at x = ±1, heart r = 0.3, sixteen electrodes); only
the rendering differs. The other four panels are unchanged in content; MATLAB's
automatic tick spacing and Netgen's re-meshing can move a few pixel edges.

## The contract

Each decomposition returns a `K × 208` matrix whose rows sum exactly to the
measured vector. Separation is then a weight vector: zero on the *separating*
row for the heart direction, above one for the lung direction. Because the
three methods share this interface, the comparison isolates the decomposition
and nothing else, and the output of a linear one-step reconstructor splits
exactly into the unprocessed difference image and the image of the removed
component (Section 2.7 of the paper; `eit_human_separating_component.m` and
the companion `p2_separating_component.m`).

The separating row is the one whose removal minimises GREIT shape deformation.
It is called separating rather than lung-dominant — the term the 2022 study
used — because Section 3.1 of the paper shows that the removed component is not
the lung component for any of the three methods.

`dagaf2.m` is a re-implementation of the published DAGAF algorithm from its
description (Lin, Tan and Tian 2022) and the reference release. Three
differences from the reference release are deliberate and are documented in the
file header: α defaults to 3 as in the paper rather than 4.0728; the
energy-ratio stop criterion follows the published equation; and the boundary
extension argument is honoured, which it is not in the reference release. A
`maskCap` argument was added, and it is required here — the published adaptive
window rule degenerates on a 208-point measurement vector, as Section 2.2 of the
paper explains.

## Human data

`montreal_data_1995` ships with EIDORS. `eit_fig11_15_repro` and
`eit_human_separating_component` load it directly; no download is needed.

## Licence

MIT — see `LICENSE`.

The code calls EIDORS rather than incorporating its source, so a permissive
licence applies here even though EIDORS itself is GPL. Anyone running these
scripts is separately bound by the EIDORS and Netgen licences, and by whatever
terms attach to the third-party EMD implementation they obtain.

## Citation

Lin Y-D, Ko Y-F (2026). *Measurement-domain component removal in thoracic
electrical impedance tomography reconstructs the removed component, not the
heart: a phantom study with a human recording.* Submitted to *Physiological
Measurement*. Code: https://github.com/yfko/EIT-DAGF. A DOI for the paper will
be added on publication; a Zenodo archive of this repository will be made at
acceptance.
