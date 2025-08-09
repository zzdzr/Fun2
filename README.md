# Fun2: Tracking the untrackable: A reinforcement learning framework for quantitative analysis of chromatin fountains/stripes

<!-- Badges (可添加 Zenodo/License/Build Status) -->
<!-- ![DOI](https://zenodo.org/badge/641802007.svg) -->

**Fun2** is a **reinforcement learning framework** designed to dynamically track **replication fountain structures**.

<div style="margin-bottom:200px; text-align:center;">
  <img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/workingModel4.svg" 
       width="1400" height="600" 
       style="border:none; display:block; margin:0 auto;">
</div>

> **Highlights**
> - Sampling box / axes control for tiled, multi-resolution scans  
> - Robust affine transforms for coordinate normalization  
> - **MCTS** (Monte Carlo Tree Search) planning for trajectory search  
> - Reproducible pipelines and exportable outputs (HDF5/CSV/figures)
---

# Table of Contents
- [Fun2](#fun2-characterize-fountain--stripe-like-chromatin-structures)
  - [Description](#description)
- [Table of Contents](#table-of-contents)
- [Getting Started / Installation](#getting-started--installation)
  - [Quick start](#quick-start)
  - [Environment notes](#environment-notes)
- [Usage](#usage)
  - [Sampling Box & Axes Configuration](#sampling-box--axes-configuration)
  - [Affine Transformation](#affine-transformation)
  - [MCTS Planning](#mcts-planning)
  - [Working Model](#working-model)
- [Examples](#examples)
  - [Minimal end-to-end run](#minimal-end-to-end-run)
  - [Batch over chromosomes](#batch-over-chromosomes)
- [Practical Considerations](#practical-considerations)
- [Understanding the Output Files](#understanding-the-output-files)
- [Repository Structure](#repository-structure)
- [License & Citation](#license--citation)

---

## Description
Replication forks in mammalian cells often progress in a **coupled** manner. Under replication stress, **type-II fountains** exhibit distinct uncoupling degrees. **Fun2** quantitatively tracks these dynamics and relates uncoupling to higher-order chromatin structures (e.g., loops), replication timing shifts, and fork-restart defects. The toolkit is intended for large-scale screening, hypothesis testing, and figure-ready output.

<div align="right">[⬆ Back to top](#table-of-contents)</div>

---

# Getting Started / Installation

### 🚀 Quick start
```bash
# 1) Clone
git clone https://github.com/zzdzr/Fun2.git
cd Fun2

# 2) Create environment (recommended)
conda env create -f environment.yml
conda activate fun2

# 3) Install in development mode
pip install -e .
```
---

## 🖼️ Details


## Sampling box and axes configuration
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/axis.svg" alt="SamplingBox" width="500" height="500" align="left"/>

We have recently published the RFpeptides protocol for using RFdiffusion to design macrocyclic peptides that bind target proteins with atomic accuracy (Rettie, Juergens, Adebomi, et al., 2025). In this section we briefly outline how to run this inference protocol. We have added two examples for running macrocycle design with the RFpeptides protocol. One for monomeric design, and one for binder design.

<br clear="all"/>

## Affine Transformation
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/affineTransform.svg" alt="Transformation" width="300" height="300" align="left"/>

We have recently published the RFpeptides protocol for using RFdiffusion to design macrocyclic peptides that bind target proteins with atomic accuracy (Rettie, Juergens, Adebomi, et al., 2025). In this section we briefly outline how to run this inference protocol. We have added two examples for running macrocycle design with the RFpeptides protocol. One for monomeric design, and one for binder design.

<br clear="all"/>

## MCTS planning
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/MCTS.svg" alt="Transformation" width="700" height="350"/>

Here is the illustration of MCTS planning.

<br clear="all"/>

## Working
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/workingModel2.svg" alt="Transformation" width="850" height="350"/>

Here is the illustration of MCTS planning.

<br clear="all"/>