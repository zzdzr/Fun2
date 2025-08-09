# Fun2 is tracking the untrackable

A **reinforcement learning framework** for quantitative analysis of **chromatin fountains/stripes**
<!-- Here is the workflow of Fun2 -->

<div style="margin-bottom:200px; text-align:center;">
  <img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/workingModel4.svg" 
       width="1400" height="600" 
       style="border:none; display:block; margin:0 auto;">
</div>

> **Highlights**
> - Developed and validated a novel method for high-precision identification and tracking of **fountain / stripe-like chromatin structures**.
> - Integrated reinforcement learning with **Monte Carlo Tree Search (MCTS) in a continuous spatial action space**, combined with geometric transformation (**affine transformation**) to jointly optimize sampling and trajectory planning while enabling cross-sample and cross-condition comparability in an automated analysis pipeline.
> - Produced reproducible outputs including visual figures, scoring tables, and configuration snapshots, enabling large-scale screening and hypothesis testing.
> - Quantitatively characterized dynamic changes between coupled and uncoupled replication fork states under various replication stress conditions.
---

# Table of Contents
- [Fun2](#fun2-is-tracking-the-untrackable)
  - [Description](#description)
- [Table of Contents](#table-of-contents)
- [Getting Started / Installation](#getting-started--installation)
  - [Quick start](#quick-start)
- [Usage](#usage)
  - [Quick usage](#quick-usage)
  - [Configuration](#configuration)
  - [Sampling Box & Axes Configuration](#sampling-box--axes-configuration)
  - [Affine Transformation](#affine-transformation)
  - [MCTS Planning](#mcts-planning)
- [Understanding the Output Files](#understanding-the-output-files)
- [Repository Structure](#repository-structure)
- [License & Citation](#license--citation)

---

## Description
Replication forks in mammalian cells often progress in a **coupled** manner. Under replication stress, **type-II fountains** exhibit distinct uncoupling degrees. **Fun2** quantitatively tracks these dynamics and relates uncoupling to higher-order chromatin structures (e.g., loops), replication timing shifts, and fork-restart defects. The toolkit is intended for large-scale screening, hypothesis testing, and figure-ready output.

[⬆️ Back to top](#table-of-contents)

---

# Getting Started / Installation

### 🚀 Quick start
```bash
# 1) Clone
git clone https://github.com/zzdzr/Fun2.git
cd Fun2

# 2) Create environment (recommended)
conda env create --name fun2 -f ./env/environment.yml
conda activate fun2

# 3) Install in development mode
pip install -e .
```
---

# 🖼️ Usage

## Quick Usage
Here is the one-step usage.
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/workingModel2.svg" alt="Transformation" width="800" height="350"/>
```bash
# 1) You just need to submit one line of command:
nohup fun2 config.yaml &
```
---

## Configuration

---

## Sampling Box & Axes Configuration
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/axis.svg" alt="SamplingBox" width="500" height="500" align="left"/>

We have recently published the RFpeptides protocol for using RFdiffusion to design macrocyclic peptides that bind target proteins with atomic accuracy (Rettie, Juergens, Adebomi, et al., 2025). In this section we briefly outline how to run this inference protocol. We have added two examples for running macrocycle design with the RFpeptides protocol. One for monomeric design, and one for binder design.

<br clear="all"/>

## Affine Transformation
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/affineTransform.svg" alt="Transformation" width="300" height="300" align="left"/>

We have recently published the RFpeptides protocol for using RFdiffusion to design macrocyclic peptides that bind target proteins with atomic accuracy (Rettie, Juergens, Adebomi, et al., 2025). In this section we briefly outline how to run this inference protocol. We have added two examples for running macrocycle design with the RFpeptides protocol. One for monomeric design, and one for binder design.

<br clear="all"/>

## MCTS planning
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/MCTS.svg" alt="Transformation" width="730" height="350"/>

Here is the illustration of MCTS planning.

<br clear="all"/>