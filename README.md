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
> - Integrated reinforcement learning with **Monte Carlo Tree Search (MCTS) in a continuous spatial action space**, combined with geometric transformation (**affine transformation**) to jointly optimize sampling and trajectory planning while enabling cross-sample and cross-condition comparability.
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
  - [One-step usage](#one-step-usage)
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

## One-step usage
Here is the one-step usage.
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/workingModel2.svg" alt="Transformation" width="800" height="350"/>
```bash
# 1) You just need to submit one line of command:
nohup fun2 config.yaml &
```

## Configuration
```yaml
# Example configuration for Fun2 (fountain detection)

## io:
  # Paths for input and output
  juicer_jar: "juicer_tools.1.9.9_jcuda.0.8.jar"
  chrom_size: "./config/genome/hg19.chrom.size"
  full_name: "Repli-HiC_test"
  working_dir: "./data/fountains/NT"
  input_hic: "./data/NT_merged.rm1k.downsample.hic"

## preprocess:
  # Script to generate cooler from Hi-C
  generate_cooler_script: "./preprocess/generate_cooler.py"
  oe_matrix_dir: "./data/fountains/NT/oe_matrix"
  cooler_output_dir: "./data/fountains/NT/oe_cooler"

  # Summit detection
  summits_output_dir: "./data/fountains/NT/summits"
  summits_threshold: 0.15
  padding: 25000    # bp

## son:
  # SoN calculation parameters
  calSoN_script: "calculate_son.py"
  SoN_params:
    - height: 100
      width: 25
      angle: 45
    - height: 200
      width: 25
      angle: 45
    - height: 300
      width: 25
      angle: 45
  resolution: 5000
  oe_normalization: "VC_SQRT"

## planner:
  # Fountain identification settings (MCTS + ES in continuous space)
  ini_height: 100
  ini_width:  25
  ini_angle:  45
  layer_height: 1
  edge_width: 8
  max_iter: 500
  exploration_constant: 10
  alpha: 0.5
  es_popsize: 8
  es_sigma: 5.0
  momentum: 0
  max_rollout_depth: 1
  eta: 1
  es_lr: 2.0
  gamma: 0.9
  action_limits: [[-2, 2], [-50, 50], [-2, 2], [0, 0]]
  angle_boundary: [40, 50]
  fountain_path: "./data/fountains/NT/Repli-HiC_test"
  n_process: 10
  seed: 1
  mode: "default"
```

## Sampling Box & Axes Configuration
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/axis.svg" alt="SamplingBox" width="500" height="500" align="left"/>

  In the contact matrix, the global coordinate system is defined by the *u*-axis and *v*-axis.  
An orientation of **0°** corresponds to a direction parallel to the *v*-axis, whereas **90°** is parallel to the *u*-axis.  
Sampling box orientations are specified with respect to this global axis system.  
  Each sampling box also possesses its own local geometric coordinates, denoted as the *x*-axis and *y*-axis, which define its intrinsic **length**, **width**, and **rotation angle**.  
This dual-axis definition allows explicit mapping between global matrix coordinates and the sampling box geometry.
  A complete sampling box is divided into multiple **layers** along its lengthwise extension.  
During analysis, layers sequentially capture signals as the box extends.  
  Each layer is further subdivided into a **central signal region** (purple) and two **edge regions** flanking it.  
The edge regions are used to estimate local background signal levels, analogous to edge detection techniques in image processing, thereby enhancing contrast between the feature of interest and its surrounding background.
<br clear="all"/>

## Affine Transformation
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/affineTransform.svg" alt="Transformation" width="300" height="300" align="left"/>

We have recently published the RFpeptides protocol for using RFdiffusion to design macrocyclic peptides that bind target proteins with atomic accuracy (Rettie, Juergens, Adebomi, et al., 2025). In this section we briefly outline how to run this inference protocol. We have added two examples for running macrocycle design with the RFpeptides protocol. One for monomeric design, and one for binder design.

<br clear="all"/>

## MCTS planning
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/MCTS.svg" alt="Transformation" width="730" height="350"/>

Here is the illustration of MCTS planning.

<br clear="all"/>