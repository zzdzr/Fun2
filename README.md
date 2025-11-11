# Fun2

A **reinforcement learning framework** for the dynamic tracing of **chromatin fountains/stripes**
<!-- Here is the workflow of Fun2 -->
> [!NOTE]
> - The project is **actively updating**, feedback and issues are welcome — I will continue to improve and update the repository regularly.
---
<div style="margin-bottom:200px; text-align:center;">
  <img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/workingModel4.svg" 
       width="1400" height="600" 
       style="border:none; display:block; margin:0 auto;">
</div>

> **Highlights**
> - Developed and validated a novel method for high-precision identification and tracking of **fountain / stripe-like chromatin structures**.
> - Integrated reinforcement learning with **Monte Carlo Tree Search (MCTS) in a continuous spatial action space**, combined with geometric transformation (**affine transformation**) to jointly optimize sampling and trajectory planning while enabling cross-sample or cross-condition comparability.
> - Quantitatively characterize the dynamic changes for the trajectory-like chromatin structures.
---

# Table of Contents
- [Fun2](#fun2)
- [Table of Contents](#table-of-contents)
- [Getting Started / Installation](#getting-started--installation)
  - [Quick start](#quick-start)
- [Usage](#usage)
  - [One-step usage](#one-step-usage)
  - [Configuration](#configuration)
  - [Sampling Box & Axes Configuration](#sampling-box--axes-configuration)
  - [Affine Transformation](#affine-transformation)
  - [Planning](#planning)
- [Understanding the Output Files](#understanding-the-output-files)
- [Repository Structure](#repository-structure)
- [License & Citation](#license--citation)

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

# Usage

## One-step usage
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/workingModel2.svg" alt="Transformation" width="800" height="350"/>

```bash
# 1) You just need to submit one line of command:
nohup fun2 --config config.yaml &
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

#### 1) Global axis and angle convention (contact matrix)
- The contact matrix uses a global coordinate system **(u, v)**.
- Orientation angles **θ** are defined **with respect to the matrix axes**:
  - **0°** → direction **parallel to the v-axis**
  - **90°** → direction **parallel to the u-axis**
- Unless stated otherwise, angles increase counter-clockwise in the *(u, v)* frame.

#### 2) Local geometry of the sampling box
- Each sampling box has its own local coordinates **(x, y)**:
  - **y-axis**: the **extension** (lengthwise) direction of the box.
  - **x-axis**: the direction **perpendicular** to the extension (across-box).
- Geometric parameters:
  - **Length (h)**: extent along the **y**-axis.
  - **Width (w)**: extent along the **x**-axis.
  - **Rotation (θ)**: box orientation expressed in the **global (u, v)** system.
- This dual-axis definition (global *(u, v)*, local *(x, y)*) permits explicit mapping between matrix coordinates and box geometry.

#### 3) Layered sampling and background estimation
- The box is partitioned into **K layers** along the **y**-axis (extension direction).
- At each layer *i*, the area is split into:
  - a **central signal region** (purple), and
  - two **edge bands** flanking it (the **edge width** is user-specified).
- The edge bands provide **local background** estimates, analogous to **edge-detection** strategies in image processing, thereby improving contrast between the trajectory signal and surrounding background.

#### 4) Notation recap
| Symbol | Meaning | Corresponding config parameter |
|:--|:--|:--|
| `(u, v)` | Global axes of the contact matrix | — |
| `(x, y)` | Local axes of sampling box (across vs. along extension) | — |
| `θ` | Orientation angle (0° ∥ v-axis; 90° ∥ u-axis) | `angle` |
| `w` | Box width along **x** | `width` |
| `h` | Box height (length) along **y** | `height` |
| `K` | Number of layers along **y** | `K = h / layer_height` |
| `edge_width` | Width of bands for estimating background signal | `edge_width` |

<br clear="all"/>



## Affine Transformation
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/affineTransform.svg" alt="Transformation" width="300" height="300" align="left"/>

In **Fun2**, an affine transformation extends the concept of a static sampling box into a **dynamic structure** capable of continuous-space transformations.  
This enables the sampling box to **trace trajectory-like patterns** in the contact matrix and quantitatively describe their properties.

Four primary forms of affine transformation are applied:

1. **Extension (`Δh`)**  
   - Adjusts the box length (*h*, along the y-axis).  
   - **−Δh** shortens the box; **+Δh** lengthens it, allowing coverage of varying segment lengths along the trajectory.

2. **Rotation (`Δθ`)**  
   - Rotates the box relative to the global *(u, v)* coordinate system.  
   - **−Δθ** rotates counter-clockwise; **+Δθ** rotates clockwise.  

3. **Expansion (`Δw`)**  
   - Modifies the box width (*w*, along the x-axis).  
   - **−Δw** narrows the box; **+Δw** widens it.

4. **Translocation (`Δx`)**  
   - Shifts the box laterally along the diagonal of contact matrix without changing its orientation or dimensions.  
   - **−Δx** shifts left; **+Δx** shifts right, enabling local repositioning while maintaining the same geometric parameters.

By combining these transformations, Fun2 optimizes box placement in a continuous spatial search space—forming the basis for trajectory detection and characterization.
#### 5) Notation recap
| Symbol / Parameter | Meaning (in this section) |
|:--|:--|
| `Δh` | Change in box length (*h*, along y-axis) |
| `Δθ` | Change in box rotation angle (global *(u, v)* frame) |
| `Δw` | Change in box width (*w*, along x-axis) |
| `Δx` | Lateral shift (translocation) along the diagonal of contact matrix |
| `action_limits` | Maximum absolute changes allowed for `[Δw, Δh, Δθ, Δx]` |
| `angle_boundary` | Angular constraints on `θ` (rotation) in degree |

> **Note:** `action_limits` directly controls how much the sampling box can change per transformation step, while `angle_boundary` defines the allowable orientation range in the global coordinate system.

<br clear="all"/>

## Planning
<p align="center">
  <img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/MCTS.svg"
       alt="MCTS planning" width="730">
</p>

Fun2 employs **Monte Carlo Tree Search (MCTS)** in a **continuous spatial search space** to optimize sampling box placement and trajectory planning.  
The MCTS process iteratively explores the search space to identify high-value geometric configurations, enabling efficient detection of trajectory-like features in contact matrices.

The planning procedure follows the standard four-stage MCTS loop, repeated for *n* iterations:

1. **Selection**  
   - Starting from the root node (*s₀*), the algorithm selects a child node using the **Upper Confidence Bound (UCB)** strategy.  
   - This balances exploration (trying new geometric configurations) and exploitation (refining known promising configurations).

2. **Expansion**  
   - A new node ($s_{t+2}$) is generated by applying a geometric transformation (affine modification) to the current sampling box state.  
   - **Progressive widening** controls the branching factor, allowing the search tree to adaptively grow in promising regions of the continuous space.

3. **Simulation**  
   - From the expanded node, a **rollout** is performed to simulate the resulting trajectory and evaluate its quality.  
   - Fun2 uses **Value Gradient Estimation** to assess the expected improvement and updates the policy (π) guiding subsequent actions.

4. **Backpropagation**  
   - The results of the simulation are propagated back through the selected path in the tree.  
   - Nodes along the trajectory update their statistics, improving the accuracy of future selection steps.

By iteratively refining box placement and orientation through MCTS, Fun2 is able to:
- Search **continuous transformation parameters** (length, width, rotation, translation).
- Adaptively optimize sampling geometry for complex trajectory patterns.
- Balance global exploration with local exploitation, ensuring both novelty and accuracy in detection.

This approach integrates **reinforcement learning principles** with **geometric search**, making the sampling process dynamic, data-driven, and optimal in high-dimensional spatial contexts.

#### 5) Notation recap
| Parameter | Meaning (in this section) |
|:--|:--|
| `max_iter` | Maximum number of MCTS iterations (*n* in the planning loop) |
| `exploration_constant` | Controls exploration vs. exploitation in the UCB selection formula (larger → more exploration) |
| `alpha` | Progressive widening coefficient, regulates branching growth in continuous space |
| `es_popsize` | Population size for Evolution Strategies (ES) used in simulation/optimization |
| `es_sigma` | Initial standard deviation of perturbations in ES search |
| `momentum` | Momentum term for policy updates during value gradient estimation |
| `max_rollout_depth` | Maximum depth for simulation rollouts from an expanded node |
| `eta` | Learning rate scaling factor for policy/value updates |
| `es_lr` | Learning rate for ES optimization of transformation parameters |
| `gamma` | Discount factor for cumulative reward in trajectory evaluation |
| `seed` | Random seed for reproducible MCTS and ES runs |

<br clear="all"/>


---
# Understanding the Output Files

After running Fun2, the output directory contains the following key files and subfolders:

| File / Folder | Description |
|:--|:--|
| `summits/` | Folder contains detected summits at the preprocessing step. |
| `results` | Tabulated metrics (e.g., fountain quality, width, angle, length, position) for each detected structure. |
| `oe_matrix/` | Serialized configuration snapshots (YAML/JSON) for reproducibility and further analysis. |
| `oe_cooler/` | Runtime logs including parameter settings, iteration summaries, and convergence diagnostics. |
### 📄 Column description of output file

| Column name         | Description | Unit / Range |
|---------------------|-------------|--------------|
| `chrom`             | Chromosome ID | — |
| `summit_start`      | Start position of the detected summit (genomic coordinate) | bp |
| `summit_end`        | End position of the detected summit (genomic coordinate) | bp |
| `identified_center` | Summit center position identified by Fun2 | bp |
| `height`            | Sampling box length/height along the y-axis | kb |
| `width`             | Sampling box width along the x-axis | kb |
| `angle`             | Sampling box rotation angle relative to the global (u, v) axes (0° ∥ v-axis; 90° ∥ u-axis) | degrees |
| `elongation_up`     | Genomic position for upstream end of elongation | bp |
| `elongation_down`   | Genomic position for downstream end of elongation  | bp |
| `width_up`          | Genomic position for upstream end of box' width | bp |
| `width_down`        | Genomic position for downstream end of box' width | bp |
| `reward`            | Reward score from MCTS planning | — (higher is better) |
| `intensity`         | Signal intensity for the sampling box | - |
| `quality`           | Quality score for identified structure | - |
| `pval_upstream`     | Raw P-value of signal significance for the sampling box relative to upstream background region | - |
| `pval_downstream`   | Raw P-value of signal significance for the sampling box relative to downstream background region | - |
| `FDR_upstream`      | Adjusted P-value for the sampling box relative to upstream background region | - |
| `FDR_downstream`    | Adjusted P-value for the sampling box relative to upstream background region | - |

> [!NOTE] 
> - Coordinates are in bp; bin size depends on the Hi-C resolution.  
> - `reward` is computed during the MCTS + ES optimization process.  
> - The definition of `angle` follows the convention in the [Sampling Box & Axes Configuration](#sampling-box--axes-configuration) section.

[⬆️ Back to top](#table-of-contents)

---
# Citation
```bibtex
@article{Zhangding2025Fun2,
  title   = {Characterizing trajectory-like chromatin architectures with Fun2},
  author  = {Zhengrong Zhangding, Xuhao Liu and Jiazhi Hu},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.02.25.640072},
  url     = {https://doi.org/10.1101/2025.02.25.640072}
}
```