# Fun2: Characterize fountain / stripe-like chromatin structures

<!-- Badges (可添加 Zenodo/License/Build Status) -->
<!-- ![DOI](https://zenodo.org/badge/641802007.svg) -->

**Fun2** is a **reinforcement learning framework** designed to dynamically track **replication fountain structures**.

<div style="margin-bottom:200px; text-align:center;">
  <img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/workingModel4.svg" 
       width="1400" height="600" 
       style="border:none; display:block; margin:0 auto;">
</div>

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
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/MCTS.svg" alt="Transformation" width="850" height="350"/>

Here is the illustration of MCTS planning.

<br clear="all"/>

## Working
<img src="https://github.com/zzdzr/Fun2/blob/main/docs/image/workingModel2.svg" alt="Transformation" width="850" height="350"/>

Here is the illustration of MCTS planning.

<br clear="all"/>

## 🚀 Installation
Clone the repository and set up the environment:

```bash
git clone https://github.com/zzdzr/Fun2.git
cd Fun2
conda env create -f environment.yml
conda activate fun2
pip install -e .
