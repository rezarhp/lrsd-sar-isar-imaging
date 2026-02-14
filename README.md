# lrsd-sar-isar-imaging
MATLAB code for fast and robust LRSD-based SAR/ISAR imaging and decomposition (arXiv:2512.10740).
# Fast and Robust LRSD-based SAR/ISAR Imaging and Decomposition (MATLAB)

This repository provides MATLAB implementations and reproducible demos for:

**Fast and Robust LRSD-based SAR/ISAR Imaging and Decomposition**  
arXiv:2512.10740  
Paper: https://arxiv.org/abs/2512.10740

The code includes two main implementations:
- **`main_conv_lrsd.m`**: a baseline / conventional LRSD update using iterative solvers (CGS).
- **`main_fast_LRSD.m`**: a faster LRSD variant using simplified closed-form updates.

---

## Repository layout

- `code/` : main scripts
- `code/utils/` : helper functions used by the scripts (patch operators, ROI extraction, sampling operators, colormap)
- `images/` : demo images used by the scripts
- `results/` : optional folder for outputs

---

## Requirements

- MATLAB R2020a+ (likely works on nearby versions)
- No special toolboxes are required beyond standard MATLAB functions used in the scripts.

---
## Citation

If you use this code, please cite:

H. R. Hashempour, M. Moradikia, H. Bastami, A. Abdelhadi and M. Soltanalian,  
"Fast and Robust LRSD-Based SAR/ISAR Imaging and Decomposition,"  
IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1–13, 2022, Art no. 5227413,  
doi: 10.1109/TGRS.2022.3172018.

cd <REPO>

