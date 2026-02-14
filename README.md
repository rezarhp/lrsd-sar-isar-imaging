# LRSD SAR/ISAR Imaging (MATLAB)

MATLAB implementation of:

**Fast and Robust LRSD-Based SAR/ISAR Imaging and Decomposition**

📄 arXiv: https://arxiv.org/abs/2512.10740  

This repository provides reproducible implementations of LRSD-based SAR/ISAR imaging algorithms, including both conventional iterative and fast variants.

---

## 📦 Contents

The repository includes two main implementations:

- **`main_conv_lrsd.m`**  
  Conventional LRSD algorithm using iterative CGS updates.

- **`main_fast_lrsd.m`**  
  Fast LRSD algorithm using simplified closed-form updates.

### Supporting functions
- `Extract_ROI.m`
- `Patch_Operator.m`
- `Inverse_Patch_Operator.m`
- `USFO.m`
- `sar_cmap.m`

### Data
Located in the `images/` folder:
- `ku6.png`
- `sub_image_from_SARPER.png`

---

## ⚙️ Requirements

- MATLAB R2020a or newer
- No additional toolboxes required

---

# Citation

If you use this code or repository in your research, please cite:

H. R. Hashempour, M. Moradikia, H. Bastami, A. Abdelhadi and M. Soltanalian,
"Fast and Robust LRSD-Based SAR/ISAR Imaging and Decomposition,"
IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1–13, 2022, Art no. 5227413.
https://doi.org/10.1109/TGRS.2022.3172018
git clone https://github.com/<YOUR_USERNAME>/lrsd-sar-isar-imaging.git
cd lrsd-sar-isar-imaging
