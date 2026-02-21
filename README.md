# Graphene-lattice-dynamics

[English](#english-version) | [中文](#中文版)

## English Version

Phonon dispersion, DOS, and specific-heat analysis for graphene and graphene nanoribbons (AGNR/ZGNR) using force-constant lattice dynamics models.

### Scripts

- `Graphene_5NNFC_2D.m`: graphene dispersion along `Gamma-M-K-Gamma`.
- `Graphene_5NNFC_3D.m`: 2D Brillouin-zone phonon surfaces; saves `data.mat` (`X`, `Y`, `W`).
- `Graphene_DOS.m`: graphene DOS + heat capacity; loads `data.mat`, saves `DOS.mat` (`frequency`, `pdf`).
- `GNR_A_DR.m`: AGNR dispersion.
- `GNR_A_DOS.m`: AGNR DOS and comparison with graphene; saves `DOS_AGNR_8_16_32.mat`.
- `GNR_Z_DR.m`: ZGNR dispersion.
- `GNR_Z_DOS.m`: ZGNR DOS and comparison with graphene; saves `DOS_ZGNR-8-16-32.mat`.

## 中文版

基于力常数晶格动力学模型，计算石墨烯及石墨烯纳米带（AGNR/ZGNR）的声子色散、DOS 与比热。

### 脚本说明

- `Graphene_5NNFC_2D.m`：石墨烯 `Gamma-M-K-Gamma` 路径色散。
- `Graphene_5NNFC_3D.m`：二维布里渊区声子频率面；保存 `data.mat`（`X`、`Y`、`W`）。
- `Graphene_DOS.m`：石墨烯 DOS 与热容；读取 `data.mat`，保存 `DOS.mat`（`frequency`、`pdf`）。
- `GNR_A_DR.m`：AGNR 色散。
- `GNR_A_DOS.m`：AGNR DOS 及与石墨烯对比；保存 `DOS_AGNR_8_16_32.mat`。
- `GNR_Z_DR.m`：ZGNR 色散。
- `GNR_Z_DOS.m`：ZGNR DOS 及与石墨烯对比；保存 `DOS_ZGNR-8-16-32.mat`。
