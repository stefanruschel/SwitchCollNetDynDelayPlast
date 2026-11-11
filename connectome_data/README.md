# Connectome data

Structural connectome of the 68 cortical regions of the Desikan–Killiany atlas, obtained from Human Connectome
Project diffusion data [1] and downsampled onto the atlas as in Ref. [2].

| File | MATLAB variable | Content |
|---|---|---|
| `sc_raw.mat` | `sc` | structural coupling weights w_ij (68 × 68, symmetric); only the first slice `sc(:,:,1)` is used |
| `D.mat` | `dist` | Euclidean distances d_ij between the regions (68 × 68, symmetric, maximum ≈ 0.3); only `dist(:,:,1)` is used |

The model code (`src/connectome_model.jl`, function `load_connectome`) rescales the weights to an average row sum
of 1 and the distances to mean 1, as stated in the Letter.

## Terms of use

These matrices are derived from WU-Minn HCP Open Access data and are redistributed under the
[WU-Minn HCP Consortium Open Access Data Use Terms](https://www.humanconnectome.org/study/hcp-young-adult/document/wu-minn-hcp-consortium-open-access-data-use-terms),
as those terms permit ("I may redistribute original WU-Minn HCP Open Access data and any derived data as long as the
data are redistributed under these same Data Use Terms"). Anyone using them agrees to those terms, in particular not to
attempt to identify or contact the subjects, to follow the rules of their institution, and to acknowledge the HCP as
below. The licence of the code in this repository does not apply to the files in this folder.

## Acknowledgement

Data were provided in part by the Human Connectome Project, WU-Minn Consortium (Principal Investigators: David Van
Essen and Kamil Ugurbil; 1U54MH091657) funded by the 16 NIH Institutes and Centers that support the NIH Blueprint
for Neuroscience Research; and by the McDonnell Center for Systems Neuroscience at Washington University.

## References

1. D. C. Van Essen *et al.* and WU-Minn HCP Consortium, NeuroImage **80**, 62 (2013).
2. M. Forrester, S. Petros, O. Cattell, Y. M. Lai, R. D. O'Dea, S. Sotiropoulos, and S. Coombes,
   PLoS Comput. Biol. **20**, e1012647 (2024).
