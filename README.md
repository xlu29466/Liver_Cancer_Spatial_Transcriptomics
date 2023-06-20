# Liver_Cancer_Spatial_Transcriptomics
We used nHDP model to learn gene expression modules (GEMs) from scRNA data, and then use the learned model to infer GEM expression level on spatial transcriptomics data. Then we tried to establish GEM-Ligand-Recptor-GEM causual relationships.

## Requirements (all lastest version):
scanpy \n
numpy
pandas
matplotlib
anndata
scipy
cv2
seaborn
sklearn
harmony-pytorch
pingouin
tqdm




## Step 1:
Download the proejct folder to local.
## Step 2:
Download raw data (HCC 1-4N, HCC 1-4L, HCC 1-4T, 12 folders in total) from "DataAccessWeb"(in this repo), put the downloaded files into a folder called "raw_data". Put "raw_data" folder in the project folder. 
## Step 3:
Download the nHDP inference folder "nHDP" from Microsoft Teams SpatialSep - Files - Spatial_Eric - "nHDP" folder. Put the folder in the project folder. 
