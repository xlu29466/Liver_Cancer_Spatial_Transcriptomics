# Liver_Cancer_Spatial_Transcriptomics
We used nHDP model to learn gene expression modules (GEMs) from scRNA data, and then use the learned model to infer GEM expression level on spatial transcriptomics data. Then we tried to establish GEM-Ligand-Recptor-GEM causual relationships.

## Requirements (lastest version):
scanpy  
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
## Step 4:
You are ready to explore the notebooks.
## Notebook Descriptions:
1. read_visium.ipynb: Basic analysis of the raw data. Batch effects removal was also explored. No nHDP results involved.
2. Prepare_for_inference.ipynb: prepare raw spatial data for nHDP inference. Can be ignored.
3. 1-4N, 1-4L, 1-4T.ipynb: Visualization of nHDP inference results on each slide.
4. Patient1-4_analysis.ipynb: Each patient's T, N, L tissues were analysed together. We used genes or GEMs to do clustering. We also developed the method to select out boundary between 2 clusters.
5. All_Patients_Analysis.ipynb: All patients' data were put together to do comprehensive GEM similarity analysis.
6. Ligand-receptor Analysis.ipynb: We screened ligand-receptor pairs that were measured in the raw data. We also developed "Ligand Diffusion" method to calculate ligand-receptor products on each slide.
7. GEM_correlation_analysis.ipynb: Since there are batch effects among slides, thus we treat each slide individually here. For each slide, we used cosine similarity to get top 30 most correlated GEM-GEM pairs. GEM pairs of the same cell type were excluded. If one GEM-GEM pair showed up in more than one slide's top 30 candidates, then this pair is called a "Pattern".
8. Align_GEM_with_LR_Version 1 and 2.ipynb: For each GEM-GEM pair on each slide, we used there average to represent the expression level of this pair. Then we apply cosine similarity again to select out top 10 most correlated ligand-receptor pairs. That is, these ligand-recptor product (one vector) matches well with GEM-GEM pair average (one vector). Then we applied detailed scatter plot and regression analysis on each of them.
9.  Global_GEM_LR_conditional_possibility.ipynb: Our model is that GEM1 -> L*R -> GEM2, the activation of GEM1 will increase the secretion of ligands in one cell, then ligands will activate more receptors on other cells' membrane, and GEM2 is the down-stream effects of the activation of the receptors. Firstly, we calculated the correlation-coefficient of GEM1 and GEM2, then we used log/lowess regression to explain the conribution of ligand-receptor to GEM2. We remove the portion of GEM2 that can be explained by ligand-receptor, the remaining was called residue. Then we calculated the correlation-coefficient of GEM1 and GEM2 residue. The expected results is that the correlation-coefficient drop to near 0. That means conditional on LR, the GEM1-GEM2 correlation was cut off.

