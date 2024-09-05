### Overlap of differentially expressed genes

maSigPro is used to find the significant genes at each time point over Day 0. This script generates an R datastructure (.rds) that stores the masigpro results. A tab seperated file with the significant genes is generated in the currene directory.

### Prerequisites


#### Packages

`maSigPro`


#### Data

(i) tab seperated normalied expression files
(ii) rds of the masigpro run, generated in the current directory for each cell type by the first script. This is then required by the next script. 

>`Rscript run_masigpro_CM.R`

>`Rscript run_masigpro_EC.R`

>`Rscript run_masigpro_FB.R`

>`Rscript run_masigpro_HC.R`

This is then used by `maSigpro_celltypes.R` to generate the upset plots that show the overlap of DE genes at each timepoint.

>`Rscript maSigpro_celltypes.R -c <CT>`

where CT is either CM, FB, HC or EC.

This generates upset plots as in Supplementary figure 3A-D.

-------------

### Overlap of differentially expressed miRNAs

maSigPro is used to find the significant miRNAs across all timepoints. The script also generates the heatmap and the upset plot as in Supplementary figure 3E-H.

### Prerequisites


#### Packages

`maSigPro`, `pheatmap`, `UpSetR` 


#### Data

(i) configuration lists out the column names that correspond to the sample names ffrom each cell type, after removal of the outliers. This is used for the maSigPro analysis.

(ii) miRNA normalized expression files in the folder `expression/<CT>/<CT>_mature_normalized_CPM.txt`


We then run the script for each cell type as:

>`Rscript miR_masigpro_CM.r`

>`Rscript miR_masigpro_remout_EC.r`

>`Rscript miR_masigpro_remout_FB.r`

>`Rscript miR_masigpro_remout_HC.r`


