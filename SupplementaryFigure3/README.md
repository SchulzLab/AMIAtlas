### Overlap of differentially expressed genes

maSigPro is used to find the significant genes at each time point over Day 0. This script generates an R datastructure (.rds) that stores the masigpro results.

>`Rscript run_masigpro_CM.R`

>`Rscript run_masigpro_EC.R`

>`Rscript run_masigpro_FB.R`

>`Rscript run_masigpro_HC.R`

This is then used by `maSigpro_celltypes.R` to generate the upset plots that show the overlap of DE genes at each timepoint.

>`Rscript maSigpro_celltypes.R -c CM`

>`Rscript maSigpro_celltypes.R -c EC`

>`Rscript maSigpro_celltypes.R -c FB`

>`Rscript maSigpro_celltypes.R -c HC`

This generates plots as in Supplementary figure 3A-D.

### Overlap of differentially expressed miRNAs

maSigPro is used to find the significant miRNAs across each timepoint.

>`Rscript miR_masigpro_CM.r`

>`Rscript miR_masigpro_remout_EC.r`

>`Rscript miR_masigpro_remout_FB.r`

>`Rscript miR_masigpro_remout_HC.r`

This generates plots as in Supplementary figure 3E-H.
