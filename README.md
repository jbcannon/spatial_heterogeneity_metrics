# spatial_heterogeneity_metrics
 R and python functions and scripts to calculate spatial heterogeneity metrics from Cannon et al. 2020
 
 Goal of these scripts are to take raster canopy cover and study boundary polygon for actual or simulation restoration treatments. The generate_gaps.R script will create polygons of gaps that meet certain canopy cover critera. This utilizes a Python script that requires installation of ArcMap with the Spatial Analyst extension activated.
 
 The spatial_heterogeneity_metrics.R script will take canopy cover raster, boundary, and gap file from generate_gaps.R script and quantify several metrics including those reported in Cannon et al. 2020. The outputs include distribution of canopy cover among various cover classes, gap size and variability charcteristics, landscape complexity (H), and relative contagion (RC2).
