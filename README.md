# Fluke-MF-FF-SED-Comparison
Data and Scripts for comparison of three diagnostic techniques for the detection and quantification of *Fasciola hepatica* and *Calicophoron daubneyi* parasite eggs. The analysis is presented in the manuscript by A. Bosco et al., "Comparison of Mini-FLOTAC, Flukefinder and Sedimentation techniques for detection and quantification of *Fasciola hepatica* and *Calicophoron daubneyi* eggs using spiked and naturally infected bovine faecal samples", 2023, under review

Description of Files:

-----
**SpikedInfection.csv** contains the data on spiked samples with a given dose (in eggs per gram, EPGs). Data is coded as follows:

Mini-FLOTAC - EPG: The eggs per gram obtained with Mini-FLOTAC.

Flukefinder - EPG: The eggs per gram obtained with Flukefinder.

Sedimentation - EPG: The eggs per gram obtained with Sedimentation.

Dose: Dosage of the spike (in EPGs).

Parasite: Parasite species, either *Fasciola hepatica* or *Calicophoron daubneyi*

-----
**FieldData.csv** contains the individual-level field data. Data is coded as follows:

flotac: The eggs per gram obtained with Mini-FLOTAC.

Parasite: Parasite species, either *Fasciola hepatica* or *Calicophoron daubneyi*.

Pool: The pool the sample will contribute to, numbered from 1 to 4.

Farm: The farm of origin, numbered from 1 to 10.

-----
**PooledSamples.csv** contains the pool-level field data. Data is coded as follows:

Flotac: The eggs per gram obtained with Mini-FLOTAC.

Parasite: Parasite species, either *Fasciola hepatica* or *Calicophoron daubneyi*.	

Pool: The pool the sample will contribute to, numbered from 1 to 4.	

Sedimentation: The eggs per gram obtained with Sedimentation.	

FlukeFinder: The eggs per gram obtained with Flukefinder.

Farm: The farm of origin, numbered from 1 to 10.

-----
**SpikedDataModel.R** contains the annotated code to run the latent class analysis model on the spiked data. This generates the expected heterogeneity (variation) with each method.

**NaturalInfectionModel.R** contains the annotated code to run the latent class analysis model on the field data. This uses both individual-level, as well as pool-level data.

**AnalysisSpikedDataModel.R** contains the annotated code to analyse the outcomes of the latent class analysis model on the spiked data. It includes the comparison of real vs simulated counts, accuracy (measured as absolute difference), expected number of negative counts and sensitivity as a function of infection intensity (in EPGs). The code for Figure 2 and Figure 3 is also included. 

**AnalysisNaturalInfectionModel.R** contains the annotated code to analyse the outcomes of the latent class analysis model on the field data. It includes the estimation of prevalence in the different farms and the code to generate Figure 4.
