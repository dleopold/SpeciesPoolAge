## Repository notes

This repository contains code associated with a plant-fungal microcosm experiment showing how species pool age affects the realtionship between regional and local biodiversity. Working title of the associated manuscript: *Greater local diversity under older species pools may arise from enhanced competitive equivalence*

**Authors:** Devin R. Leopold & Tadashi Fukami

*******

* All code for bioinformatic processing and analyses is available in the `code/` folder and dependencies are specified in the Makefile in the project directory.  

* Raw Illumina data is archived in the NCBI SRA under [BioProject PRJNA613615](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA613615). In order to run the bioinformatic processing scripts, the raw data should be downloaded to a subfolder `data/MiSeq/` in the project directory.  

* This repository includes the output of the bioinformatic processing, and all associated meta data, which allows the analyses to be recreated without recreating the bioinformatics workflow.  

* The primary analyses associated with this project can be previewed here:  
[Diversity analyses](https://htmlpreview.github.io/?https://github.com/dleopold/SpeciesPoolAge/master/output/html/DiversityAnalyses.html)  
[Compositional analyses](https://htmlpreview.github.io/?https://github.com/dleopold/SpeciesPoolAge/master/output/html/CompositionalAnalyses.html)

* In addition to the main project code, additional R scripts used during project planning can be found in folder, `data/ProjectPlanning`. Code showing the process used to design the mock community control samples can be seen [here](https://github.com/dleopold/SpeciesPoolAge/blob/master/data/ProjectPlanning/MocComDesign.R). Details of the random sampling of isolates to create the experimental species pools can be found [here](https://htmlpreview.github.io/?https://github.com/dleopold/SpeciesPoolAge/master/data/ProjectPlanning/PoolDesign.html). 