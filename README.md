## Repository notes

This repository contains code associated with a plant-fungal microcosm experiment showing how species pool age affects the realtionship between regional and local biodiversity: *Greater local diversity under older species pools may arise from enhanced competitive equivalence*, ***Ecology Letters*** (accepted).

**Authors:** Devin R. Leopold & Tadashi Fukami

**Abstract:** Ecological communities typically contain more species when located within geologically older regions. This pattern is traditionally attributed to the long-term accumulation of species in the regional species pool, with local species interactions playing a minor role. We provide evidence suggesting a more important role of local species interactions than generally assumed. We assembled 320 communities of root-associated fungi under 80 species pools, varying species-pool richness and the mean age of the sites from which the fungi were collected across a 4-myr soil chronosequence. We found that local diversity increased more with increasing species-pool richness when species were from older sites. We also found that older species pools had lower functional and phylogenetic diversity, indicating greater competitive equivalence among species. Our results suggest that older regions have higher local richness not simply because older pools are more speciose, but also because species have evolved traits that allow them to locally co-occur.

*******

* All code for bioinformatic processing and analyses is available in the `code/` folder and dependencies are specified in the Makefile in the project directory.  

* Raw Illumina data is archived in the NCBI SRA under [BioProject PRJNA613615](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA613615). In order to run the bioinformatic processing scripts, the raw data should be downloaded to a subfolder `data/MiSeq/` in the project directory.  

* A preview of the bioinformatic processing of the raw Illumina data can be viewed [here](https://htmlpreview.github.io/?https://github.com/dleopold/SpeciesPoolAge/master/output/html/MiSeqProcessing.html).

* This repository includes the output of the bioinformatic processing, and all associated metadata, which allows the analyses to be recreated without recreating the bioinformatics workflow.  

* The primary analyses associated with this project can be previewed here:  
[Diversity analyses](https://htmlpreview.github.io/?https://github.com/dleopold/SpeciesPoolAge/master/output/html/DiversityAnalyses.html)  
[Compositional analyses](https://htmlpreview.github.io/?https://github.com/dleopold/SpeciesPoolAge/master/output/html/CompositionalAnalyses.html)

* In addition to the main project code, additional R scripts used during project planning can be found in folder, `data/ProjectPlanning`. Code showing the process used to design the mock community control samples can be seen [here](https://github.com/dleopold/SpeciesPoolAge/blob/master/data/ProjectPlanning/MocComDesign.R). Details of the random sampling of isolates to create the experimental species pools can be found [here](https://htmlpreview.github.io/?https://github.com/dleopold/SpeciesPoolAge/master/data/ProjectPlanning/PoolDesign.html). 