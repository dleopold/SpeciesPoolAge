all: analyses
	
clean:
	rm -r output
	$(info you have a clean slate)

### Create output directories ###
housekeeping:
	mkdir -p output/MiSeq
	mkdir -p output/html
	mkdir -p output/csv
	mkdir -p output/rds
	mkdir -p output/figs

### Process raw MiSeq data ###
miseq: output/html/MiSeqProcessing.html
output/html/MiSeqProcessing.html output/rds/phy.rds output/csv/pools.dat.csv: code/MiSeqProcessing.R | housekeeping
	Rscript -e 'rmarkdown::render("$<", output_dir = "output/html", knit_root_dir = "$(CURDIR)")'

### Process raw phenotypic microarray data ###
biolog: output/csv/biolog.csv
output/csv/biolog.csv: code/biolog.R | housekeeping
	Rscript $< 

### Analyses ###
analyses: diversity composition supplemental

diversity: output/html/DiversityAnalyses.html
output/html/DiversityAnalyses.html: code/DiversityAnalyses.R output/rds/phy.rds | housekeeping
	Rscript -e 'rmarkdown::render("$<", output_dir = "output/html", knit_root_dir = "$(CURDIR)")'
	
composition: output/html/CompositionalAnalyses.html
output/html/CompositionalAnalyses.html: code/CompositionalAnalyses.R output/rds/phy.rds output/csv/biolog.csv | housekeeping
	Rscript -e 'rmarkdown::render("$<", output_dir = "output/html", knit_root_dir = "$(CURDIR)")'
	
supplemental: output/figs/FigS1.jpg output/figs/FigS2.jpg output/figs/FigS3.jpg
output/figs/FigS1.jpg: code/phylotree.R | housekeeping
	Rscript $<
output/figs/FigS3.jpg: code/RelativeAbundanceOfTaxa.R | housekeeping
	Rscript $<
output/figs/FigS2.jpg: code/plantData.R output/csv/pools.dat.csv | housekeeping
	Rscript $<
	
