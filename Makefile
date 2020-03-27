all: miseq biolog
	
clean:
	rm -r output
	$(info you have a clean slate)

housekeeping:
	mkdir -p output/html
	mkdir -p output/csv
	mkdir -p output/rds
	mkdir -p output/figs

miseq: output/html/MiSeqProcessing.html
output/html/MiSeqProcessing.html: code/MiSeqProcessing.R | housekeeping
	Rscript -e 'rmarkdown::render("$<", output_dir = "output/html", knit_root_dir = "$(CURDIR)")'

biolog: output/csv/biolog.csv
output/csv/biolog.csv: code/biolog.R | housekeeping
	Rscript $< 

