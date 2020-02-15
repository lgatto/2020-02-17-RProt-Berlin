
%.md: %.Rmd
	/opt/Rpatched/lib/R/bin/Rscript -e 'require("knitr"); knit("$^")'

