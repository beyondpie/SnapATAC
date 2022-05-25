doc:
	Rscript -e "devtools::document(pkg = '.')"

install:
	Rscript -e 'remotes::install_github("beyondpie/SnapATAC", ref =
	"szu", upgrade = "never")'
