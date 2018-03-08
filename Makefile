R?=R
PACKAGE=otoclass

all: doc build check test install

doc:
	@echo "\033[0;32mUpdating documentation\033[0;0m"
	rm -f ${PACKAGE}/src/*.so
	$(R) -q -e 'roxygen2::roxygenize("${PACKAGE}")'

build: doc
	@echo "\033[0;32mBuilding package\033[0;0m"
	$(R) CMD build ${PACKAGE}

check: doc build
	@echo "\033[0;32mChecking package as cran\033[0;0m"
	$(eval VERSION = $(shell Rscript -e "l<-readLines(\"${PACKAGE}/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))"))
	$(R) CMD check --as-cran ${PACKAGE}_${VERSION}.tar.gz

install: doc build
	@echo "\033[0;32mInstalling package\033[0;0m"
	$(eval VERSION = $(shell Rscript -e "l<-readLines(\"${PACKAGE}/DESCRIPTION\");cat(gsub(\"Version: \",\"\",l[grepl(\"Version: \",l)]))"))
	@echo "Version: ${VERSION}"
	$(R) CMD INSTALL ${PACKAGE}_${VERSION}.tar.gz

install_dependencies:
	@echo "\033[0;32mInstalling package dependencies\033[0;0m"
	@echo "source('https://raw.githubusercontent.com/calbertsen/caMisc/master/R/build_from_github.R'); \
	installDependencies('${PACKAGE}/DESCRIPTION',dependencies=c(\"Depends\",\"LinkingTo\",\"Imports\",\"Suggests\",\"Enhances\"))" | $(R) -q --slave

test:
	@echo "\033[0;31mNothing yet\033[0;0m"

clean:
	@echo "\033[0;32mCleaning directory\033[0;0m"
	git clean -f -d

uninstall:
	@echo "\033[0;32mUninstalling package\033[0;0m"
	$(R) CMD REMOVE ${PACKAGE}
