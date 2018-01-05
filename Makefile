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

test:
	@echo "\033[0;31mNothing yet\033[0;0m"

clean:
	@echo "\033[0;32mCleaning directory\033[0;0m"
	git clean -f -d

uninstall:
	@echo "\033[0;32mUninstalling package\033[0;0m"
	$(R) CMD REMOVE ${PACKAGE}
