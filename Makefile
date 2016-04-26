R=R
PACKAGE= growmod

VERSION := $(shell sed -n '/^Version: /s///p' growmod/DESCRIPTION)

TARBALL := $(PACKAGE)_$(VERSION).tar.gz
ZIPFILE := =$(PACKAGE)_$(VERSION).zip

CPP_SRC := $(PACKAGE)/src/*.cpp

all:
	make doc-update
	make build-package
	make install
	make pdf


build-package: $(TARBALL)
$(TARBALL): $(PACKAGE)/NAMESPACE $(CPP_SRC)
	$(R) CMD build --resave-data=no $(PACKAGE)

install: $(TARBALL)
	$(R) CMD INSTALL --preclean $<
	@touch $@

quick-install: $(PACKAGE)/src/growmod.so
	$(R) CMD INSTALL $(PACKAGE)

$(PACKAGE)/src/growmod.so: $(PACKAGE)/src/growmod.cpp
	cd $(PACKAGE)/src; echo "library(TMB); compile('growmod.cpp')" | R --slave
