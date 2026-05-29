
.PHONY: all clean

SOURCE=theory.sinp.msu.ru/~pukhov/SPECTRA
data = AtProduction-AntiP.dat  AtProduction-Gamma.dat  AtProduction-Nuel.dat  AtProduction-Numu.dat  AtProduction-Nuta.dat  AtProduction-Positrons.dat 

ifeq ($(MAKECMDGOALS),clean)
clean: 
	rm -f $(data)
endif

all: 
ifneq (, $(shell command -v curl 2> /dev/null))
	curl https://$(SOURCE)/CosmiXs.tgz  -O -L
else 
ifneq (, $(shell command -v wget 2> /dev/null))
	wget https://$(SOURCE)/CosmiXs.tgz
else 
	$(error "Neither wget nor curl are available, please install wget or curl or change SSUSY.make accordingly.")
endif
endif
	tar -xvzf CosmiXs.tgz
	rm CosmiXs.tgz
