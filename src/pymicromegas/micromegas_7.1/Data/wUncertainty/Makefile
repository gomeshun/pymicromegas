
.PHONY: all clean

SOURCE=theory.sinp.msu.ru/~pukhov/SPECTRA
data = AtProduction-AntiP.dat  AtProduction-Ga.dat  AtProduction-Positrons.dat 

ifeq ($(MAKECMDGOALS),clean)
clean: 
	rm -f $(data)
endif

all: 
ifneq (, $(shell command -v curl 2> /dev/null))
	curl https://$(SOURCE)/Pythia8_U.tgz  -O -L
else 
ifneq (, $(shell command -v wget 2> /dev/null))
	wget https://$(SOURCE)/Pythia8_U.tgz
else 
	$(error "Neither wget nor curl are available, please install wget or curl or change SSUSY.make accordingly.")
endif
endif
	tar -xvzf Pythia8_U.tgz
	rm Pythia8_U.tgz
