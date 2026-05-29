
.PHONY: all clean

SOURCE=theory.sinp.msu.ru/~pukhov/SPECTRA
data = AtProduction_antideuterons.dat  AtProduction_antiprotons.dat  AtProduction_gammas.dat  AtProduction_neutrinos_e.dat  AtProduction_neutrinos_mu.dat  AtProduction_neutrinos_tau.dat  AtProduction_positrons.dat

ifeq ($(MAKECMDGOALS),clean)
clean: 
	rm -f $(data)
endif

all: 
ifneq (, $(shell command -v curl 2> /dev/null))
	curl https://$(SOURCE)/PPPC.tgz  -O -L
else 
ifneq (, $(shell command -v wget 2> /dev/null))
	wget https://$(SOURCE)/PPPC.tgz
else 
	$(error "Neither wget nor curl are available, please install wget or curl or change SSUSY.make accordingly.")
endif
endif
	tar -xvzf PPPC.tgz
	rm PPPC.tgz
