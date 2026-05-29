.PHONY: clean init

ifeq ($(MAKECMDGOALS),init)
init:  
	$(MAKE) -C sources
else 

ifeq ($(MAKECMDGOALS),clean)
clean:
	$(MAKE) -C sources clean
	$(MAKE) -C main clean
	-find . -name '*.DS_Store' -print | xargs rm -rf
	-find . -name '*.dSYM' -print | xargs rm -rf
	-find . -name '._*' -print | xargs rm -rf
else

all:
	$(MAKE) -C main

endif

endif
