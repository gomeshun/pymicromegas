
.PHONY: all clean flags


all:include/microPath.h
	$(MAKE) -C CalcHEP_src  MICROMEGAS=MICROMEGAS
	$(MAKE) -C sources

include/microPath.h:
	printf '%s\n' '#include <stdlib.h>' '#define micrO (getenv("PYMICROMEGAS_MICROPATH") ? getenv("PYMICROMEGAS_MICROPATH") : "$(CURDIR)")' > include/microPath.h
        
   
clean:  
	rm -f include/microPath.h
	./clean
	rm -rf */*.dSYM 

flags: 
	$(MAKE) -C CalcHEP_src flags