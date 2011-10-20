###################################################################

all: inchworm chrysalis Meryl 

inchworm:
	cd Inchworm && ./configure --prefix=`pwd` && make install

chrysalis:
	cd Chrysalis && $(MAKE) UNSUPPORTED=yes

Meryl:
	cd trinity-plugins/kmer && $(MAKE) install

clean:
	cd Inchworm && make clean
	cd Chrysalis && $(MAKE) clean UNSUPPORTED=yes
	cd trinity-plugins/kmer && $(MAKE) clean
	find trinity-plugins/kmer -name "*C.d" -exec rm -f {} \;

###################################################################


