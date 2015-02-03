#
# Examples Makefile.
#
#                  M. Kirsanov 07.04.2006
#                     Modified 18.11.2006
#                     26.03.2008 CLHEP dependency removed

SHELL = /bin/sh
CXX = g++
# CXX_COMMON=-ansi -pedantic -W -Wall -Wshadow -fPIC

# Set PYTHIA8, FASTJET, HEPMC yourself
PYTHIA8 = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/pythia8/185
# PYTHIA8 = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/pythia8/185
# FASTJET = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/fastjet/3.1.0
# FASTJET = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/fastjet/3.0.3-cms
FASTJET = /users/ra12451/fastjet-install
HEPMC=/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/hepmc/2.06.07-cms3
# HEPMC=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/hepmc/2.06.07-cms3
# CXX_COMMON = -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_COMMON:=-I$(PYTHIA8)/include $(CXX_COMMON) -Wl,-rpath $(PYTHIA8)/lib -ldl

 -include $(PYTHIA8)/config.mk
# -include config.mk

# Location of directories.
PYTHIA8_INC=$(PYTHIA8)/include
PYTHIA8_LIB=$(PYTHIA8)/lib
PYTHIA8_LIBARCH=$(PYTHIA8)/lib/archive
PYTHIA8_BIN=$(PYTHIA8)/bin

HEPMC_INC=$(HEPMC)/include
HEPMC_LIB=$(HEPMC)/lib

# There is no default behaviour, so remind user.
# Run using e.g. make bbtata
.cc:

	$(CXX) -I$(PYTHIA8_INC) $@.cc  \
	`$(FASTJET)/bin/fastjet-config --cxxflags  --plugins` \
	-I$(HEPMC_INC) \
	$(CXXFLAGSNO) -o $@.exe \
	-L$(HEPMC_LIB) -Wl,-rpath $(HEPMC_LIB) -lHepMC \
	-L$(PYTHIA8_LIB) -L$(PYTHIA8)/lib/archive -lpythia8 -llhapdfdummy \
	-lpythia8tohepmc  -Wl,-rpath $(HEPMC_LIB) -lHepMC \
	-L$(FASTJET)/lib \
	`$(FASTJET)/bin/fastjet-config --libs --plugins` \
	-pthread -std=c++11 -O2 $(CXX_COMMON)
# 	# -llhapdfdummy \

# main41:
	# For pythia 8.2X
	# $(CXX) $@.cc $(PYTHIA8)/lib/libpythia8.so -o $@ -I$(HEPMC_INC) $(CXX_COMMON) \
	# -L$(HEPMC_LIB) -Wl,-rpath $(HEPMC_LIB) -lHepMC

	# For pythia 8.1X
	# $(CXX) $@.cc -o $@ -I$(HEPMC_INC) $(CXX_COMMON) \
	# -L$(HEPMC_LIB) \
	# -L$(PYTHIA8)/lib/archive -lpythia8 -llhapdfdummy  \
	# -lpythia8tohepmc \
	# -Wl,-rpath $(HEPMC_LIB) -lHepMC

#	$(CXX) $(CXXFLAGS) -I$(PYTHIA8_INC) $@.cc -o $@.exe  \
#	-L$(PYTHIA8_LIBARCH) -lpythia8 -llhapdfdummy
#	mv $(PYTHIA8_BIN)/$@.exe $@.exe


# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	# rm -rf $(PYTHIA8_BIN)
	rm -rf *.exe
	rm -f *~;
	# rm -f \#*; rm -f core*
	rm main41