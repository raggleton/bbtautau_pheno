#
# Examples Makefile.
#
#                  M. Kirsanov 07.04.2006
#                     Modified 18.11.2006
#                     26.03.2008 CLHEP dependency removed

SHELL = /bin/sh
CXX = g++

# Set PYTHIA8LOCATION, FASTJETLOCATION, HEPMC yourself
PYTHIA8LOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/pythia8/185
# PYTHIA8LOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/pythia8/185

# FASTJETLOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/fastjet/3.1.0
FASTJETLOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/fastjet/3.0.3-cms
# FASTJETLOCATION = /users/ra12451/fastjet-install

HEPMC = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/hepmc/2.06.07-cms3
# HEPMC=/cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/hepmc/2.06.07-cms3

# CXX_COMMON = -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_COMMON:=-I$(PYTHIA8LOCATION)/include $(CXX_COMMON) -Wl,-rpath $(PYTHIA8LOCATION)/lib -ldl

 -include $(PYTHIA8LOCATION)/config.mk
# -include config.mk

# Location of directories.
INCDIR=include
SRCDIR=src
LIBDIR=lib
LIBDIRARCH=lib/archive
BINDIR=bin

PYTHIA8LOCATION_INC=$(PYTHIA8LOCATION)/include
PYTHIA8LOCATION_LIB=$(PYTHIA8LOCATION)/lib
PYTHIA8LOCATION_LIBARCH=$(PYTHIA8LOCATION)/lib/archive
PYTHIA8LOCATION_BIN=$(PYTHIA8LOCATION)/bin

HEPMC_INC=$(HEPMC)/include
HEPMC_LIB=$(HEPMC)/lib

# There is no default behaviour, so remind user.
# Run using e.g. make bbtata
.cc:
	# $(CXX) -I$(PYTHIA8LOCATION)/$(INCDIR) $@.cc  \
	# `$(FASTJETLOCATION)/bin/fastjet-config --cxxflags  --plugins` \
	# -I$(HEPMC_INC) \
	# $(CXXFLAGSNO) -o $@.exe \
	# -L$(HEPMC_LIB) -Wl,-rpath $(HEPMC_LIB) -lHepMC \
	# -L$(PYTHIA8LOCATION_LIB) -L$(PYTHIA8LOCATION)/lib/archive -lpythia8 -llhapdfdummy \
	# -lpythia8tohepmc  -Wl,-rpath $(HEPMC_LIB) -lHepMC \
	# -L$(FASTJETLOCATION)/lib \
	# `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins` \
	# -pthread -std=c++11 -O2 $(CXX_COMMON)
# 	# -llhapdfdummy \

    # Note: $(CXXFLAGS) is after Fastjet flags as Fastjet includes
    #       optimisation/debug flags which may be unwanted (e.g. -g -O2)
	echo $(HEPMCLOCATION)
	$(CXX) -I$(PYTHIA8LOCATION)/$(INCDIR) $@.cc \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	$(CXXFLAGS) -Wno-shadow \
	-I$(PYTHIA8LOCATION)/$(INCDIR) -I$(HEPMC)/include \
	-o $@.exe \
	-L$(PYTHIA8LOCATION)/$(LIBDIR) -L$(PYTHIA8LOCATION)/$(LIBDIRARCH) -lpythia8 -llhapdfdummy \
	-lpythia8tohepmc \
	-L$(HEPMC)/lib -Wl,-rpath $(HEPMC_LIB) -lHepMC \
	-L$(FASTJETLOCATION)/lib \
	`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins` \
	-pthread -std=c++11
	# -L$(LHAPDFLOCATION) $(LHAPDFLIBNAME) \
	# -L$(LHAPDFLOCATION)/lib \


# main41:
	# For pythia 8.2X
	# $(CXX) $@.cc $(PYTHIA8LOCATION)/lib/libpythia8.so -o $@ -I$(HEPMC_INC) $(CXX_COMMON) \
	# -L$(HEPMC_LIB) -Wl,-rpath $(HEPMC_LIB) -lHepMC

	# For pythia 8.1X
	# $(CXX) $@.cc -o $@ -I$(HEPMC_INC) $(CXX_COMMON) \
	# -L$(HEPMC_LIB) \
	# -L$(PYTHIA8LOCATION)/lib/archive -lpythia8 -llhapdfdummy  \
	# -lpythia8tohepmc \
	# -Wl,-rpath $(HEPMC_LIB) -lHepMC

#	$(CXX) $(CXXFLAGS) -I$(PYTHIA8LOCATION_INC) $@.cc -o $@.exe  \
#	-L$(PYTHIA8LOCATION_LIBARCH) -lpythia8 -llhapdfdummy
#	mv $(PYTHIA8LOCATION_BIN)/$@.exe $@.exe


# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	# rm -rf $(PYTHIA8LOCATION_BIN)
	rm -rf *.exe
	rm -f *~;
	# rm -f \#*; rm -f core*
	rm main41