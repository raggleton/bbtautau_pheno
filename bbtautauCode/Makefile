#
# Examples Makefile.
#
#                  M. Kirsanov 07.04.2006
#                     Modified 18.11.2006
#                     26.03.2008 CLHEP dependency removed

SHELL = /bin/sh
CXX = g++

# Location of directories.
INCDIR=include
SRCDIR=src
LIBDIR=lib
LIBDIRARCH=lib/archive
BINDIR=bin

# Set PYTHIA8LOCATION, FASTJETLOCATION, HEPMC yourself
# PYTHIA8LOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/pythia8/185
PYTHIA8LOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/pythia8/201
PYTHIA8LOCATION_INC=$(PYTHIA8LOCATION)/include
PYTHIA8LOCATION_LIB=$(PYTHIA8LOCATION)/lib
PYTHIA8LOCATION_LIBARCH=$(PYTHIA8LOCATION)/lib/archive
PYTHIA8LOCATION_BIN=$(PYTHIA8LOCATION)/bin

# FASTJETLOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/fastjet/3.1.0
# FASTJETLOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/fastjet/3.0.3-cms
FASTJETLOCATION = /users/ra12451/fastjet-install
FASTJETLOCATION_LIB = $(FASTJETLOCATION)/lib

# Don't call it HEPMCLOCATION beacuse it defaults to /build/degano not cvmsfs for some stupid reason!!
# HEPMC = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/hepmc/2.06.07-cms3
HEPMC = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/hepmc/2.06.07-cms3
HEPMCLOCATION_INC=$(HEPMC)/include
HEPMCLOCATION_LIB=$(HEPMC)/lib

# For boost libraries
BOOSTLOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/boost/1.57.0/
BOOSTLOCATION_INC = -I$(BOOSTLOCATION)/include
BOOSTLOCATION_LIB = -L$(BOOSTLOCATION)/lib -lboost_program_options -Wl,-rpath $(BOOSTLOCATION)/lib

# CXX_COMMON = -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_COMMON:=-I$(PYTHIA8LOCATION)/include $(CXX_COMMON) -Wl,-rpath $(PYTHIA8LOCATION)/lib -ldl

 -include $(PYTHIA8LOCATION)/config.mk
# -include config.mk



# There is no default behaviour, so remind user.
# Run using e.g. make bbtata
.cc:
	# $(CXX) -I$(PYTHIA8LOCATION_INC) $@.cc  \
	# `$(FASTJETLOCATION)/bin/fastjet-config --cxxflags  --plugins` \
	# -I$(HEPMCLOCATION_INC) \
	# $(CXXFLAGSNO) -o $@.exe \
	# -L$(HEPMCLOCATION_LIB) -Wl,-rpath $(HEPMCLOCATION_LIB) -lHepMC \
	# -L$(PYTHIA8LOCATION_LIB) -L$(PYTHIA8LOCATION_LIBARCH) -lpythia8 -llhapdfdummy \
	# -lpythia8tohepmc  -Wl,-rpath $(HEPMCLOCATION_LIB) -lHepMC \
	# -L$(FASTJETLOCATION_LIB) \
	# `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins` \
	# -pthread -std=c++11 -O2 $(CXX_COMMON)
# 	# -llhapdfdummy \

	# For Pythia 8.1X
    # Note: $(CXXFLAGS) is after Fastjet flags as Fastjet includes
    #       optimisation/debug flags which may be unwanted (e.g. -g -O2)
	# $(CXX) -I$(PYTHIA8LOCATION_INC) $@.cc \
	# `$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	# $(CXXFLAGS) -Wno-shadow \
	# -I$(PYTHIA8LOCATION_INC) -I$(HEPMCLOCATION_INC) \
	# -o $@.exe \
	# -L$(PYTHIA8LOCATION)/lib -L$(PYTHIA8LOCATION_LIBARCH) -lpythia8 -llhapdfdummy \
	# -lpythia8tohepmc \
	# -L$(HEPMCLOCATION_LIB) -Wl,-rpath $(HEPMCLOCATION_LIB) -lHepMC \
	# -L$(FASTJETLOCATION_LIB) \
	# `$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins` \
	# -pthread -std=c++11
	# -L$(LHAPDFLOCATION) $(LHAPDFLIBNAME) \
	# -L$(LHAPDFLOCATION)/lib \

	# For pythia 8.2X
	# make sure you do cmsenv in CMSSW_7_3_0/src
	$(CXX) -I$(PYTHIA8LOCATION_INC) $@.cc \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	$(CXXFLAGS) -Wno-shadow \
	-I$(PYTHIA8LOCATION)/include -I$(HEPMCLOCATION_INC) \
	$(BOOSTLOCATION_INC) \
	-o $@.exe \
	-L$(PYTHIA8LOCATION)/lib -L$(PYTHIA8LOCATION_LIBARCH) -lpythia8 \
	-L$(HEPMCLOCATION_LIB) -Wl,-rpath $(HEPMCLOCATION_LIB) -lHepMC \
	-L$(FASTJETLOCATION_LIB) \
	`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins` \
	$(BOOSTLOCATION_LIB) \
	-pthread -std=c++11 -O2

main41:
	# For pythia 8.2X
	# $(CXX) $@.cc $(PYTHIA8LOCATION)/lib/libpythia8.so -o $@ -I$(HEPMCLOCATION_INC) $(CXX_COMMON) \
	# -L$(HEPMCLOCATION_LIB) -Wl,-rpath $(HEPMCLOCATION_LIB) -lHepMC

	# For pythia 8.1X
	$(CXX) $@.cc -o $@ -I$(HEPMCLOCATION_INC) $(CXX_COMMON) \
	-L$(HEPMCLOCATION_LIB) \
	-L$(PYTHIA8LOCATION)/lib/archive -lpythia8 -llhapdfdummy  \
	-lpythia8tohepmc \
	-Wl,-rpath $(HEPMCLOCATION_LIB) -lHepMC

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