SHELL = /bin/sh
CXX = g++

# Location of directories.
INCDIR = include
SRCDIR = src
LIBDIR = lib
LIBDIRARCH = lib/archive
BINDIR = bin

# Set PYTHIA8LOCATION, FASTJETLOCATION, HEPMC yourself
# PYTHIA8LOCATION = /cvmfs/sft.cern.ch/lcg/external/MCGenerators_lcgcmt67b/pythia8/205/x86_64-slc6-gcc48-opt
# PYTHIA8LOCATION = /cvmfs/sft.cern.ch/lcg/external/MCGenerators_lcgcmt67c/pythia8/201/x86_64-slc6-gcc48-opt
# PYTHIA8LOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/pythia8/205
PYTHIA8LOCATION = /users/ra12451/pythia8205
PYTHIA8LOCATION_INC = $(PYTHIA8LOCATION)/$(INCDIR)
PYTHIA8LOCATION_LIB = $(PYTHIA8LOCATION)/$(LIBDIR)

# FASTJETLOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/fastjet/3.0.3-cms
FASTJETLOCATION = /users/ra12451/fastjet-install

# Don't call it HEPMCLOCATION beacuse it defaults to /build/degano not cvmsfs for some stupid reason!!
# HEPMC = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/hepmc/2.06.07-cms3
# HEPMC = /cvmfs/sft.cern.ch/lcg/external/HepMC/2.06.08/x86_64-slc6-gcc48-opt
HEPMC = /users/ra12451/HepMC/install
HEPMCLOCATION_INC = $(HEPMC)/$(INCDIR)
HEPMCLOCATION_LIB = $(HEPMC)/$(LIBDIR)

# For boost libraries
# /cvmfs/sft.cern.ch/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/
# BOOSTLOCATION = /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/boost/1.57.0/
# BOOSTLOCATION = /cvmfs/sft.cern.ch/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/
# BOOSTLOCATION_INC = /usr/include/boost/
BOOSTLOCATION_INC = /cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/boost/1.57.0-cms/include/
BOOSTLOCATION_LIB = /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_4_ROOT5/external/slc6_amd64_gcc491/lib
# BOOSTLOCATION_LIB = /usr/lib64/
BOOST_LIB = -L$(BOOSTLOCATION_LIB) -lboost_program_options -lboost_iostreams -Wl,-rpath $(BOOSTLOCATION_LIB)

# For ROOT libraries
ROOTDIR = /cvmfs/cms.cern.ch/slc6_amd64_gcc491/lcg/root/5.34.22-ilphmn/
ROOT_INCLUDE = $(ROOTDIR)/include/
ROOT_BIN = $(ROOTDIR)/bin/
ROOT_LIB = $(ROOTDIR)/lib/

# CXX_COMMON = -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_COMMON = -std=c++11 -Wno-shadow -pthread -D_USE_XOPEN2K8 -W -Wall

# Make using e.g. make bbtautau
# For pythia 8.2X
.cc:
	$(CXX) $@.cc $(PYTHIA8LOCATION_LIB)/libpythia8.a -o $@.exe \
	-w -I$(ROOT_INCLUDE) \
	-I$(BOOSTLOCATION_INC) \
	-I$(HEPMCLOCATION_INC) \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags --plugins` \
	-I$(PYTHIA8LOCATION_INC) \
	$(CXX_COMMON) \
	-I./ -L./ -Wl,-rpath ./ -lz -Wl,-rpath $(PYTHIA8LOCATION_LIB) -ldl \
	$(BOOST_LIB) \
	-L$(HEPMCLOCATION_LIB) -Wl,-rpath $(HEPMCLOCATION_LIB) -lHepMC \
	`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins` \
	-Wl,-rpath $(ROOT_LIB) `$(ROOT_BIN)root-config --glibs` \


# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	rm *.exe
