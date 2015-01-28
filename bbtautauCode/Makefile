#
# Examples Makefile.
#
#                  M. Kirsanov 07.04.2006
#                     Modified 18.11.2006
#                     26.03.2008 CLHEP dependency removed

SHELL = /bin/sh
CXX = g++
# Set p8 and FASTJETLOCATION yourself
p8 = /Users/robina/pythia8/pythia8205/
FASTJETLOCATION = /Users/robina/fastjet-install/

 -include $(p8)/config.mk
# -include config.mk

# Location of directories.
INCDIR=$(p8)/include
LIBDIR=$(p8)/lib
LIBDIRARCH=$(p8)/lib/archive
BINDIR=$(p8)/bin

# There is no default behaviour, so remind user.
# Run using e.g. make bbtata
.cc:

	$(CXX) -I$(INCDIR) $@.cc  \
	`$(FASTJETLOCATION)/bin/fastjet-config --cxxflags  --plugins` \
	$(CXXFLAGSNO) -o $@.exe \
	-L$(LIBDIR) -lpythia8 \
	-L$(FASTJETLOCATION)/lib \
	`$(FASTJETLOCATION)/bin/fastjet-config --libs --plugins` \
	-pthread -std=c++11
	# -llhapdfdummy \


#	$(CXX) $(CXXFLAGS) -I$(INCDIR) $@.cc -o $@.exe  \
#	-L$(LIBDIRARCH) -lpythia8 -llhapdfdummy
#	mv $(BINDIR)/$@.exe $@.exe


# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	# rm -rf $(BINDIR)
	rm -rf *.exe
	rm -f *~; rm -f \#*; rm -f core*
