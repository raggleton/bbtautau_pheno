For pheno work on bbtautau final state

# Installation
(In this order)

0. __ROOT__

Easiest way is to use MacPorts:

    ```
    sudo port install root5 +cocoa +gcc49 -gcc48 +graphviz +gsl +http +minuit2 +opengl +python27 +roofit +soversion +ssl +tmva +xml
    ```
(need gcc49 & python 27 first)

1. __FASTJET__

    http://fastjet.fr/quickstart.html
    Download & unpack.
    Then in resultant folder:

    ```
    ./configure --prefix=$PWD/../fastjet-install
    make -j4
    make check
    make install
    ```
FastJet now installed in `fastjet-install/`

2. __HEPMC__

(For PYTHIA8)
http://lcgapp.cern.ch/project/simu/HepMC/
Instructions online don't work. Get Mac version?

3. __LHAPDF6__

(For PYTHIA8)
https://www.hepforge.org/downloads/lhapdf
Download & unpack.
Requires BOOST. (install either via MacPorts or manually)
Then in the resultant folder:
    ```
    ./configure --prefix=$PWD/../LHAPDF6-install
    make -j4
    make install
    ```
LHAPDF now installed in `LHAPDF6-install/`. Now need to download the PDF sets you want and put them in `LHAPDF6-install/share/LHAPDF`. Typically: CT10 (CT10nlo?), cteq6l1, MSTW2008.

4. __PYTHIA8__

http://home.thep.lu.se/~torbjorn/Pythia.html
Download & unpack. In resultant folder, do:
    ```
    ./configure --with-fastjet3=<path to>/fastjet-install/ --with-hepmc2=<path to>/HepMC/x86_64-mac106-gcc42-opt/ --with-lhapdf6=<path to>/LHAPDF/LHAPDF-6.1.5-install/ --with-root-bin=/opt/local/libexec/root5/bin/ --with-root-lib=/opt/local/libexec/root5/lib/root/ --with-root-include=/opt/local/libexec/root5/include/root/ --with-boost
    make -j4
    ```
Ugly ROOT options because MacPorts or ROOT puts things in e.g. `include/root/TH1D.h`, not `include/TH1D.h`. May need to specify boost library location.

5. __DELPHES 3__

    ```
    git clone https://github.com/delphes/delphes.git
    cd delphes
    make -j4
    ```
I put mine in $HOME, but can go wherever.

6. __MadGraph5_aMC@NLO__

https://launchpad.net/mg5amcnlo
Download & unpack. No need to compile or anything.
To use FastJet inside:
    ```
    set fastjet /Users/robina/fastjet/fastjet-3.1.1/fastjet-config
    ```
Note, NOT the fastjet-install dir.
7.