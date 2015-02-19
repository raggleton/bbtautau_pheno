#include "PythiaProgramOpts.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


using namespace boost::algorithm;
using namespace std;
using namespace Pythia8;


/**
 * @brief Get current time & date
 * @return std::string with time & date
 */
std::string getCurrentTime() {
    time_t now = time(0);
    char* dt = ctime(&now);
    std::string str1 = std::string(dt);
    trim(str1);
    return str1;
}

/**
 * @brief Class to run Pythia simulation of gg -> h(125) -> 2a, a->2b or 2tau
 */
class Simul {

public:
    Simul(PythiaProgramOpts opts);
    ~Simul();
    void run(int nEvnt);
private:
    // Interface for conversion from Pythia8::Event to HepMC event.
    HepMC::Pythia8ToHepMC ToHepMC;
    HepMC::IO_GenEvent * ascii_io;

    Pythia pythia;
    PythiaProgramOpts opts_;
    // Text file to write progress - handy for monitoring during jobs
    ofstream myfile;
};


Simul::Simul(PythiaProgramOpts opts):
    opts_(opts)
{
    if (opts_.writeToHEPMC()) {
        ascii_io = new HepMC::IO_GenEvent(opts_.filename(), std::ios::out);
    }

    // Setup file to write progress
    std::string ext = ".hepmc";
    std::string stem = opts_.filename().substr(0, opts_.filename().size() - ext.size());
    myfile.open(stem + "_progress.txt");

    ///////////////////////
    // SETUP PYTHIA HERE //
    ///////////////////////
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + boost::lexical_cast<std::string>(opts_.seed()));
    pythia.readString("Next:numberShowEvent = 00");

    /**
     * Setup pythia to produce gg->h(125)->AA,
     * and A->bb or A->tautau
     */
    pythia.readString("HiggsSM:gg2H = on");
    pythia.readString("25:m0 = 125.");
    for (int i = 0; i < 76; i++) {
        // turn off all h(125) decay modes
        stringstream SSt;
        SSt << i;
        pythia.readString("25:" + SSt.str() + ":onMode = off");
    }
    pythia.readString("25:addChannel = 1 .05 100 36 36");
    pythia.readString("36:m0 = " + boost::lexical_cast<std::string>(opts_.mass()));
    pythia.readString("36:mMin = 9.5");
    pythia.readString("36:mWidth = 0.1");
    pythia.readString("36:oneChannel = 1 0.5 100 5 -5");
    pythia.readString("36:addChannel = 1 0.5 100 15 -15");
    pythia.readString("36:1:bRatio = 0.5");
    pythia.readString("36:2:bRatio = 0.5");
    pythia.readString("Beams:eCM = 13000.");
    // pythia.readString("PartonLevel:all = off");
    // pythia.readString("PartonLevel:ISR = on");
    // pythia.readString("PartonLevel:FSR = on");
    // pythia.readString("PartonLevel:all = off");
    pythia.init();

}


Simul::~Simul() {
}

void Simul::run(int nEvnt) {
    int outputEvery = 50;
    for (int iEvent = 0; iEvent < nEvnt; ++iEvent) {
        // output progress info
        if (iEvent % outputEvery == 0) {
            cout << "iEvent: " << iEvent << " - " << getCurrentTime() << endl;
            myfile << "iEvent: " << iEvent << " - " << getCurrentTime() << endl;
        }

        // Generate event safely
        if (!pythia.next()) {
            break;
        }

        // Output to screen if wanted
        if (iEvent < 2 && opts_.printEvent()) {
            pythia.info.list();
            pythia.event.list();
            pythia.process.list();
        }

        // Construct new empty HepMC event and fill it.
        // Write the HepMC event to file. Done with it.
        if (opts_.writeToHEPMC()) {
            HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
            ToHepMC.fill_next_event(pythia, hepmcevt);
            *ascii_io << hepmcevt;
            delete hepmcevt;
        }

    }

    myfile.close();
    pythia.stat();
}
