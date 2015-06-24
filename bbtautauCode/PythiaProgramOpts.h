#ifndef PYTHIAPROGRAMOPTS_H
#define PYTHIAPROGRAMOPTS_H

#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

using std::cout;
using std::endl;

namespace po = boost::program_options;

class PythiaProgramOpts
{
    private:
        bool printEvent_;
        bool writeToHEPMC_;
        bool writeToLHE_;
        int nEvents_;
        std::string filenameHEPMC_;
        std::string filenameLHE_;
        double mass_;
        bool bbDecay_, tautauDecay_;
        bool verbose_;
        int seed_;

    public:
        // constructor, parses input
        PythiaProgramOpts(int argc, char* argv[]):
            printEvent_(false),
            writeToHEPMC_(false),
            writeToLHE_(false),
            nEvents_(1),
            filenameHEPMC_(""),
            filenameLHE_(""),
            mass_(15.),
            bbDecay_(true),
            tautauDecay_(true),
            verbose_(false),
            seed_(0)
        {
            po::options_description desc("\nProduces gg -> h(125) -> AA, with " \
                "A->2b or A-> 2 tau.\nAllowed options:");
            desc.add_options()
                ("help,h", "Produce help message")
                ("printEvent", po::bool_switch(&printEvent_)->default_value(false),
                    "Prints complete event listing of first event to screen")
                ("hepmc", po::bool_switch(&writeToHEPMC_)->default_value(false),
                    "write events to file in HepMC format")
                ("lhe", po::bool_switch(&writeToLHE_)->default_value(false),
                    "write events to file in LHE format")
                ("number,n", po::value<int>(&nEvents_)->default_value(1),
                    "Number of events to run over [default = 1]. " \
                    "If writeHLT enabled, counts # events passing HLT. " \
                    "Otherwise, counts # events with 2+ muons.")
                ("nameHEPMC", po::value<std::string>(&filenameHEPMC_),
                    "Filename for output HepMC filename. " \
                    "If you don't provide a value but enable --hepmc, " \
                    "the default filename will be ma1_<mass>_<seed>.hepmc")
                ("nameLHE", po::value<std::string>(&filenameLHE_),
                    "Filename for output LHE filename. " \
                    "If you don't provide a value but enable --lhe, " \
                    "the default filename will be ma1_<mass>_<seed>.lhe")
                ("mass", po::value<double>(&mass_)->default_value(15),
                    "Mass of a1 boson in GeV")
                ("bb", po::bool_switch(&bbDecay_)->default_value(false),
                    "Turn on a1 -> bb decay mode")
                ("tautau", po::bool_switch(&tautauDecay_)->default_value(false),
                    "Turn on a1 -> tautau decay mode")
                ("seed", po::value<int>(&seed_)->default_value(0),
                    "Seed for random number generator. 0 = uses time. " \
                    "WARNING: DON'T USE 0 FOR BATCH SYSTEM. " \
                    "Get simultaneous start = same seed = same events. " \
                    "Set seed explicitly instead (e.g. file number).")
                ("verbose,v", po::bool_switch(&verbose_)->default_value(false),
                    "Output debugging statements")
            ;

            po::variables_map vm;
            try {
                po::store(po::parse_command_line(argc, argv, desc), vm);
            } catch (boost::program_options::invalid_option_value e) {
                cout << "Invalid option value: " << e.what() << endl;
                cout << desc << endl;
                cout << "Exiting" << endl;
                exit(1);
            } catch (boost::program_options::unknown_option e) {
                cout << "Unrecognised option: " << e.what() << endl;
                cout << desc << endl;
                cout << "Exiting" << endl;
                exit(1);
            }

            po::notify(vm);

            if (vm.count("help")) {
                cout << desc << endl;
                exit(1);
            }

            // Setup filenames
            if (filenameHEPMC_ == "") {
                filenameHEPMC_ = "ma1_" + boost::lexical_cast<std::string>(mass_) +
                                 "_" + boost::lexical_cast<std::string>(seed_) + ".hepmc";
            }
            if (filenameLHE_ == "") {
                filenameLHE_ = "ma1_" + boost::lexical_cast<std::string>(mass_) +
                               "_" + boost::lexical_cast<std::string>(seed_) + ".lhe";
            }

            // Check if there's already an extension on filename, if not add one
            std::string filenameHEPMClower(filenameHEPMC_);
            boost::algorithm::to_lower(filenameHEPMClower);
            if(!boost::algorithm::ends_with(filenameHEPMClower, ".hepmc")) {
                filenameHEPMC_ += ".hepmc";
            }

            // Check if there's already an extension on filename, if not add one
            std::string filenameLHElower(filenameLHE_);
            boost::algorithm::to_lower(filenameLHElower);
            if(!boost::algorithm::ends_with(filenameLHElower, ".lhe")) {
                filenameLHE_ += ".lhe";
            }

            // check we have a decay channel
            if(!(bbDecay_ || tautauDecay_)) {
                cout << "No decay channel specified" << endl;
                exit(1);
            }

        } // end of constructor

        // Getters
        bool printEvent() { return printEvent_; }
        bool writeToHEPMC() { return writeToHEPMC_; }
        bool writeToLHE() { return writeToLHE_; }
        int nEvents() { return nEvents_; }
        std::string filenameHEPMC() { return filenameHEPMC_; }
        std::string filenameLHE() { return filenameLHE_; }
        double mass() { return mass_; }
        bool bbDecay() { return bbDecay_; }
        bool tautauDecay() { return tautauDecay_; }
        bool verbose() { return verbose_; }
        int seed() { return seed_; }

        // This should really be in a separate .cc file...
        void printProgramOptions() {
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "PYTHIA PROGRAM OPTIONS" << endl;
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            if (writeToHEPMC_)
                cout << "Writing events to hepmc file " << filenameHEPMC_ << endl;
            if (writeToLHE_)
                cout << "Writing events to lhe file " << filenameLHE_ << endl;
            cout << "Doing " << nEvents_ << " events" << endl;
            cout << "Random seed: " << seed_ << endl;
            cout << "Mass of a1: " << mass_ << endl;
            cout << "Decays:" << endl;
            if (bbDecay_)
                cout << "a1 -> bb" << endl;
            if (tautauDecay_)
                cout << "a1 -> tautau" << endl;
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        }

};

#endif
