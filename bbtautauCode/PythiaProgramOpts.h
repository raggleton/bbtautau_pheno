#ifndef PYTHIAPROGRAMOPTS_H
#define PYTHIAPROGRAMOPTS_H

#include <iostream>

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
        int nEvents_;
        std::string filename_;
        double mass_;
        bool verbose_;
        int seed_;

    public:
        // constructor, parses input
        PythiaProgramOpts(int argc, char* argv[]):
            printEvent_(false),
            writeToHEPMC_(false),
            nEvents_(1),
            filename_(""),
            mass_(15.),
            verbose_(false),
            seed_(0)
        {
            po::options_description desc("\nProduces gg -> h(125) -> AA, with " \
                "A->2b or A-> 2 tau.\nAllowed options:");
            desc.add_options()
                ("help,h", "Produce help message")
                ("printEvent", po::bool_switch(&printEvent_)->default_value(false),
                    "Prints complete event listing of first event to screen")
                ("write", po::bool_switch(&writeToHEPMC_)->default_value(false),
                    "write events to file in HepMC format")
                ("number,n", po::value<int>(&nEvents_)->default_value(1),
                    "Number of events to run over [default = 1]. " \
                    "If writeHLT enabled, counts # events passing HLT. " \
                    "Otherwise, counts # events with 2+ muons.")
                ("name", po::value<std::string>(&filename_),
                    "Filename for output HepMC filename. " \
                    "If you don't provide a value but enable -write, " \
                    "the default filename will be ma1_<mass>_<seed>.hepmc")
                ("mass", po::value<double>(&mass_)->default_value(15),
                    "Mass of a1 boson in GeV")
                ("seed", po::value<int>(&seed_),
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

            // Setup filename
            if (writeToHEPMC_ && filename_ == "") {
                filename_ = "ma1_" + boost::lexical_cast<std::string>(mass_) +
                            "_" + boost::lexical_cast<std::string>(seed_) + ".hepmc";
            }

            // Check if there's already an extension on filename, if not add one
            if(!boost::algorithm::ends_with(filename_, ".hepmc")) {
                filename_ += ".hepmc";
            }


        } // end of constructor

        // Getters
        bool printEvent() { return printEvent_; }
        bool writeToHEPMC() { return writeToHEPMC_; }
        int nEvents() { return nEvents_; }
        std::string filename() { return filename_; }
        double mass() { return mass_; }
        bool verbose() { return verbose_; }
        int seed() { return seed_; }

        // This should really be in a separate .cc file...
        void printProgramOptions() {
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "PYTHIA PROGRAM OPTIONS" << endl;
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

            if (printEvent_)
                cout << "Outputting first event" << endl;
            if (writeToHEPMC_)
                cout << "Writing events to hepmc file " << filename_ << endl;
            cout << "Doing " << nEvents_ << " events" << endl;
            cout << "Random seed: " << seed_ << endl;
            cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        }

};

#endif
