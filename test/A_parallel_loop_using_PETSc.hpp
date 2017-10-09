#ifndef TESTAPARALLELLOOP_HPP_
#define TESTAPARALLELLOOP_HPP_

// Includes first
#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <fstream>
#include <string>

// Standard Chaste classes
#include "OutputFileHandler.hpp"

#include "PetscTools.hpp"

// Should always be last
#include "PetscSetupAndFinalize.hpp"

using namespace std;


class TestAParallelLoop : public CxxTest::TestSuite
{
public:

    void TestApdPredictionsForCrumbData() throw (Exception)
    {
        // Set up an output directory - must be called collectively by every process.
        boost::shared_ptr<OutputFileHandler> p_base_handler(new OutputFileHandler("TestLoops", true)); // true wipes the folder!
        PetscTools::Barrier("Output folder created"); // Shouldn't be needed but seemed to avoid an error once!
        std::cout << "\n\nOutput folder created\n\n";

        // Everyone does their own thing now...
        PetscTools::IsolateProcesses(true);

        // Loop over each compound
        for (unsigned loop_idx = 0; loop_idx<100u; loop_idx++)
        {
            // If we are running in parallel share out the drugs between processes.
            if (loop_idx % PetscTools::GetNumProcs() != PetscTools::GetMyRank())
            {
                // Let another processor do this loop
                continue;
            }

            std::stringstream this_loop_output_file;
            this_loop_output_file << "loop " << loop_idx;
            out_stream p_file = p_base_handler->OpenOutputFile(this_loop_output_file.str());
            *p_file << "Process " << PetscTools::GetMyRank() << "wrote the data file for loop " << loop_idx << std::endl;





        double toSet_gNa_, toSet_gCaL_, toSet_gK_;
	string inFilePath = "../out/conductanceDataTest.txt";
        ifstream infile(inFilePath.c_str());
        if (infile.is_open())
        {
        while (infile >> toSet_gNa_ >> toSet_gCaL_ >> toSet_gK_) { *p_file << toSet_gNa_ << toSet_gCaL_ << toSet_gK_ << std::endl;}
        infile.close();
        }
        else std::cout << "Unable to open file.\n";




            p_file->close();
            std::cout << "Process " << PetscTools::GetMyRank() << " did loop index " << loop_idx << std::endl;
        }

        PetscTools::IsolateProcesses(false);
        std::cout << "Run complete." << std::endl;
    }
};

#endif // TESTAPARALLELLOOP_HPP_
