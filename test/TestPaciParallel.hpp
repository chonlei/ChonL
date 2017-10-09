/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*
 *
 *  = Running simulation of Paci-ventricular model IN PARALLEL =
 *
 *  Here we run simulation of Paci-ventricular model using different parameters.
 *
 *  Created by Chon Lei
 *  Last edited Apr 2017
 *
 */

#ifndef TESTPACIPARALLEL_HPP_
#define TESTPACIPARALLEL_HPP_

/*
 * == Include headers ==
 *
 */
#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "paci_hyttinen_aaltosetala_severi_ventricularVersionCvode.hpp"
/* This test is always run sequentially (never in parallel)*/
//#include "FakePetscSetup.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

// Standard Chaste classes
#include "OutputFileHandler.hpp"

#include "PetscTools.hpp"

// Should always be last
#include "PetscSetupAndFinalize.hpp"

/* Include this section if to_string is not working in current compiler*/
namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

using namespace std;




/*
 * Define the test class, which inherited from {{{CxxTest::TestSuite}}}
 * as usual, and the (public) test method
 *
 */
class TestParameterSweepParallel : public CxxTest::TestSuite
{
public:
    void TestPaciSimulation() throw(Exception)
    {
#ifdef CHASTE_CVODE
        
        // Set up an output directory - must be called collectively by every process.
        boost::shared_ptr<OutputFileHandler> p_base_handler(new OutputFileHandler("out_PaciSweep", true)); // true wipes the folder!
        PetscTools::Barrier("Output folder created"); // Shouldn't be needed but seemed to avoid an error once!
        std::cout << "\n\nOutput folder created\n\n";


        // Pre-read state variables starting from Steady State 
        ifstream sssv_infile("../sweep_input/PaciSSStateVar1Hz_2");
        double sssv_temp;
        std::vector<double> readSSStateVariable;
        while (sssv_infile >> sssv_temp ) {
            readSSStateVariable.push_back(sssv_temp);
        }
        sssv_infile.close();


        // Everyone does their own thing now
        PetscTools::IsolateProcesses(true);
        



        /*
         * == Define a CVODE model ==
         *
         */
        boost::shared_ptr<RegularStimulus> p_stimulus;
        // The parameters are magnitude, duration, period, and start time of stimulus.
        //boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(0, 2.0, 1000.0, 500));
	//boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-25.5, 2.0, 1000.0, 0));
        // Setup a CVODE model that has empty solver (which requires stimulus)
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new Cellpaci_hyttinen_aaltosetala_severi_ventricularVersionFromCellMLCvode(p_solver, p_stimulus));

        // Use default stimulus if p_stimulus not set
        // To check default stimulus: {{{p_model->HasCellMLDefaultStimulus()}}}
        boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
        double stimPeriod = 905.1;
        p_regular_stim->SetPeriod(stimPeriod); // 1000.0 1375.08 1176.23 932.9 905.1/735.9 


        /*
         * Some optional useful methods:
         * {{{double max_timestep = p_regular_stim->GetDuration();}}}
         * {{{p_model->SetMaxSteps(1e5);}}}
         * {{{p_model->SetTolerances(1e-6,1e-8);}}} // default being (1e-5,1e-7)
         * {{{p_model->ForceUseOfNumericalJacobian();}}}
         *
         */
        p_model->SetTolerances(1e-7,1e-9);


        /*
         * == Get Model Parameters ==
         * File .txt should have format of 3 columns in the order of
         * g_Na_	g_CaL_	g_K_sth_
         *
         */
        string inFilePath = "../sweep_input/conductanceDataTest.txt";
        ifstream infile(inFilePath.c_str());
        if (infile.is_open())
        {
        /** Simulation counter */
        int loop_idx, numSim = 0;
        /** Input of model parameters from infile */
        double toSet_gNa_, toSet_gCaL_, toSet_gto_, toSet_gKr_, toSet_gKs_, toSet_gK1_, toSet_gNaCa_, toSet_gNaK_, toSet_gf_, toSet_gCap_;
        // Repeat simulation for each set of model parameters

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// MAKE loop_idx in the inputfile


        while (infile >> loop_idx >> toSet_gNa_ >> toSet_gCaL_ >> toSet_gto_ >> toSet_gKr_ >> toSet_gKs_ >> toSet_gK1_ >> toSet_gNaCa_ >> toSet_gNaK_ >> toSet_gf_ >> toSet_gCap_) {
        	//cout << "Test: " << toSet_gNa_ << toSet_gCaL_ << toSet_gK_ << ".\n";
        	//printf("Testing ratio of original gNa_ %f; gCaL_ %f; gK_ %f. \n",toSet_gNa_,toSet_gCaL_,toSet_gK_);

            // If we are running in parallel share out the drugs between processes.
            if (loop_idx % PetscTools::GetNumProcs() != PetscTools::GetMyRank())
            {
                // Let another processor do this loop
                continue;
            }

            /* == Re-set Model State Variable == 
             * Always starts with Original Model Steady State - potentially faster to get to steady state
             *
             */
            p_model->SetStateVariables(readSSStateVariable);
	    //p_model->SetStateVariables(p_model->GetInitialConditions());


            /*
             * == Change Model Parameters ==
             * Change Parameter(s) and examine the impact on APD.
             *
             * Any parameters that are labelled in the cell model can be changed.
             * Instructions for parameters annotation can be found at [wiki:ChasteGuides/CodeGenerationFromCellML]
             *
             */
            // g_Na_
            p_model->SetParameter("membrane_fast_sodium_current_conductance", toSet_gNa_*3671.2302);

            // g_CaL_
            p_model->SetParameter("membrane_L_type_calcium_current_conductance", toSet_gCaL_*8.635702e-5);

            // g_K1_
            p_model->SetParameter("membrane_inward_rectifier_potassium_current_conductance", toSet_gK1_*28.1492);
            // g_Kr_
            p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", toSet_gKr_*29.8667);
            // g_Ks_
            p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", toSet_gKs_*2.041);
	    // g_to_
            p_model->SetParameter("membrane_transient_outward_current_conductance", toSet_gto_*29.9038);
	    // g_NaCa_
	    p_model->SetParameter("membrane_sodium_calcium_exchanger_current_conductance", toSet_gNaCa_*4900);

	    // g_NaK_
	    p_model->SetParameter("membrane_sodium_potassium_pump_current_permeability", toSet_gNaK_*1.841424);
	    // g_f_
	    p_model->SetParameter("membrane_hyperpolarisation_activated_funny_current_potassium_component_conductance", toSet_gf_*30.10312);
	    // g_Cap_
	    p_model->SetParameter("membrane_calcium_pump_current_conductance", toSet_gCap_*0.4125);


            // Time it
            time_t timer_start = time(0);

            /*
             * == Run model to steady state ==
             *
             * TO BE UPDATED FOR SELF-STIMULUS MODEL
             *
             */
	    //boost::shared_ptr<RegularStimulus> p_reg_stim =
            //             boost::static_pointer_cast<RegularStimulus>(p_model->GetStimulusFunction());
	    //const double pacing_cycle_length = p_reg_stim->GetPeriod(); //ms
	    //double maximum_time_step = p_reg_stim->GetDuration();  // ms
	    //p_model->Solve(0, 1000*1000, 0.01);



            SteadyStateRunner steady_runner(p_model,true);
            steady_runner.SetMaxNumPaces(100000u);  // Default 1e5
            bool result_ss=true;
            //result_ss = steady_runner.RunToSteadyState();
	    steady_runner.RunToSteadyState();
	    //OdeSolution sol = p_model->Compute(0, 10000, 0.1);


            /*
             * == Get detail for paces of interest ==
             *
             */
            double max_timestep = 0.01;
            p_model->SetMaxTimestep(max_timestep);

            double sampling_timestep = max_timestep*10.0;
            double start_time = 0.0;
            double end_time = stimPeriod; //1000.0;
            OdeSolution solution = p_model->Compute(start_time, end_time, sampling_timestep);

            // Time it till here
            time_t timer_end = time(0);
            double time_taken = difftime(timer_end,timer_start);
            //cout << "Time taken: " << time_taken << " s.\n";

            /*
             * == Write the data out to a file ==
             *
             */
            // This might export unnecessary information.
            //solution.WriteToFile("TestCvodeCells","Paci2011Cvode","ms");

            // Get index of voltage stored
            unsigned voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
            std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);
            std::vector<double> time = solution.rGetTimes();

            // Export whatever is needed

/**          == OLD METHOD ==
            ofstream outfile;
            if (result_ss) {
                string outFilePath = "../out/outputPaciSimulation" + patch::to_string(loop_idx) + ".txt";
            } else {
                string outFilePath = "../out/outputPaciSimulation" + patch::to_string(loop_idx) + "_nonSS" + ".txt";
            }
            outfile.open (outFilePath.c_str());
            if (outfile.is_open())
            {
            for (unsigned i=0; i<voltages.size(); i++) {
                outfile << time[i] << "\t" << voltages[i] << "\n";
            }
            outfile.close();
            }
            else cout << "Unable to open output file.\n";
*/
            std::stringstream this_loop_output_file;
            // If model is in steady state.
            if (result_ss) {
                this_loop_output_file << "outputPaciSimulation" << loop_idx << ".txt";
            } else {
                this_loop_output_file << "outputPaciSimulation" << loop_idx << "_nonSS.txt";
                //std::vector<double> GetStdVecStateVariables()
                //void SetStateVariables(const std::vector<double>& rVariables);
            }
            out_stream p_file = p_base_handler->OpenOutputFile(this_loop_output_file.str());
            for (unsigned i=0; i<voltages.size(); i++) {
                *p_file << time[i] << "\t" << voltages[i] << "\n";
            }
            p_file->close();


            /*
             * == Calculating APD and Upstroke Velocity ==
             *
             * Calculate APD and upstroke velocity using {{{CellProperties}}}
             */
            //CellProperties cell_props(voltages, time);

            //double apd = cell_props.GetLastActionPotentialDuration(90);
            //double upstroke_velocity = cell_props.GetLastMaxUpstrokeVelocity();

            // Check that the values are equal to the ones we expect
            //TS_ASSERT_DELTA(apd, 212.41, 1e-2);
            //TS_ASSERT_DELTA(upstroke_velocity, 338, 1.25);

	    


            numSim++;
            std::cout << "Process " << PetscTools::GetMyRank() << " ran loop index " << loop_idx << "; of " << numSim << "th simulations; taken" << time_taken << "s." << std::endl;
        }
        infile.close();
        }
        else std::cout << "Unable to open file.\n";




        PetscTools::IsolateProcesses(false);
        std::cout << "Run complete." << std::endl;
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

#endif /*TESTPACIPARALLEL_HPP_ */
