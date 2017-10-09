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

#ifndef ICAL_IV_CURVE_CLAMP_HPP_
#define ICAL_IV_CURVE_CLAMP_HPP_

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
// Use modified version which allows clamping ion concentration
#include "paci_hyttinen_aaltosetala_severi_ventricularVersionAllowClampCvode.hpp"
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


using namespace std;


/*
 * Define the test class, which inherited from {{{CxxTest::TestSuite}}}
 * as usual, and the (public) test method
 *
 */
class TestParameterSweepParallel : public CxxTest::TestSuite
{
public:
    void TestPaciICaL_IV() throw(Exception)
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



        // here should match experimental setup
	double holding_potential = -80;
	double steady_state_time = 100000;
        double start_time = 0.0;
        double end_time = 200;  // duration [ms] for the test_potentials

	std::vector<double> test_potentials;
	for (int i=0; i<25; i++) {
		test_potentials.push_back(0.01 + 0.01*(-6000.0 + 500*i));
	}


        /*
         * == Parallel Computing ==
         * Everyone does their own thing now
         *
         */
        PetscTools::IsolateProcesses(true);
        

        /*
         * == Define a CVODE model ==
         *
         */
        //boost::shared_ptr<RegularStimulus> p_stimulus;
        // The parameters are magnitude, duration, period, and start time of stimulus.
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(0, 2.0, 1000.0, 500));
        // Setup a CVODE model that has empty solver (which requires stimulus)
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        // use model as <class> for boost::shared_ptr, st p_model can access functions defined in the model, such as SetNaiDerivativeToZero() and SetCaiDerivativeToZero().
        boost::shared_ptr<Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode> p_model(new Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode(p_solver, p_stimulus));

        // set model to imitate experimental condition
        /************************************ Set experimental used clamped concentration values and temperature here ************************************/
        p_model->SetParameter("temperature",/**__value_from_exp__*/298.15);
        p_model->SetParameter("cytosolic_potassium_concentration",/**__value_from_exp__*/140.0);
        p_model->SetParameter("extracellular_calcium_concentration",/**__value_from_exp__*/2.0);
        p_model->SetParameter("extracellular_potassium_concentration",/**__value_from_exp__*/5.0);
        p_model->SetParameter("extracellular_sodium_concentration",/**__value_from_exp__*/140.0);
        p_model->SetStateVariable("cytosolic_sodium_concentration",/**__value_from_exp__*/10.0);
        p_model->SetStateVariable("cytosolic_calcium_concentration",/**__value_from_exp__*/0.1 * 1.95916e-05);  // put it lower? e.g. 0.1 * 1.95916e-05
        p_model->SetStateVariable("JSR_calcium_concentration",/**__value_from_exp__*/0.295465);
	p_model->SetVoltageDerivativeToZero(true);
	p_model->SetNaiDerivativeToZero(true);
	p_model->SetCaiDerivativeToZero(true);

        // Use default stimulus if p_stimulus not set
        // To check default stimulus: {{{p_model->HasCellMLDefaultStimulus()}}}
        //boost::shared_ptr<RegularStimulus> p_regular_stim = p_model->UseCellMLDefaultStimulus();
        //p_regular_stim->SetPeriod(1000.0);

        /*
         * Some optional useful methods:
         * {{{double max_timestep = p_regular_stim->GetDuration();}}}
         * {{{p_model->SetMaxSteps(1e5);}}}
         * {{{p_model->SetTolerances(1e-6,1e-8);}}} // default being (1e-5,1e-7)
         * {{{p_model->ForceUseOfNumericalJacobian();}}}
         *
         */
        p_model->SetTolerances(1e-7,1e-9);



        // OdeSolution::CalculateDerivedQuantitiesAndParameters() only accept this form
        Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode paci_cvode_system(p_solver,p_stimulus);
        // repeat what is in the p_model
        /************************************ Set experimental used clamped concentration values and temperature so as here ************************************/
        paci_cvode_system.SetParameter("temperature",/**__value_from_exp__*/298.15);
        paci_cvode_system.SetParameter("cytosolic_potassium_concentration",/**__value_from_exp__*/140.0);
        paci_cvode_system.SetParameter("extracellular_calcium_concentration",/**__value_from_exp__*/2.0);
        paci_cvode_system.SetParameter("extracellular_potassium_concentration",/**__value_from_exp__*/5.0);
        paci_cvode_system.SetParameter("extracellular_sodium_concentration",/**__value_from_exp__*/140.0);
        paci_cvode_system.SetStateVariable("cytosolic_sodium_concentration",/**__value_from_exp__*/10.0);
        paci_cvode_system.SetStateVariable("cytosolic_calcium_concentration",/**__value_from_exp__*/0.1 * 1.95916e-05);  // put it lower? e.g. 0.1 * 1.95916e-05
        paci_cvode_system.SetStateVariable("JSR_calcium_concentration",/**__value_from_exp__*/0.295465);
        paci_cvode_system.SetVoltageDerivativeToZero(true);
	paci_cvode_system.SetNaiDerivativeToZero(true);
	paci_cvode_system.SetCaiDerivativeToZero(true);


        /*
         * == Get Model Parameters ==
         * File .txt should have format of 3 columns in the order of
         * idx	g_Na_	g_CaL_	g_K_sth_
         *
         */
        string inFilePath = "../sweep_input/conductanceDataTest.txt";
        ifstream infile(inFilePath.c_str());

        if (infile.is_open())
        {
        /** Simulation counter */
        int loop_idx, numSim = 0;
        /** Input of model parameters from infile */
        double toSet_gNa_, toSet_gCaL_, toSet_gK_;
        // Repeat simulation for each set of model parameters

        /***************************************** Main loop for parameter sweep *****************************************/

        while (infile >> loop_idx >> toSet_gNa_ >> toSet_gCaL_ >> toSet_gK_) {
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
            p_model->SetParameter("membrane_inward_rectifier_potassium_current_conductance", 28.1492);
            // g_Kr_
            p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", toSet_gK_*29.8667);
            // g_Ks_
            p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", 2.041);


            // pass parameters from p_model to paci_cvode_system
            // g_Na_
            paci_cvode_system.SetParameter("membrane_fast_sodium_current_conductance", p_model->GetParameter("membrane_fast_sodium_current_conductance"));
            // g_CaL_
            paci_cvode_system.SetParameter("membrane_L_type_calcium_current_conductance", p_model->GetParameter("membrane_L_type_calcium_current_conductance"));
            // g_K1_
            paci_cvode_system.SetParameter("membrane_inward_rectifier_potassium_current_conductance", p_model->GetParameter("membrane_inward_rectifier_potassium_current_conductance"));
            // g_Kr_
            paci_cvode_system.SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", p_model->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance"));
            // g_Ks_
            paci_cvode_system.SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", p_model->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance"));



            /*
             * == Run model to steady state ==
             *
             * TO BE UPDATED FOR SELF-STIMULUS MODEL
             *
             */
            //SteadyStateRunner steady_runner(p_model,true);
            //steady_runner.SetMaxNumPaces(50000u);  // Default 1e5
            //bool result_ss;
            //result_ss = steady_runner.RunToSteadyState();
            //steady_runner.RunToSteadyState();  // If not checking Steady State.
            // Get number of cell model 'paces/beats' that had to be evaluated.
            //{{{steady_runner.GetNumEvaluations();}}}

            // Detect for steady state alternans
            //{{{SteadyStateRunner steady_runner(p_model, true);}}}

            // Check that the model has NOT reached steady state
            //TS_ASSERT_EQUALS(result_ss,true);
	    //



	    // set holding potential
	    p_model->SetVoltage(holding_potential);
	    //p_model->SetMinimalReset(true); // in ssRunner, needed?
	    // solve here, i.e. run for steady_state_time long to steady it at the holding potential
	    p_model->Compute(0, steady_state_time);


	    std::vector<double> holding_state;
	    // save the steady state; reset the model to this steady state later
	    CopyToStdVector(p_model->rGetStateVariables(),holding_state);

	    // set fine simulation
	    double max_timestep = 0.01;
	    p_model->SetMaxTimestep(max_timestep);

	    double sampling_timestep = max_timestep*10.0;

	    // run from the steady state, one per test potential defined in the test_potentials.
	    for (unsigned i=0; i<test_potentials.size(); i++) {

		    // set to steady state
		    p_model->SetStateVariables(holding_state);
		    // set to test_potentials
		    p_model->SetVoltage(test_potentials[i]);

		    OdeSolution solution = p_model->Compute(start_time, end_time, sampling_timestep);

                    // done simulation;  get output for post-processing
                    paci_cvode_system.SetVoltage(test_potentials[i]);  // make sure we are at the correctly clamped potential - very important
                    solution.CalculateDerivedQuantitiesAndParameters(&paci_cvode_system);

                    // double check they are clamped
		    unsigned voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
                    unsigned Cai_index = p_model->GetSystemInformation()->GetStateVariableIndex("cytosolic_calcium_concentration");

                    // should be clamped
		    std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);
		    std::vector<double> Cai = solution.GetVariableAtIndex(Cai_index);
                    // interested quantity (I_CaL current)
                    std::vector<double> iCaL = solution.GetAnyVariable("membrane_L_type_calcium_current");
		    std::vector<double> time = solution.rGetTimes();


		    /** Manually calculate the ICaL - just for verification - ?? does not match model calculation. From here... *****
                    unsigned gated_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_L_type_calcium_current_d_gate");
                    unsigned gatef1_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_L_type_calcium_current_f_gate");
                    unsigned gatef2_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_L_type_calcium_current_f2_gate");
                    unsigned gatefCa_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_L_type_calcium_current_fCa_gate");
		    double ICaL_P1 = p_model->GetParameter(1);
		    double ICaL_P4 = p_model->GetParameter(4);
		    double ICaL_P17 = p_model->GetParameter(17);
		    // I_CaL
		    std::vector<double> gated = solution.GetVariableAtIndex(gated_index);
		    std::vector<double> gatef1 = solution.GetVariableAtIndex(gatef1_index);
		    std::vector<double> gatef2 = solution.GetVariableAtIndex(gatef2_index);
		    std::vector<double> gatefCa = solution.GetVariableAtIndex(gatefCa_index);
                    // from model code
                    double i_CaL = ((((ICaL_P4 * (4.0 * (voltages[i] * 9309421124.3716221))) / (8.3144720000000003 * ICaL_P17)) * ((Cai[i] * exp((2.0 * (voltages[i] * 96485.341499999995)) / (8.3144720000000003 * ICaL_P17))) - (0.34100000000000003 * ICaL_P1))) / (exp((2.0 * (voltages[i] * 96485.341499999995)) / (8.3144720000000003 * ICaL_P17)) - 1.0)) * (gated[i] * (gatef1[i] * (gatef2[i] * gatefCa[i]))); // A_per_F
		    ***** ...to here */


		    // export result here...
                    std::stringstream this_loop_output_file;
                    this_loop_output_file << "outPaciICaL" << test_potentials[i] << ".txt";
                    cout << this_loop_output_file.str() << " outputing format: time \\t voltages \\t iCaL \\t Cai" << endl;

                    out_stream p_file = p_base_handler->OpenOutputFile(this_loop_output_file.str());
                    for (unsigned i=0; i<voltages.size(); i++) {

			    *p_file << time[i] << "\t" << voltages[i] << "\t" << iCaL[i]<< "\t" << Cai[i] << "\n";
                    }
                    p_file->close();

	    }



            numSim++;
            std::cout << "Process " << PetscTools::GetMyRank() << " ran loop index " << "N/A" << "; of " << "N/A" << "th simulations; taken" << "N/A" << "s." << std::endl;
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

#endif /*ICAL_IV_CURVE_CLAMP_HPP_ */
