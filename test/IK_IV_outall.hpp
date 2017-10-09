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

#ifndef IK_IV_CURVE_CLAMP_HPP_
#define IK_IV_CURVE_CLAMP_HPP_

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
    void TestPaciIK_IV() throw(Exception)
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
	double holding_potential = -80.01;
	double steady_state_time = 100000;
        double final_potential = -80.01;
	double holdoutput_time = -100;
        double start_time = 0.0;
        double end_time = 500;  // duration [ms] for the test_potentials
        double tail_time = 1000;  // duration for final_potential = (tail_time - end_time) [ms]

	std::vector<double> test_potentials;
	for (int i=0; i<20; i++) {
		test_potentials.push_back(0.01 + 0.01*(-4000.0 + 500*i));
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
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(0, 2.0, 1000.0, 500)); // i.e. no stimulus
        // Setup a CVODE model that has empty solver (which requires stimulus)
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        // use model as <class> for boost::shared_ptr, st p_model can access functions defined in the model, such as SetNaiDerivativeToZero() and SetCaiDerivativeToZero().
        boost::shared_ptr<Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode> p_model(new Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode(p_solver, p_stimulus));

        // set model to imitate experimental condition
        /************************************ Set experimental used clamped concentration values and temperature here ************************************/
        p_model->SetParameter("temperature",/**__value_from_exp__*/298.15);
        p_model->SetParameter("cytosolic_potassium_concentration",/**__value_from_exp__*/110.0);
        p_model->SetParameter("extracellular_calcium_concentration",/**__value_from_exp__*/1.2);
        p_model->SetParameter("extracellular_potassium_concentration",/**__value_from_exp__*/4.0);
        p_model->SetParameter("extracellular_sodium_concentration",/**__value_from_exp__*/150.0);
        p_model->SetStateVariable("cytosolic_sodium_concentration",/**__value_from_exp__*/10.0);
        p_model->SetStateVariable("cytosolic_calcium_concentration",/**__value_from_exp__*/0.1 * 1.95916e-05);  // put it lower? e.g. 0.1 * 1.95916e-05
        //p_model->SetStateVariable("JSR_calcium_concentration",/**__value_from_exp__*/);
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
        p_model->SetTolerances(1e-9,1e-11);



        // OdeSolution::CalculateDerivedQuantitiesAndParameters() only accept this form
        Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode paci_cvode_system(p_solver,p_stimulus);
        // repeat what is in the p_model
        /************************************ Set experimental used clamped concentration values and temperature so as here ************************************/
        paci_cvode_system.SetParameter("temperature",/**__value_from_exp__*/298.15);
        paci_cvode_system.SetParameter("cytosolic_potassium_concentration",/**__value_from_exp__*/110.0);
        paci_cvode_system.SetParameter("extracellular_calcium_concentration",/**__value_from_exp__*/1.2);
        paci_cvode_system.SetParameter("extracellular_potassium_concentration",/**__value_from_exp__*/4.0);
        paci_cvode_system.SetParameter("extracellular_sodium_concentration",/**__value_from_exp__*/150.0);
        paci_cvode_system.SetStateVariable("cytosolic_sodium_concentration",/**__value_from_exp__*/10.0);
        paci_cvode_system.SetStateVariable("cytosolic_calcium_concentration",/**__value_from_exp__*/0.1 * 1.95916e-05);  // put it lower? e.g. 0.1 * 1.95916e-05
        //paci_cvode_system.SetStateVariable("JSR_calcium_concentration",/**__value_from_exp__*/);
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
        double toSet_gNa_, toSet_gCaL_, toSet_gto_, toSet_gKr_, toSet_gKs_;
        // Repeat simulation for each set of model parameters

        /***************************************** Main loop for parameter sweep *****************************************/


        while (infile >> loop_idx >> toSet_gNa_ >> toSet_gCaL_ >> toSet_gto_ >> toSet_gKr_ >> toSet_gKs_) {
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
            p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", toSet_gKr_*29.8667);
            // g_Ks_
            p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", toSet_gKs_*2.041);
	    // g_to_
            p_model->SetParameter("membrane_transient_outward_current_conductance", toSet_gto_*29.9038);


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
            // g_to_
            paci_cvode_system.SetParameter("membrane_transient_outward_current_conductance", p_model->GetParameter("membrane_transient_outward_current_conductance"));



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








            OdeSolution solutionHold = p_model->Compute(holdoutput_time, start_time, 0.1);

                    // done simulation;  get output for post-processing
                    paci_cvode_system.SetVoltage(holding_potential);  // make sure we are at the correctly clamped potential - very important
                    solutionHold.CalculateDerivedQuantitiesAndParameters(&paci_cvode_system);









	    std::vector<double> holding_state;
	    // save the steady state; reset the model to this steady state later
	    CopyToStdVector(p_model->rGetStateVariables(),holding_state);

	    // set fine simulation
	    double max_timestep = 0.01;
	    p_model->SetMaxTimestep(max_timestep);

	    double sampling_timestep = max_timestep*10.0;

	    // run from the steady state, one per test potential defined in the test_potentials.
	    for (unsigned i=0; i<test_potentials.size(); i++) {
                    // Simulation 1: set clamp voltage at test_potentials
		    // set to steady state
		    p_model->SetStateVariables(holding_state);
		    // set to test_potentials
		    p_model->SetVoltage(test_potentials[i]);
                    //p_model->SetVoltageDerivativeToZero(true);

		    OdeSolution solutionStep = p_model->Compute(start_time, end_time, sampling_timestep);

                    // Simulation 2: set clamp voltage at final_potential
                    // in case it doesn't work...
                    std::vector<double> step_state;
                    // save the steady state; reset the model to this steady state later
                    CopyToStdVector(p_model->rGetStateVariables(),step_state);
		    //p_model->SetStateVariables(step_state);
                    // set to final_potential
		    p_model->SetVoltage(final_potential);
                    //p_model->SetVoltageDerivativeToZero(true);

		    OdeSolution solutionTail = p_model->Compute(end_time, tail_time, sampling_timestep);


                    // done simulation;  get output for post-processing
                    paci_cvode_system.SetStateVariables(holding_state);
                    paci_cvode_system.SetVoltage(test_potentials[i]);  // make sure we are at the correctly clamped potential - very important
                    solutionStep.CalculateDerivedQuantitiesAndParameters(&paci_cvode_system);
                    paci_cvode_system.SetStateVariables(step_state);
                    paci_cvode_system.SetVoltage(final_potential);  // make sure we are at the correctly clamped potential - very important
                    solutionTail.CalculateDerivedQuantitiesAndParameters(&paci_cvode_system);










                    // interested quantities (I_K currents)
                    std::vector<double> iK1Step = solutionStep.GetAnyVariable("membrane_inward_rectifier_potassium_current");
                    std::vector<double> iKrStep = solutionStep.GetAnyVariable("membrane_rapid_delayed_rectifier_potassium_current");
                    std::vector<double> iKsStep = solutionStep.GetAnyVariable("membrane_slow_delayed_rectifier_potassium_current");
                    std::vector<double> iKtoStep = solutionStep.GetAnyVariable("membrane_transient_outward_current");

                    std::vector<double> iCapStep = solutionStep.GetAnyVariable("membrane_calcium_pump_current");
                    std::vector<double> ifStep = solutionStep.GetAnyVariable("membrane_hyperpolarisation_activated_funny_current");
                    std::vector<double> iNaCaStep = solutionStep.GetAnyVariable("membrane_sodium_calcium_exchanger_current");
                    std::vector<double> iNaKStep = solutionStep.GetAnyVariable("membrane_sodium_potassium_pump_current");

                    std::vector<double> iCaLStep = solutionStep.GetAnyVariable("membrane_L_type_calcium_current");
                    std::vector<double> iNaStep = solutionStep.GetAnyVariable("membrane_fast_sodium_current");


		    std::vector<double> timeStep = solutionStep.rGetTimes();





/**
		    std::vector<double> tempVector;
                    // should be clamped
		    std::vector<double> voltagesStep = solutionStep.GetVariableAtIndex(voltage_index);
		    std::vector<double> CaiStep = solutionStep.GetVariableAtIndex(Cai_index);
		    std::vector<double> NaiStep = solutionStep.GetVariableAtIndex(Nai_index);
		    std::vector<double> voltagesTail = solutionTail.GetVariableAtIndex(voltage_index);
		    std::vector<double> CaiTail = solutionTail.GetVariableAtIndex(Cai_index);
		    std::vector<double> NaiTail = solutionTail.GetVariableAtIndex(Nai_index);
                    // combine step and tail simulations for clamped quantities
                    voltagesStep.insert(voltagesStep.end(), voltagesTail.begin(), voltagesTail.end());
                    std::vector<double>().swap(voltagesTail);
                    CaiStep.insert(CaiStep.end(), CaiTail.begin(), CaiTail.end());
                    std::vector<double>().swap(CaiTail);
                    NaiStep.insert(NaiStep.end(), NaiTail.begin(), NaiTail.end());
                    std::vector<double>().swap(NaiTail);
                    // interested quantities (I_K currents)
                    std::vector<double> iK1Step = solutionStep.GetAnyVariable("membrane_inward_rectifier_potassium_current");
                    std::vector<double> iKrStep = solutionStep.GetAnyVariable("membrane_rapid_delayed_rectifier_potassium_current");
                    std::vector<double> iKsStep = solutionStep.GetAnyVariable("membrane_slow_delayed_rectifier_potassium_current");
                    std::vector<double> iK1Tail = solutionTail.GetAnyVariable("membrane_inward_rectifier_potassium_current");
                    std::vector<double> iKrTail = solutionTail.GetAnyVariable("membrane_rapid_delayed_rectifier_potassium_current");
                    std::vector<double> iKsTail = solutionTail.GetAnyVariable("membrane_slow_delayed_rectifier_potassium_current");
		    std::vector<double> timeStep = solutionStep.rGetTimes();
		    std::vector<double> timeTail = solutionTail.rGetTimes();
                    // combine step and tail simulations for interested quantities
                    iK1Step.insert(iK1Step.end(), iK1Tail.begin(), iK1Tail.end());
                    std::vector<double>().swap(iK1Tail);
                    iKrStep.insert(iKrStep.end(), iKrTail.begin(), iKrTail.end());
                    std::vector<double>().swap(iKrTail);
                    iKsStep.insert(iKsStep.end(), iKsTail.begin(), iKsTail.end());
                    std::vector<double>().swap(iKsTail);
                    timeStep.insert(timeStep.end(), timeTail.begin(), timeTail.end());
                    std::vector<double>().swap(timeTail);

**/
		    // export result here...
                    std::stringstream this_loop_output_file;
                    this_loop_output_file << "outPaciICaL" << test_potentials[i] << ".txt";
                    cout << this_loop_output_file.str() << " outputing format: time \\t \\t iK1 \\t iKr \\t iKs \\t iKto \\t iCap \\t if \\t iNaCa \\t iNaK \\t iCaL \\t iNa " << endl;

                    out_stream p_file = p_base_handler->OpenOutputFile(this_loop_output_file.str());
                    for (unsigned i=0; i<timeStep.size(); i++) {
                            // SomethingStep but it is actually combined result from step and tail simulations
			    *p_file << timeStep[i] << "\t" << iK1Step[i] << "\t" << iKrStep[i]<< "\t" << iKsStep[i] << "\t" << iKtoStep[i] << "\t" << iCapStep[i]<< "\t" << ifStep[i] << "\t" << iNaCaStep[i] << "\t" << iNaKStep[i] << "\t" << iCaLStep[i]<< "\t" << iNaStep[i] << "\n";
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

#endif /*IK_IV_CURVE_CLAMP_HPP_ */
