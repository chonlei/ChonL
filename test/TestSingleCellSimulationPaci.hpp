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
 *  = Running simulation of Paci-ventricular model =
 *
 *  Here we run simulation of Paci-ventricular model using different parameters.
 *
 *  Created by Chon Lei
 *  Last edited Apr 2017
 *
 */

#ifndef TESTSINGLECELLSIMULATIONTUTORIAL_HPP_
#define TESTSINGLECELLSIMULATIONTUTORIAL_HPP_

/*
 * == Include headers ==
 *
 */
#include <cxxtest/TestSuite.h>
#include "CellProperties.hpp"
#include "SteadyStateRunner.hpp"
#include "AbstractCvodeCell.hpp"
#include "RegularStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "paci_hyttinen_aaltosetala_severi_ventricularVersionCvode.hpp"
/* This test is always run sequentially (never in parallel)*/
#include "FakePetscSetup.hpp"
#include <iostream>
#include <fstream>
#include <string>

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
class TestSingleCellParameterSweep : public CxxTest::TestSuite
{
public:
    void TestPaciSimulation() throw(Exception)
    {
#ifdef CHASTE_CVODE
        /*
         * == Define a CVODE model ==
         *
         */
        //boost::shared_ptr<RegularStimulus> p_stimulus;
        // The parameters are magnitude, duration, period, and start time of stimulus.
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(0, 2.0, 1000.0, 500));
        // Setup a CVODE model that has empty solver (which requires stimulus)
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        boost::shared_ptr<AbstractCvodeCell> p_model(new Cellpaci_hyttinen_aaltosetala_severi_ventricularVersionFromCellMLCvode(p_solver, p_stimulus));

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

        /*
         * == Get Model Parameters ==
         * File .txt should have format of 3 columns in the order of
         * g_Na_	g_CaL_	g_K_sth_
         *
         */
        string inFilePath = "../out/conductanceDataTest.txt";
        ifstream infile(inFilePath.c_str());
        if (infile.is_open())
        {
        /** Simulation counter */
        int numSim = 0;
        /** Input of model parameters from infile */
        double toSet_gNa_, toSet_gCaL_, toSet_gK_;
        // Repeat simulation for each set of model parameters
        while (infile >> toSet_gNa_ >> toSet_gCaL_ >> toSet_gK_) {
        	cout << "Test: " << toSet_gNa_ << toSet_gCaL_ << toSet_gK_ << ".\n";

            /*
             * == Change Model Parameters ==
             * Change Parameter(s) and examine the impact on APD.
             *
             * Any parameters that are labelled in the cell model can be changed.
             * Instructions for parameters annotation can be found at [wiki:ChasteGuides/CodeGenerationFromCellML]
             *
             */

            // g_Na_
            p_model->SetParameter("membrane_fast_sodium_current_conductance", 3671.2302);

            // g_CaL_
            p_model->SetParameter("membrane_L_type_calcium_current_conductance", 8.635702e-5);

            // g_K1_
            p_model->SetParameter("membrane_inward_rectifier_potassium_current_conductance", 28.1492);
            // g_Kr_
            p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", 29.8667);
            // g_Ks_
            p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", 2.041);


            /*
             * == Run model to steady state ==
             *
             * TO BE UPDATED FOR SELF-STIMULUS MODEL
             *
             */

            SteadyStateRunner steady_runner(p_model);
            steady_runner.SetMaxNumPaces(1000u);  // Default 1e5
            //bool result;
            //result = steady_runner.RunToSteadyState();
            steady_runner.RunToSteadyState();

            // Detect for steady state alternans
            //{{{SteadyStateRunner steady_runner(p_model, true);}}}

            // Check that the model has NOT reached steady state
            //TS_ASSERT_EQUALS(result,false);


            /*
             * == Get detail for paces of interest ==
             *
             */
            double max_timestep = 0.1;
            p_model->SetMaxTimestep(max_timestep);

            double sampling_timestep = max_timestep;
            double start_time = 0.0;
            double end_time = 3500.0;
            OdeSolution solution = p_model->Compute(start_time, end_time, sampling_timestep);

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
            ofstream outfile;
            string outFilePath = "../out/outputPaciSimulation" + patch::to_string(numSim) + ".txt";
            outfile.open (outFilePath.c_str());
            if (outfile.is_open())
            {
            for (unsigned i=0; i<voltages.size(); i++) {
                outfile << time[i] << "\t" << voltages[i] << "\n";
            }
            outfile.close();
            }
            else cout << "Unable to open output file.\n";

            /*
             * == Calculating APD and Upstroke Velocity ==
             *
             * Calculate APD and upstroke velocity using {{{CellProperties}}}
             */
            CellProperties cell_props(voltages, time);

            //double apd = cell_props.GetLastActionPotentialDuration(90);
            //double upstroke_velocity = cell_props.GetLastMaxUpstrokeVelocity();

            // Check that the values are equal to the ones we expect
            //TS_ASSERT_DELTA(apd, 212.41, 1e-2);
            //TS_ASSERT_DELTA(upstroke_velocity, 338, 1.25);

            numSim++;
        }
        infile.close();
        }
        else std::cout << "Unable to open file.\n";
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

#endif /*TESTANOTHERBIDOMAINSIMULATIONTUTORIAL_HPP_*/
