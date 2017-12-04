/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef TESTCIPAOHARAMODEL_HPP_
#define TESTCIPAOHARAMODEL_HPP_

// Standard libraries - might not need
#include <fstream>
#include <iostream>
#include <string>

#include <cxxtest/TestSuite.h>
#include "AbstractCvodeCell.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "CellProperties.hpp"
#include "FileFinder.hpp"
#include "RegularStimulus.hpp"
#include "SteadyStateRunner.hpp"
#include "ohara_rudy_2011_endo_dyHergCvode.hpp"

/* This test is always run sequentially (never in parallel)*/
#include "FakePetscSetup.hpp"

class TestCipaOHaraModelSolver : public CxxTest::TestSuite
{
public:
    void TestCipaVersionOfOharaModelSolver()
    {

#ifdef CHASTE_CVODE

        // Use CiPA stimulus (magnitudeOfStimulus, duration, period, startTime, stopTime = DBL_MAX)
        boost::shared_ptr<RegularStimulus> p_stimulus(new RegularStimulus(-80, 0.5, 2000.0, 0.0)); // Not default but match what CiPA use
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
        // Use CellML version of the CiPA model
        boost::shared_ptr<AbstractCvodeCell> p_model(new Cellohara_rudy_2011_endo_dyHergFromCellMLCvode(p_solver, p_stimulus));

        /*
		 * == Define Variable Order ==
		 *
		 */
        FileFinder variable_order("projects/ChonL/CiPA/Chaste_CiPA_var_order.txt", RelativeTo::ChasteSourceRoot);
        std::ifstream Chaste_CiPA_var_file(variable_order.GetAbsolutePath());
        unsigned tmp_unsigned;
        std::map<unsigned, unsigned> map_Chaste_CiPA_var;
        unsigned var_i = 0;
        while (Chaste_CiPA_var_file >> tmp_unsigned)
        {
            map_Chaste_CiPA_var.insert(std::make_pair(var_i, tmp_unsigned));
            var_i++;
        }
        Chaste_CiPA_var_file.close();

        // Print out what we can see
        std::cout << "Number of Chaste (visible) parameters: " << p_model->GetNumberOfParameters() << "\n";
        std::cout << "Number of variables: " << p_model->GetNumberOfStateVariables() << "\n";
        // Match the number of variables
        TS_ASSERT_EQUALS(var_i, p_model->GetNumberOfStateVariables());

        // Get the name of the state variables (just in case)
        const std::vector<std::string> YName = p_model->rGetStateVariableNames();

        /*
		 * ====== Control Test ======
		 *
		 * Test the model using default setting and control condition.
		 *
		 */
        {
            /* == Read Model State Variable ==
			 * Try to set initial values to match CiPA one
			 *
			 */
            FileFinder state_vars_filename("projects/ChonL/CiPA/Chaste_newordherg_states_CL2000.txt", RelativeTo::ChasteSourceRoot);
            std::ifstream CiPA_stateVariables_file(state_vars_filename.GetAbsolutePath());
            double tmp_double;
            std::vector<double> CiPA_stateVariables;
            while (CiPA_stateVariables_file >> tmp_double)
            {
                CiPA_stateVariables.push_back(tmp_double);
            }
            CiPA_stateVariables_file.close();

            TS_ASSERT_EQUALS(p_model->GetNumberOfStateVariables(), CiPA_stateVariables.size());
            TS_ASSERT_EQUALS(CiPA_stateVariables.size(), map_Chaste_CiPA_var.size());

            // Change var order
            std::vector<double> CiPA_stateVariables_Chaste_order;
            for (unsigned i = 0; i < CiPA_stateVariables.size(); i++)
            {
                CiPA_stateVariables_Chaste_order.push_back(CiPA_stateVariables[map_Chaste_CiPA_var[i]]);
            }

            /* == Re-set Model State Variable ==
			 * Start with the state variables that CiPA use
			 *
			 */
            p_model->SetStateVariables(CiPA_stateVariables_Chaste_order);
            // Try different tolerances
            p_model->SetTolerances(1e-15, 1e-15);
	    p_model->SetMaxSteps(1e8);

            /*
			 * == Getting detail for paces of interest ==
			 * Now we solve for the number of paces we are interested in.
			 *
			 */
            double max_timestep = 5e-1;
            p_model->SetMaxTimestep(max_timestep);

            // Set sampling step size as the one CiPA use
            double sampling_timestep = 1.0;
            double start_time = 0.0;
            double end_time = 2000.0;
            OdeSolution solution = p_model->Compute(start_time, end_time, sampling_timestep);

            // Write the data out to a file.
            solution.WriteToFile("TestCipaOHaraModel", "OHaraDyHergCvode", "ms");

            unsigned voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
            std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);

	    /*
	     		 * == Get Reference Solution ==
			 */
	    std::vector<double> state_variable_ref;
	    for (unsigned i = 0; i < YName.size(); i++)
	    {
		unsigned state_variable_index = p_model->GetSystemInformation()->GetStateVariableIndex(YName[i]);
		std::vector<double> state_variable_sol = solution.GetVariableAtIndex(state_variable_index);
		for (unsigned j = 0; j < state_variable_sol.size(); j++)
		{
		    state_variable_ref.push_back(state_variable_sol[j]);
		}
	    }

	    /*
	     		 * == Compare Different Tol ==
			 */
            FileFinder tol_test_filename("projects/ChonL/plot_output", RelativeTo::ChasteSourceRoot);
            std::ofstream tol_test_file(tol_test_filename.GetAbsolutePath()+"solver_tol_test.txt");
	    tol_test_file << "# tol    Chaste_err    CiPA_err \n";
	    for (unsigned toli = 0; toli < 10; toli++)
	    {
                /* == Re-set Model State Variable ==
		             * Start with the state variables that CiPA use
			     *
			     */
                p_model->SetStateVariables(CiPA_stateVariables_Chaste_order);
                // Try different tolerances
                p_model->SetTolerances(1e-6*pow(1e-1,toli), 1e-6*pow(1e-1,toli)); 
	        p_model->SetMaxSteps(1e8); // Use default

		// Load pre-saved different tol solution from CiPA
                FileFinder CiPA_out_tol_filename("projects/ChonL/CiPA/fixed/control_" + std::to_string(toli) + "/Chaste_out.txt", RelativeTo::ChasteSourceRoot);
                std::ifstream CiPA_out_tol_file(CiPA_out_tol_filename.GetAbsolutePath());
                std::vector< std::vector<double> > test_CiPA_out_tmp(voltages.size());
                for (unsigned i = 0; i < voltages.size(); i++)
                {
                    test_CiPA_out_tmp[i].resize(YName.size());
                    for (unsigned j = 0; j < YName.size(); j++)
                    {
			CiPA_out_tol_file >> test_CiPA_out_tmp[i][j];
                    }
                }
                CiPA_out_tol_file.close();
                std::vector<double> test_CiPA_out;
                for (unsigned i = 0; i < YName.size(); i++)
                {
                    for (unsigned j = 0; j < voltages.size(); j++)
                    {
			CiPA_out_tol_file >> tmp_double;
                        test_CiPA_out.push_back(test_CiPA_out_tmp[j][map_Chaste_CiPA_var[i]]);
                    }
                }

		// Run again in Chaste
            	OdeSolution test_solution = p_model->Compute(start_time, end_time, sampling_timestep);
		std::vector<double> state_variable_test_tol;
		for (unsigned i = 0; i < YName.size(); i++)
		{
		    unsigned state_variable_index = p_model->GetSystemInformation()->GetStateVariableIndex(YName[i]);
		    std::vector<double> state_variable_sol = test_solution.GetVariableAtIndex(state_variable_index);
		    for (unsigned j = 0; j < state_variable_sol.size(); j++)
		    {
			state_variable_test_tol.push_back(state_variable_sol[j]);
		    }
		}

		double err_Chaste = 0;
		double err_CiPA = 0;
		for (unsigned i = 0; i < state_variable_ref.size(); i++)
		{
		    err_Chaste += abs(state_variable_ref[i] - state_variable_test_tol[i]);
		    err_CiPA += abs(state_variable_ref[i] - test_CiPA_out[i]);
		}

		std::cout << "At " << toli << ", error = (Chaste) " << err_Chaste << ", (CiPA) " << err_CiPA << "\n";
		tol_test_file << 1e-6*pow(1e-1,toli) << "    " << err_Chaste << "    " << err_CiPA << "\n";
		
	    } // Done tol test
	    tol_test_file.close();

        } // End of ====== Control Test ======

        std::cout << "End of control test" << std::endl;


#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

#endif /*TESTCIPAOHARAMODEL_HPP_*/
