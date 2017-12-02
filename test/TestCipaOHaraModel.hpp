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

class TestCipaOHaraModel : public CxxTest::TestSuite
{
public:
    void TestCipaVersionOfOharaModel()
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
		 * ====== Control test ======
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
            p_model->SetTolerances(1e-6, 1e-8);

            /*
			 * == Quick Check dy (RHS) At t=0 ==
			 */
            double getTime = 0;
            N_Vector Y = p_model->GetStateVariables();
            N_Vector dY = N_VNew_Serial(YName.size());
            p_model->EvaluateYDerivatives(getTime, Y, dY);
            // Read CiPA dy at t=0 results
            FileFinder CiPA_derivatives_from_R_filename("projects/ChonL/CiPA/fixed/control/Chaste_derivs_t0.txt", RelativeTo::ChasteSourceRoot);
            std::ifstream CiPA_dy_t0_file(CiPA_derivatives_from_R_filename.GetAbsolutePath());
            std::vector<double> CiPA_dy_t0;
            while (CiPA_dy_t0_file >> tmp_double)
            {
                CiPA_dy_t0.push_back(tmp_double);
            }
            CiPA_dy_t0_file.close();
            for (unsigned i = 0; i < YName.size(); i++)
            {
                TS_ASSERT_DELTA(NV_Ith_S(dY, i), CiPA_dy_t0[map_Chaste_CiPA_var[i]], 1e-10); // This can be quite accurate
            }
            DeleteVector(dY); // Tidy up memory

            // Print out current for checking
            /*
			Cellohara_rudy_2011_endo_dyHergFromCellMLCvode oharady_cvode_system(p_solver, p_stimulus);
			OdeSolution temp_sol = p_model->Compute(0, 1, 0.01);
			temp_sol.CalculateDerivedQuantitiesAndParameters(&oharady_cvode_system);
			std::vector<double> iNa = temp_sol.GetAnyVariable("membrane_fast_sodium_current");
			std::vector<double> iNaL = temp_sol.GetAnyVariable("membrane_persistent_sodium_current");
			std::vector<double> ito = temp_sol.GetAnyVariable("membrane_transient_outward_current");
			std::vector<double> iCaL = temp_sol.GetAnyVariable("membrane_L_type_calcium_current");
			std::vector<double> iKr = temp_sol.GetAnyVariable("membrane_rapid_delayed_rectifier_potassium_current");
			std::vector<double> iKs = temp_sol.GetAnyVariable("membrane_slow_delayed_rectifier_potassium_current");
			std::vector<double> iK1 = temp_sol.GetAnyVariable("membrane_inward_rectifier_potassium_current");
			std::cout << "iNa: " << iNa[0] << "\n";
			std::cout << "iNaL: " << iNaL[0] << "\n";
			std::cout << "ito: " << ito[0] << "\n";
			std::cout << "ICaL: " << iCaL[0] << "\n";
			std::cout << "iKr: " << iKr[0] << "\n";
			std::cout << "iKs: " << iKs[0] << "\n";
			std::cout << "iK1: " << iK1[0] << "\n";
			p_model->SetStateVariables(CiPA_stateVariables_Chaste_order);
			*/

            /* == Double check the state variable ==
			 * To make we loaded the state variables properly.
			 *
			 */
            for (unsigned i = 0u; i < CiPA_stateVariables.size(); i++)
            {
                TS_ASSERT_DELTA(CiPA_stateVariables_Chaste_order[i], p_model->GetStdVecStateVariables()[i], 1e-6);
            }

            /*
			 * == Getting detail for paces of interest ==
			 * Now we solve for the number of paces we are interested in.
			 *
			 */
            double max_timestep = 0.1;
            p_model->SetMaxTimestep(max_timestep);

            // Set sampling step size as the one CiPA use
            double sampling_timestep = 1.0;
            double start_time = 0.0;
            double end_time = 2000.0;
            OdeSolution solution = p_model->Compute(start_time, end_time, sampling_timestep);

            // Write the data out to a file.
            solution.WriteToFile("TestCipaOHaraModel", "OHaraDyHergCvode", "ms");

            /*
			 * == Calculating APD and Upstroke Velocity ==
			 */
            unsigned voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
            std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);
            CellProperties cell_props(voltages, solution.rGetTimes());

            double apd = cell_props.GetLastActionPotentialDuration(90);

            TS_ASSERT_DELTA(apd, 290.186, 1e-2);
            std::cout << "APD90: " << apd << "\n";

            /*
			 * == Compare Voltages with CiPA Output ==
			 */
            FileFinder CiPA_v_filename("projects/ChonL/CiPA/fixed/control/Chaste_v.txt", RelativeTo::ChasteSourceRoot);
            std::ifstream CiPA_v_file(CiPA_v_filename.GetAbsolutePath());
            std::vector<double> CiPA_voltages;
            while (CiPA_v_file >> tmp_double)
            {
                CiPA_voltages.push_back(tmp_double);
            }
            CiPA_v_file.close();
            for (unsigned i = 0; i < voltages.size(); i++)
            {
                TS_ASSERT_DELTA(voltages[i], CiPA_voltages[i], 1e-2);
            }

            /*
			 * == Compare All Other Variables ==
			 */
            FileFinder CiPA_out_filename("projects/ChonL/CiPA/fixed/control/Chaste_out.txt", RelativeTo::ChasteSourceRoot);
            std::ifstream CiPA_out_file(CiPA_out_filename.GetAbsolutePath());
            std::vector<std::vector<double> > CiPA_out(voltages.size());
            for (unsigned i = 0; i < voltages.size(); i++)
            {
                CiPA_out[i].resize(YName.size());
                for (unsigned j = 0; j < YName.size(); j++)
                {
                    CiPA_out_file >> CiPA_out[i][j];
                }
            }
            CiPA_out_file.close();
            for (unsigned i = 1; i < YName.size(); i++)
            {
                unsigned state_variable_index = p_model->GetSystemInformation()->GetStateVariableIndex(YName[i]);
                std::vector<double> state_variable_sol = solution.GetVariableAtIndex(state_variable_index);
                for (unsigned j = 0; j < state_variable_sol.size(); j++)
                {
                    TS_ASSERT_DELTA(state_variable_sol[j], CiPA_out[j][map_Chaste_CiPA_var[i]], 1e-4);
                }
            }

            /*
			 * == Check dy At All Times ==
			 * The most important check!! Compare the RHS.
			 *
			 */
            FileFinder CiPA_dy_filename("projects/ChonL/CiPA/fixed/control/Chaste_derivs.txt", RelativeTo::ChasteSourceRoot);
            std::ifstream CiPA_dy_file(CiPA_dy_filename.GetAbsolutePath());
            std::vector<std::vector<double> > CiPA_dy(voltages.size());
            for (unsigned i = 0; i < voltages.size(); i++)
            {
                CiPA_dy[i].resize(YName.size());
                for (unsigned j = 0; j < YName.size(); j++)
                {
                    CiPA_dy_file >> CiPA_dy[i][j];
                }
            }
            CiPA_dy_file.close();
            double computeTime = 0;
            N_Vector CiPA_y = N_VNew_Serial(YName.size());
            N_Vector Chaste_dy = N_VNew_Serial(YName.size());
            for (unsigned i = 0; i < voltages.size() - 1; i++)
            { // The last one it start pacing again (not in CiPA one)

                for (unsigned j = 0; j < YName.size(); j++)
                {
                    // CiPA_dy corresponds to the RHS output at the state of CiPA_out
                    NV_Ith_S(CiPA_y, j) = CiPA_out[i][map_Chaste_CiPA_var[j]];
                }
                p_model->EvaluateYDerivatives(computeTime, CiPA_y, Chaste_dy);
                for (unsigned j = 0; j < YName.size(); j++)
                {
                    TS_ASSERT_DELTA(NV_Ith_S(Chaste_dy, j), CiPA_dy[i][map_Chaste_CiPA_var[j]], 1e-10);
                }
                computeTime++;
            }
            DeleteVector(Y);
            DeleteVector(CiPA_y);
            DeleteVector(Chaste_dy);
        } // End of ====== Control Test ======

        std::cout << "End of control test" << std::endl;

        /*
		 * ====== Drug test ======
		 *
		 * Test the model using drug (dofetilide) setting.
		 *
		 */
        {
            /*
			 * == Set Drug Condition ==
			 * Set both conductance and dynamic hERG parameters
			 *
			 */
            double tmp_double;
            std::string dummyLine;
            // For drug
            FileFinder CiPA_drug_parameters_filename("projects/ChonL/CiPA/fixed/dofetilide/Chaste_drug_pars.txt", RelativeTo::ChasteSourceRoot);
            std::ifstream CiPA_drug_parameters_file(CiPA_drug_parameters_filename.GetAbsolutePath());
            std::getline(CiPA_drug_parameters_file, dummyLine);
            std::vector<double> CiPA_drug_parameters;
            while (CiPA_drug_parameters_file >> tmp_double)
            {
                CiPA_drug_parameters.push_back(tmp_double);
            }
            CiPA_drug_parameters_file.close();
            // For rescaling opt scaling factors back to default ones
            FileFinder CiPA_opt_parameters_filename("projects/ChonL/CiPA/Chaste_newordherg_opt_scaling.txt", RelativeTo::ChasteSourceRoot);
            std::ifstream CiPA_opt_parameters_file(CiPA_opt_parameters_filename.GetAbsolutePath());
            std::getline(CiPA_opt_parameters_file, dummyLine);
            std::vector<double> CiPA_opt_parameters;
            while (CiPA_opt_parameters_file >> tmp_double)
            {
                CiPA_opt_parameters.push_back(tmp_double);
            }
            CiPA_opt_parameters_file.close();

            // Set conductance
            // GKrfc
            tmp_double = p_model->GetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance");
            tmp_double = tmp_double * CiPA_opt_parameters[0] / CiPA_drug_parameters[0];
            p_model->SetParameter("membrane_rapid_delayed_rectifier_potassium_current_conductance", tmp_double);
            // GNaLfc
            tmp_double = p_model->GetParameter("membrane_persistent_sodium_current_conductance");
            tmp_double = tmp_double * CiPA_opt_parameters[1] / CiPA_drug_parameters[1];
            p_model->SetParameter("membrane_persistent_sodium_current_conductance", tmp_double);
            // GNafc
            tmp_double = p_model->GetParameter("membrane_fast_sodium_current_conductance");
            tmp_double = tmp_double * CiPA_opt_parameters[2] / CiPA_drug_parameters[2];
            p_model->SetParameter("membrane_fast_sodium_current_conductance", tmp_double);
            // GKsfc
            tmp_double = p_model->GetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance");
            tmp_double = tmp_double * CiPA_opt_parameters[3] / CiPA_drug_parameters[3];
            p_model->SetParameter("membrane_slow_delayed_rectifier_potassium_current_conductance", tmp_double);
            // GK1fc
            tmp_double = p_model->GetParameter("membrane_inward_rectifier_potassium_current_conductance");
            tmp_double = tmp_double * CiPA_opt_parameters[4] / CiPA_drug_parameters[4];
            p_model->SetParameter("membrane_inward_rectifier_potassium_current_conductance", tmp_double);
            // PCafc
            tmp_double = p_model->GetParameter("membrane_L_type_calcium_current_conductance");
            tmp_double = tmp_double * CiPA_opt_parameters[5] / CiPA_drug_parameters[5];
            p_model->SetParameter("membrane_L_type_calcium_current_conductance", tmp_double);
            // Gtofc
            tmp_double = p_model->GetParameter("membrane_transient_outward_current_conductance");
            tmp_double = tmp_double * CiPA_opt_parameters[6] / CiPA_drug_parameters[6];
            p_model->SetParameter("membrane_transient_outward_current_conductance", tmp_double);
            // Set dynamic hERG parameters
            // Kmax
            p_model->SetParameter("Dynamic_hERG_Kmax", CiPA_drug_parameters[7]);
            // Ku
            p_model->SetParameter("Dynamic_hERG_Ku", CiPA_drug_parameters[8]);
            // n
            p_model->SetParameter("Dynamic_hERG_n", CiPA_drug_parameters[9]);
            // halfmax
            p_model->SetParameter("Dynamic_hERG_halfmax", CiPA_drug_parameters[10]);
            // Kt
            p_model->SetParameter("Dynamic_hERG_Kt", CiPA_drug_parameters[11]);
            // Vhalf
            p_model->SetParameter("Dynamic_hERG_Vhalf", CiPA_drug_parameters[12]);

            /* == Read Model State Variable ==
			 * Try to set initial values to match CiPA one
			 *
			 */
            FileFinder CiPA_stateVariables_filename("projects/ChonL/CiPA/fixed/dofetilide/Chaste_newordherg_drug_states_CL2000.txt", RelativeTo::ChasteSourceRoot);
            std::ifstream CiPA_stateVariables_file(CiPA_stateVariables_filename.GetAbsolutePath());
            std::vector<double> CiPA_stateVariables;
            while (CiPA_stateVariables_file >> tmp_double)
            {
                CiPA_stateVariables.push_back(tmp_double);
            }
            CiPA_stateVariables_file.close();
            // Change var order
            std::vector<double> CiPA_stateVariables_Chaste_order;
            for (unsigned i = 0; i < CiPA_stateVariables.size(); i++)
            {
                CiPA_stateVariables_Chaste_order.push_back(CiPA_stateVariables[map_Chaste_CiPA_var[i]]);
            }
            /* == Re-set Model State Variable ==
			 * Start with the state variablesthat CiPA use
			 *
			 */
            p_model->SetStateVariables(CiPA_stateVariables_Chaste_order);
            // Try different tolerances
            p_model->SetTolerances(1e-6, 1e-8);

            /* == Double check the state variable ==
			 * To make we loaded the state variables properly.
			 *
			 */
            for (unsigned i = 0u; i < CiPA_stateVariables.size(); i++)
            {
                TS_ASSERT_DELTA(CiPA_stateVariables_Chaste_order[i], p_model->GetStdVecStateVariables()[i], 1e-6);
            }

            /* == Run Longer ==
			 * If interested in how it looks like for drug effect.
			 * Remember to comment out the tests below, should only pass for the first beat comparison.
			 */
            /*
			SteadyStateRunner steady_runner(p_model,true);
			steady_runner.SetMaxNumPaces(1000u);
			steady_runner.RunToSteadyState();
			*/

            /* == Getting detail for paces of interest ==
			 * Now we solve for the number of paces we are interested in.
			 *
			 */
            double max_timestep = 0.1;
            p_model->SetMaxTimestep(max_timestep);

            // Set sampling step size as the one CiPA use
            double sampling_timestep = 1.0;
            double start_time = 0.0;
            double end_time = 2000.0;
            OdeSolution solution = p_model->Compute(start_time, end_time, sampling_timestep);

            // Write the data out to a file.
            solution.WriteToFile("TestCipaOHaraModel", "OHaraDyHergCvode_drug", "ms", 1, false); //Not clean the folder

            /*
			 * == Calculating APD and Upstroke Velocity ==
			 */
            unsigned voltage_index = p_model->GetSystemInformation()->GetStateVariableIndex("membrane_voltage");
            std::vector<double> voltages = solution.GetVariableAtIndex(voltage_index);
            CellProperties cell_props(voltages, solution.rGetTimes());

            double apd = cell_props.GetLastActionPotentialDuration(90);

            //TS_ASSERT_DELTA(apd, , 1e-2);
            std::cout << "APD90_drug: " << apd << "\n"; // This is just the first beat!

            /*
			 * == Compare All Variables ==
			 */
            FileFinder CiPA_out_filename("projects/ChonL/CiPA/fixed/dofetilide/Chaste_out.txt", RelativeTo::ChasteSourceRoot);
            std::ifstream CiPA_out_file(CiPA_out_filename.GetAbsolutePath());
            std::vector<std::vector<double> > CiPA_out(voltages.size());
            for (unsigned i = 0; i < voltages.size(); i++)
            {
                CiPA_out[i].resize(YName.size());
                for (unsigned j = 0; j < YName.size(); j++)
                {
                    CiPA_out_file >> CiPA_out[i][j];
                }
            }
            CiPA_out_file.close();
            for (unsigned i = 1; i < YName.size(); i++)
            {
                unsigned state_variable_index = p_model->GetSystemInformation()->GetStateVariableIndex(YName[i]);
                std::vector<double> state_variable_sol = solution.GetVariableAtIndex(state_variable_index);
                for (unsigned j = 0; j < state_variable_sol.size(); j++)
                {
                    TS_ASSERT_DELTA(state_variable_sol[j], CiPA_out[j][map_Chaste_CiPA_var[i]], 1e-4);
                }
            }

            /*
			 * == Check dy At All Times ==
			 * The most important check!! Compare the RHS.
			 */
            FileFinder CiPA_dy_filename("projects/ChonL/CiPA/fixed/dofetilide/Chaste_derivs.txt", RelativeTo::ChasteSourceRoot);
            std::ifstream CiPA_dy_file(CiPA_dy_filename.GetAbsolutePath());
            std::vector<std::vector<double> > CiPA_dy(voltages.size());
            for (unsigned i = 0; i < voltages.size(); i++)
            {
                CiPA_dy[i].resize(YName.size());
                for (unsigned j = 0; j < YName.size(); j++)
                {
                    CiPA_dy_file >> CiPA_dy[i][j];
                }
            }
            CiPA_dy_file.close();
            double computeTime = 0;

            N_Vector CiPA_y = N_VNew_Serial(YName.size());
            N_Vector Chaste_dy = N_VNew_Serial(YName.size());

            for (unsigned i = 0; i < voltages.size() - 1; i++)
            { // The last one it start pacing again (not in CiPA one)
                for (unsigned j = 0; j < YName.size(); j++)
                {
                    // CiPA_dy corresponds to the RHS output at the state of CiPA_out
                    NV_Ith_S(CiPA_y, j) = CiPA_out[i][map_Chaste_CiPA_var[j]];
                }
                p_model->EvaluateYDerivatives(computeTime, CiPA_y, Chaste_dy);
                for (unsigned j = 0; j < YName.size(); j++)
                {
                    TS_ASSERT_DELTA(NV_Ith_S(Chaste_dy, j), CiPA_dy[i][map_Chaste_CiPA_var[j]], 1e-7);
                }
                computeTime++;
            }

            DeleteVector(CiPA_y);
            DeleteVector(Chaste_dy);

        } // End of ====== Drug test ======

#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }
};

#endif /*TESTCIPAOHARAMODEL_HPP_*/
