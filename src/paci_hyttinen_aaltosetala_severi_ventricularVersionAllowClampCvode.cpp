#ifdef CHASTE_CVODE
//! @file
//! 
//! This source file was generated from CellML and modified by ChonL on Wed May 05 2017
//! 
//! Model: paci_hyttinen_aaltosetala_severi_ventricularVersion
//! 


#include "paci_hyttinen_aaltosetala_severi_ventricularVersionAllowClampCvode.hpp"
#include <cmath>
#include <cassert>
#include <memory>
#include "Exception.hpp"
#include "OdeSystemInformation.hpp"
#include "RegularStimulus.hpp"
#include "HeartConfig.hpp"
#include "IsNan.hpp"
#include "MathsCustomFunctions.hpp"


    // 
    // allow clamp ion concentration
    // 

    void Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode::SetNaiDerivativeToZero(bool clamp)
    {
        // reset solver was in AbstractCvodeCell::SetVoltageDerivativeToZero(); not sure if needed
        if (clamp != mSetNaiDerivativeToZero)
        {
            //ResetSolver();
        }

        mSetNaiDerivativeToZero = clamp;
        if (clamp)
        {
            mFixedNai = GetAnyVariable("cytosolic_sodium_concentration");
        }
    }

    void Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode::SetCaiDerivativeToZero(bool clamp)
    {
        // reset solver was in AbstractCvodeCell::SetVoltageDerivativeToZero(); not sure if needed
        if (clamp != mSetCaiDerivativeToZero)
        {
            //ResetSolver();
        }

        mSetCaiDerivativeToZero = clamp;
        if (clamp)
        {
            mFixedCai = GetAnyVariable("cytosolic_calcium_concentration");
        }
    }






    boost::shared_ptr<RegularStimulus> Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode::UseCellMLDefaultStimulus()
    {
        // Use the default stimulus specified by CellML metadata
        const double var_chaste_interface__stim_mode__i_stim_Start = 0.0; // millisecond
        const double var_chaste_interface__stim_mode__i_stim_Period = 1000.0; // millisecond
        const double var_chaste_interface__stim_mode__i_stim_PulseDuration = 5.0; // millisecond
        const double var_chaste_interface__stim_mode__i_stim_Amplitude = (5.4999999999999996e-10 / NV_Ith_S(mParameters, 8)) * HeartConfig::Instance()->GetCapacitance(); // uA_per_cm2
        boost::shared_ptr<RegularStimulus> p_cellml_stim(new RegularStimulus(
                -fabs(var_chaste_interface__stim_mode__i_stim_Amplitude),
                var_chaste_interface__stim_mode__i_stim_PulseDuration,
                var_chaste_interface__stim_mode__i_stim_Period,
                var_chaste_interface__stim_mode__i_stim_Start
                ));
        mpIntracellularStimulus = p_cellml_stim;
        return p_cellml_stim;
    }
    
    double Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode::GetIntracellularCalciumConcentration()
    {
        return NV_Ith_S(mStateVariables, 15);
    }
    
    Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode::Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver /* unused; should be empty */, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
        : AbstractCvodeCell(
                pOdeSolver,
                18,
                0,
                pIntracellularStimulus)
    {
        // Time units: millisecond
        // 
        this->mpSystemInfo = OdeSystemInformation<Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode>::Instance();
        Init();

        // We have a default stimulus specified in the CellML file metadata
        this->mHasDefaultStimulusFromCellML = true;
        NV_Ith_S(this->mParameters, 0) = 150; // (c,model_parameters__Ki) [millimolar]
        NV_Ith_S(this->mParameters, 1) = 1.8; // (c,model_parameters__Cao) [millimolar]
        NV_Ith_S(this->mParameters, 2) = 5.4; // (c,model_parameters__Ko) [millimolar]
        NV_Ith_S(this->mParameters, 3) = 151; // (c,model_parameters__Nao) [millimolar]
        NV_Ith_S(this->mParameters, 4) = 8.635702e-5; // (c,i_CaL__g_CaL) [metre_cube_per_F_per_s]
        NV_Ith_S(this->mParameters, 5) = 0.69264; // (c,i_b_Ca__g_b_Ca) [S_per_F]
        NV_Ith_S(this->mParameters, 6) = 0.9; // (c,i_b_Na__g_b_Na) [S_per_F]
        NV_Ith_S(this->mParameters, 7) = 0.4125; // (c,i_PCa__g_PCa) [A_per_F]
        NV_Ith_S(this->mParameters, 8) = 9.87109e-11; // (c,model_parameters__Cm) [farad]
        NV_Ith_S(this->mParameters, 9) = 3671.2302; // (c,i_Na__g_Na) [S_per_F]
        NV_Ith_S(this->mParameters, 10) = 30.10312; // (c,i_f__g_f) [S_per_F]
        NV_Ith_S(this->mParameters, 11) = 28.1492; // (c,i_K1__g_K1) [S_per_F]
        NV_Ith_S(this->mParameters, 12) = 29.8667; // (c,i_Kr__g_Kr) [S_per_F]
        NV_Ith_S(this->mParameters, 13) = 2.041; // (c,i_Ks__g_Ks) [S_per_F]
        NV_Ith_S(this->mParameters, 14) = 4900; // (c,i_NaCa__kNaCa) [A_per_F]
        NV_Ith_S(this->mParameters, 15) = 1.841424; // (c,i_NaK__PNaK) [A_per_F]
        NV_Ith_S(this->mParameters, 16) = 29.9038; // (c,i_to__g_to) [S_per_F]
        NV_Ith_S(this->mParameters, 17) = 310; // (c,model_parameters__T) [kelvin]
        mUseAnalyticJacobian = false;
        mHasAnalyticJacobian = false;
        mSetNaiDerivativeToZero = false;
        mSetCaiDerivativeToZero = false;
    }
    
    Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode::~Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode()
    {
    }
    
    double Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode::GetIIonic(const std::vector<double>* pStateVariables)
    {
        // For state variable interpolation (SVI) we read in interpolated state variables,
        // otherwise for ionic current interpolation (ICI) we use the state variables of this model (node).
        N_Vector rY;
        bool made_new_cvode_vector = false;
        if (!pStateVariables)
        {
            rY = rGetStateVariables();
        }
        else
        {
            made_new_cvode_vector = true;
            rY = MakeNVector(*pStateVariables);
        }
        
        double var_chaste_interface__Membrane__Vm = (mSetVoltageDerivativeToZero ? this->mFixedVoltage : NV_Ith_S(rY, 0));
        // Units: millivolt; Initial value: -74.3340057624
        double var_chaste_interface__i_Na_m_gate__m = NV_Ith_S(rY, 1);
        // Units: dimensionless; Initial value: 0.102953468725004
        double var_chaste_interface__i_Na_h_gate__h = NV_Ith_S(rY, 2);
        // Units: dimensionless; Initial value: 0.786926637881461
        double var_chaste_interface__i_Na_j_gate__j = NV_Ith_S(rY, 3);
        // Units: dimensionless; Initial value: 0.253943221774722
        double var_chaste_interface__i_CaL_d_gate__d = NV_Ith_S(rY, 4);
        // Units: dimensionless; Initial value: 8.96088425225182e-5
        double var_chaste_interface__i_CaL_f1_gate__f1 = NV_Ith_S(rY, 5);
        // Units: dimensionless; Initial value: 0.970411811263976
        double var_chaste_interface__i_CaL_f2_gate__f2 = NV_Ith_S(rY, 6);
        // Units: dimensionless; Initial value: 0.999965815466749
        double var_chaste_interface__i_CaL_fCa_gate__fCa = NV_Ith_S(rY, 7);
        // Units: dimensionless; Initial value: 0.998925296531804
        double var_chaste_interface__i_Kr_Xr1_gate__Xr1 = NV_Ith_S(rY, 8);
        // Units: dimensionless; Initial value: 0.00778547011240132
        double var_chaste_interface__i_Kr_Xr2_gate__Xr2 = NV_Ith_S(rY, 9);
        // Units: dimensionless; Initial value: 0.432162576531617
        double var_chaste_interface__i_Ks_Xs_gate__Xs = NV_Ith_S(rY, 10);
        // Units: dimensionless; Initial value: 0.0322944866983666
        double var_chaste_interface__i_f_Xf_gate__Xf = NV_Ith_S(rY, 11);
        // Units: dimensionless; Initial value: 0.100615100568753
        double var_chaste_interface__i_to_q_gate__q = NV_Ith_S(rY, 12);
        // Units: dimensionless; Initial value: 0.839295925773219
        double var_chaste_interface__i_to_r_gate__r = NV_Ith_S(rY, 13);
        // Units: dimensionless; Initial value: 0.00573289893326379
        //double var_chaste_interface__sodium_dynamics__Nai = NV_Ith_S(rY, 14);
        double var_chaste_interface__sodium_dynamics__Nai = (mSetNaiDerivativeToZero ? this->mFixedNai : NV_Ith_S(rY, 14));
        // Units: millimolar; Initial value: 10.9248496211574
        //double var_chaste_interface__calcium_dynamics__Cai = NV_Ith_S(rY, 15);
        double var_chaste_interface__calcium_dynamics__Cai = (mSetCaiDerivativeToZero ? this->mFixedCai : NV_Ith_S(rY, 15));
        // Units: millimolar; Initial value: 1.80773974140477e-5
        
        const double var_i_CaL__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_CaL__i_CaL = ((((NV_Ith_S(mParameters, 4) * (4.0 * (var_i_CaL__Vm * 9309421124.3716221))) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * ((var_chaste_interface__calcium_dynamics__Cai * exp((2.0 * (var_i_CaL__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17)))) - (0.34100000000000003 * NV_Ith_S(mParameters, 1)))) / (exp((2.0 * (var_i_CaL__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) - 1.0)) * (var_chaste_interface__i_CaL_d_gate__d * (var_chaste_interface__i_CaL_f1_gate__f1 * (var_chaste_interface__i_CaL_f2_gate__f2 * var_chaste_interface__i_CaL_fCa_gate__fCa))); // A_per_F
        const double var_i_K1__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_electric_potentials__E_K = ((8.3144720000000003 * NV_Ith_S(mParameters, 17)) * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 2) / NV_Ith_S(mParameters, 0)); // volt
        const double var_i_K1__alpha_K1 = 3.9100000000000001 / (1.0 + exp(0.59419999999999995 * (((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0)) - 200.0))); // dimensionless
        const double var_i_K1__i_K1 = NV_Ith_S(mParameters, 11) * ((var_i_K1__alpha_K1 / (var_i_K1__alpha_K1 + ((( -1.5089999999999999 * exp(0.00020000000000000001 * (((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0)) + 100.0))) + exp(0.58860000000000001 * (((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0)) - 10.0))) / (1.0 + exp(0.45469999999999999 * ((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0))))))) * ((var_i_K1__Vm - var_electric_potentials__E_K) * sqrt(NV_Ith_S(mParameters, 2) * 0.18518518518518517))); // A_per_F
        const double var_i_f__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_f__i_f = NV_Ith_S(mParameters, 10) * (var_chaste_interface__i_f_Xf_gate__Xf * (var_i_f__Vm -  -0.017000000000000001)); // A_per_F
        const double var_electric_potentials__E_Na = ((8.3144720000000003 * NV_Ith_S(mParameters, 17)) * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 3) / var_chaste_interface__sodium_dynamics__Nai); // volt
        const double var_i_Na__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_Na__i_Na = 1.0 * (NV_Ith_S(mParameters, 9) * (pow(var_chaste_interface__i_Na_m_gate__m, 3.0) * (var_chaste_interface__i_Na_h_gate__h * (var_chaste_interface__i_Na_j_gate__j * (var_i_Na__Vm - var_electric_potentials__E_Na))))); // A_per_F
        const double var_i_Kr__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_Kr__i_Kr = 1.0 * (NV_Ith_S(mParameters, 12) * ((var_i_Kr__Vm - var_electric_potentials__E_K) * (var_chaste_interface__i_Kr_Xr1_gate__Xr1 * (var_chaste_interface__i_Kr_Xr2_gate__Xr2 * sqrt(NV_Ith_S(mParameters, 2) * 0.18518518518518517))))); // A_per_F
        const double var_i_Ks__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_Ks__i_Ks = 1.0 * (NV_Ith_S(mParameters, 13) * ((var_i_Ks__Vm - (((8.3144720000000003 * NV_Ith_S(mParameters, 17)) * 1.0364268649036186e-05) * log((NV_Ith_S(mParameters, 2) + (0.029999999999999999 * NV_Ith_S(mParameters, 3))) / (NV_Ith_S(mParameters, 0) + (0.029999999999999999 * var_chaste_interface__sodium_dynamics__Nai))))) * (pow(var_chaste_interface__i_Ks_Xs_gate__Xs, 2.0) * (1.0 + (0.59999999999999998 / (1.0 + pow(3.8000000000000002e-05 / var_chaste_interface__calcium_dynamics__Cai, 1.3999999999999999))))))); // A_per_F
        const double var_i_to__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_to__i_to = NV_Ith_S(mParameters, 16) * ((var_i_to__Vm - var_electric_potentials__E_K) * (var_chaste_interface__i_to_q_gate__q * var_chaste_interface__i_to_r_gate__r)); // A_per_F
        const double var_i_PCa__i_PCa = (NV_Ith_S(mParameters, 7) * var_chaste_interface__calcium_dynamics__Cai) / (var_chaste_interface__calcium_dynamics__Cai + 0.00050000000000000001); // A_per_F
        const double var_i_NaK__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_NaK__i_NaK = ((((NV_Ith_S(mParameters, 15) * NV_Ith_S(mParameters, 2)) / (NV_Ith_S(mParameters, 2) + 1.0)) * var_chaste_interface__sodium_dynamics__Nai) / (var_chaste_interface__sodium_dynamics__Nai + 40.0)) / (1.0 + ((0.1245 * exp(( -0.10000000000000001 * (var_i_NaK__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17)))) + (0.035299999999999998 * exp(((-var_i_NaK__Vm) * 96485.341499999995) / (8.3144720000000003 * NV_Ith_S(mParameters, 17)))))); // A_per_F
        const double var_i_NaCa__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_NaCa__i_NaCa = (NV_Ith_S(mParameters, 14) * ((exp((0.34999999999999998 * (var_i_NaCa__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * (pow(var_chaste_interface__sodium_dynamics__Nai, 3.0) * NV_Ith_S(mParameters, 1))) - (exp(( -0.65000000000000002 * (var_i_NaCa__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * (pow(NV_Ith_S(mParameters, 3), 3.0) * (var_chaste_interface__calcium_dynamics__Cai * 2.8571431999999999))))) / ((669921.875 + pow(NV_Ith_S(mParameters, 3), 3.0)) * ((1.3799999999999999 + NV_Ith_S(mParameters, 1)) * (1.0 + (0.10000000000000001 * exp(( -0.65000000000000002 * (var_i_NaCa__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))))))); // A_per_F
        const double var_i_b_Ca__i_b_Ca = NV_Ith_S(mParameters, 5) * ((0.001 * var_chaste_interface__Membrane__Vm) - (((0.5 * (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 1) / var_chaste_interface__calcium_dynamics__Cai))); // A_per_F
        const double var_i_b_Na__i_b_Na = NV_Ith_S(mParameters, 6) * ((0.001 * var_chaste_interface__Membrane__Vm) - var_electric_potentials__E_Na); // A_per_F
        const double var_chaste_interface__i_NaK__i_NaK = var_i_NaK__i_NaK; // A_per_F
        const double var_chaste_interface__i_K1__i_K1 = var_i_K1__i_K1; // A_per_F
        const double var_chaste_interface__i_f__i_f = var_i_f__i_f; // A_per_F
        const double var_chaste_interface__i_Na__i_Na = var_i_Na__i_Na; // A_per_F
        const double var_chaste_interface__i_to__i_to = var_i_to__i_to; // A_per_F
        const double var_chaste_interface__i_Kr__i_Kr = var_i_Kr__i_Kr; // A_per_F
        const double var_chaste_interface__i_Ks__i_Ks = var_i_Ks__i_Ks; // A_per_F
        const double var_chaste_interface__i_CaL__i_CaL = var_i_CaL__i_CaL; // A_per_F
        const double var_chaste_interface__i_PCa__i_PCa = var_i_PCa__i_PCa; // A_per_F
        const double var_chaste_interface__i_NaCa__i_NaCa = var_i_NaCa__i_NaCa; // A_per_F
        const double var_chaste_interface__i_ionic = (var_chaste_interface__i_K1__i_K1 + var_chaste_interface__i_to__i_to + var_chaste_interface__i_Kr__i_Kr + var_chaste_interface__i_Ks__i_Ks + var_chaste_interface__i_CaL__i_CaL + var_chaste_interface__i_NaK__i_NaK + var_chaste_interface__i_Na__i_Na + var_chaste_interface__i_NaCa__i_NaCa + var_chaste_interface__i_PCa__i_PCa + var_chaste_interface__i_f__i_f + var_i_b_Na__i_b_Na + var_i_b_Ca__i_b_Ca) * HeartConfig::Instance()->GetCapacitance(); // uA_per_cm2
        
        const double i_ionic = var_chaste_interface__i_ionic;
        if (made_new_cvode_vector)
        {
            DeleteVector(rY);
        }
        EXCEPT_IF_NOT(!std::isnan(i_ionic));
        return i_ionic;
    }
    
    void Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode::EvaluateYDerivatives(double var_chaste_interface__environment__time, const N_Vector rY, N_Vector rDY)
    {
        // Inputs:
        // Time units: millisecond
        double var_chaste_interface__Membrane__Vm = (mSetVoltageDerivativeToZero ? this->mFixedVoltage : NV_Ith_S(rY, 0));
        // Units: millivolt; Initial value: -74.3340057624
        double var_chaste_interface__i_Na_m_gate__m = NV_Ith_S(rY, 1);
        // Units: dimensionless; Initial value: 0.102953468725004
        double var_chaste_interface__i_Na_h_gate__h = NV_Ith_S(rY, 2);
        // Units: dimensionless; Initial value: 0.786926637881461
        double var_chaste_interface__i_Na_j_gate__j = NV_Ith_S(rY, 3);
        // Units: dimensionless; Initial value: 0.253943221774722
        double var_chaste_interface__i_CaL_d_gate__d = NV_Ith_S(rY, 4);
        // Units: dimensionless; Initial value: 8.96088425225182e-5
        double var_chaste_interface__i_CaL_f1_gate__f1 = NV_Ith_S(rY, 5);
        // Units: dimensionless; Initial value: 0.970411811263976
        double var_chaste_interface__i_CaL_f2_gate__f2 = NV_Ith_S(rY, 6);
        // Units: dimensionless; Initial value: 0.999965815466749
        double var_chaste_interface__i_CaL_fCa_gate__fCa = NV_Ith_S(rY, 7);
        // Units: dimensionless; Initial value: 0.998925296531804
        double var_chaste_interface__i_Kr_Xr1_gate__Xr1 = NV_Ith_S(rY, 8);
        // Units: dimensionless; Initial value: 0.00778547011240132
        double var_chaste_interface__i_Kr_Xr2_gate__Xr2 = NV_Ith_S(rY, 9);
        // Units: dimensionless; Initial value: 0.432162576531617
        double var_chaste_interface__i_Ks_Xs_gate__Xs = NV_Ith_S(rY, 10);
        // Units: dimensionless; Initial value: 0.0322944866983666
        double var_chaste_interface__i_f_Xf_gate__Xf = NV_Ith_S(rY, 11);
        // Units: dimensionless; Initial value: 0.100615100568753
        double var_chaste_interface__i_to_q_gate__q = NV_Ith_S(rY, 12);
        // Units: dimensionless; Initial value: 0.839295925773219
        double var_chaste_interface__i_to_r_gate__r = NV_Ith_S(rY, 13);
        // Units: dimensionless; Initial value: 0.00573289893326379
        //double var_chaste_interface__sodium_dynamics__Nai = NV_Ith_S(rY, 14);
        double var_chaste_interface__sodium_dynamics__Nai = (mSetNaiDerivativeToZero ? this->mFixedNai : NV_Ith_S(rY, 14));
        // Units: millimolar; Initial value: 10.9248496211574
        //double var_chaste_interface__calcium_dynamics__Cai = NV_Ith_S(rY, 15);
        double var_chaste_interface__calcium_dynamics__Cai = (mSetCaiDerivativeToZero ? this->mFixedCai : NV_Ith_S(rY, 15));
        // Units: millimolar; Initial value: 1.80773974140477e-5
        double var_chaste_interface__calcium_dynamics__Ca_SR = NV_Ith_S(rY, 16);
        // Units: millimolar; Initial value: -0.2734234751931
        double var_chaste_interface__calcium_dynamics__g = NV_Ith_S(rY, 17);
        // Units: dimensionless; Initial value: 0.999999981028517
        
        
        // Mathematics
        double d_dt_chaste_interface__Membrane__Vm;
        const double var_i_CaL__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_CaL__i_CaL = ((((NV_Ith_S(mParameters, 4) * (4.0 * (var_i_CaL__Vm * 9309421124.3716221))) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * ((var_chaste_interface__calcium_dynamics__Cai * exp((2.0 * (var_i_CaL__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17)))) - (0.34100000000000003 * NV_Ith_S(mParameters, 1)))) / (exp((2.0 * (var_i_CaL__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) - 1.0)) * (var_chaste_interface__i_CaL_d_gate__d * (var_chaste_interface__i_CaL_f1_gate__f1 * (var_chaste_interface__i_CaL_f2_gate__f2 * var_chaste_interface__i_CaL_fCa_gate__fCa))); // A_per_F
        const double var_i_f__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_electric_potentials__E_Na = ((8.3144720000000003 * NV_Ith_S(mParameters, 17)) * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 3) / var_chaste_interface__sodium_dynamics__Nai); // volt
        const double var_i_Na__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_Na__i_Na = 1.0 * (NV_Ith_S(mParameters, 9) * (pow(var_chaste_interface__i_Na_m_gate__m, 3.0) * (var_chaste_interface__i_Na_h_gate__h * (var_chaste_interface__i_Na_j_gate__j * (var_i_Na__Vm - var_electric_potentials__E_Na))))); // A_per_F
        const double var_i_Kr__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_Ks__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_to__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_PCa__i_PCa = (NV_Ith_S(mParameters, 7) * var_chaste_interface__calcium_dynamics__Cai) / (var_chaste_interface__calcium_dynamics__Cai + 0.00050000000000000001); // A_per_F
        const double var_i_NaK__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_NaK__i_NaK = ((((NV_Ith_S(mParameters, 15) * NV_Ith_S(mParameters, 2)) / (NV_Ith_S(mParameters, 2) + 1.0)) * var_chaste_interface__sodium_dynamics__Nai) / (var_chaste_interface__sodium_dynamics__Nai + 40.0)) / (1.0 + ((0.1245 * exp(( -0.10000000000000001 * (var_i_NaK__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17)))) + (0.035299999999999998 * exp(((-var_i_NaK__Vm) * 96485.341499999995) / (8.3144720000000003 * NV_Ith_S(mParameters, 17)))))); // A_per_F
        const double var_i_NaCa__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_NaCa__i_NaCa = (NV_Ith_S(mParameters, 14) * ((exp((0.34999999999999998 * (var_i_NaCa__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * (pow(var_chaste_interface__sodium_dynamics__Nai, 3.0) * NV_Ith_S(mParameters, 1))) - (exp(( -0.65000000000000002 * (var_i_NaCa__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * (pow(NV_Ith_S(mParameters, 3), 3.0) * (var_chaste_interface__calcium_dynamics__Cai * 2.8571431999999999))))) / ((669921.875 + pow(NV_Ith_S(mParameters, 3), 3.0)) * ((1.3799999999999999 + NV_Ith_S(mParameters, 1)) * (1.0 + (0.10000000000000001 * exp(( -0.65000000000000002 * (var_i_NaCa__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))))))); // A_per_F
        const double var_i_b_Ca__i_b_Ca = NV_Ith_S(mParameters, 5) * ((0.001 * var_chaste_interface__Membrane__Vm) - (((0.5 * (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 1) / var_chaste_interface__calcium_dynamics__Cai))); // A_per_F
        const double var_i_b_Na__i_b_Na = NV_Ith_S(mParameters, 6) * ((0.001 * var_chaste_interface__Membrane__Vm) - var_electric_potentials__E_Na); // A_per_F
        const double var_i_CaL_f1_gate__f1_inf = 1.0 / (1.0 + exp(((var_i_CaL__Vm * 1000.0) + 26.0) * 0.33333333333333331)); // dimensionless
        const double var_i_CaL_fCa_gate__fCa_inf = ((1.0 / (1.0 + pow(var_chaste_interface__calcium_dynamics__Cai * 1666.6666666666667, 8.0))) + ((0.10000000000000001 / (1.0 + exp((var_chaste_interface__calcium_dynamics__Cai - 0.00089999999999999998) * 10000.0))) + (0.29999999999999999 / (1.0 + exp((var_chaste_interface__calcium_dynamics__Cai - 0.00075000000000000002) * 1250.0))))) * 0.76010945576162958; // dimensionless
        const double var_calcium_dynamics__g_inf = (var_chaste_interface__calcium_dynamics__Cai <= 0.00035) ? (1.0 / (1.0 + pow(var_chaste_interface__calcium_dynamics__Cai * 2857.1428571428573, 6.0))) : (1.0 / (1.0 + pow(var_chaste_interface__calcium_dynamics__Cai * 2857.1428571428573, 16.0))); // dimensionless
        const double var_calcium_dynamics__i_rel = (8.2319999999999993 + ((16.463999999999999 * pow(var_chaste_interface__calcium_dynamics__Ca_SR, 2.0)) / (0.0625 + pow(var_chaste_interface__calcium_dynamics__Ca_SR, 2.0)))) * (var_chaste_interface__i_CaL_d_gate__d * (var_chaste_interface__calcium_dynamics__g * 0.041099999999999998)); // millimolar_per_second
        const double var_calcium_dynamics__i_up = 0.56064000000000003 / (1.0 + (6.2499999999999997e-08 / pow(var_chaste_interface__calcium_dynamics__Cai, 2.0))); // millimolar_per_second
        const double var_calcium_dynamics__i_leak = (var_chaste_interface__calcium_dynamics__Ca_SR - var_chaste_interface__calcium_dynamics__Cai) * 0.00044443999999999999; // millimolar_per_second
        const double d_dt_chaste_interface__i_Na_m_gate__m = 0.001 * (((1.0 / pow(1.0 + exp((((-var_i_Na__Vm) * 1000.0) - 34.100000000000001) * 0.16949152542372881), 0.33333333333333331)) - var_chaste_interface__i_Na_m_gate__m) / ((1.0 * ((1.0 / (1.0 + exp((((-var_i_Na__Vm) * 1000.0) - 60.0) * 0.20000000000000001))) * ((0.10000000000000001 / (1.0 + exp(((var_i_Na__Vm * 1000.0) + 35.0) * 0.20000000000000001))) + (0.10000000000000001 / (1.0 + exp(((var_i_Na__Vm * 1000.0) - 50.0) * 0.0050000000000000001)))))) * 0.001)); // 'per millisecond'
        const double d_dt_chaste_interface__i_Na_h_gate__h = 0.001 * (((1.0 / sqrt(1.0 + exp(((var_i_Na__Vm * 1000.0) + 72.099999999999994) * 0.17543859649122806))) - var_chaste_interface__i_Na_h_gate__h) / ((var_i_Na__Vm <  -0.040000000000000001) ? (1.5 / ((((var_i_Na__Vm <  -0.040000000000000001) ? (0.057000000000000002 * exp((-((var_i_Na__Vm * 1000.0) + 80.0)) * 0.14705882352941177)) : 0.0) + ((var_i_Na__Vm <  -0.040000000000000001) ? ((2.7000000000000002 * exp(0.079000000000000001 * (var_i_Na__Vm * 1000.0))) + (3.1000000000000001 * (100000.0 * exp(0.34849999999999998 * (var_i_Na__Vm * 1000.0))))) : (0.77000000000000002 / (0.13 * (1.0 + exp(((var_i_Na__Vm * 1000.0) + 10.66) *  -0.0900900900900901)))))) * 1000.0)) : 0.002542)); // 'per millisecond'
        const double d_dt_chaste_interface__i_Na_j_gate__j = 0.001 * (((1.0 / sqrt(1.0 + exp(((var_i_Na__Vm * 1000.0) + 72.099999999999994) * 0.17543859649122806))) - var_chaste_interface__i_Na_j_gate__j) / (7.0 / ((((var_i_Na__Vm <  -0.040000000000000001) ? (((( -25428.0 * exp(0.24440000000000001 * (var_i_Na__Vm * 1000.0))) - (6.9480000000000004 * (9.9999999999999995e-07 * exp( -0.043909999999999998 * (var_i_Na__Vm * 1000.0))))) * ((var_i_Na__Vm * 1000.0) + 37.780000000000001)) / (1.0 + exp(0.311 * ((var_i_Na__Vm * 1000.0) + 79.230000000000004)))) : 0.0) + ((var_i_Na__Vm <  -0.040000000000000001) ? ((0.024240000000000001 * exp( -0.01052 * (var_i_Na__Vm * 1000.0))) / (1.0 + exp( -0.13780000000000001 * ((var_i_Na__Vm * 1000.0) + 40.140000000000001)))) : ((0.59999999999999998 * exp(0.057000000000000002 * (var_i_Na__Vm * 1000.0))) / (1.0 + exp( -0.10000000000000001 * ((var_i_Na__Vm * 1000.0) + 32.0)))))) * 1000.0))); // 'per millisecond'
        const double d_dt_chaste_interface__i_CaL_d_gate__d = 0.001 * (((1.0 / (1.0 + exp((-((var_i_CaL__Vm * 1000.0) + 9.0999999999999996)) * 0.14285714285714285))) - var_chaste_interface__i_CaL_d_gate__d) / (((((0.25 + (1.3999999999999999 / (1.0 + exp((((-var_i_CaL__Vm) * 1000.0) - 35.0) * 0.076923076923076927)))) * (1.3999999999999999 / (1.0 + exp(((var_i_CaL__Vm * 1000.0) + 5.0) * 0.20000000000000001)))) + (1.0 / (1.0 + exp((((-var_i_CaL__Vm) * 1000.0) + 50.0) * 0.050000000000000003)))) * 1.0) * 0.001)); // 'per millisecond'
        const double d_dt_chaste_interface__i_CaL_f1_gate__f1 = 0.001 * ((var_i_CaL_f1_gate__f1_inf - var_chaste_interface__i_CaL_f1_gate__f1) / (((20.0 + ((1102.5 * exp(-pow(pow((var_i_CaL__Vm * 1000.0) + 27.0, 2.0) * 0.066666666666666666, 2.0))) + ((200.0 / (1.0 + exp((13.0 - (var_i_CaL__Vm * 1000.0)) * 0.10000000000000001))) + (180.0 / (1.0 + exp((30.0 + (var_i_CaL__Vm * 1000.0)) * 0.10000000000000001)))))) * (((var_i_CaL_f1_gate__f1_inf - var_chaste_interface__i_CaL_f1_gate__f1) > 0.0) ? (1.0 + (1433.0 * (var_chaste_interface__calcium_dynamics__Cai - 4.9999999999999996e-05))) : 1.0)) * 0.001)); // 'per millisecond'
        const double d_dt_chaste_interface__i_CaL_f2_gate__f2 = 0.001 * (((0.33000000000000002 + (0.67000000000000004 / (1.0 + exp(((var_i_CaL__Vm * 1000.0) + 35.0) * 0.25)))) - var_chaste_interface__i_CaL_f2_gate__f2) / ((((600.0 * exp((-pow((var_i_CaL__Vm * 1000.0) + 25.0, 2.0)) * 0.0058823529411764705)) + ((31.0 / (1.0 + exp((25.0 - (var_i_CaL__Vm * 1000.0)) * 0.10000000000000001))) + (16.0 / (1.0 + exp((30.0 + (var_i_CaL__Vm * 1000.0)) * 0.10000000000000001))))) * 1.0) * 0.001)); // 'per millisecond'
        const double d_dt_chaste_interface__i_CaL_fCa_gate__fCa = 0.001 * (((((var_i_CaL__Vm >  -0.059999999999999998) && (var_i_CaL_fCa_gate__fCa_inf > var_chaste_interface__i_CaL_fCa_gate__fCa)) ? 0.0 : 1.0) * (var_i_CaL_fCa_gate__fCa_inf - var_chaste_interface__i_CaL_fCa_gate__fCa)) * 500.0); // 'per millisecond'
        const double d_dt_chaste_interface__i_Kr_Xr1_gate__Xr1 = 0.001 * (((1.0 / (1.0 + exp(((1000.0 * (((( -8.3144720000000003 * NV_Ith_S(mParameters, 17)) * 4.5062037604505159e-06) * log(pow(1.0 + (NV_Ith_S(mParameters, 1) * 0.38461538461538458), 4.0) / (0.025000000000000001 * pow(1.0 + (NV_Ith_S(mParameters, 1) * 1.7241379310344829), 4.0)))) - 0.019)) - (var_i_Kr__Vm * 1000.0)) * 0.2040816326530612))) - var_chaste_interface__i_Kr_Xr1_gate__Xr1) / ((1.0 * ((450.0 / (1.0 + exp(( -45.0 - (var_i_Kr__Vm * 1000.0)) * 0.10000000000000001))) * (6.0 / (1.0 + exp((30.0 + (var_i_Kr__Vm * 1000.0)) * 0.086956521739130432))))) * 0.001)); // 'per millisecond'
        const double d_dt_chaste_interface__i_Kr_Xr2_gate__Xr2 = 0.001 * (((1.0 / (1.0 + exp(((var_i_Kr__Vm * 1000.0) + 88.0) * 0.02))) - var_chaste_interface__i_Kr_Xr2_gate__Xr2) / ((1.0 * ((3.0 / (1.0 + exp(( -60.0 - (var_i_Kr__Vm * 1000.0)) * 0.050000000000000003))) * (1.1200000000000001 / (1.0 + exp(( -60.0 + (var_i_Kr__Vm * 1000.0)) * 0.050000000000000003))))) * 0.001)); // 'per millisecond'
        const double d_dt_chaste_interface__i_Ks_Xs_gate__Xs = 0.001 * (((1.0 / (1.0 + exp((((-var_i_Ks__Vm) * 1000.0) - 20.0) * 0.0625))) - var_chaste_interface__i_Ks_Xs_gate__Xs) / ((1.0 * ((1100.0 / sqrt(1.0 + exp(( -10.0 - (var_i_Ks__Vm * 1000.0)) * 0.16666666666666666))) * (1.0 / (1.0 + exp(( -60.0 + (var_i_Ks__Vm * 1000.0)) * 0.050000000000000003))))) * 0.001)); // 'per millisecond'
        const double d_dt_chaste_interface__i_f_Xf_gate__Xf = 0.001 * (((1.0 / (1.0 + exp(((var_i_f__Vm * 1000.0) + 77.849999999999994) * 0.20000000000000001))) - var_chaste_interface__i_f_Xf_gate__Xf) / ((1900.0 / (1.0 + exp(((var_i_f__Vm * 1000.0) + 15.0) * 0.10000000000000001))) * 0.001)); // 'per millisecond'
        const double d_dt_chaste_interface__i_to_q_gate__q = 0.001 * (((1.0 / (1.0 + exp(((var_i_to__Vm * 1000.0) + 53.0) * 0.076923076923076927))) - var_chaste_interface__i_to_q_gate__q) / ((6.0599999999999996 + (39.101999999999997 / ((0.56999999999999995 * exp( -0.080000000000000002 * ((var_i_to__Vm * 1000.0) + 44.0))) + (0.065000000000000002 * exp(0.10000000000000001 * ((var_i_to__Vm * 1000.0) + 45.93)))))) * 0.001)); // 'per millisecond'
        const double d_dt_chaste_interface__i_to_r_gate__r = 0.001 * (((1.0 / (1.0 + exp((-((var_i_to__Vm * 1000.0) - 22.300000000000001)) * 0.053333333333333337))) - var_chaste_interface__i_to_r_gate__r) / ((2.75352 + (14.40516 / ((1.0369999999999999 * exp(0.089999999999999997 * ((var_i_to__Vm * 1000.0) + 30.609999999999999))) + (0.36899999999999999 * exp( -0.12 * ((var_i_to__Vm * 1000.0) + 23.84)))))) * 0.001)); // 'per millisecond'
        //const double d_dt_chaste_interface__sodium_dynamics__Nai = 0.001 * (((-NV_Ith_S(mParameters, 8)) * (var_i_Na__i_Na + (var_i_b_Na__i_b_Na + ((3.0 * var_i_NaK__i_NaK) + (3.0 * var_i_NaCa__i_NaCa))))) * 1177757801.0268393); // 'millimole per litre per millisecond'

        double d_dt_chaste_interface__sodium_dynamics__Nai;
        if (mSetNaiDerivativeToZero)
        {
            d_dt_chaste_interface__sodium_dynamics__Nai = 0.0;
        }
        else
        {
            d_dt_chaste_interface__sodium_dynamics__Nai = 0.001 * (((-NV_Ith_S(mParameters, 8)) * (var_i_Na__i_Na + (var_i_b_Na__i_b_Na + ((3.0 * var_i_NaK__i_NaK) + (3.0 * var_i_NaCa__i_NaCa))))) * 1177757801.0268393); // 'millimole per litre per millisecond'
        }

        //const double d_dt_chaste_interface__calcium_dynamics__Cai = 0.001 * ((1.0 / (1.0 + (0.00025000000000000001 / pow(var_chaste_interface__calcium_dynamics__Cai + 0.001, 2.0)))) * (((var_calcium_dynamics__i_leak - var_calcium_dynamics__i_up) + var_calcium_dynamics__i_rel) - ((((var_i_CaL__i_CaL + (var_i_b_Ca__i_b_Ca + var_i_PCa__i_PCa)) - (2.0 * var_i_NaCa__i_NaCa)) * NV_Ith_S(mParameters, 8)) * 588878900.51341963))); // 'millimole per litre per millisecond'

        double d_dt_chaste_interface__calcium_dynamics__Cai;
        if (mSetCaiDerivativeToZero)
        {
            d_dt_chaste_interface__calcium_dynamics__Cai = 0.0;
        }
        else
        {
            d_dt_chaste_interface__calcium_dynamics__Cai = 0.001 * ((1.0 / (1.0 + (0.00025000000000000001 / pow(var_chaste_interface__calcium_dynamics__Cai + 0.001, 2.0)))) * (((var_calcium_dynamics__i_leak - var_calcium_dynamics__i_up) + var_calcium_dynamics__i_rel) - ((((var_i_CaL__i_CaL + (var_i_b_Ca__i_b_Ca + var_i_PCa__i_PCa)) - (2.0 * var_i_NaCa__i_NaCa)) * NV_Ith_S(mParameters, 8)) * 588878900.51341963))); // 'millimole per litre per millisecond'
        }

        const double d_dt_chaste_interface__calcium_dynamics__Ca_SR = 0.001 * ((((1.0 / (1.0 + (3.0 / pow(var_chaste_interface__calcium_dynamics__Ca_SR + 0.29999999999999999, 2.0)))) * 8800.0) * 0.0017131207921470542) * (var_calcium_dynamics__i_up - (var_calcium_dynamics__i_rel + var_calcium_dynamics__i_leak))); // 'millimole per litre per millisecond'
        const double d_dt_chaste_interface__calcium_dynamics__g = 0.001 * (((((var_calcium_dynamics__g_inf > var_chaste_interface__calcium_dynamics__g) && ((0.001 * var_chaste_interface__Membrane__Vm) >  -0.059999999999999998)) ? 0.0 : 1.0) * (var_calcium_dynamics__g_inf - var_chaste_interface__calcium_dynamics__g)) * 500.0); // 'per millisecond'
        
        if (mSetVoltageDerivativeToZero)
        {
            d_dt_chaste_interface__Membrane__Vm = 0.0;
        }
        else
        {
            const double var_i_K1__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
            const double var_electric_potentials__E_K = ((8.3144720000000003 * NV_Ith_S(mParameters, 17)) * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 2) / NV_Ith_S(mParameters, 0)); // volt
            const double var_i_K1__alpha_K1 = 3.9100000000000001 / (1.0 + exp(0.59419999999999995 * (((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0)) - 200.0))); // dimensionless
            const double var_i_K1__i_K1 = NV_Ith_S(mParameters, 11) * ((var_i_K1__alpha_K1 / (var_i_K1__alpha_K1 + ((( -1.5089999999999999 * exp(0.00020000000000000001 * (((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0)) + 100.0))) + exp(0.58860000000000001 * (((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0)) - 10.0))) / (1.0 + exp(0.45469999999999999 * ((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0))))))) * ((var_i_K1__Vm - var_electric_potentials__E_K) * sqrt(NV_Ith_S(mParameters, 2) * 0.18518518518518517))); // A_per_F
            const double var_i_f__i_f = NV_Ith_S(mParameters, 10) * (var_chaste_interface__i_f_Xf_gate__Xf * (var_i_f__Vm -  -0.017000000000000001)); // A_per_F
            const double var_i_Kr__i_Kr = 1.0 * (NV_Ith_S(mParameters, 12) * ((var_i_Kr__Vm - var_electric_potentials__E_K) * (var_chaste_interface__i_Kr_Xr1_gate__Xr1 * (var_chaste_interface__i_Kr_Xr2_gate__Xr2 * sqrt(NV_Ith_S(mParameters, 2) * 0.18518518518518517))))); // A_per_F
            const double var_i_Ks__i_Ks = 1.0 * (NV_Ith_S(mParameters, 13) * ((var_i_Ks__Vm - (((8.3144720000000003 * NV_Ith_S(mParameters, 17)) * 1.0364268649036186e-05) * log((NV_Ith_S(mParameters, 2) + (0.029999999999999999 * NV_Ith_S(mParameters, 3))) / (NV_Ith_S(mParameters, 0) + (0.029999999999999999 * var_chaste_interface__sodium_dynamics__Nai))))) * (pow(var_chaste_interface__i_Ks_Xs_gate__Xs, 2.0) * (1.0 + (0.59999999999999998 / (1.0 + pow(3.8000000000000002e-05 / var_chaste_interface__calcium_dynamics__Cai, 1.3999999999999999))))))); // A_per_F
            const double var_i_to__i_to = NV_Ith_S(mParameters, 16) * ((var_i_to__Vm - var_electric_potentials__E_K) * (var_chaste_interface__i_to_q_gate__q * var_chaste_interface__i_to_r_gate__r)); // A_per_F
            const double var_chaste_interface__stim_mode__i_stim = -GetIntracellularAreaStimulus(var_chaste_interface__environment__time);
            d_dt_chaste_interface__Membrane__Vm = -((var_i_K1__i_K1 + var_i_to__i_to + var_i_Kr__i_Kr + var_i_Ks__i_Ks + var_i_CaL__i_CaL + var_i_NaK__i_NaK + var_i_Na__i_Na + var_i_NaCa__i_NaCa + var_i_PCa__i_PCa + var_i_f__i_f + var_i_b_Na__i_b_Na + var_i_b_Ca__i_b_Ca) - (var_chaste_interface__stim_mode__i_stim / HeartConfig::Instance()->GetCapacitance())); // 'millivolt per millisecond'
        }
        
        NV_Ith_S(rDY, 0) = d_dt_chaste_interface__Membrane__Vm;
        NV_Ith_S(rDY, 1) = d_dt_chaste_interface__i_Na_m_gate__m;
        NV_Ith_S(rDY, 2) = d_dt_chaste_interface__i_Na_h_gate__h;
        NV_Ith_S(rDY, 3) = d_dt_chaste_interface__i_Na_j_gate__j;
        NV_Ith_S(rDY, 4) = d_dt_chaste_interface__i_CaL_d_gate__d;
        NV_Ith_S(rDY, 5) = d_dt_chaste_interface__i_CaL_f1_gate__f1;
        NV_Ith_S(rDY, 6) = d_dt_chaste_interface__i_CaL_f2_gate__f2;
        NV_Ith_S(rDY, 7) = d_dt_chaste_interface__i_CaL_fCa_gate__fCa;
        NV_Ith_S(rDY, 8) = d_dt_chaste_interface__i_Kr_Xr1_gate__Xr1;
        NV_Ith_S(rDY, 9) = d_dt_chaste_interface__i_Kr_Xr2_gate__Xr2;
        NV_Ith_S(rDY, 10) = d_dt_chaste_interface__i_Ks_Xs_gate__Xs;
        NV_Ith_S(rDY, 11) = d_dt_chaste_interface__i_f_Xf_gate__Xf;
        NV_Ith_S(rDY, 12) = d_dt_chaste_interface__i_to_q_gate__q;
        NV_Ith_S(rDY, 13) = d_dt_chaste_interface__i_to_r_gate__r;
        NV_Ith_S(rDY, 14) = d_dt_chaste_interface__sodium_dynamics__Nai;
        NV_Ith_S(rDY, 15) = d_dt_chaste_interface__calcium_dynamics__Cai;
        NV_Ith_S(rDY, 16) = d_dt_chaste_interface__calcium_dynamics__Ca_SR;
        NV_Ith_S(rDY, 17) = d_dt_chaste_interface__calcium_dynamics__g;
    }


/********************************************************************************** Turn off Analytic Jacobian **********************************************************************************
    
    void Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode::EvaluateAnalyticJacobian(long int N, double var_chaste_interface__environment__time, N_Vector rY, N_Vector rDY, CHASTE_CVODE_DENSE_MATRIX rJacobian, N_Vector rTmp1, N_Vector rTmp2, N_Vector rTmp3)
    {
        double var_chaste_interface__Membrane__Vm = (mSetVoltageDerivativeToZero ? this->mFixedVoltage : NV_Ith_S(rY, 0));
        // Units: millivolt; Initial value: -74.3340057624
        double var_chaste_interface__i_Na_m_gate__m = NV_Ith_S(rY, 1);
        // Units: dimensionless; Initial value: 0.102953468725004
        double var_chaste_interface__i_Na_h_gate__h = NV_Ith_S(rY, 2);
        // Units: dimensionless; Initial value: 0.786926637881461
        double var_chaste_interface__i_Na_j_gate__j = NV_Ith_S(rY, 3);
        // Units: dimensionless; Initial value: 0.253943221774722
        double var_chaste_interface__i_CaL_d_gate__d = NV_Ith_S(rY, 4);
        // Units: dimensionless; Initial value: 8.96088425225182e-5
        double var_chaste_interface__i_CaL_f1_gate__f1 = NV_Ith_S(rY, 5);
        // Units: dimensionless; Initial value: 0.970411811263976
        double var_chaste_interface__i_CaL_f2_gate__f2 = NV_Ith_S(rY, 6);
        // Units: dimensionless; Initial value: 0.999965815466749
        double var_chaste_interface__i_CaL_fCa_gate__fCa = NV_Ith_S(rY, 7);
        // Units: dimensionless; Initial value: 0.998925296531804
        double var_chaste_interface__i_Kr_Xr1_gate__Xr1 = NV_Ith_S(rY, 8);
        // Units: dimensionless; Initial value: 0.00778547011240132
        double var_chaste_interface__i_Kr_Xr2_gate__Xr2 = NV_Ith_S(rY, 9);
        // Units: dimensionless; Initial value: 0.432162576531617
        double var_chaste_interface__i_Ks_Xs_gate__Xs = NV_Ith_S(rY, 10);
        // Units: dimensionless; Initial value: 0.0322944866983666
        double var_chaste_interface__i_f_Xf_gate__Xf = NV_Ith_S(rY, 11);
        // Units: dimensionless; Initial value: 0.100615100568753
        double var_chaste_interface__i_to_q_gate__q = NV_Ith_S(rY, 12);
        // Units: dimensionless; Initial value: 0.839295925773219
        double var_chaste_interface__i_to_r_gate__r = NV_Ith_S(rY, 13);
        // Units: dimensionless; Initial value: 0.00573289893326379
        //double var_chaste_interface__sodium_dynamics__Nai = NV_Ith_S(rY, 14);
        double var_chaste_interface__sodium_dynamics__Nai = (mSetNaiDerivativeToZero ? this->mFixedNai : NV_Ith_S(rY, 14));
        // Units: millimolar; Initial value: 10.9248496211574
        //double var_chaste_interface__calcium_dynamics__Cai = NV_Ith_S(rY, 15);
        double var_chaste_interface__calcium_dynamics__Cai = (mSetCaiDerivativeToZero ? this->mFixedCai : NV_Ith_S(rY, 15));
        // Units: millimolar; Initial value: 1.80773974140477e-5
        double var_chaste_interface__calcium_dynamics__Ca_SR = NV_Ith_S(rY, 16);
        // Units: millimolar; Initial value: -0.2734234751931
        double var_chaste_interface__calcium_dynamics__g = NV_Ith_S(rY, 17);
        // Units: dimensionless; Initial value: 0.999999981028517
        
        const double var_chaste_interface__environment__fake_dt = 0.001; // second
        const double var_Membrane__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_t2 = 8.3144720000000003 * NV_Ith_S(mParameters, 17); // dimensionless
        const double var_t8 = (var_t2 * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 2) / NV_Ith_S(mParameters, 0)); // dimensionless
        const double var_t11 = exp(((594.20000000000005 * var_Membrane__Vm) - (594.20000000000005 * var_t8)) - 118.84); // dimensionless
        const double var_t12 = 1.0 + var_t11; // dimensionless
        const double var_t14 = 1.0 / pow(var_t12, 2.0); // dimensionless
        const double var_t16 = 1.0 / var_t12; // dimensionless
        const double var_t21 = exp(((0.20000000000000001 * var_Membrane__Vm) - (0.20000000000000001 * var_t8)) + 0.02); // dimensionless
        const double var_t26 = exp(((588.60000000000002 * var_Membrane__Vm) - (588.60000000000002 * var_t8)) - 5.8860000000000001); // dimensionless
        const double var_t27 = ( -1.5089999999999999 * var_t21) + var_t26; // dimensionless
        const double var_t31 = exp((454.69999999999999 * var_Membrane__Vm) - (454.69999999999999 * var_t8)); // dimensionless
        const double var_t32 = 1.0 + var_t31; // dimensionless
        const double var_t33 = 1.0 / var_t32; // dimensionless
        const double var_t35 = (3.9100000000000001 * var_t16) + (var_t27 * var_t33); // dimensionless
        const double var_t36 = 1.0 / var_t35; // dimensionless
        const double var_t38 = var_Membrane__Vm - var_t8; // dimensionless
        const double var_t39 = pow(NV_Ith_S(mParameters, 2), 0.5); // dimensionless
        const double var_t40 = var_t38 * var_t39; // dimensionless
        const double var_t44 = NV_Ith_S(mParameters, 11) * var_t16; // dimensionless
        const double var_t71 = 1.0 * NV_Ith_S(mParameters, 12); // dimensionless
        const double var_t81 = 1.0 * NV_Ith_S(mParameters, 13); // dimensionless
        const double var_t82 = pow(var_chaste_interface__i_Ks_Xs_gate__Xs, 2.0); // dimensionless
        const double var_t83 = 1.0 / var_chaste_interface__calcium_dynamics__Cai; // dimensionless
        const double var_t86 = 1.0 + (6.481821026e-07 * pow(var_t83, 1.3999999999999999)); // dimensionless
        const double var_t89 = 1.0 + (0.59999999999999998 / var_t86); // dimensionless
        const double var_t95 = 1.0 / NV_Ith_S(mParameters, 17); // dimensionless
        const double var_t96 = 0.12027221933034352 * var_t95; // dimensionless
        const double var_t98 = (var_Membrane__Vm * 96485.341499999995) * var_t96; // dimensionless
        const double var_t100 = exp(2.0 * var_t98); // dimensionless
        const double var_t103 = (var_chaste_interface__calcium_dynamics__Cai * var_t100) - (0.34100000000000003 * NV_Ith_S(mParameters, 1)); // dimensionless
        const double var_t104 = var_t96 * var_t103; // dimensionless
        const double var_t106 = var_t100 - 1.0; // dimensionless
        const double var_t107 = 1.0 / var_t106; // dimensionless
        const double var_t109 = var_chaste_interface__i_CaL_f1_gate__f1 * var_chaste_interface__i_CaL_f2_gate__f2; // dimensionless
        const double var_t110 = var_t109 * var_chaste_interface__i_CaL_fCa_gate__fCa; // dimensionless
        const double var_t113 = (((((4.0 * NV_Ith_S(mParameters, 4)) * 9309421124.3716221) * var_t104) * var_t107) * var_chaste_interface__i_CaL_d_gate__d) * var_t110; // dimensionless
        const double var_t114 = NV_Ith_S(mParameters, 4) * var_Membrane__Vm; // dimensionless
        const double var_t116 = (var_t114 * 9309421124.3716221) * 96485.341499999995; // dimensionless
        const double var_t121 = 0.014465406742646259 / pow(NV_Ith_S(mParameters, 17), 2.0); // dimensionless
        const double var_t126 = ((var_t100 * var_t107) * var_chaste_interface__i_CaL_d_gate__d) * var_t110; // dimensionless
        const double var_t128 = (((8.0 * var_t116) * var_t121) * var_chaste_interface__calcium_dynamics__Cai) * var_t126; // dimensionless
        const double var_t135 = var_chaste_interface__i_CaL_f2_gate__f2 * var_chaste_interface__i_CaL_fCa_gate__fCa; // dimensionless
        const double var_t139 = (((((((8.0 * var_t116) * var_t121) * var_t103) / pow(var_t106, 2.0)) * var_chaste_interface__i_CaL_d_gate__d) * var_chaste_interface__i_CaL_f1_gate__f1) * var_t135) * var_t100; // dimensionless
        const double var_t140 = NV_Ith_S(mParameters, 15) * NV_Ith_S(mParameters, 2); // dimensionless
        const double var_t142 = 1.0 / (NV_Ith_S(mParameters, 2) + 1.0); // dimensionless
        const double var_t143 = var_t140 * var_t142; // dimensionless
        const double var_t144 = var_chaste_interface__sodium_dynamics__Nai + 40.0; // dimensionless
        const double var_t145 = 1.0 / var_t144; // dimensionless
        const double var_t148 = exp( -0.10000000000000001 * var_t98); // dimensionless
        const double var_t150 = exp(-var_t98); // dimensionless
        const double var_t152 = (1.0 + (0.1245 * var_t148)) + (0.035299999999999998 * var_t150); // dimensionless
        const double var_t165 = (((var_t143 * var_chaste_interface__sodium_dynamics__Nai) * var_t145) / pow(var_t152, 2.0)) * ((( -144.47610163038613 * var_t95) * var_t148) - ((409.63906727330362 * var_t95) * var_t150)); // dimensionless
        const double var_t170 = 1.0 * NV_Ith_S(mParameters, 9); // dimensionless
        const double var_t171 = pow(var_chaste_interface__i_Na_m_gate__m, 2.0); // dimensionless
        const double var_t172 = var_t171 * var_chaste_interface__i_Na_m_gate__m; // dimensionless
        const double var_t173 = var_t172 * var_chaste_interface__i_Na_h_gate__h; // dimensionless
        const double var_t175 = (var_t170 * var_t173) * var_chaste_interface__i_Na_j_gate__j; // dimensionless
        const double var_t179 = 11604.506155051095 * var_t95; // dimensionless
        const double var_t181 = exp((0.34999999999999998 * var_Membrane__Vm) * var_t179); // dimensionless
        const double var_t183 = pow(var_chaste_interface__sodium_dynamics__Nai, 2.0); // dimensionless
        const double var_t184 = var_t183 * var_chaste_interface__sodium_dynamics__Nai; // dimensionless
        const double var_t193 = exp(( -0.65000000000000002 * var_Membrane__Vm) * var_t179); // dimensionless
        const double var_t195 = pow(NV_Ith_S(mParameters, 3), 2.0) * NV_Ith_S(mParameters, 3); // dimensionless
        const double var_t196 = var_t193 * var_t195; // dimensionless
        const double var_t198 = (var_t196 * var_chaste_interface__calcium_dynamics__Cai) * 2.8571431999999999; // dimensionless
        const double var_t205 = 1.0 / (669921.875 + var_t195); // dimensionless
        const double var_t207 = 1.0 / (1.3799999999999999 + NV_Ith_S(mParameters, 1)); // dimensionless
        const double var_t208 = var_t205 * var_t207; // dimensionless
        const double var_t210 = 1.0 + (0.10000000000000001 * var_t193); // dimensionless
        const double var_t211 = 1.0 / var_t210; // dimensionless
        const double var_t212 = var_t208 * var_t211; // dimensionless
        const double var_t213 = (NV_Ith_S(mParameters, 14) * (((((4061.5771542678826 * var_t95) * var_t181) * var_t184) * NV_Ith_S(mParameters, 1)) - (( -62715.471975 * var_t96) * var_t198))) * var_t212; // dimensionless
        const double var_t217 = NV_Ith_S(mParameters, 14) * (((var_t181 * var_t184) * NV_Ith_S(mParameters, 1)) - var_t198); // dimensionless
        const double var_t226 = ((((((var_t217 * var_t208) / pow(var_t210, 2.0)) * 0.10000000000000001) *  -0.65000000000000002) * 96485.341499999995) * var_t96) * var_t193; // dimensionless
        const double var_t230 = var_chaste_interface__i_Na_h_gate__h * var_chaste_interface__i_Na_j_gate__j; // dimensionless
        const double var_t231 = 1.0 / var_chaste_interface__sodium_dynamics__Nai; // dimensionless
        const double var_t236 = var_Membrane__Vm - ((var_t2 * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 3) * var_t231)); // dimensionless
        const double var_t247 = ((var_t114 * 9309421124.3716221) * 0.12027221933034352) * var_t95; // dimensionless
        const double var_t248 = var_t103 * var_t107; // dimensionless
        const double var_t257 = var_chaste_interface__i_CaL_d_gate__d * var_chaste_interface__i_CaL_f1_gate__f1; // dimensionless
        const double var_t278 = 1.0 / (NV_Ith_S(mParameters, 0) + (0.029999999999999999 * var_chaste_interface__sodium_dynamics__Nai)); // dimensionless
        const double var_t283 = var_Membrane__Vm - ((var_t2 * 1.0364268649036186e-05) * log((NV_Ith_S(mParameters, 2) + (0.029999999999999999 * NV_Ith_S(mParameters, 3))) * var_t278)); // dimensionless
        const double var_t290 = NV_Ith_S(mParameters, 16) * var_t38; // dimensionless
        const double var_t300 = 1.0 / var_t152; // dimensionless
        const double var_t302 = ((var_t140 * var_t142) * var_t145) * var_t300; // dimensionless
        const double var_t307 = ((var_t143 * var_chaste_interface__sodium_dynamics__Nai) / pow(var_t144, 2.0)) * var_t300; // dimensionless
        const double var_t310 = NV_Ith_S(mParameters, 17) * 1.0364268649036186e-05; // dimensionless
        const double var_t311 = var_t310 * var_t231; // dimensionless
        const double var_t313 = (((var_t170 * var_t173) * var_chaste_interface__i_Na_j_gate__j) * 8.3144720000000003) * var_t311; // dimensionless
        const double var_t317 = var_t207 * var_t211; // dimensionless
        const double var_t319 = ((((NV_Ith_S(mParameters, 14) * var_t181) * var_t183) * NV_Ith_S(mParameters, 1)) * var_t205) * var_t317; // dimensionless
        const double var_t322 = (NV_Ith_S(mParameters, 6) * 8.3144720000000003) * var_t311; // dimensionless
        const double var_t329 = pow(var_chaste_interface__calcium_dynamics__Cai, 2.0); // dimensionless
        const double var_t330 = 1.0 / var_t329; // dimensionless
        const double var_t336 = (4.0 * var_t247) * var_t126; // dimensionless
        const double var_t341 = ((((NV_Ith_S(mParameters, 14) * var_t193) * var_t195) * 2.8571431999999999) * var_t205) * var_t317; // dimensionless
        const double var_t342 = var_chaste_interface__calcium_dynamics__Cai + 0.00050000000000000001; // dimensionless
        const double var_t343 = 1.0 / var_t342; // dimensionless
        const double var_t344 = NV_Ith_S(mParameters, 7) * var_t343; // dimensionless
        const double var_t345 = NV_Ith_S(mParameters, 7) * var_chaste_interface__calcium_dynamics__Cai; // dimensionless
        const double var_t348 = var_t345 / pow(var_t342, 2.0); // dimensionless
        const double var_t352 = (((0.5 * NV_Ith_S(mParameters, 5)) * 8.3144720000000003) * var_t310) * var_t83; // dimensionless
        const double var_t356 = exp(( -169.4915254 * var_Membrane__Vm) - 5.7796610160000004); // dimensionless
        const double var_t357 = 1.0 + var_t356; // dimensionless
        const double var_t358 = pow(var_t357, 0.33333333333333331); // dimensionless
        const double var_t362 = 200.0 * var_Membrane__Vm; // dimensionless
        const double var_t364 = exp((-var_t362) - 12.0); // dimensionless
        const double var_t365 = 1.0 + var_t364; // dimensionless
        const double var_t367 = exp(var_t362 + 7.0); // dimensionless
        const double var_t368 = 1.0 + var_t367; // dimensionless
        const double var_t373 = exp((5.0 * var_Membrane__Vm) - 0.25); // dimensionless
        const double var_t374 = 1.0 + var_t373; // dimensionless
        const double var_t377 = (0.10000000000000001 / var_t368) + (0.10000000000000001 / var_t374); // dimensionless
        const double var_t378 = 1.0 / var_t377; // dimensionless
        const double var_t379 = var_t365 * var_t378; // dimensionless
        const double var_t383 = (1.0 / var_t358) - var_chaste_interface__i_Na_m_gate__m; // dimensionless
        const double var_t406 = exp((175.43859649999999 * var_Membrane__Vm) + 12.64912281); // dimensionless
        const double var_t407 = 1.0 + var_t406; // dimensionless
        const double var_t408 = pow(var_t407, 0.5); // dimensionless
        const double var_t411 = ((1.0 / var_t408) / var_t407) * var_t406; // dimensionless
        const double var_t412 = var_Membrane__Vm <  -0.040000000000000001; // dimensionless
        const double var_t415 = exp(( -147.05882349999999 * var_Membrane__Vm) - 11.764705879999999); // dimensionless
        const double var_t422 = exp(348.5 * var_Membrane__Vm); // dimensionless
        const double var_t436 = var_t412 ? (0.0015 / ((var_t412 ? (0.057000000000000002 * var_t415) : 0.0) + (var_t412 ? ((2.7000000000000002 * exp(79.0 * var_Membrane__Vm)) + (310000.0 * var_t422)) : (0.77000000000000002 / (0.13 + (0.13 * exp(( -90.090090090000004 * var_Membrane__Vm) - 0.96036036039999995))))))) : 0.002542; // dimensionless
        const double var_t437 = 1.0 / var_t436; // dimensionless
        const double var_t440 = 1.0 / var_t408; // dimensionless
        const double var_t445 = var_Membrane__Vm <  -0.040000000000000001; // dimensionless
        const double var_t448 = exp(79.0 * var_Membrane__Vm); // dimensionless
        const double var_t466 = exp(244.40000000000001 * var_Membrane__Vm); // dimensionless
        const double var_t469 = exp( -43.909999999999997 * var_Membrane__Vm); // dimensionless
        const double var_t470 = 6.9480000000000002e-06 * var_t469; // dimensionless
        const double var_t472 = 1000.0 * var_Membrane__Vm; // dimensionless
        const double var_t484 = 137.80000000000001 * var_Membrane__Vm; // dimensionless
        const double var_t493 = 100.0 * var_Membrane__Vm; // dimensionless
        const double var_t501 = (var_t412 ? (((( -25428.0 * var_t466) - var_t470) * (var_t472 + 37.780000000000001)) / (1.0 + exp((311.0 * var_Membrane__Vm) + 24.640529999999998))) : 0.0) + (var_t412 ? ((0.024240000000000001 * exp( -10.52 * var_Membrane__Vm)) / (1.0 + exp((-var_t484) - 5.5312919999999997))) : ((0.59999999999999998 * exp(57.0 * var_Membrane__Vm)) / (1.0 + exp((-var_t493) - 3.2000000000000002)))); // dimensionless
        const double var_t509 = (1000.0 * var_Membrane__Vm) + 37.780000000000001; // dimensionless
        const double var_t513 = exp((311.0 * var_Membrane__Vm) + 24.640530009999999); // dimensionless
        const double var_t514 = 1.0 + var_t513; // dimensionless
        const double var_t515 = 1.0 / var_t514; // dimensionless
        const double var_t518 = ( -25428.0 * var_t466) - var_t470; // dimensionless
        const double var_t530 = exp( -10.52 * var_Membrane__Vm); // dimensionless
        const double var_t532 = exp((-var_t484) - 5.5312919950000001); // dimensionless
        const double var_t533 = 1.0 + var_t532; // dimensionless
        const double var_t544 = exp(57.0 * var_Membrane__Vm); // dimensionless
        const double var_t547 = exp(( -100.0 * var_Membrane__Vm) - 3.2000000000000002); // dimensionless
        const double var_t548 = 1.0 + var_t547; // dimensionless
        const double var_t565 = exp(( -142.85714285714286 * var_Membrane__Vm) - 1.3); // dimensionless
        const double var_t566 = 1.0 + var_t565; // dimensionless
        const double var_t570 = 76.92307692307692 * var_Membrane__Vm; // dimensionless
        const double var_t572 = exp((-var_t570) - 2.6923076923076925); // dimensionless
        const double var_t573 = 1.0 + var_t572; // dimensionless
        const double var_t576 = 0.25 + (1.3999999999999999 / var_t573); // dimensionless
        const double var_t578 = exp(var_t362 + 1.0); // dimensionless
        const double var_t579 = 1.0 + var_t578; // dimensionless
        const double var_t580 = 1.0 / var_t579; // dimensionless
        const double var_t583 = 50.0 * var_Membrane__Vm; // dimensionless
        const double var_t585 = exp((-var_t583) + 2.5); // dimensionless
        const double var_t586 = 1.0 + var_t585; // dimensionless
        const double var_t589 = ((0.0014 * var_t576) * var_t580) + (0.001 / var_t586); // dimensionless
        const double var_t590 = 1.0 / var_t589; // dimensionless
        const double var_t617 = exp((333.33333333333331 * var_Membrane__Vm) + 8.6666666666666661); // dimensionless
        const double var_t618 = 1.0 + var_t617; // dimensionless
        const double var_t622 = var_t472 + 27.0; // dimensionless
        const double var_t623 = pow(var_t622, 2.0); // dimensionless
        const double var_t626 = exp( -0.0044444444444444444 * pow(var_t623, 2.0)); // dimensionless
        const double var_t628 = 100.0 * var_Membrane__Vm; // dimensionless
        const double var_t630 = exp(1.3 - var_t628); // dimensionless
        const double var_t631 = 1.0 + var_t630; // dimensionless
        const double var_t635 = exp(3.0 + var_t628); // dimensionless
        const double var_t636 = 1.0 + var_t635; // dimensionless
        const double var_t637 = 1.0 / var_t636; // dimensionless
        const double var_t639 = ((20.0 + (1102.5 * var_t626)) + (200.0 / var_t631)) + (180.0 * var_t637); // dimensionless
        const double var_t640 = 1.0 / var_t639; // dimensionless
        const double var_t642 = (1.0 / var_t618) - var_chaste_interface__i_CaL_f1_gate__f1; // dimensionless
        const double var_t643 = 0.0 < var_t642; // dimensionless
        const double var_t646 = var_t643 ? (0.92835000000000001 + (1433.0 * var_chaste_interface__calcium_dynamics__Cai)) : 1.0; // dimensionless
        const double var_t647 = 1.0 / var_t646; // dimensionless
        const double var_t648 = var_t640 * var_t647; // dimensionless
        const double var_t663 = (1.0 / pow(var_t636, 2.0)) * var_t635; // dimensionless
        const double var_t670 = var_t642 * var_t640; // dimensionless
        const double var_t672 = 1.0 / pow(var_t646, 2.0); // dimensionless
        const double var_t689 = exp((250.0 * var_Membrane__Vm) + 8.75); // dimensionless
        const double var_t690 = 1.0 + var_t689; // dimensionless
        const double var_t697 = exp( -0.0058823529411764705 * pow(var_t472 + 25.0, 2.0)); // dimensionless
        const double var_t700 = exp(2.5 - var_t628); // dimensionless
        const double var_t701 = 1.0 + var_t700; // dimensionless
        const double var_t705 = ((0.59999999999999998 * var_t697) + (0.031 / var_t701)) + (0.016 * var_t637); // dimensionless
        const double var_t706 = 1.0 / var_t705; // dimensionless
        const double var_t727 =  -0.059999999999999998 < var_Membrane__Vm; // dimensionless
        const double var_t728 = pow(var_t329, 2.0); // dimensionless
        const double var_t729 = pow(var_t728, 2.0); // dimensionless
        const double var_t731 = 1.0 + (5.9537418169999996e+25 * var_t729); // dimensionless
        const double var_t733 = 0.7601094558 / var_t731; // dimensionless
        const double var_t736 = exp((10000.0 * var_chaste_interface__calcium_dynamics__Cai) - 9.0); // dimensionless
        const double var_t737 = 1.0 + var_t736; // dimensionless
        const double var_t739 = 0.076010945580000003 / var_t737; // dimensionless
        const double var_t742 = exp((1250.0 * var_chaste_interface__calcium_dynamics__Cai) - 0.9375); // dimensionless
        const double var_t743 = 1.0 + var_t742; // dimensionless
        const double var_t745 = 0.22803283669999999 / var_t743; // dimensionless
        const double var_t748 = var_t727 && (var_chaste_interface__i_CaL_fCa_gate__fCa < ((var_t733 + var_t739) + var_t745)); // dimensionless
        const double var_t753 = ((var_t748 ? 0.0 : 0.0) * (((var_t733 + var_t739) + var_t745) - var_chaste_interface__i_CaL_fCa_gate__fCa)) * 500.0; // dimensionless
        const double var_t754 = var_t748 ? 0.0 : 1.0; // dimensionless
        const double var_t759 = var_t329 * var_chaste_interface__calcium_dynamics__Cai; // dimensionless
        const double var_t760 = var_t728 * var_t759; // dimensionless
        const double var_t795 = exp(((((( -204.0816327 * var_t2) * 1.0364268649036186e-05) * 0.43478260869565222) * log((pow(pow(1.0 + (0.3846153846 * NV_Ith_S(mParameters, 1)), 2.0), 2.0) * 40.0) / pow(pow(1.0 + (1.724137931 * NV_Ith_S(mParameters, 1)), 2.0), 2.0))) - 3.8775510209999999) - (204.0816327 * var_Membrane__Vm)); // dimensionless
        const double var_t796 = 1.0 + var_t795; // dimensionless
        const double var_t801 = exp( -4.5 - var_t628); // dimensionless
        const double var_t802 = 1.0 + var_t801; // dimensionless
        const double var_t805 = exp(2.6086956520000002 + (86.956521739999999 * var_Membrane__Vm)); // dimensionless
        const double var_t806 = 1.0 + var_t805; // dimensionless
        const double var_t807 = var_t802 * var_t806; // dimensionless
        const double var_t811 = (1.0 / var_t796) - var_chaste_interface__i_Kr_Xr1_gate__Xr1; // dimensionless
        const double var_t822 = exp((20.0 * var_Membrane__Vm) + 1.76); // dimensionless
        const double var_t823 = 1.0 + var_t822; // dimensionless
        const double var_t828 = exp( -3.0 - var_t583); // dimensionless
        const double var_t829 = 1.0 + var_t828; // dimensionless
        const double var_t831 = exp( -3.0 + var_t583); // dimensionless
        const double var_t832 = 1.0 + var_t831; // dimensionless
        const double var_t833 = var_t829 * var_t832; // dimensionless
        const double var_t837 = (1.0 / var_t823) - var_chaste_interface__i_Kr_Xr2_gate__Xr2; // dimensionless
        const double var_t848 = exp(( -62.5 * var_Membrane__Vm) - 1.25); // dimensionless
        const double var_t849 = 1.0 + var_t848; // dimensionless
        const double var_t855 = exp( -1.6666666666666667 - (166.66666666666666 * var_Membrane__Vm)); // dimensionless
        const double var_t857 = pow(1.0 + var_t855, 0.5); // dimensionless
        const double var_t858 = var_t857 * var_t832; // dimensionless
        const double var_t862 = (1.0 / var_t849) - var_chaste_interface__i_Ks_Xs_gate__Xs; // dimensionless
        const double var_t874 = exp(var_t362 + 15.57); // dimensionless
        const double var_t875 = 1.0 + var_t874; // dimensionless
        const double var_t880 = exp(var_t628 + 1.5); // dimensionless
        const double var_t881 = 1.0 + var_t880; // dimensionless
        const double var_t890 = exp(var_t570 + 4.0769230769230766); // dimensionless
        const double var_t891 = 1.0 + var_t890; // dimensionless
        const double var_t897 = exp(( -80.0 * var_Membrane__Vm) - 3.52); // dimensionless
        const double var_t900 = exp(var_t493 + 4.593); // dimensionless
        const double var_t902 = (0.56999999999999995 * var_t897) + (0.065000000000000002 * var_t900); // dimensionless
        const double var_t905 = 0.0060600000000000003 + (0.039101999999999998 / var_t902); // dimensionless
        const double var_t906 = 1.0 / var_t905; // dimensionless
        const double var_t925 = exp(( -53.333333330000002 * var_Membrane__Vm) + 1.189333333); // dimensionless
        const double var_t926 = 1.0 + var_t925; // dimensionless
        const double var_t932 = exp((90.0 * var_Membrane__Vm) + 2.7549000000000001); // dimensionless
        const double var_t936 = exp(( -120.0 * var_Membrane__Vm) - 2.8607999999999998); // dimensionless
        const double var_t938 = (1.0369999999999999 * var_t932) + (0.36899999999999999 * var_t936); // dimensionless
        const double var_t941 = 0.0027535200000000002 + (0.01440516 / var_t938); // dimensionless
        const double var_t942 = 1.0 / var_t941; // dimensionless
        const double var_t968 = NV_Ith_S(mParameters, 8) * 1.0; // dimensionless
        const double var_t977 = (var_t968 * NV_Ith_S(mParameters, 9)) * var_t172; // dimensionless
        const double var_t1002 = var_chaste_interface__calcium_dynamics__Cai + 0.001; // dimensionless
        const double var_t1003 = pow(var_t1002, 2.0); // dimensionless
        const double var_t1006 = 1.0 + (0.00025000000000000001 / var_t1003); // dimensionless
        const double var_t1007 = 1.0 / var_t1006; // dimensionless
        const double var_t1012 = NV_Ith_S(mParameters, 8) * 0.00011363636363636364; // dimensionless
        const double var_t1013 = var_t1012 * 1.0364268649036186e-05; // dimensionless
        const double var_t1016 = pow(var_chaste_interface__calcium_dynamics__Ca_SR, 2.0); // dimensionless
        const double var_t1019 = 0.0625 + var_t1016; // dimensionless
        const double var_t1020 = 1.0 / var_t1019; // dimensionless
        const double var_t1022 = 8.2319999999999993 + ((16.463999999999999 * var_t1016) * var_t1020); // dimensionless
        const double var_t1038 = ((var_t1007 * NV_Ith_S(mParameters, 4)) * var_Membrane__Vm) * var_t179; // dimensionless
        const double var_t1039 = var_t248 * var_chaste_interface__i_CaL_d_gate__d; // dimensionless
        const double var_t1063 = (var_chaste_interface__calcium_dynamics__Ca_SR - var_chaste_interface__calcium_dynamics__Cai) * 0.00044443999999999999; // dimensionless
        const double var_t1066 = 1.0 + (6.2499999999999997e-08 * var_t330); // dimensionless
        const double var_t1068 = 0.56064000000000003 / var_t1066; // dimensionless
        const double var_t1070 = (var_t1022 * var_chaste_interface__i_CaL_d_gate__d) * var_chaste_interface__calcium_dynamics__g; // dimensionless
        const double var_t1102 = ((1.1212800000000001 / pow(var_t1066, 2.0)) * 6.2499999999999997e-08) / var_t759; // dimensionless
        const double var_t1122 = 0.00044443999999999999 + (((0.082199999999999995 * (((16.463999999999999 * var_chaste_interface__calcium_dynamics__Ca_SR) * var_t1020) - (((16.463999999999999 * var_t1016) * var_chaste_interface__calcium_dynamics__Ca_SR) / pow(var_t1019, 2.0)))) * var_chaste_interface__i_CaL_d_gate__d) * var_chaste_interface__calcium_dynamics__g); // dimensionless
        const double var_t1128 = var_chaste_interface__calcium_dynamics__Ca_SR + 0.29999999999999999; // dimensionless
        const double var_t1129 = pow(var_t1128, 2.0); // dimensionless
        const double var_t1132 = 1.0 + (3.0 / var_t1129); // dimensionless
        const double var_t1134 = (1.0 / var_t1132) * 8800.0; // dimensionless
        const double var_t1136 = 0.0017131207921470542 * var_t1022; // dimensionless
        const double var_t1164 = (5.43991024e+20 * var_t728) * var_t329; // dimensionless
        const double var_t1168 = 1.9720198860000001e+55 * pow(var_t729, 2.0); // dimensionless
        const double var_t1171 = (var_chaste_interface__calcium_dynamics__Cai <= 0.00035) ? (1.0 / (1.0 + var_t1164)) : (1.0 / (1.0 + var_t1168)); // dimensionless
        const double var_t1173 = (var_chaste_interface__calcium_dynamics__g < var_t1171) && var_t727; // dimensionless
        const double var_t1178 = ((var_t1173 ? 0.0 : 0.0) * (var_t1171 - var_chaste_interface__calcium_dynamics__g)) * 500.0; // dimensionless
        const double var_t1179 = var_t1173 ? 0.0 : 1.0; // dimensionless
        
        // Matrix entries
        DENSE_ELEM(rJacobian, 0, 0) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : ((((((((((((((((((((999.79860140000005 * NV_Ith_S(mParameters, 11)) * var_t14) * var_t36) * var_t40) * var_t11) + ((((1.6825960980000001 * var_t44) / pow(var_t35, 2.0)) * var_t40) * (((( -2323.3220000000001 * var_t14) * var_t11) + ((( -0.30180000000000001 * var_t21) + (588.60000000000002 * var_t26)) * var_t33)) - (((454.69999999999999 * var_t27) / pow(var_t32, 2.0)) * var_t31)))) - (((1.6825960980000001 * var_t44) * var_t36) * var_t39)) - ((NV_Ith_S(mParameters, 16) * var_chaste_interface__i_to_q_gate__q) * var_chaste_interface__i_to_r_gate__r)) - ((((0.43033148290000001 * var_t71) * var_chaste_interface__i_Kr_Xr1_gate__Xr1) * var_chaste_interface__i_Kr_Xr2_gate__Xr2) * var_t39)) - ((var_t81 * var_t82) * var_t89)) - var_t113) - var_t128) + var_t139) + var_t165) - var_t175) - var_t213) + var_t226) - (NV_Ith_S(mParameters, 10) * var_chaste_interface__i_f_Xf_gate__Xf)) - NV_Ith_S(mParameters, 6)) - NV_Ith_S(mParameters, 5)));
        DENSE_ELEM(rJacobian, 1, 0) = var_chaste_interface__environment__fake_dt * (0.001 * ((((((56497.175130000003 / var_t358) / var_t357) * var_t356) * var_t379) - (((200000.0 * var_t383) * var_t364) * var_t378)) - ((((1000.0 * var_t383) * var_t365) / pow(var_t377, 2.0)) * ((( -20.0 / pow(var_t368, 2.0)) * var_t367) - ((0.5 / pow(var_t374, 2.0)) * var_t373)))));
        DENSE_ELEM(rJacobian, 2, 0) = var_chaste_interface__environment__fake_dt * (0.001 * ((( -87.719298249999994 * var_t411) * var_t437) - (((var_t440 - var_chaste_interface__i_Na_h_gate__h) / pow(var_t436, 2.0)) * (var_t445 ? (( -1.5 / pow(((57.0 * var_t415) + (2700.0 * var_t448)) + (310000000.0 * var_t422), 2.0)) * ((( -8382.3529409999992 * var_t415) + (213300.0 * var_t448)) + (108035000000.0 * var_t422))) : ( -0.040000000000000001 <= var_Membrane__Vm) ? 0.0 : NAN))));
        DENSE_ELEM(rJacobian, 3, 0) = var_chaste_interface__environment__fake_dt * (0.001 * ((( -12531.328320000001 * var_t411) * var_t501) + ((142.85714285714286 * (var_t440 - var_chaste_interface__i_Na_j_gate__j)) * ((var_t445 ? (((((( -6214603.2000000002 * var_t466) + (0.00030508667999999998 * var_t469)) * var_t509) * var_t515) + ((1000.0 * var_t518) * var_t515)) - ((((311.0 * var_t518) * var_t509) / pow(var_t514, 2.0)) * var_t513)) : ( -0.040000000000000001 <= var_Membrane__Vm) ? 0.0 : NAN) + ((var_Membrane__Vm <=  -0.040000000000000001) ? ((( -0.25500479999999998 * var_t530) / var_t533) + (((3.3402720000000001 * var_t530) / pow(var_t533, 2.0)) * var_t532)) : ( -0.040000000000000001 < var_Membrane__Vm) ? (((34.200000000000003 * var_t544) / var_t548) + (((60.0 * var_t544) / pow(var_t548, 2.0)) * var_t547)) : NAN)))));
        DENSE_ELEM(rJacobian, 4, 0) = var_chaste_interface__environment__fake_dt * (0.001 * ((((142.85714285714286 / pow(var_t566, 2.0)) * var_t565) * var_t590) - ((((1.0 / var_t566) - var_chaste_interface__i_CaL_d_gate__d) / pow(var_t589, 2.0)) * (((((0.15076923079999999 / pow(var_t573, 2.0)) * var_t572) * var_t580) - (((0.28000000000000003 * var_t576) / pow(var_t579, 2.0)) * var_t578)) + ((0.050000000000000003 / pow(var_t586, 2.0)) * var_t585)))));
        DENSE_ELEM(rJacobian, 5, 0) = var_chaste_interface__environment__fake_dt * (0.001 * (((( -333333.33333333331 / pow(var_t618, 2.0)) * var_t617) * var_t648) - ((((1000.0 * var_t642) / pow(var_t639, 2.0)) * var_t647) * ((((( -19600.0 * var_t623) * var_t622) * var_t626) + ((20000.0 / pow(var_t631, 2.0)) * var_t630)) - (18000.0 * var_t663)))));
        DENSE_ELEM(rJacobian, 6, 0) = var_chaste_interface__environment__fake_dt * (0.001 * (((( -167.5 / pow(var_t690, 2.0)) * var_t689) * var_t706) - ((((0.33000000000000002 + (0.67000000000000004 / var_t690)) - var_chaste_interface__i_CaL_f2_gate__f2) / pow(var_t705, 2.0)) * ((((0.59999999999999998 * (( -11764.705882352941 * var_Membrane__Vm) - 294.11764705882354)) * var_t697) + ((3.1000000000000001 / pow(var_t701, 2.0)) * var_t700)) - (1.6000000000000001 * var_t663)))));
        DENSE_ELEM(rJacobian, 7, 0) = var_chaste_interface__environment__fake_dt * (0.001 * var_t753);
        DENSE_ELEM(rJacobian, 8, 0) = var_chaste_interface__environment__fake_dt * (0.001 * (((((75.585789890000001 / pow(var_t796, 2.0)) * var_t795) * var_t807) - (((37.037037037037038 * var_t811) * var_t801) * var_t806)) + (((32.20611916 * var_t811) * var_t802) * var_t805)));
        DENSE_ELEM(rJacobian, 9, 0) = var_chaste_interface__environment__fake_dt * (0.001 * ((((( -5952.3809520000004 / pow(var_t823, 2.0)) * var_t822) * var_t833) - (((14880.952380000001 * var_t837) * var_t828) * var_t832)) + (((14880.952380000001 * var_t837) * var_t829) * var_t831)));
        DENSE_ELEM(rJacobian, 10, 0) = var_chaste_interface__environment__fake_dt * (0.001 * (((((56.81818181818182 / pow(var_t849, 2.0)) * var_t848) * var_t858) - ((((75.757575757575751 * var_t862) / var_t857) * var_t832) * var_t855)) + (((45.454545454545453 * var_t862) * var_t857) * var_t831)));
        DENSE_ELEM(rJacobian, 11, 0) = var_chaste_interface__environment__fake_dt * (0.001 * (((( -105.26315789473684 / pow(var_t875, 2.0)) * var_t874) * var_t881) + ((52.631578947368418 * ((1.0 / var_t875) - var_chaste_interface__i_f_Xf_gate__Xf)) * var_t880)));
        DENSE_ELEM(rJacobian, 12, 0) = var_chaste_interface__environment__fake_dt * (0.001 * (((( -76.92307692307692 / pow(var_t891, 2.0)) * var_t890) * var_t906) + ((((0.039101999999999998 * ((1.0 / var_t891) - var_chaste_interface__i_to_q_gate__q)) / pow(var_t905, 2.0)) / pow(var_t902, 2.0)) * (( -45.600000000000001 * var_t897) + (6.5 * var_t900)))));
        DENSE_ELEM(rJacobian, 13, 0) = var_chaste_interface__environment__fake_dt * (0.001 * ((((53.333333330000002 / pow(var_t926, 2.0)) * var_t925) * var_t942) + ((((0.01440516 * ((1.0 / var_t926) - var_chaste_interface__i_to_r_gate__r)) / pow(var_t941, 2.0)) / pow(var_t938, 2.0)) * ((93.329999999999998 * var_t932) - (44.280000000000001 * var_t936)))));
        DENSE_ELEM(rJacobian, 14, 0) = var_chaste_interface__environment__fake_dt * (0.001 * ((( -1e+18 * NV_Ith_S(mParameters, 8)) * ((((var_t175 + NV_Ith_S(mParameters, 6)) - (3.0 * var_t165)) + (3.0 * var_t213)) - (3.0 * var_t226))) * 1.1777578010268393e-09));
        DENSE_ELEM(rJacobian, 15, 0) = var_chaste_interface__environment__fake_dt * (0.001 * ((( -5e+17 * var_t1007) * (((((var_t113 + var_t128) - var_t139) + NV_Ith_S(mParameters, 5)) - (2.0 * var_t213)) + (2.0 * var_t226))) * var_t1013));
        DENSE_ELEM(rJacobian, 17, 0) = var_chaste_interface__environment__fake_dt * (0.001 * var_t1178);
        DENSE_ELEM(rJacobian, 0, 1) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * (((( -3.0 * var_t170) * var_t171) * var_t230) * var_t236)));
        DENSE_ELEM(rJacobian, 1, 1) = var_chaste_interface__environment__fake_dt * ( -1000.0 * var_t379);
        DENSE_ELEM(rJacobian, 14, 1) = var_chaste_interface__environment__fake_dt * ((((((( -3e+18 * var_t968) * NV_Ith_S(mParameters, 9)) * var_t171) * var_t230) * var_t236) * 1.0364268649036186e-05) * 0.00011363636363636364);
        DENSE_ELEM(rJacobian, 0, 2) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * ((((-var_t170) * var_t172) * var_chaste_interface__i_Na_j_gate__j) * var_t236)));
        DENSE_ELEM(rJacobian, 2, 2) = var_chaste_interface__environment__fake_dt * (-var_t437);
        DENSE_ELEM(rJacobian, 14, 2) = var_chaste_interface__environment__fake_dt * (((( -1e+18 * var_t977) * var_chaste_interface__i_Na_j_gate__j) * var_t236) * 1.1777578010268393e-09);
        DENSE_ELEM(rJacobian, 0, 3) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * (((-var_t170) * var_t173) * var_t236)));
        DENSE_ELEM(rJacobian, 3, 3) = var_chaste_interface__environment__fake_dt * ( -142.85714285714286 * var_t501);
        DENSE_ELEM(rJacobian, 14, 3) = var_chaste_interface__environment__fake_dt * (((( -1e+18 * var_t977) * var_chaste_interface__i_Na_h_gate__h) * var_t236) * 1.1777578010268393e-09);
        DENSE_ELEM(rJacobian, 0, 4) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * ((( -4.0 * var_t247) * var_t248) * var_t110)));
        DENSE_ELEM(rJacobian, 4, 4) = var_chaste_interface__environment__fake_dt * (-var_t590);
        DENSE_ELEM(rJacobian, 15, 4) = var_chaste_interface__environment__fake_dt * (var_t1007 * (((0.041099999999999998 * var_t1022) * var_chaste_interface__calcium_dynamics__g) - (((((((((2e+18 * var_t114) * 96485.341499999995) * var_t104) * var_t107) * var_chaste_interface__i_CaL_f1_gate__f1) * var_chaste_interface__i_CaL_f2_gate__f2) * var_chaste_interface__i_CaL_fCa_gate__fCa) * NV_Ith_S(mParameters, 8)) * 0.00011363636363636364)));
        DENSE_ELEM(rJacobian, 16, 4) = var_chaste_interface__environment__fake_dt * ((( -0.041099999999999998 * var_t1134) * var_t1136) * var_chaste_interface__calcium_dynamics__g);
        DENSE_ELEM(rJacobian, 0, 5) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * ((((( -4.0 * var_t247) * var_t248) * var_chaste_interface__i_CaL_d_gate__d) * var_chaste_interface__i_CaL_f2_gate__f2) * var_chaste_interface__i_CaL_fCa_gate__fCa)));
        DENSE_ELEM(rJacobian, 5, 5) = var_chaste_interface__environment__fake_dt * (( -1000.0 * var_t648) - (((1000.0 * var_t670) * var_t672) * 0.0));
        DENSE_ELEM(rJacobian, 15, 5) = var_chaste_interface__environment__fake_dt * (((( -2e+18 * var_t1038) * var_t1039) * var_t135) * var_t1012);
        DENSE_ELEM(rJacobian, 0, 6) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * (((( -4.0 * var_t247) * var_t248) * var_t257) * var_chaste_interface__i_CaL_fCa_gate__fCa)));
        DENSE_ELEM(rJacobian, 6, 6) = var_chaste_interface__environment__fake_dt * (-var_t706);
        DENSE_ELEM(rJacobian, 15, 6) = var_chaste_interface__environment__fake_dt * ((((( -2e+18 * var_t1038) * var_t1039) * var_chaste_interface__i_CaL_f1_gate__f1) * var_chaste_interface__i_CaL_fCa_gate__fCa) * var_t1012);
        DENSE_ELEM(rJacobian, 0, 7) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * (((( -4.0 * var_t247) * var_t248) * var_t257) * var_chaste_interface__i_CaL_f2_gate__f2)));
        DENSE_ELEM(rJacobian, 7, 7) = var_chaste_interface__environment__fake_dt * (var_t753 - (var_t754 * 500.0));
        DENSE_ELEM(rJacobian, 15, 7) = var_chaste_interface__environment__fake_dt * (((( -2e+18 * var_t1038) * var_t1039) * var_t109) * var_t1012);
        DENSE_ELEM(rJacobian, 0, 8) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * (((( -0.43033148290000001 * var_t71) * var_t38) * var_chaste_interface__i_Kr_Xr2_gate__Xr2) * var_t39)));
        DENSE_ELEM(rJacobian, 8, 8) = var_chaste_interface__environment__fake_dt * ( -0.37037037037037035 * var_t807);
        DENSE_ELEM(rJacobian, 0, 9) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * (((( -0.43033148290000001 * var_t71) * var_t38) * var_chaste_interface__i_Kr_Xr1_gate__Xr1) * var_t39)));
        DENSE_ELEM(rJacobian, 9, 9) = var_chaste_interface__environment__fake_dt * ( -297.61904759999999 * var_t833);
        DENSE_ELEM(rJacobian, 0, 10) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * (((( -2.0 * var_t81) * var_t283) * var_chaste_interface__i_Ks_Xs_gate__Xs) * var_t89)));
        DENSE_ELEM(rJacobian, 10, 10) = var_chaste_interface__environment__fake_dt * ( -0.90909090909090906 * var_t858);
        DENSE_ELEM(rJacobian, 0, 11) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * ((-NV_Ith_S(mParameters, 10)) * (var_Membrane__Vm -  -0.017000000000000001))));
        DENSE_ELEM(rJacobian, 11, 11) = var_chaste_interface__environment__fake_dt * ( -0.52631578947368418 * var_t881);
        DENSE_ELEM(rJacobian, 0, 12) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * ((-var_t290) * var_chaste_interface__i_to_r_gate__r)));
        DENSE_ELEM(rJacobian, 12, 12) = var_chaste_interface__environment__fake_dt * (-var_t906);
        DENSE_ELEM(rJacobian, 0, 13) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * ((-var_t290) * var_chaste_interface__i_to_q_gate__q)));
        DENSE_ELEM(rJacobian, 13, 13) = var_chaste_interface__environment__fake_dt * (-var_t942);
        DENSE_ELEM(rJacobian, 0, 14) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * ((((((((((((-var_t81) * var_t2) * 1.0364268649036186e-05) * var_t278) * 0.029999999999999999) * var_t82) * var_t89) - var_t302) + var_t307) - var_t313) - (3.0 * var_t319)) - var_t322)));
        DENSE_ELEM(rJacobian, 14, 14) = var_chaste_interface__environment__fake_dt * ((( -1e+18 * NV_Ith_S(mParameters, 8)) * ((((var_t313 + var_t322) + (3.0 * var_t302)) - (3.0 * var_t307)) + (9.0 * var_t319))) * 1.1777578010268393e-09);
        DENSE_ELEM(rJacobian, 15, 14) = var_chaste_interface__environment__fake_dt * (((((((3e+18 * var_t1007) * NV_Ith_S(mParameters, 14)) * var_t181) * var_t183) * NV_Ith_S(mParameters, 1)) * var_t212) * var_t1013);
        DENSE_ELEM(rJacobian, 0, 15) = var_chaste_interface__environment__fake_dt * (mSetVoltageDerivativeToZero ? 0.0 : (1000.0 * ((((((((((( -5.4447296619999995e-07 * var_t81) * var_t283) * var_t82) / pow(var_t86, 2.0)) * pow(var_t83, 0.40000000000000002)) * var_t330) - var_t336) + var_t341) - var_t344) + var_t348) - var_t352)));
        DENSE_ELEM(rJacobian, 5, 15) = var_chaste_interface__environment__fake_dt * ((( -1000.0 * var_t670) * var_t672) * (var_t643 ? 1433.0 : 0.0));
        DENSE_ELEM(rJacobian, 7, 15) = var_chaste_interface__environment__fake_dt * (var_t753 + ((var_t754 * (((( -3.620396362e+26 / pow(var_t731, 2.0)) * var_t760) - ((760.10945579999998 / pow(var_t737, 2.0)) * var_t736)) - ((285.04104589999997 / pow(var_t743, 2.0)) * var_t742))) * 500.0));
        DENSE_ELEM(rJacobian, 14, 15) = var_chaste_interface__environment__fake_dt * ((((((((3e+18 * NV_Ith_S(mParameters, 8)) * NV_Ith_S(mParameters, 14)) * var_t196) * 2.8571431999999999) * var_t208) * var_t211) * 1.0364268649036186e-05) * 0.00011363636363636364);
        DENSE_ELEM(rJacobian, 15, 15) = var_chaste_interface__environment__fake_dt * ((((((2.0 / pow(var_t1006, 2.0)) * (((var_t1063 - var_t1068) + (0.041099999999999998 * var_t1070)) - (((5e+17 * ((((((4.0 * var_t247) * var_t1039) * var_t110) + (NV_Ith_S(mParameters, 5) * (var_Membrane__Vm - (((0.5 * var_t2) * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 1) * var_t83))))) + (var_t345 * var_t343)) - ((2.0 * var_t217) * var_t212))) * NV_Ith_S(mParameters, 8)) * 1.1777578010268393e-09))) * 0.00025000000000000001) / var_t1003) / var_t1002) + (var_t1007 * (( -0.00044443999999999999 - var_t1102) - (((5e+17 * ((((var_t336 + var_t352) + var_t344) - var_t348) + (2.0 * var_t341))) * NV_Ith_S(mParameters, 8)) * 1.1777578010268393e-09))));
        DENSE_ELEM(rJacobian, 16, 15) = var_chaste_interface__environment__fake_dt * ((var_t1134 * 0.0017131207921470542) * (var_t1102 + 0.00044443999999999999));
        DENSE_ELEM(rJacobian, 17, 15) = var_chaste_interface__environment__fake_dt * (var_t1178 + ((var_t1179 * ((var_chaste_interface__calcium_dynamics__Cai <= 0.00035) ? ((( -3.2639461439999999e+21 / pow(1.0 + var_t1164, 2.0)) * var_t728) * var_chaste_interface__calcium_dynamics__Cai) : (0.00035 < var_chaste_interface__calcium_dynamics__Cai) ? ((( -3.1552318180000002e+56 / pow(1.0 + var_t1168, 2.0)) * var_t729) * var_t760) : NAN)) * 500.0));
        DENSE_ELEM(rJacobian, 15, 16) = var_chaste_interface__environment__fake_dt * (var_t1007 * var_t1122);
        DENSE_ELEM(rJacobian, 16, 16) = var_chaste_interface__environment__fake_dt * (((((((((2.0 / pow(var_t1132, 2.0)) * 8800.0) * 0.0017131207921470542) * ((var_t1068 - (0.041099999999999998 * var_t1070)) - var_t1063)) * 10.0) * 0.29999999999999999) / var_t1129) / var_t1128) - ((var_t1134 * 0.0017131207921470542) * var_t1122));
        DENSE_ELEM(rJacobian, 15, 17) = var_chaste_interface__environment__fake_dt * (((0.041099999999999998 * var_t1007) * var_t1022) * var_chaste_interface__i_CaL_d_gate__d);
        DENSE_ELEM(rJacobian, 16, 17) = var_chaste_interface__environment__fake_dt * ((( -0.041099999999999998 * var_t1134) * var_t1136) * var_chaste_interface__i_CaL_d_gate__d);
        DENSE_ELEM(rJacobian, 17, 17) = var_chaste_interface__environment__fake_dt * (var_t1178 - (var_t1179 * 500.0));
    }

*********************************************************************************************************************************************************************************/
    
    N_Vector Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode::ComputeDerivedQuantities(double var_chaste_interface__environment__time, const N_Vector & rY)
    {
        // Inputs:
        // Time units: millisecond
        double var_chaste_interface__Membrane__Vm = (mSetVoltageDerivativeToZero ? this->mFixedVoltage : NV_Ith_S(rY, 0));
        // Units: millivolt; Initial value: -74.3340057624
        double var_chaste_interface__i_Na_m_gate__m = NV_Ith_S(rY, 1);
        // Units: dimensionless; Initial value: 0.102953468725004
        double var_chaste_interface__i_Na_h_gate__h = NV_Ith_S(rY, 2);
        // Units: dimensionless; Initial value: 0.786926637881461
        double var_chaste_interface__i_Na_j_gate__j = NV_Ith_S(rY, 3);
        // Units: dimensionless; Initial value: 0.253943221774722
        double var_chaste_interface__i_CaL_d_gate__d = NV_Ith_S(rY, 4);
        // Units: dimensionless; Initial value: 8.96088425225182e-5
        double var_chaste_interface__i_CaL_f1_gate__f1 = NV_Ith_S(rY, 5);
        // Units: dimensionless; Initial value: 0.970411811263976
        double var_chaste_interface__i_CaL_f2_gate__f2 = NV_Ith_S(rY, 6);
        // Units: dimensionless; Initial value: 0.999965815466749
        double var_chaste_interface__i_CaL_fCa_gate__fCa = NV_Ith_S(rY, 7);
        // Units: dimensionless; Initial value: 0.998925296531804
        double var_chaste_interface__i_Kr_Xr1_gate__Xr1 = NV_Ith_S(rY, 8);
        // Units: dimensionless; Initial value: 0.00778547011240132
        double var_chaste_interface__i_Kr_Xr2_gate__Xr2 = NV_Ith_S(rY, 9);
        // Units: dimensionless; Initial value: 0.432162576531617
        double var_chaste_interface__i_Ks_Xs_gate__Xs = NV_Ith_S(rY, 10);
        // Units: dimensionless; Initial value: 0.0322944866983666
        double var_chaste_interface__i_f_Xf_gate__Xf = NV_Ith_S(rY, 11);
        // Units: dimensionless; Initial value: 0.100615100568753
        double var_chaste_interface__i_to_q_gate__q = NV_Ith_S(rY, 12);
        // Units: dimensionless; Initial value: 0.839295925773219
        double var_chaste_interface__i_to_r_gate__r = NV_Ith_S(rY, 13);
        // Units: dimensionless; Initial value: 0.00573289893326379
        //double var_chaste_interface__sodium_dynamics__Nai = NV_Ith_S(rY, 14);
        double var_chaste_interface__sodium_dynamics__Nai = (mSetNaiDerivativeToZero ? this->mFixedNai : NV_Ith_S(rY, 14));
        // Units: millimolar; Initial value: 10.9248496211574
        //double var_chaste_interface__calcium_dynamics__Cai = NV_Ith_S(rY, 15);
        double var_chaste_interface__calcium_dynamics__Cai = (mSetCaiDerivativeToZero ? this->mFixedCai : NV_Ith_S(rY, 15));
        // Units: millimolar; Initial value: 1.80773974140477e-5
        
        
        // Mathematics
        const double var_i_CaL__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_CaL__i_CaL = ((((NV_Ith_S(mParameters, 4) * (4.0 * (var_i_CaL__Vm * 9309421124.3716221))) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * ((var_chaste_interface__calcium_dynamics__Cai * exp((2.0 * (var_i_CaL__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17)))) - (0.34100000000000003 * NV_Ith_S(mParameters, 1)))) / (exp((2.0 * (var_i_CaL__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) - 1.0)) * (var_chaste_interface__i_CaL_d_gate__d * (var_chaste_interface__i_CaL_f1_gate__f1 * (var_chaste_interface__i_CaL_f2_gate__f2 * var_chaste_interface__i_CaL_fCa_gate__fCa))); // A_per_F
        const double var_i_K1__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_electric_potentials__E_K = ((8.3144720000000003 * NV_Ith_S(mParameters, 17)) * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 2) / NV_Ith_S(mParameters, 0)); // volt
        const double var_i_K1__alpha_K1 = 3.9100000000000001 / (1.0 + exp(0.59419999999999995 * (((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0)) - 200.0))); // dimensionless
        const double var_i_K1__i_K1 = NV_Ith_S(mParameters, 11) * ((var_i_K1__alpha_K1 / (var_i_K1__alpha_K1 + ((( -1.5089999999999999 * exp(0.00020000000000000001 * (((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0)) + 100.0))) + exp(0.58860000000000001 * (((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0)) - 10.0))) / (1.0 + exp(0.45469999999999999 * ((var_i_K1__Vm * 1000.0) - (var_electric_potentials__E_K * 1000.0))))))) * ((var_i_K1__Vm - var_electric_potentials__E_K) * sqrt(NV_Ith_S(mParameters, 2) * 0.18518518518518517))); // A_per_F
        const double var_i_f__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_f__i_f = NV_Ith_S(mParameters, 10) * (var_chaste_interface__i_f_Xf_gate__Xf * (var_i_f__Vm -  -0.017000000000000001)); // A_per_F
        const double var_electric_potentials__E_Na = ((8.3144720000000003 * NV_Ith_S(mParameters, 17)) * 1.0364268649036186e-05) * log(NV_Ith_S(mParameters, 3) / var_chaste_interface__sodium_dynamics__Nai); // volt
        const double var_i_Na__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_Na__i_Na = 1.0 * (NV_Ith_S(mParameters, 9) * (pow(var_chaste_interface__i_Na_m_gate__m, 3.0) * (var_chaste_interface__i_Na_h_gate__h * (var_chaste_interface__i_Na_j_gate__j * (var_i_Na__Vm - var_electric_potentials__E_Na))))); // A_per_F
        const double var_i_Kr__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_Kr__i_Kr = 1.0 * (NV_Ith_S(mParameters, 12) * ((var_i_Kr__Vm - var_electric_potentials__E_K) * (var_chaste_interface__i_Kr_Xr1_gate__Xr1 * (var_chaste_interface__i_Kr_Xr2_gate__Xr2 * sqrt(NV_Ith_S(mParameters, 2) * 0.18518518518518517))))); // A_per_F
        const double var_i_Ks__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_Ks__i_Ks = 1.0 * (NV_Ith_S(mParameters, 13) * ((var_i_Ks__Vm - (((8.3144720000000003 * NV_Ith_S(mParameters, 17)) * 1.0364268649036186e-05) * log((NV_Ith_S(mParameters, 2) + (0.029999999999999999 * NV_Ith_S(mParameters, 3))) / (NV_Ith_S(mParameters, 0) + (0.029999999999999999 * var_chaste_interface__sodium_dynamics__Nai))))) * (pow(var_chaste_interface__i_Ks_Xs_gate__Xs, 2.0) * (1.0 + (0.59999999999999998 / (1.0 + pow(3.8000000000000002e-05 / var_chaste_interface__calcium_dynamics__Cai, 1.3999999999999999))))))); // A_per_F
        const double var_i_to__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_to__i_to = NV_Ith_S(mParameters, 16) * ((var_i_to__Vm - var_electric_potentials__E_K) * (var_chaste_interface__i_to_q_gate__q * var_chaste_interface__i_to_r_gate__r)); // A_per_F
        const double var_i_PCa__i_PCa = (NV_Ith_S(mParameters, 7) * var_chaste_interface__calcium_dynamics__Cai) / (var_chaste_interface__calcium_dynamics__Cai + 0.00050000000000000001); // A_per_F
        const double var_i_NaK__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_NaK__i_NaK = ((((NV_Ith_S(mParameters, 15) * NV_Ith_S(mParameters, 2)) / (NV_Ith_S(mParameters, 2) + 1.0)) * var_chaste_interface__sodium_dynamics__Nai) / (var_chaste_interface__sodium_dynamics__Nai + 40.0)) / (1.0 + ((0.1245 * exp(( -0.10000000000000001 * (var_i_NaK__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17)))) + (0.035299999999999998 * exp(((-var_i_NaK__Vm) * 96485.341499999995) / (8.3144720000000003 * NV_Ith_S(mParameters, 17)))))); // A_per_F
        const double var_i_NaCa__Vm = 0.001 * var_chaste_interface__Membrane__Vm; // volt
        const double var_i_NaCa__i_NaCa = (NV_Ith_S(mParameters, 14) * ((exp((0.34999999999999998 * (var_i_NaCa__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * (pow(var_chaste_interface__sodium_dynamics__Nai, 3.0) * NV_Ith_S(mParameters, 1))) - (exp(( -0.65000000000000002 * (var_i_NaCa__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))) * (pow(NV_Ith_S(mParameters, 3), 3.0) * (var_chaste_interface__calcium_dynamics__Cai * 2.8571431999999999))))) / ((669921.875 + pow(NV_Ith_S(mParameters, 3), 3.0)) * ((1.3799999999999999 + NV_Ith_S(mParameters, 1)) * (1.0 + (0.10000000000000001 * exp(( -0.65000000000000002 * (var_i_NaCa__Vm * 96485.341499999995)) / (8.3144720000000003 * NV_Ith_S(mParameters, 17))))))); // A_per_F
        const double var_chaste_interface__stim_mode__i_stim = -GetIntracellularAreaStimulus(var_chaste_interface__environment__time);
        const double var_chaste_interface__i_NaK__i_NaK = var_i_NaK__i_NaK; // A_per_F
        const double var_chaste_interface__i_K1__i_K1 = var_i_K1__i_K1; // A_per_F
        const double var_chaste_interface__i_f__i_f = var_i_f__i_f; // A_per_F
        const double var_chaste_interface__i_Na__i_Na = var_i_Na__i_Na; // A_per_F
        const double var_chaste_interface__i_to__i_to = var_i_to__i_to; // A_per_F
        const double var_chaste_interface__i_Kr__i_Kr = var_i_Kr__i_Kr; // A_per_F
        const double var_chaste_interface__i_Ks__i_Ks = var_i_Ks__i_Ks; // A_per_F
        const double var_chaste_interface__i_CaL__i_CaL = var_i_CaL__i_CaL; // A_per_F
        const double var_chaste_interface__i_PCa__i_PCa = var_i_PCa__i_PCa; // A_per_F
        const double var_chaste_interface__i_NaCa__i_NaCa = var_i_NaCa__i_NaCa; // A_per_F
        
        N_Vector dqs = N_VNew_Serial(12);
        NV_Ith_S(dqs, 0) = var_chaste_interface__i_CaL__i_CaL;
        NV_Ith_S(dqs, 1) = var_chaste_interface__i_PCa__i_PCa;
        NV_Ith_S(dqs, 2) = var_chaste_interface__i_Na__i_Na;
        NV_Ith_S(dqs, 3) = var_chaste_interface__i_f__i_f;
        NV_Ith_S(dqs, 4) = var_chaste_interface__i_K1__i_K1;
        NV_Ith_S(dqs, 5) = var_chaste_interface__i_Kr__i_Kr;
        NV_Ith_S(dqs, 6) = var_chaste_interface__i_Ks__i_Ks;
        NV_Ith_S(dqs, 7) = var_chaste_interface__i_NaCa__i_NaCa;
        NV_Ith_S(dqs, 8) = var_chaste_interface__i_NaK__i_NaK;
        NV_Ith_S(dqs, 9) = var_chaste_interface__stim_mode__i_stim;
        NV_Ith_S(dqs, 10) = var_chaste_interface__i_to__i_to;
        NV_Ith_S(dqs, 11) = var_chaste_interface__environment__time;
        return dqs;
    }
    
template<>
void OdeSystemInformation<Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode>::Initialise(void)
{
    this->mSystemName = "paci_hyttinen_aaltosetala_severi_ventricularVersion";
    this->mFreeVariableName = "time";
    this->mFreeVariableUnits = "millisecond";
    
    this->mVariableNames.push_back("membrane_voltage");
    this->mVariableUnits.push_back("millivolt");
    this->mInitialConditions.push_back(-74.3340057624);

    this->mVariableNames.push_back("membrane_fast_sodium_current_m_gate");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.102953468725004);

    this->mVariableNames.push_back("membrane_fast_sodium_current_h_gate");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.786926637881461);

    this->mVariableNames.push_back("membrane_fast_sodium_current_j_gate");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.253943221774722);

    this->mVariableNames.push_back("membrane_L_type_calcium_current_d_gate");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(8.96088425225182e-5);

    this->mVariableNames.push_back("membrane_L_type_calcium_current_f_gate");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.970411811263976);

    this->mVariableNames.push_back("membrane_L_type_calcium_current_f2_gate");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.999965815466749);

    this->mVariableNames.push_back("membrane_L_type_calcium_current_fCa_gate");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.998925296531804);

    this->mVariableNames.push_back("i_Kr_Xr1_gate__Xr1");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00778547011240132);

    this->mVariableNames.push_back("i_Kr_Xr2_gate__Xr2");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.432162576531617);

    this->mVariableNames.push_back("i_Ks_Xs_gate__Xs");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0322944866983666);

    this->mVariableNames.push_back("membrane_hyperpolarisation_activated_funny_current_single_gate");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.100615100568753);

    this->mVariableNames.push_back("i_to_q_gate__q");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.839295925773219);

    this->mVariableNames.push_back("i_to_r_gate__r");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.00573289893326379);

    this->mVariableNames.push_back("cytosolic_sodium_concentration");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(10.9248496211574);

    this->mVariableNames.push_back("cytosolic_calcium_concentration");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(1.80773974140477e-5);

    this->mVariableNames.push_back("JSR_calcium_concentration");
    this->mVariableUnits.push_back("millimolar");
    this->mInitialConditions.push_back(-0.2734234751931);

    this->mVariableNames.push_back("calcium_dynamics__g");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.999999981028517);

    this->mParameterNames.push_back("cytosolic_potassium_concentration");
    this->mParameterUnits.push_back("millimolar");
    
    this->mParameterNames.push_back("extracellular_calcium_concentration");
    this->mParameterUnits.push_back("millimolar");
    
    this->mParameterNames.push_back("extracellular_potassium_concentration");
    this->mParameterUnits.push_back("millimolar");
    
    this->mParameterNames.push_back("extracellular_sodium_concentration");
    this->mParameterUnits.push_back("millimolar");
    
    this->mParameterNames.push_back("membrane_L_type_calcium_current_conductance");
    this->mParameterUnits.push_back("metre_cube_per_F_per_s");
    
    this->mParameterNames.push_back("membrane_background_calcium_current_conductance");
    this->mParameterUnits.push_back("S_per_F");
    
    this->mParameterNames.push_back("membrane_background_sodium_current_conductance");
    this->mParameterUnits.push_back("S_per_F");
    
    this->mParameterNames.push_back("membrane_calcium_pump_current_conductance");
    this->mParameterUnits.push_back("A_per_F");
    
    this->mParameterNames.push_back("membrane_capacitance");
    this->mParameterUnits.push_back("farad");
    
    this->mParameterNames.push_back("membrane_fast_sodium_current_conductance");
    this->mParameterUnits.push_back("S_per_F");
    
    this->mParameterNames.push_back("membrane_hyperpolarisation_activated_funny_current_potassium_component_conductance");
    this->mParameterUnits.push_back("S_per_F");
    
    this->mParameterNames.push_back("membrane_inward_rectifier_potassium_current_conductance");
    this->mParameterUnits.push_back("S_per_F");
    
    this->mParameterNames.push_back("membrane_rapid_delayed_rectifier_potassium_current_conductance");
    this->mParameterUnits.push_back("S_per_F");
    
    this->mParameterNames.push_back("membrane_slow_delayed_rectifier_potassium_current_conductance");
    this->mParameterUnits.push_back("S_per_F");
    
    this->mParameterNames.push_back("membrane_sodium_calcium_exchanger_current_conductance");
    this->mParameterUnits.push_back("A_per_F");
    
    this->mParameterNames.push_back("membrane_sodium_potassium_pump_current_permeability");
    this->mParameterUnits.push_back("A_per_F");
    
    this->mParameterNames.push_back("membrane_transient_outward_current_conductance");
    this->mParameterUnits.push_back("S_per_F");
    
    this->mParameterNames.push_back("temperature");
    this->mParameterUnits.push_back("kelvin");
    
    this->mDerivedQuantityNames.push_back("membrane_L_type_calcium_current");
    this->mDerivedQuantityUnits.push_back("A_per_F");
    
    this->mDerivedQuantityNames.push_back("membrane_calcium_pump_current");
    this->mDerivedQuantityUnits.push_back("A_per_F");
    
    this->mDerivedQuantityNames.push_back("membrane_fast_sodium_current");
    this->mDerivedQuantityUnits.push_back("A_per_F");
    
    this->mDerivedQuantityNames.push_back("membrane_hyperpolarisation_activated_funny_current");
    this->mDerivedQuantityUnits.push_back("A_per_F");
    
    this->mDerivedQuantityNames.push_back("membrane_inward_rectifier_potassium_current");
    this->mDerivedQuantityUnits.push_back("A_per_F");
    
    this->mDerivedQuantityNames.push_back("membrane_rapid_delayed_rectifier_potassium_current");
    this->mDerivedQuantityUnits.push_back("A_per_F");
    
    this->mDerivedQuantityNames.push_back("membrane_slow_delayed_rectifier_potassium_current");
    this->mDerivedQuantityUnits.push_back("A_per_F");
    
    this->mDerivedQuantityNames.push_back("membrane_sodium_calcium_exchanger_current");
    this->mDerivedQuantityUnits.push_back("A_per_F");
    
    this->mDerivedQuantityNames.push_back("membrane_sodium_potassium_pump_current");
    this->mDerivedQuantityUnits.push_back("A_per_F");
    
    this->mDerivedQuantityNames.push_back("membrane_stimulus_current");
    this->mDerivedQuantityUnits.push_back("uA_per_cm2");
    
    this->mDerivedQuantityNames.push_back("membrane_transient_outward_current");
    this->mDerivedQuantityUnits.push_back("A_per_F");
    
    this->mDerivedQuantityNames.push_back("time");
    this->mDerivedQuantityUnits.push_back("millisecond");
    
    this->mInitialised = true;
}


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode)
#endif // CHASTE_CVODE
