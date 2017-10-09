#ifdef CHASTE_CVODE
#ifndef CELLPACI_HYTTINEN_AALTOSETALA_SEVERI_VENTRICULARVERSION_ALLOWCLAMP_FROMCELLMLCVODE_HPP_
#define CELLPACI_HYTTINEN_AALTOSETALA_SEVERI_VENTRICULARVERSION_ALLOWCLAMP_FROMCELLMLCVODE_HPP_

//! @file
//! 
//! This source file was generated from CellML and modified by ChonL on Wed May 05 2017
//! 
//! Model: paci_hyttinen_aaltosetala_severi_ventricularVersion
//! 

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCvodeCell.hpp"
#include "AbstractStimulusFunction.hpp"

class Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode : public AbstractCvodeCell
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCvodeCell >(*this);
    }
    

    // 
    // allow clamp variables
    // 

private:
    bool mSetNaiDerivativeToZero;
    bool mSetCaiDerivativeToZero;
    double mFixedNai;
    double mFixedCai;
    //bool mSetCaSRDerivativeToZero;


    // 
    // Settable parameters and readable variables
    // 
    
public:
    boost::shared_ptr<RegularStimulus> UseCellMLDefaultStimulus();
    double GetIntracellularCalciumConcentration();
    Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver /* unused; should be empty */, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
    ~Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode();
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);
    void EvaluateYDerivatives(double var_chaste_interface__environment__time, const N_Vector rY, N_Vector rDY);
    //void EvaluateAnalyticJacobian(long int N, double var_chaste_interface__environment__time, N_Vector rY, N_Vector rDY, CHASTE_CVODE_DENSE_MATRIX rJacobian, N_Vector rTmp1, N_Vector rTmp2, N_Vector rTmp3);
    N_Vector ComputeDerivedQuantities(double var_chaste_interface__environment__time, const N_Vector & rY);

    // 
    // allow clamp variables
    // 
    void SetNaiDerivativeToZero(bool clamp);
    void SetCaiDerivativeToZero(bool clamp);

};


// Needs to be included last
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }
        
        template<class Archive>
        inline void load_construct_data(
            Archive & ar, Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)Cellpaci_hyttinen_aaltosetala_severi_ventricularVersion_allowClamp_fromCellMLCvode(p_solver, p_stimulus);
        }
        
    }
    
}

#endif // CELLPACI_HYTTINEN_AALTOSETALA_SEVERI_VENTRICULARVERSION_ALLOWCLAMP_FROMCELLMLCVODE_HPP_
#endif // CHASTE_CVODE
