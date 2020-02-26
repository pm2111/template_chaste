#ifdef CHASTE_CVODE
#ifndef CELLTORORD_FKATP_ENDO_MIDFROMCELLMLCVODEOPT_HPP_
#define CELLTORORD_FKATP_ENDO_MIDFROMCELLMLCVODEOPT_HPP_

//! @file
//! 
//! This source file was generated from CellML.
//! 
//! Model: Tomek_model13mid
//! 
//! Processed by pycml - CellML Tools in Python
//!     (translators: , pycml: , optimize: )
//! on Thu Feb 13 07:34:59 2020
//! 
//! <autogenerated>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCvodeCell.hpp"
#include "AbstractStimulusFunction.hpp"

class CellToRORd_fkatp_endo_midFromCellMLCvodeOpt : public AbstractCvodeCell
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCvodeCell >(*this);
    }
    
    // 
    // Settable parameters and readable variables
    // 
    
public:
    boost::shared_ptr<RegularStimulus> UseCellMLDefaultStimulus();
    double GetIntracellularCalciumConcentration();
    CellToRORd_fkatp_endo_midFromCellMLCvodeOpt(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver /* unused; should be empty */, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
    ~CellToRORd_fkatp_endo_midFromCellMLCvodeOpt();
    AbstractLookupTableCollection* GetLookupTableCollection();
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);
    void EvaluateYDerivatives(double var_chaste_interface__environment__time, const N_Vector rY, N_Vector rDY);
};


// Needs to be included last
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellToRORd_fkatp_endo_midFromCellMLCvodeOpt)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const CellToRORd_fkatp_endo_midFromCellMLCvodeOpt * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }
        
        template<class Archive>
        inline void load_construct_data(
            Archive & ar, CellToRORd_fkatp_endo_midFromCellMLCvodeOpt * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)CellToRORd_fkatp_endo_midFromCellMLCvodeOpt(p_solver, p_stimulus);
        }
        
    }
    
}

#endif // CELLTORORD_FKATP_ENDO_MIDFROMCELLMLCVODEOPT_HPP_
#endif // CHASTE_CVODE
