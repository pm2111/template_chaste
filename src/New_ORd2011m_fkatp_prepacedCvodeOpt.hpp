#ifdef CHASTE_CVODE
#ifndef CELLNEW_ORD2011M_FKATP_PREPACEDFROMCELLMLCVODEOPT_HPP_
#define CELLNEW_ORD2011M_FKATP_PREPACEDFROMCELLMLCVODEOPT_HPP_

//! @file
//! 
//! This source file was generated from CellML.
//! 
//! Model: ohara_rudy_2011
//! 
//! Processed by pycml - CellML Tools in Python
//!     (translators: , pycml: , optimize: )
//! on Fri Dec 21 16:57:41 2018
//! 
//! <autogenerated>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCvodeCell.hpp"
#include "AbstractStimulusFunction.hpp"

class CellNew_ORd2011m_fkatp_prepacedFromCellMLCvodeOpt : public AbstractCvodeCell
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
    CellNew_ORd2011m_fkatp_prepacedFromCellMLCvodeOpt(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver /* unused; should be empty */, boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);
    ~CellNew_ORd2011m_fkatp_prepacedFromCellMLCvodeOpt();
    AbstractLookupTableCollection* GetLookupTableCollection();
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);
    void EvaluateYDerivatives(double var_chaste_interface__environment__time, const N_Vector rY, N_Vector rDY);
};


// Needs to be included last
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellNew_ORd2011m_fkatp_prepacedFromCellMLCvodeOpt)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const CellNew_ORd2011m_fkatp_prepacedFromCellMLCvodeOpt * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }
        
        template<class Archive>
        inline void load_construct_data(
            Archive & ar, CellNew_ORd2011m_fkatp_prepacedFromCellMLCvodeOpt * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)CellNew_ORd2011m_fkatp_prepacedFromCellMLCvodeOpt(p_solver, p_stimulus);
        }
        
    }
    
}

#endif // CELLNEW_ORD2011M_FKATP_PREPACEDFROMCELLMLCVODEOPT_HPP_
#endif // CHASTE_CVODE
