#pragma once
#ifndef EA_HEADERS_EXECUTEVARIABLES_H_
#define EA_HEADERS_EXECUTEVARIABLES_H_

#include "variables.h"
/// <summary>
/// For each thread of EA instance, create a CC to prevent from using same global variables.
/// </summary>
namespace evolutionalgorithm {
    class ExecutiveVariables : public Variables{
    public:
        ExecutiveVariables();
        ExecutiveVariables(const Variables& v) : Variables(v) {}; // copy
        virtual ~ExecutiveVariables();
        void UpdateDE_CR(double);
        double GetDE_CR();
        void UpdateDE_NP(size_t);
        size_t GetDE_NP();
        void UpdateDE_NextGenNP(size_t);
        size_t GetDE_NextGenNP();
        void UpdateDE_Dim(unsigned);
        unsigned GetDE_Dim();
        void UpdateDE_F(double);
        double GetDE_F();
        void UpdateDE_K(double);
        double GetDE_K();
        void UpdateOut_B(bool);
        bool GetOut_B();
        size_t GetDE_PMin(); // Gets at least how many individuals in np should be choosed as pbest
        double GetDE_PMaxScale(); // Gets max proportion of np which can be choosed when executing DE/cur-to-pbest/n/bin



        size_t GetArchiveSize();
        void UpdateDE_P(double);
        double GetDE_P();
        double* GetBoundMin();
        double* GetBoundMax();
        double GetBoundMin(unsigned);
        void SetBoundMin(size_t, double);
        double GetBoundMax(unsigned);
        void SetBoundMax(size_t, double);
        size_t GetDE_NumOfStrategies();
        unsigned GetDE_MaxFE();
        void SetDE_MaxFE(size_t);
        unsigned GetDE_MaxGen();
        double GetDE_Optimal_();
        size_t GetRun();
        std::map<corematerial::MutantAndCrossover, size_t>* GetDE_StrategiesIndexMap();

        // SaDE
        size_t GetSaDE_LP();
        double GetSaDE_Epsilion();
        // END SaDE

        // SHADE
        size_t GetSHADE_H();
        // END SHADE

        // PAEHN
        size_t GetPAEHN_H();
        size_t GetPAEHN_SubpopulationSize();
        size_t GetPAEHN_StagnationUpper();
        // END PAEHN

        corematerial::Sample GetDE_SamplePattern();
        corematerial::InitialPopulation GetDE_InitialPopulationPattern();
        corematerial::Archive GetDE_ArchivePattern();
        corematerial::Memory GetDE_Memory_Pattern();
        corematerial::MemoryCR GetDE_MemoryCR_Pattern();
        corematerial::MemoryCR_Setting GetDE_MemoryCR_SettingPattern();
        corematerial::MemoryF GetDE_MemoryF_Pattern();
        corematerial::MutantAndCrossover GetDE_StrategiesPattern();
        corematerial::OutBoundHandler GetDE_OutBoundHandlerPattern();
        corematerial::InitialPopulation GetInitialPopulationPattern();
        corematerial::StopCriterion GetStopCriterionPattern();
        corematerial::StrategiesSelection GetDE_StrategiesSelectionPattern();
        corematerial::Selection GetDE_SelectionPattern();
        corematerial::ParametersControl GetParametersControlPattern();
        corematerial::ParametersControlMethod GetParametersControlMethodPattern();
        corematerial::Subpopulation GetSubpopulationPattern();
        corematerial::PAEHN GetPAEHN_Pattern();
    protected:
    private:
    };
}
#endif // EA_HEADERS_EXECUTEVARIABLES_H_