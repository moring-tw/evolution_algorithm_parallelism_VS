#pragma once
#ifndef EA_HEADERS_DIFFERENTIAL_EVOLUTION_VARIANTS_PAEHN_H_
#define EA_HEADERS_DIFFERENTIAL_EVOLUTION_VATIANTS_PAEHN_H_


#include <iostream>
#include <fstream>

#include "basic.h"
#include "../behavior.h"
#include "../localvariables.h"
#include "../../cells.h"
#include "../../publicfunction.h"
#include "../../subpopulation.h"
namespace evolutionalgorithm{
	namespace differentialevolution {
		class PAEHN : public Basic {
		public:
			PAEHN();
			~PAEHN();
			void NextRun();

			//ini
			void initialPAEHN_Archive();

			//evolution
			void RandomRi();
			void SelectSHADE_SubpopCR();
			void SelectSHADE_SubpopF();
			void SubpopulationMutationAndCrossover();

			//selection
			void CalculateSuccessDeltaFitness();
			void SuccessSelectionSHADE_SubpopMemoryCR();
			void SuccessSelectionSHADE_SubpopMemoryF();
			void SuccessSelectionCrowdingArchive(double* dimensions, const double fitness, bool random_promising_archive = false);
			void FormSubpopulation(bool enable_recursive_clearing = true);
			void FormNewPopulationNearArchiveIndividuals(bool reinitializing_crowding_promising_archive = true);
			void SubpopulationSelection(bool random_promising_archive = false);
			bool AreAllSubpopStagnation();
			
			
			void NextGen();

			//other
			void StagnationSubpop(const std::vector<paehn::Subpopulation>::iterator&);

			void ShowAllSubpopBest(std::vector<std::vector<double>>&);
			bool IsArchiveRepeat(const size_t base_index);
			static inline double ProduceCR(double mean) {
				double cr = PublicFunction::NormalDistributionRandom(mean, 0.1);
				if (cr > 1) cr = 1;
				else if (cr < 0) cr = 0;
				return cr;
			}

			static inline double ProduceF(double mean) {
				double f = PublicFunction::CauchyDistributionRandom(mean, 0.1);
				if (f > 1) f = 1;
				while (f <= 0) f = PublicFunction::CauchyDistributionRandom(mean, 0.1);
				return f;
			}
		protected:
		private:
			size_t FindNearestIndividualInArchive(double*); // returns a nearest Individual in Population'
			void CleanInTheRadius(size_t base_individual,double distance);
			void InitializeOfstream();
			void OutputPopulation();
			void OutputArchive();
			void OutputSubpop();
			void OutputBeforeAfterPopulation();
			void OutputBeforAfterArchive();
			std::vector<ofstream> file_streams_;
			size_t subpop_output_counter = 0;
			size_t no_convergence_counter_ = 0;
		};
	} // namespace differentialevolution
} // namespace evolutionalgorithm
#endif // !EA_HEADERS_DIFFERENTIAL_EVOLUTION_VARIANTS_PAEHN_H_
