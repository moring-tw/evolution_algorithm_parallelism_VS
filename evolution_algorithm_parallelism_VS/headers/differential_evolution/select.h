#pragma once
#ifndef EA_HEADERS_DIFFERENTIAL_EVOLUTION_SELECT_H_
#define EA_HEADERS_DIFFERENTIAL_EVOLUTION_SELECT_H_

#include "behavior.h"

namespace evolutionalgorithm {
	namespace differentialevolution {
		/// <summary>
		///		checks some variables are valid before invoking corresponding algorithms.
		///		This class is about selecting parameters and recording parameters.
		/// </summary>
		class Select : public Behavior
		{
		public:
			Select();
			~Select();
			void SelectDE();
			void AfterSelectDE();
			void BeforeSelectSaDE_StrategiesSelection();
			void AfterSelectSaDE_SuccessStrategiesSelection();
			void AfterSelectSaDE_FailStrategiesSelection();
			void AfterSelectSaDE_SuccessMemoryCR();
			
			void AfterSelectSHADE_SuccessMemoryCR();
			void AfterSelectSHADE_SuccessMemoryF();
			void AfterSelectSHADE_SuccessArchive();

			void AfterSelectPAEHN_SuccessArchive();
			void AfterSelectPAEHN_FormSubpopulation(bool enable_recursive_clearing = true, bool random_new_population = false, bool reinitializing_crowding_promising_archive = true);
		private:
			static inline void SelectionError(std::string func_name, std::string error_variable) {
				std::cout << "Cannnot execute\"" << func_name << "\". Please check: " << error_variable << std::endl;
				exit(0);
			}
		};
	}
}
#endif // !EA_HEADERS_DIFFERENTIAL_EVOLUTION_SELECT_H_
