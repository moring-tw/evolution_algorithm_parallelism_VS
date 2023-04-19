#pragma once
#ifndef EA_HEADERS_DIFFERENTIAL_EVOLUTION_PROCESS_CONTROLLER_H_
#define EA_HEADERS_DIFFERENTIAL_EVOLUTION_PROCESS_CONTROLLER_H_

#include <Windows.h>

#include <iostream>
#include <string>

#include "../../benchmark/inc/Benchmarks.h"
#include "../../benchmark/inc/MMOP/cec2013.h"
#include "../executivevariables.h"
#include "../publicfunction.h"
#include "strategies.h"
#include "memorize.h"
#include "evolve.h"
#include "select.h"
#include "sample.h"
#include "parameter_control.h"
#include "localvariables.h"
#include "variants/de.h"

//
// Here is the core Class who controls the algorithm flow of which algorithms/features should be executed.
// For the maintainability, it's extremely NOT recommended to control when to using different algorithms/features in other class.
// However, don't invoke the algorithm/feature directly, use the functions in the classes who extend from the behavior class.
// Those class will check every relative parameters to keep the program from crash. Then they will invoke the algorithms/features.
//
namespace evolutionalgorithm {
	namespace differentialevolution {
		class ProcessController {
		public:
			ProcessController();
			virtual ~ProcessController();
			void NextRun();
			void Execute();
			void SetExecutiveVariables(ExecutiveVariables&);
			void SetOutputFileStream(const string&, const string&, const string&);
			void SetBenchmark(Benchmarks*);
			double GetBest();
			std::vector<size_t>& GetHowManyOptima();
		protected:

		private:
			void SetBehaviors(Behavior&);
			void SetVariants(Basic&);
			void InitialVariants(Basic&);
			void BeforeInitial();
			void Initial();
			void AfterInitial();
			void BeforeEvolution();
			void Evolution();
			void AfterEvolution();
			void BeforeSelection();
			void Selection();
			void AfterSelection();
			void NextGen();
			bool ReachStopCriterion();
			void RecycleMemory();
			bool is_initialed_ = false;
			bool debug_ = false;
			
			// behaviors
			Sample sample_;
			Memorize memory_;
			ParameterControl parameter_control_;
			Evolve evolve_;
			Select select_;
			// local variables
			LocalVariables lv_;
			// strategies
			Strategies ss_;
			// variants
			DE de_;
			SaDE sade_;
			SHADE shade_;
			Random random_;
			PAEHN paehn_;

			static inline void parameter_control_error(std::string which_parameter) {
				std::cout << "Cannot control parameter "<< which_parameter<<" without specified method."<<std::endl;
				exit(0);
			}
		};
	} // namespace differentialevolution
} // namespace evolutionalgorithm
#endif // !EA_HEADERS_DIFFERENTIAL_EVOLUTION_PROCESS_CONTROLLER_H_
