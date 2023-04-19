#include "../../headers/differential_evolution/localvariables.h"

namespace evolutionalgorithm {
	namespace differentialevolution {
		LocalVariables::LocalVariables() : file_stream_("123")
		{
			sade_CRm_ = nullptr;
			sade_fail_memory_ = nullptr;
			sade_fail_memory_sum_ = nullptr;
			sade_memory_CR_ = nullptr;
			sade_success_memory_ = nullptr;
			sade_success_memory_sum_ = nullptr;
			strategies_probability_ = nullptr;
			shade_memory_CR_ = nullptr;
			shade_memory_F_ = nullptr;
		}

		LocalVariables::~LocalVariables()
		{
		}

		void LocalVariables::NextRun() {
			LocalVariables::population_.SetBestFitness(DBL_MAX);
			LocalVariables::population_.best_individual = SIZE_MAX;
			LocalVariables::now_ffe_ = 0;
			LocalVariables::now_gen_ = 0;
		}
	}

}