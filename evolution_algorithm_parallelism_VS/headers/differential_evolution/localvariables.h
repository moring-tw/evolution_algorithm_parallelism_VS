#pragma once
#ifndef EA_HEADERS_DIFFERENTIAL_EVOLUTION_LOCAL_VARIABLES_H_
#define EA_HEADERS_DIFFERENTIAL_EVOLUTION_LOCAL_VARIABLES_H_

#include <fstream>

#include "../cells.h"
#include "../population.h"
#include "../subpopulation.h"
#include "../../benchmark/inc/Benchmarks.h"
#include "../executivevariables.h"

namespace evolutionalgorithm {
	namespace differentialevolution {
		/// <summary>
		/// Every variables will be changed through the program should be stored here.
		/// 
		/// </summary>
		class LocalVariables
		{
		public:
			LocalVariables();
			~LocalVariables();
			void NextRun();

		private:
			// widely used
			Population population_;
			Population next_generation_;
			Benchmarks* fp_ = nullptr;
			ExecutiveVariables ev_;
			double* fitness_ = nullptr;
			//double now_best_ = DBL_MAX;
			//size_t now_best_individual_ = SIZE_MAX;
			unsigned now_ffe_ = 0;
			unsigned now_gen_ = 0;
			//size_t* np_rank = nullptr;
			std::map<size_t, Individual> wait_free_map_of_individual_dimension_;
			//std::map<size_t, double*> wait_free_map_of_individual_dimension_;
			std::map<size_t, double> wait_free_map_of_individual_fitness_;

			double best_fitness_ever_ = DBL_MAX;
			double worst_fitness_ever_ = DBL_MIN;
			ofstream file_stream_;
			string file_path_ = "";
			string filename_ = "";
			// END widely used

			bool is_initialed_ = false;

			// SaDE
			size_t** sade_success_memory_ = nullptr,
				* sade_success_memory_sum_ = nullptr,
				** sade_fail_memory_ = nullptr,
				* sade_fail_memory_sum_ = nullptr;
			double* strategies_probability_ = nullptr, * sade_CRm_ = nullptr;
			std::vector<double>*** sade_memory_CR_ = nullptr;
			// END SaDE

			// SHADE
			Population archive_;
			std::vector<double> archive_fitness_;
			size_t shade_insert_h_ = 0; // , shade_archive_now_size_ = 0;
			size_t * shade_ri_ = nullptr;
			double* shade_memory_CR_ = nullptr, * shade_memory_F_ = nullptr; // fixed H
			std::vector<double> delta_fitness_; // clear each generation
			std::vector<double> shade_success_CR_, shade_success_F_; // clear each generation
			// END SHADE

			// PAEHN
			Population paehn_archive_;
			std::vector<double> paehn_archive_fitness_;

				//subpopulation parameters
				size_t paehn_subpopulation_size_;
				std::vector<paehn::Subpopulation> paehn_subpopulations_;

					// SHADE parameters
					size_t paehn_insert_H_ = 0;
					bool paehn_is_subpop_fitness_copied_ = false;
					bool paehn_is_subpop_delta_fitness_calculated_ = false;
					bool paehn_is_ri_produced_ = false;
					// end SHADE parameters 

				bool paehn_is_this_gen_create_subpopulations_ = false; // using when debug output
				size_t paehn_form_subpopulations_gen_counter_ = 0; // using when debug output
				std::vector<size_t> total_found_optima_; // using when final output
				// end subpopulation parameters
			// END PAEHN

			friend class Random;
			friend class DE;
			friend class SaDE;
			friend class SHADE;
			friend class PAEHN;
			friend class ProcessController;
			friend class Memorize;
			friend class Evolve;
			friend class Select;
			friend class Basic;
			friend class ParameterControl;
			friend class Strategies;
		};
	}
}



#endif // !EA_HEADERS_DIFFERENTIAL_EVOLUTION_LOCAL_VARIABLES_H_
