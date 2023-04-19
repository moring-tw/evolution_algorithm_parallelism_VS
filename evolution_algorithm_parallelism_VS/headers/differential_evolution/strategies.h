#pragma once
#ifndef EA_HEADERS_STRATEGIES_H_
#define EA_HEADERS_STRATEGIES_H_

#include <cmath>

#include "../publicfunction.h"
#include "../population.h"
#include "../executivevariables.h"
#include "localvariables.h"

namespace evolutionalgorithm {
	namespace differentialevolution {
		/// <summary>
		///		This class is about implementing evolving strategies.
		/// </summary>
		class Strategies {
		public:
			Strategies();
			virtual ~Strategies();
			void DE_Rand(Population&, const size_t, Individual&, const size_t number_of_vector_differences = 1);
			//double* DE_Rand(Population&, const size_t, const size_t, const double, const double, const size_t number_of_vector_differences = 1);
			void DE_Best(Population&, const size_t, Individual&, const size_t number_of_vector_differences = 1);
			//double* DE_Best(Population&, const size_t, const size_t, const size_t, const double, const double, const size_t number_of_vector_differences = 1);
			void DE_Cur_To_Rand(Population&, const size_t, Individual&, const size_t number_of_vector_differences = 1);
			//double* DE_Cur_To_Rand(Population&, const size_t, const size_t, const double, const double, const size_t number_of_vector_differences = 1);
			void DE_Cur_To_Best(Population&, const size_t, Individual&, const size_t number_of_vector_differences = 1);
			//double* DE_Cur_To_Best(Population&, const size_t, const size_t, const size_t, const double, const double, const size_t number_of_vector_differences = 1);
			void DE_Cur_To_PBest(Population&, const size_t, Population&, Individual&, const size_t number_of_vector_differences = 1);
			//double* DE_Cur_To_PBest(Population&, const size_t, const size_t, const size_t*, Population&, const size_t, const double, const double, const size_t number_of_vector_differences = 1);

			void DE_Cur_To_PBest(Population&, const size_t, std::vector<Individual>&, const size_t, Individual&, const size_t number_of_vector_differences = 1);

			void SetExecutiveVariables(ExecutiveVariables&);
			void SetLocalVariables(LocalVariables&);
			bool IsOutBound(double, unsigned);
			void OutBoundHandler(Population&, unsigned, double*, unsigned);
			void OutBoundHandler(Individual&, double*, unsigned);
		protected:
			std::vector<std::pair<unsigned, unsigned>> VectorDifferencesPair(unsigned*, unsigned,const size_t, const size_t);
			// archive may be used by different algorithms, using pointer to remain flexibility.
			Population* archive_ = nullptr;
			size_t* archive_now_size_ = nullptr;
			bool archive_is_set_ = false;

			static inline void CheckArchiveIsSet(size_t archive_size, bool is_set) {
				if (archive_size > 0) {
					if (!is_set) {
						std::cout << "Before using DE/Current-to-pBest, Archive must be set if setting of Archive size in config file is larger than 0." << std::endl;
						exit(0);
					}
				}
			}

			static inline double* CreateTrailIndividual(const size_t dim) {
				if (dim == 0) {
					std::cout << "Unable to malloc trail_individual when DE_Rand(). Dim value:" << dim << std::endl;
					exit(0);
				}
				double* trail_individual = (double*)malloc(sizeof(double) * dim);
				if (trail_individual == nullptr || dim == 0) {
					std::cout << "Unable to malloc trail_individual when DE_Rand(). Dim value:" << dim << std::endl;
					exit(0);
				}
				return trail_individual;
			}
		private:
			ExecutiveVariables* ev_ = nullptr;
			LocalVariables* lv_ = nullptr;
			void ExecutiveVariablesCheck();
		};
	} // namespace differentialevolution
} // namespace evolutionalgorithm

#endif // EA_HEADERS_STRATEGIES_H_