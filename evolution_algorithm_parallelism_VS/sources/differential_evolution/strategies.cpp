#include "../../headers/differential_evolution/strategies.h"

namespace evolutionalgorithm {
	namespace differentialevolution {
		Strategies::Strategies() {

		}

		Strategies::~Strategies() {
			Strategies::archive_ = nullptr;
		}

		void Strategies::SetExecutiveVariables(ExecutiveVariables& ev) {
			Strategies::ev_ = &ev;
		}

		void Strategies::SetLocalVariables(LocalVariables& lv) {
			Strategies::lv_ = &lv;
		}

		void Strategies::ExecutiveVariablesCheck() {
			if (Strategies::ev_ == nullptr) {
				std::cout << "Before using any strategy, ev_(ExecutiveVariables) must be set already." << std::endl;
				exit(0);
			}
		}

		/// <summary>
		/// DE/Rand/x/bin
		/// </summary>
		/// <param name="population"></param>
		/// <param name="np"></param>
		/// <param name="base_individual"></param>
		/// <param name="CR"></param>
		/// <param name="F"></param>
		/// <param name="number_of_vector_differences"></param>
		/// <returns></returns>
		void Strategies::DE_Rand(Population& population, const size_t base_individual, Individual& next_individual, const size_t number_of_vector_differences) {
			ExecutiveVariablesCheck();
			const size_t np = population.GetNP();
			unsigned* individuals = PublicFunction::ShufflePopulation(np, base_individual, number_of_vector_differences * 2 + 1);
			const size_t dim = Strategies::ev_->GetDE_Dim();
			//double* trail_individual = CreateTrailIndividual(dim);
			unsigned rnbr = PublicFunction::RandomInt(0, dim - 1);

			for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
				if (PublicFunction::RandomDouble(0, 1) < next_individual.GetCR() || dim_counter == rnbr) {
					next_individual[dim_counter] = population[individuals[0]][dim_counter];
					for (unsigned vector_diff_counter = 0; vector_diff_counter < number_of_vector_differences; ++vector_diff_counter) {
						next_individual[dim_counter] += next_individual.GetF() * (population[individuals[vector_diff_counter * 2 + 1]][dim_counter] - population[individuals[vector_diff_counter * 2 + 2]][dim_counter]);
					}
					if (IsOutBound(next_individual[dim_counter], dim_counter)) {
						OutBoundHandler(population, base_individual, &next_individual[0], dim_counter);
					}
				}
				else {
					next_individual[dim_counter] = population[base_individual][dim_counter];
				}
			}
			free(individuals);
		}
		/// <summary>
		/// DE/Best/x/bin
		/// </summary>
		/// <param name="population"></param>
		/// <param name="np"></param>
		/// <param name="base_individual"></param>
		/// <param name="best_individual"></param>
		/// <param name="CR"></param>
		/// <param name="F"></param>
		/// <param name="number_of_vector_differences"></param>
		/// <returns></returns>
		void Strategies::DE_Best(Population& population, const size_t base_individual, Individual& next_individual, const size_t number_of_vector_differences) {
			ExecutiveVariablesCheck();
			const size_t np = population.GetNP();
			const size_t best_individual = population.best_individual;
			unsigned* individuals = PublicFunction::ShufflePopulation(np, base_individual, number_of_vector_differences * 2);
			const size_t dim = Strategies::ev_->GetDE_Dim();
			//double* trail_individual = CreateTrailIndividual(dim);
			unsigned rnbr = PublicFunction::RandomInt(0, dim - 1);

			for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
				if (PublicFunction::RandomDouble(0, 1) < next_individual.GetCR() || dim_counter == rnbr) {
					next_individual[dim_counter] = population[best_individual][dim_counter];
					for (unsigned vector_diff_counter = 0; vector_diff_counter < number_of_vector_differences; ++vector_diff_counter) {
						next_individual[dim_counter] += next_individual.GetF() * (population[individuals[vector_diff_counter * 2]][dim_counter] - population[individuals[vector_diff_counter * 2 + 1]][dim_counter]);
					}
					if (IsOutBound(next_individual[dim_counter], dim_counter)) {
						OutBoundHandler(population, base_individual, &next_individual[0], dim_counter);
					}
				}
				else {
					next_individual[dim_counter] = population[base_individual][dim_counter];
				}
			}
			free(individuals);
		}
		/// <summary>
		/// 
		/// </summary>
		/// <param name="population"></param>
		/// <param name="np"></param>
		/// <param name="base_individual"></param>
		/// <param name="F"></param>
		/// <param name="K"></param>
		/// <param name="number_of_vector_differences"></param>
		/// <returns></returns>
		void Strategies::DE_Cur_To_Rand(Population& population, const size_t base_individual, Individual& next_individual, const size_t number_of_vector_differences) {
			ExecutiveVariablesCheck();
			const size_t np = population.GetNP();
			unsigned* individuals = PublicFunction::ShufflePopulation(np, base_individual, number_of_vector_differences * 2 + 1);
			const size_t dim = Strategies::ev_->GetDE_Dim();
			//double* trail_individual = CreateTrailIndividual(dim);

			for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
				next_individual[dim_counter] = population[base_individual][dim_counter] + next_individual.GetK() * (population[individuals[0]][dim_counter] - population[base_individual][dim_counter]);
				for (unsigned vector_diff_counter = 0; vector_diff_counter < number_of_vector_differences; ++vector_diff_counter) {
					next_individual[dim_counter] += next_individual.GetF() * (population[individuals[vector_diff_counter * 2 + 1]][dim_counter] - population[individuals[vector_diff_counter * 2 + 2]][dim_counter]);
				}
				if (IsOutBound(next_individual[dim_counter], dim_counter)) {
					OutBoundHandler(population, base_individual, &next_individual[0], dim_counter);
				}
			}
			free(individuals);
		}
		/// <summary>
		/// 
		/// </summary>
		/// <param name="population"></param>
		/// <param name="np"></param>
		/// <param name="base_individual"></param>
		/// <param name="best_individual"></param>
		/// <param name="CR"></param>
		/// <param name="F"></param>
		/// <param name="number_of_vector_differences"></param>
		/// <returns></returns>
		void Strategies::DE_Cur_To_Best(Population& population, const size_t base_individual, Individual& next_individual, const size_t number_of_vector_differences) {
			ExecutiveVariablesCheck();
			const size_t np = population.GetNP();
			const size_t best_individual = population.best_individual;
			unsigned* individuals = PublicFunction::ShufflePopulation(np, base_individual, number_of_vector_differences * 2);
			const size_t dim = Strategies::ev_->GetDE_Dim();
			//double* trail_individual = CreateTrailIndividual(dim);
			unsigned rnbr = PublicFunction::RandomInt(0, dim - 1);

			for (unsigned dim_counter = 0; dim_counter < dim; ++dim_counter) {
				if (PublicFunction::RandomDouble(0, 1) < next_individual.GetCR() || dim_counter == rnbr) {
					next_individual[dim_counter] = population[base_individual][dim_counter] + next_individual.GetF() * (population[best_individual][dim_counter] - population[base_individual][dim_counter]);
					for (unsigned vector_diff_counter = 0; vector_diff_counter < number_of_vector_differences; ++vector_diff_counter) {
						next_individual[dim_counter] += next_individual.GetF() * (population[individuals[vector_diff_counter * 2]][dim_counter] - population[individuals[vector_diff_counter * 2 + 1]][dim_counter]);
					}
					if (IsOutBound(next_individual[dim_counter], dim_counter)) {
						OutBoundHandler(population, base_individual, &next_individual[0], dim_counter);
					}
				}
				else {
					next_individual[dim_counter] = population[base_individual][dim_counter];
				}
			}
			free(individuals);
		}
		/// <summary>
		/// 
		/// </summary>
		/// <param name="population"></param>
		/// <param name="np"></param>
		/// <param name="base_individual"></param>
		/// <param name="sorted_individuals"></param>
		/// <param name="archive"></param>
		/// <param name="archive_size"></param>
		/// <param name="CR"></param>
		/// <param name="F"></param>
		/// <param name="number_of_vector_differences"></param>
		/// <returns></returns>
		void Strategies::DE_Cur_To_PBest(Population& population, const size_t base_individual, Population& archive, Individual& next_individual, const size_t number_of_vector_differences) {
			ExecutiveVariablesCheck();
			//
			// If trail individual is successfully selected but it had been regenerated because of out-of-bound. Should this success of CR and F be added into success-history?
			// This Strategy implementation(DE/Current-to-pBest) is follow SHADE paper, but this "ADDITIONAL" consideration of out-of-bound is not in general SHADE.
			// BTW, it is using kHalfBoundDistance when any dimension is outside the boundaries in original SHADE.
			//
			Strategies::ev_->UpdateOut_B(false);
			const size_t np = population.GetNP();
			unsigned total_size = np + archive.GetNowSize();
			unsigned* individuals = PublicFunction::ShufflePopulation(total_size, base_individual, total_size - 1);
			std::vector<std::pair<unsigned, unsigned>> difference_individual_pairs = VectorDifferencesPair(individuals, number_of_vector_differences, np, archive.GetNowSize());
			const size_t dim = Strategies::ev_->GetDE_Dim();
			//double* trail_individual = CreateTrailIndividual(dim);
			unsigned rnbr = PublicFunction::RandomInt(0, dim - 1);
			const double P = next_individual.GetP();
			//unsigned pbest_individual = sorted_individuals[PublicFunction::RandomInt(0, static_cast<int>(std::round(np * P)))];
			size_t pbest_individual = population.np_rank_[PublicFunction::RandomInt(0, static_cast<int>(std::round(np * P)))];
			double difference = 0;

			for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
				if (PublicFunction::RandomDouble(0, 1) < next_individual.GetCR() || dim_counter == rnbr) {
					next_individual[dim_counter] = population[base_individual][dim_counter] + next_individual.GetF() * (population[pbest_individual][dim_counter] - population[base_individual][dim_counter]);
					for (unsigned vector_differences_counter = 0; vector_differences_counter < number_of_vector_differences; ++vector_differences_counter) {
						difference = population[difference_individual_pairs[vector_differences_counter].first][dim_counter];
						// The function VectorDifferenecesPair() guarantees that the first value of the pair in the vector is selected from population.
						if (difference_individual_pairs[vector_differences_counter].second >= np) {
							// It means that the individual should be taken from archive if the second value is larger than or equal to NP
							difference -= archive[difference_individual_pairs[vector_differences_counter].second - np][dim_counter];
						}
						else {
							// or from population.
							difference -= population[difference_individual_pairs[vector_differences_counter].second][dim_counter];
						}
						next_individual[dim_counter] += next_individual.GetF() * difference;
					}
					if (IsOutBound(next_individual[dim_counter], dim_counter)) {
						ev_->UpdateOut_B(true);
						OutBoundHandler(population, base_individual, &next_individual[0], dim_counter);
					}
				}
				else {
					next_individual[dim_counter] = population[base_individual][dim_counter];
				}
			}
			free(individuals);
		}

		void Strategies::DE_Cur_To_PBest(Population& population, const size_t base_individual, std::vector<Individual>& archive, const size_t archive_size, Individual& next_individual, const size_t number_of_vector_differences) {
			ExecutiveVariablesCheck();
			//
			// If trail individual is successfully selected but it had been regenerated because of out-of-bound. Should this success of CR and F be added into success-history?
			// This Strategy implementation(DE/Current-to-pBest) is follow SHADE paper, but this "ADDITIONAL" consideration of out-of-bound is not in general SHADE.
			// BTW, it is using kHalfBoundDistance when any dimension is outside the boundaries in original SHADE.
			//
			Strategies::ev_->UpdateOut_B(false);
			const size_t np = population.GetNP();
			unsigned total_size = np + archive.size();
			unsigned* individuals = PublicFunction::ShufflePopulation(total_size, base_individual, total_size - 1);
			std::vector<std::pair<unsigned, unsigned>> difference_individual_pairs = VectorDifferencesPair(individuals, number_of_vector_differences, np, archive.size());
			const size_t dim = Strategies::ev_->GetDE_Dim();
			//double* trail_individual = CreateTrailIndividual(dim);
			unsigned rnbr = PublicFunction::RandomInt(0, dim - 1);
			//const double P = next_individual.GetP();
			const double P = 0.5; // if archive size = 0, and P = 0, this function is equal to cur-to-best
			//unsigned pbest_individual = sorted_individuals[PublicFunction::RandomInt(0, static_cast<int>(std::round(np * P)))];
			size_t pbest_individual = population.np_rank_[PublicFunction::RandomInt(0, static_cast<int>(std::round(np * P)))];
			double difference = 0;

			for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
				if (PublicFunction::RandomDouble(0, 1) < next_individual.GetCR() || dim_counter == rnbr) {
					next_individual[dim_counter] = population[base_individual][dim_counter] + next_individual.GetF() * (population[pbest_individual][dim_counter] - population[base_individual][dim_counter]);
					for (unsigned vector_differences_counter = 0; vector_differences_counter < number_of_vector_differences; ++vector_differences_counter) {
						difference = population[difference_individual_pairs[vector_differences_counter].first][dim_counter];
						//difference = population[pbest_individual][dim_counter];//------------------------------------
						// The function VectorDifferenecesPair() guarantees that the first value of the pair in the vector is selected from population.
						if (difference_individual_pairs[vector_differences_counter].second >= np) {
							// It means that the individual should be taken from archive if the second value is larger than or equal to NP
							difference -= archive[difference_individual_pairs[vector_differences_counter].second - np][dim_counter];
						}
						else {
							// or from population.
							difference -= population[difference_individual_pairs[vector_differences_counter].second][dim_counter];
						}
						next_individual[dim_counter] += next_individual.GetF() * difference;
					}
					if (IsOutBound(next_individual[dim_counter], dim_counter)) {
						ev_->UpdateOut_B(true);
						OutBoundHandler(population, base_individual, &next_individual[0], dim_counter);
					}
				}
				else {
					next_individual[dim_counter] = population[base_individual][dim_counter];
				}
			}
			free(individuals);
		}

		bool Strategies::IsOutBound(double value, unsigned base_dim) {
			if (value < Strategies::ev_->GetBoundMin(base_dim) || value > Strategies::ev_->GetBoundMax(base_dim))
				return true;
			return false;
		}

		void Strategies::OutBoundHandler(Population& population, unsigned base_individual, double* trail_individual, unsigned base_dim) {
			corematerial::OutBoundHandler obh_pattern = Strategies::ev_->GetDE_OutBoundHandlerPattern();
			const double bound_min = Strategies::ev_->GetBoundMin(base_dim), bound_max = Strategies::ev_->GetBoundMax(base_dim);
			if (obh_pattern == corematerial::OutBoundHandler::kBackToBound) {
				if (trail_individual[base_dim] < bound_min) trail_individual[base_dim] = bound_min;
				else if (trail_individual[base_dim] > bound_max) trail_individual[base_dim] = bound_max;
			} else if (obh_pattern == corematerial::OutBoundHandler::kHalfBoundDistance) {
				if (trail_individual[base_dim] < bound_min)
					// Middle point of the "BASE" individual's dimension to the lower bound.
					trail_individual[base_dim] = (population[base_individual][base_dim] + bound_min) / 2;
				else if (trail_individual[base_dim] > bound_max)
					// Middle point of the upper bound to the "BASE" individual's dimension.
					trail_individual[base_dim] = (bound_max + population[base_individual][base_dim]) / 2;
			} else if (obh_pattern == corematerial::OutBoundHandler::kRerandomDimension) {
				trail_individual[base_dim] = PublicFunction::RandomDouble(bound_min, bound_max);
			} else if (obh_pattern == corematerial::OutBoundHandler::kRerandomIndividual) {
				for (size_t dim_counter = 0; dim_counter < Strategies::ev_->GetDE_Dim(); ++dim_counter) {
					trail_individual[dim_counter] = PublicFunction::RandomDouble(Strategies::ev_->GetBoundMin(dim_counter), Strategies::ev_->GetBoundMax(dim_counter));
				}
			}
		}

		void Strategies::OutBoundHandler(Individual& individual, double* trail_individual, unsigned base_dim) {
			corematerial::OutBoundHandler obh_pattern = Strategies::ev_->GetDE_OutBoundHandlerPattern();
			const double bound_min = Strategies::ev_->GetBoundMin(base_dim), bound_max = Strategies::ev_->GetBoundMax(base_dim);
			if (obh_pattern == corematerial::OutBoundHandler::kBackToBound) {
				if (trail_individual[base_dim] < bound_min) trail_individual[base_dim] = bound_min;
				else if (trail_individual[base_dim] > bound_max) trail_individual[base_dim] = bound_max;
			}
			else if (obh_pattern == corematerial::OutBoundHandler::kHalfBoundDistance) {
				if (trail_individual[base_dim] < bound_min)
					// Middle point of the "BASE" individual's dimension to the lower bound.
					trail_individual[base_dim] = (individual[base_dim] + bound_min) / 2;
				else if (trail_individual[base_dim] > bound_max)
					// Middle point of the upper bound to the "BASE" individual's dimension.
					trail_individual[base_dim] = (bound_max + individual[base_dim]) / 2;
			}
			else if (obh_pattern == corematerial::OutBoundHandler::kRerandomDimension) {
				trail_individual[base_dim] = PublicFunction::RandomDouble(bound_min, bound_max);
			}
			else if (obh_pattern == corematerial::OutBoundHandler::kRerandomIndividual) {
				for (size_t dim_counter = 0; dim_counter < Strategies::ev_->GetDE_Dim(); ++dim_counter) {
					trail_individual[dim_counter] = PublicFunction::RandomDouble(Strategies::ev_->GetBoundMin(dim_counter), Strategies::ev_->GetBoundMax(dim_counter));
				}
			}
		}

		std::vector<std::pair<unsigned, unsigned>> Strategies::VectorDifferencesPair(unsigned* individuals, unsigned number_of_vector_differences,const size_t NP, const size_t archive_size) {
			//
			// Only the second value of the pair in the vector is randomly select from the union of the population and the archive.
			// Otherwise is merely selected from population.
			//
			//const size_t NP = Strategies::ev_->GetDE_NP();
			//const size_t archive_size = lv_->shade_archive_now_size_;
			std::vector<std::pair<unsigned, unsigned>> difference_individual_pairs;
			for (unsigned vector_differences_counter = 0; vector_differences_counter < number_of_vector_differences; ++vector_differences_counter) {
				unsigned individual_first;
				unsigned index = 0;
				while (index < NP + archive_size) {
					if (individuals[index] < NP && individuals[index] != UINT_MAX) {
						// Only choose from population.
						individual_first = individuals[index];
						individuals[index] = UINT_MAX;
						break;
					}
					if (index == NP + archive_size - 1) {
						std::cout << "Can't find any individuals[index] falling within the interval [0:" << Strategies::ev_->GetDE_NP() << "]." << std::endl;
						std::cout << "In strategies.cpp, function Strategies::VectorDifferencesPair(unsigned*,unsigned)." << std::endl;
						exit(0);
					}
					++index;
				}
				index = 0;
				while (index < NP + archive_size) {
					// choose from the union of the population and the archive.
					if (individuals[index] != UINT_MAX) {
						difference_individual_pairs.push_back(std::make_pair(individual_first, individuals[index]));
						individuals[index] = UINT_MAX;
						break;
					}
					++index;
				}
			}
			return difference_individual_pairs;
		}
	} // differentialevolution
} // evolutionalgorithm