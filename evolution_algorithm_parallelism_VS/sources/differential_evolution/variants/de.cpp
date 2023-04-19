#pragma warning(disable:6386; disable:6385)

#include "../../../headers/differential_evolution/variants/de.h"

namespace evolutionalgorithm {
	namespace differentialevolution {
		DE::DE()
		{
		}

		DE::~DE()
		{
		}

		void DE::InitialNP_Rank() {
			Basic::lv_->population_.InitialNpRank();
		}

		void DE::InitialPopulation() {
			Basic::CheckLV();
			Basic::lv_->population_.Initial(Basic::lv_->ev_.GetDE_NP(), Basic::lv_->ev_.GetDE_Dim());
			for (size_t np_counter = 0; np_counter < Basic::lv_->ev_.GetDE_NP(); ++np_counter) {
				for (size_t dim_counter = 0; dim_counter < Basic::lv_->ev_.GetDE_Dim(); ++dim_counter) {
					Basic::lv_->population_[np_counter][dim_counter] = PublicFunction::RandomDouble(Basic::lv_->ev_.GetBoundMin(dim_counter), Basic::lv_->ev_.GetBoundMax(dim_counter));
				}
			}
		}

		void DE::InitialNextGeneration() {
			Basic::lv_->next_generation_.Initial(Basic::lv_->ev_.GetDE_NextGenNP(), Basic::lv_->ev_.GetDE_Dim());
		}

		void DE::InitialFitness() {
			Basic::CheckFP();
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			if (Basic::lv_->is_initialed_ || Basic::lv_->fitness_ != nullptr) {
				free(Basic::lv_->fitness_);
				Basic::lv_->fitness_ = nullptr;
			}
			//if (np == SIZE_MAX) exit(0);
			Basic::lv_->fitness_ = (double*)malloc(sizeof(double) * np);
			if (Basic::lv_->fitness_ == nullptr || np == 0) {
				exit(0);
			}
			for (size_t np_counter = 0; np_counter < np; ++np_counter) {
				Basic::lv_->fitness_[np_counter] = Basic::lv_->fp_->compute(Basic::lv_->population_[np_counter].GetDimensions());
				if (Basic::lv_->fitness_[np_counter] < Basic::lv_->population_.GetBestFitness()) {
					Basic::lv_->population_.SetBestFitness(Basic::lv_->fitness_[np_counter]);
					Basic::lv_->population_.best_individual = np_counter;
					Basic::lv_->best_fitness_ever_ = Basic::lv_->fitness_[np_counter];
				}
				if (Basic::lv_->fitness_[np_counter] > Basic::lv_->worst_fitness_ever_) {
					Basic::lv_->worst_fitness_ever_ = Basic::lv_->fitness_[np_counter];
				}
			}
			Basic::lv_->now_ffe_ += np;
			Basic::lv_->now_gen_ = 1;
			if (!Basic::lv_->is_initialed_) Basic::lv_->is_initialed_ = true;
		}

		void DE::DemallocFitness() {
			free(Basic::lv_->fitness_);
			Basic::lv_->fitness_ = nullptr;
		}

		void DE::InitialUseStrategy() {
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			for (size_t next_gen_np_counter = 0; next_gen_np_counter < np; ++next_gen_np_counter) {
				Basic::lv_->next_generation_[next_gen_np_counter].use_strategy_ = Basic::lv_->ev_.GetDE_StrategiesPattern();
				Basic::lv_->next_generation_[next_gen_np_counter].use_strategy_index_ = 0;
			}
		}

		void DE::InitialUseCR() {
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			for (size_t next_gen_np_counter = 0; next_gen_np_counter < np; ++next_gen_np_counter) {
				Basic::lv_->next_generation_[next_gen_np_counter].SetCR(Basic::lv_->ev_.GetDE_CR());
			}
		}

		void DE::InitialUseF() {
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			for (size_t next_gen_np_counter = 0; next_gen_np_counter < np; ++next_gen_np_counter) {
				Basic::lv_->next_generation_[next_gen_np_counter].SetF(Basic::lv_->ev_.GetDE_F());
			}
		}

		void DE::InitialUseK() {
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			for (size_t next_gen_np_counter = 0; next_gen_np_counter < np; ++next_gen_np_counter) {
				Basic::lv_->next_generation_[next_gen_np_counter].SetK(Basic::lv_->ev_.GetDE_K());
			}
		}

		void DE::InitialUseP() {
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			for (size_t next_gen_np_counter = 0; next_gen_np_counter < np; ++next_gen_np_counter) {
				Basic::lv_->next_generation_[next_gen_np_counter].SetP(Basic::lv_->ev_.GetDE_P());
			}
		}

		void DE::SortNP_Rank() {
			//PublicFunction::SortPopulation(Basic::lv_->fitness_, &(Basic::lv_->population_.np_rank_[0]), Basic::lv_->ev_.GetDE_NP());
			PublicFunction::SortPopulation(&Basic::lv_->fitness_[0], &(Basic::lv_->population_.np_rank_[0]), Basic::lv_->ev_.GetDE_NP());
		}

		void DE::Evolution() {
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			const size_t next_gen_np = Basic::lv_->ev_.GetDE_NextGenNP();
			for (size_t next_gen_np_counter = 0; next_gen_np_counter < next_gen_np; ++next_gen_np_counter) {
				if ((Basic::lv_->next_generation_[next_gen_np_counter].use_strategy_ & corematerial::MutantAndCrossover::kCurToPBest) != corematerial::MutantAndCrossover::kCurToPBest) {
					Basic::StrategiesExecuter(
						Basic::lv_->population_,
						Basic::lv_->next_generation_[next_gen_np_counter],
						next_gen_np_counter,
						Basic::lv_->next_generation_[next_gen_np_counter].use_strategy_
					);
				}
				else {
					if ((Basic::lv_->ev_.GetDE_Memory_Pattern() & corematerial::Memory::kArchive) == corematerial::Memory::kArchive) {
						switch (Basic::lv_->ev_.GetDE_ArchivePattern()) {
						case(corematerial::Archive::kSHADE_Archive):
							Basic::StrategiesExecuter(
								Basic::lv_->population_,
								Basic::lv_->next_generation_[next_gen_np_counter],
								Basic::lv_->archive_,
								next_gen_np_counter,
								Basic::lv_->next_generation_[next_gen_np_counter].use_strategy_
							);
							break;
						case(corematerial::Archive::kPAEHN_Archive):
							Basic::StrategiesExecuter(
								Basic::lv_->population_,
								Basic::lv_->next_generation_[next_gen_np_counter],
								Basic::lv_->paehn_archive_,
								next_gen_np_counter,
								Basic::lv_->next_generation_[next_gen_np_counter].use_strategy_
							);
							break;
						default:
							printf("Unknown archive pattern.");
							exit(0);
						}
					}
				}
			}
		}

		void DE::Selection() {
			//
			// In this selection method, np must equal to next_gen_np.
			//
			double temp_fitness = 0.0;
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			for (size_t np_counter = 0; np_counter < np; ++np_counter) {
				Basic::lv_->next_generation_[np_counter].InitialComparisonResult();
				temp_fitness = Basic::lv_->fp_->compute(Basic::lv_->next_generation_[np_counter].GetDimensions());
				Basic::lv_->population_[np_counter].SetFitnessDifference(temp_fitness - Basic::lv_->fitness_[np_counter]);
				//Don't forget to set the index of compared parent.
				Basic::lv_->next_generation_[np_counter].SetComparedParent(np_counter); // Set which one parent is compared with this offspring.
				if (temp_fitness <= Basic::lv_->fitness_[np_counter]) {
					// if temp_fitness is better or equal to parent

					if (temp_fitness == Basic::lv_->fitness_[np_counter]) {
						Basic::lv_->next_generation_[np_counter].is_equal_ = true;
					}
					else {
						Basic::lv_->next_generation_[np_counter].is_less_ = true;
					}

					// 
					// Originally, population[np_counter] should be deallocated here, but the conception of "archive" needs to store the dimensions of parent which is worse than its offspring.
					// So, stores parent's dimensions into the map with its index of population as the key value.
					// Once the dimensions are stored into the archive, the second value in the map must be set to nullptr, or it will be wrongly deallocated even archive still needs that.
					// And who takes the dimensions has the responsibility to deallocated it when not using anymore.
					//

					//if (Basic::lv_->wait_free_map_of_individual_dimension_[np_counter] != nullptr) {
					//	//free(Basic::lv_->wait_free_map_of_individual_dimension_[np_counter]);
					//	Basic::lv_->wait_free_map_of_individual_dimension_[np_counter] = nullptr;
					//}

					// Records the dimension and fitness of the parents which were worse than the offspring.
					// They will be moved into archive or deallocated later. It's depended on the features/algorithms.
					if (Basic::lv_->wait_free_map_of_individual_dimension_.find(np_counter) == Basic::lv_->wait_free_map_of_individual_dimension_.end()) {
						//if np_counter doesn't exist in wait_free_map_of_individual_dimension_
						Basic::lv_->wait_free_map_of_individual_dimension_.emplace(np_counter, Basic::lv_->ev_.GetDE_Dim());
					}
					Basic::lv_->wait_free_map_of_individual_dimension_[np_counter] = Basic::lv_->population_[np_counter].GetDimensions();
					Basic::lv_->wait_free_map_of_individual_fitness_[np_counter] = Basic::lv_->fitness_[np_counter];

					// Moves the offsping into population.
					Basic::lv_->fitness_[np_counter] = temp_fitness;
					Basic::lv_->population_[np_counter] = Basic::lv_->next_generation_[np_counter].GetDimensions();
					Basic::lv_->next_generation_[np_counter] = nullptr; // Important !! or cause memory leak.

					// records if it become best individual
					if (temp_fitness < Basic::lv_->population_.GetBestFitness()) {
						Basic::lv_->population_.SetBestFitness(temp_fitness);
						Basic::lv_->population_.best_individual = np_counter;
					}
				}
				else {
					// if temp_fitness is worse than parent
					Basic::lv_->next_generation_[np_counter].is_larger_ = true;
				}
			}
			Basic::lv_->now_ffe_ += np;
		}

		void DE::AfterSelection() {
			// clear next generation 
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			for (size_t np_counter = 0; np_counter < np; ++np_counter) {
				if (Basic::lv_->next_generation_[np_counter].GetDimensions() != nullptr) {
					//free(Basic::lv_->next_generation_[np_counter].GetDimensions());
					Basic::lv_->next_generation_[np_counter] = nullptr;
				}
			}
			++Basic::lv_->now_gen_;

			// clear useless individual dimensions
			std::map<size_t, Individual>::iterator it = Basic::lv_->wait_free_map_of_individual_dimension_.begin();
			for (; it != Basic::lv_->wait_free_map_of_individual_dimension_.end(); ++it) {
				if (it->second.GetDimensions() != nullptr) {
					//free(it->second);
					Basic::lv_->wait_free_map_of_individual_dimension_[it->first] = nullptr;
				}
			}
		}

		void DE::NextGen() {
			//++Basic::lv_->now_gen_;
		}
	} // namespace evolutionalgorithm
} // namespace differentialevolution