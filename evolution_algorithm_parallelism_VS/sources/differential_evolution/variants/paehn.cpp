#include "../../../headers/differential_evolution/variants/paehn.h"
constexpr auto kOutput = false;
constexpr auto kAccuracyLevel = 0.00001;
namespace evolutionalgorithm {
	namespace differentialevolution {
		PAEHN::PAEHN() {

		}

		PAEHN::~PAEHN() {

		}

		void PAEHN::NextRun() {
			Basic::lv_->paehn_subpopulations_.clear();
			Basic::lv_->paehn_archive_fitness_.clear();
			Basic::lv_->paehn_is_ri_produced_ = false;
			Basic::lv_->paehn_is_subpop_delta_fitness_calculated_ = false;
			Basic::lv_->paehn_is_subpop_fitness_copied_ = false;
			Basic::lv_->paehn_insert_H_ = 0;
			PAEHN::subpop_output_counter = 0;
			Basic::lv_->paehn_form_subpopulations_gen_counter_ = 0;
			Basic::lv_->paehn_is_this_gen_create_subpopulations_ = false;
			for (auto fs = PAEHN::file_streams_.begin(); fs != PAEHN::file_streams_.end(); ++fs) {
				fs->close();
			}
			PAEHN::file_streams_.clear();
			PAEHN::InitializeOfstream();
		}

		void PAEHN::initialPAEHN_Archive() {
			//std::cout << "initial paehn archive!!! " << Basic::lv_->ev_.GetArchiveSize() << std::endl;
			size_t np = Basic::lv_->ev_.GetDE_NP();
			size_t dim = Basic::lv_->ev_.GetDE_Dim();
			size_t archive_size = Basic::lv_->ev_.GetArchiveSize();

			Basic::lv_->paehn_archive_.Initial(archive_size, dim);
			for (size_t np_counter = 0; np_counter < np && np_counter < archive_size; ++np_counter) {
				Basic::lv_->paehn_archive_[np_counter].CopyDimension(Basic::lv_->population_[np_counter].GetDimensions());
				Basic::lv_->paehn_archive_fitness_.push_back(Basic::lv_->fitness_[np_counter]);
				Basic::lv_->paehn_archive_.AddNowSize();
			}
			while (Basic::lv_->paehn_archive_.GetNowSize() < archive_size) {
				for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter)
					Basic::lv_->paehn_archive_[Basic::lv_->paehn_archive_.GetNowSize()][dim_counter] = PublicFunction::RandomDouble(Basic::lv_->ev_.GetBoundMin(dim_counter), Basic::lv_->ev_.GetBoundMax(dim_counter));
				Basic::lv_->paehn_archive_fitness_.push_back(Basic::lv_->fp_->compute(Basic::lv_->paehn_archive_[Basic::lv_->paehn_archive_.GetNowSize()].GetDimensions()));
				++Basic::lv_->now_ffe_;
				Basic::lv_->paehn_archive_.AddNowSize();
			}

			PAEHN::OutputPopulation();
			PAEHN::OutputArchive();
		}

		void PAEHN::SuccessSelectionCrowdingArchive(double* dimensions, const double fitness, bool random_promising_archive) {
			const size_t archive_size = Basic::lv_->ev_.GetArchiveSize();
			size_t replace_index = 0;
			if (Basic::lv_->paehn_archive_.GetNowSize() < archive_size) {
				replace_index = Basic::lv_->paehn_archive_.GetNowSize();
				Basic::lv_->paehn_archive_[replace_index] = &dimensions[0];
				Basic::lv_->paehn_archive_fitness_.push_back(fitness);
				Basic::lv_->paehn_archive_.AddNowSize();
			}
			else {
				if (random_promising_archive)
					replace_index = PublicFunction::RandomInt(0, archive_size - 1);
				else
					replace_index = FindNearestIndividualInArchive(&dimensions[0]);
				if (fitness < Basic::lv_->paehn_archive_fitness_[replace_index]) {
					Basic::lv_->paehn_archive_[replace_index] = &dimensions[0];
					Basic::lv_->paehn_archive_fitness_[replace_index] = fitness;
				}
			}
		}

		size_t PAEHN::FindNearestIndividualInArchive(double* individual) {
			double min_distance = DBL_MAX;
			size_t min_distance_individual_index = SIZE_MAX;
			double current_distance = 0;
			const size_t archive_size = Basic::lv_->ev_.GetArchiveSize();
			for (size_t archive_counter = 0; archive_counter < archive_size; ++archive_counter) {
				current_distance = PublicFunction::CalculateEuclideanDistance(individual, Basic::lv_->paehn_archive_[archive_counter].GetDimensions(), Basic::lv_->ev_.GetDE_Dim());
				if (current_distance < min_distance) {
					min_distance = current_distance;
					min_distance_individual_index = archive_counter;
				}
			}
			if (min_distance_individual_index == SIZE_MAX)
				std::cout << "???";
			return min_distance_individual_index;
		}

		void PAEHN::SelectSHADE_SubpopCR() {
			if (Basic::lv_->paehn_subpopulations_.size() == 0) return;
			PAEHN::RandomRi();
			const size_t next_gen_subpop_np = Basic::lv_->ev_.GetPAEHN_SubpopulationSize();
			auto now_subpop = Basic::lv_->paehn_subpopulations_.begin();

			while (now_subpop != Basic::lv_->paehn_subpopulations_.end()) {
				if (!Basic::lv_->paehn_is_subpop_fitness_copied_)
					(*now_subpop).copied_fitness_ = (*now_subpop).fitness_; // copy all
				for (size_t subpop_np_counter = 0; subpop_np_counter < next_gen_subpop_np; ++subpop_np_counter) {
					(*now_subpop).children_[subpop_np_counter].SetCR(PAEHN::ProduceCR((*now_subpop).memory_CR_[(*now_subpop).random_ri_[subpop_np_counter]]));
				}
				++now_subpop;
			}
			Basic::lv_->paehn_is_subpop_fitness_copied_ = true;
		}

		void PAEHN::SelectSHADE_SubpopF() {
			if (Basic::lv_->paehn_subpopulations_.size() == 0) return;
			PAEHN::RandomRi();
			const size_t next_gen_subpop_np = Basic::lv_->ev_.GetPAEHN_SubpopulationSize();
			auto now_subpop = Basic::lv_->paehn_subpopulations_.begin();
			while (now_subpop != Basic::lv_->paehn_subpopulations_.end()) {
				if (!Basic::lv_->paehn_is_subpop_fitness_copied_)
					(*now_subpop).copied_fitness_ = (*now_subpop).fitness_; // copy all
				for (size_t subpop_np_counter = 0; subpop_np_counter < next_gen_subpop_np; ++subpop_np_counter) {
					(*now_subpop).children_[subpop_np_counter].SetF(PAEHN::ProduceF((*now_subpop).memory_F_[(*now_subpop).random_ri_[subpop_np_counter]]));
				}
				++now_subpop;
			}
			Basic::lv_->paehn_is_subpop_fitness_copied_ = true;
		}

		void PAEHN::SubpopulationMutationAndCrossover() {
			if (Basic::lv_->paehn_subpopulations_.size() == 0) return;
			auto now_subpop = Basic::lv_->paehn_subpopulations_.begin();
			auto child = (*now_subpop).children_.begin();
			size_t base_individual = 0;
			while (now_subpop != Basic::lv_->paehn_subpopulations_.end()) {
				if (now_subpop->is_stagnation) {
					++now_subpop;
					continue;
				}
				child = (*now_subpop).children_.begin();
				base_individual = 0;

				PublicFunction::SortPopulation(now_subpop->fitness_, now_subpop->np_rank_, Basic::lv_->ev_.GetPAEHN_SubpopulationSize());
				while (child != (*now_subpop).children_.end()) {
					size_t individual_counter = 0;
					Basic::ss_->DE_Cur_To_PBest((*now_subpop), base_individual, (*now_subpop).archive_, (*now_subpop).archive_.size(), (*child));
					/*child->SetCR(0.5);
					child->SetF(0.3);
					Basic::ss_->DE_Cur_To_Best((*now_subpop), base_individual, (*child));*/
					++base_individual;
					++child;
				}
				++now_subpop;
			}
		}

		void PAEHN::FormSubpopulation(bool enable_recursive_clearing) {
			const size_t subpopulation_size = Basic::lv_->ev_.GetPAEHN_SubpopulationSize(); // should be move to localvariable
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			const size_t H = Basic::lv_->ev_.GetPAEHN_H();
			size_t form_how_many_subpops_counter = 0;
			std::vector<size_t> top_individual_indice;

			PublicFunction::SortPopulation(Basic::lv_->fitness_, &Basic::lv_->population_.np_rank_[0], np);
			++Basic::lv_->paehn_form_subpopulations_gen_counter_;
			for (size_t np_counter = 0; np_counter < np; ++np_counter) { // only top_percent individuals will be choosen as speciation seeds
				if (Basic::lv_->population_[Basic::lv_->population_.np_rank_[np_counter]].is_choosen_)
					continue;
				Basic::lv_->paehn_subpopulations_.emplace_back();
				paehn::Subpopulation& new_subpopulation = Basic::lv_->paehn_subpopulations_.back();
				new_subpopulation.Initial(subpopulation_size, dim);
				new_subpopulation.children_.Initial(subpopulation_size, dim);
				new_subpopulation.Initial_Memory_CR(H, 0.5);
				new_subpopulation.Initial_Memory_F(H, 0.5);
				// the individual hasn't been choosen by speciation.
				std::vector<double> all_dis;
				std::vector<size_t> dis_rank;
				std::vector<size_t> dis_indices;
				double now_dis = 0.0;
				for (size_t remain_np_counter = np_counter, index_counter = 0; remain_np_counter < np; ++remain_np_counter) {
					if (Basic::lv_->population_[Basic::lv_->population_.np_rank_[remain_np_counter]].is_choosen_)
						continue;
					now_dis = PublicFunction::CalculateEuclideanDistance(
						Basic::lv_->population_[Basic::lv_->population_.np_rank_[np_counter]].GetDimensions(),
						Basic::lv_->population_[Basic::lv_->population_.np_rank_[remain_np_counter]].GetDimensions(),
						dim
					);
					all_dis.push_back(now_dis);
					dis_rank.push_back(Basic::lv_->population_.np_rank_[remain_np_counter]);
					dis_indices.push_back(index_counter++);

				}
				if (all_dis.size() < subpopulation_size) {
					// no enough individuals to form a new subpopulation.
					Basic::lv_->paehn_subpopulations_.pop_back();
					break;
				}
				PublicFunction::SortPopulation(all_dis, dis_indices, all_dis.size());
				new_subpopulation.SetBestFitness(DBL_MAX);
				// initial centroid
				new_subpopulation.centroid.reserve(dim);
				for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
					new_subpopulation.centroid.push_back(0.0);
				}

				// copy top subpop_np individuals into subpop and calculate centroid
				// and create new individual into population
				for (size_t top_individuals_counter = 0; top_individuals_counter < subpopulation_size; ++top_individuals_counter) {
					new_subpopulation[top_individuals_counter] = Basic::lv_->population_[dis_rank[dis_indices[top_individuals_counter]]].GetDimensions(); // copy dimensions
					Basic::lv_->population_[dis_rank[dis_indices[top_individuals_counter]]].is_choosen_ = true;
					new_subpopulation.fitness_.push_back(Basic::lv_->fitness_[dis_rank[dis_indices[top_individuals_counter]]]); // copy fitness
					if (Basic::lv_->fitness_[dis_rank[dis_indices[top_individuals_counter]]] < new_subpopulation.GetBestFitness()) {
						new_subpopulation.SetBestFitness(Basic::lv_->fitness_[dis_rank[dis_indices[top_individuals_counter]]]);
						new_subpopulation.best_individual = top_individuals_counter;
					}
					for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) { // calculate the centroid
						new_subpopulation.centroid[dim_counter] += (new_subpopulation[top_individuals_counter].GetDimensions()[dim_counter] / subpopulation_size);
					}
				}
				// calculate all the distance between the individuals in the subpop and subpop centroid
				double max_dis = 0, dis = 0, min_dis = DBL_MAX;
				for (size_t subpop_counter = 0; subpop_counter < subpopulation_size; ++subpop_counter) {
					dis = PublicFunction::CalculateEuclideanDistance(new_subpopulation[subpop_counter].GetDimensions(), new_subpopulation.centroid, dim);
					if (dis > max_dis) max_dis = dis;
					if (dis < min_dis) min_dis = dis;
				}
				new_subpopulation.extended_distance = max_dis;
				new_subpopulation.min_internal_distance = min_dis;
				if (enable_recursive_clearing) {
					for (size_t subpop_counter = 0; subpop_counter < subpopulation_size; ++subpop_counter) {
						PAEHN::CleanInTheRadius(dis_rank[dis_indices[subpop_counter]], max_dis);
					}
				}
				//new_subpopulation.furthest_distance = DBL_MIN;
				double temp_dis;
				for (size_t np_counter = 0; np_counter < np; ++np_counter) {
					if (Basic::lv_->population_[np_counter].is_choosen2_) {
						temp_dis = PublicFunction::CalculateEuclideanDistance(Basic::lv_->population_[np_counter].GetDimensions(), new_subpopulation.centroid, dim);
						if (new_subpopulation.extended_distance < temp_dis) {
							new_subpopulation.extended_distance = temp_dis;
						}
						Basic::lv_->population_[np_counter].is_choosen2_ = false;
					}
				}
				++form_how_many_subpops_counter;
			}

			for (size_t np_counter = 0; np_counter < np; ++np_counter) {
				Basic::lv_->population_[np_counter].is_choosen_ = false;
			}

			PAEHN::OutputSubpop();
		}

		void PAEHN::FormNewPopulationNearArchiveIndividuals(bool reinitializing_crowding_promising_archive) {
			if (Basic::lv_->paehn_subpopulations_.size() == 0) {
				return;
			}

			PAEHN::OutputArchive();
			const size_t archive_np = Basic::lv_->paehn_archive_.GetNowSize();
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			const size_t subpop_np = Basic::lv_->ev_.GetPAEHN_SubpopulationSize();
			std::vector<size_t> reinitial_archive_index;

			PublicFunction::SortPopulation(Basic::lv_->paehn_archive_fitness_, Basic::lv_->paehn_archive_.np_rank_, archive_np);
			double form_radius = 0.0;
			for (auto subpop = Basic::lv_->paehn_subpopulations_.begin(); subpop != Basic::lv_->paehn_subpopulations_.end(); ++subpop) {
				form_radius += subpop->extended_distance;
			}
			form_radius /= Basic::lv_->paehn_subpopulations_.size();
			const double worst_fitness_sharing = Basic::lv_->paehn_archive_fitness_[Basic::lv_->paehn_archive_.np_rank_[archive_np - 1]];
			double total = 0.0;
			for (size_t archive_np_counter = 0; archive_np_counter < archive_np; ++archive_np_counter) {
				total += worst_fitness_sharing - Basic::lv_->paehn_archive_fitness_[archive_np_counter];
			}
			std::vector<double> choosing_probability;
			for (size_t archive_np_counter = 0; archive_np_counter < archive_np; ++archive_np_counter) {
				choosing_probability.push_back((worst_fitness_sharing - Basic::lv_->paehn_archive_fitness_[archive_np_counter]) / total);
				if (archive_np_counter != 0) {
					choosing_probability[archive_np_counter] += choosing_probability[archive_np_counter - 1]; // accumulate
				}
			}
			double now_rand;
			for (size_t np_counter = 0; np_counter < np; np_counter += subpop_np) {
				now_rand = PublicFunction::RandomDouble(0, 1);
				double form_proportion = 0.0;
				bool is_positive = 0;
				for (size_t archive_np_counter = 0; archive_np_counter < archive_np; ++archive_np_counter) {
					if (now_rand > choosing_probability[archive_np_counter])
						continue;
					Basic::lv_->population_[np_counter] = Basic::lv_->paehn_archive_[archive_np_counter];
					Basic::lv_->fitness_[np_counter] = Basic::lv_->paehn_archive_fitness_[archive_np_counter];
					for (size_t new_individuals_counter = 1; new_individuals_counter < subpop_np; ++new_individuals_counter) {
						// create individual near the archive individual
						form_proportion = PublicFunction::RandomDouble(0, 1);
						std::vector<double> position_vector;
						total = 0;
						for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
							position_vector.push_back(PublicFunction::RandomDouble(0, 1));
							total += fabs(position_vector[dim_counter]);
						}
						for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
							position_vector[dim_counter] /= total;
							Basic::lv_->population_[np_counter + new_individuals_counter][dim_counter] = Basic::lv_->paehn_archive_[archive_np_counter][dim_counter];
							is_positive = PublicFunction::RandomInt(0, 1) == 1;
							if (is_positive) {
								Basic::lv_->population_[np_counter + new_individuals_counter][dim_counter] += sqrt(pow(form_radius * form_proportion, 2) * position_vector[dim_counter]);
							}
							else {
								Basic::lv_->population_[np_counter + new_individuals_counter][dim_counter] -= sqrt(pow(form_radius * form_proportion, 2) * position_vector[dim_counter]);
							}
							if (ss_->IsOutBound(Basic::lv_->population_[np_counter + new_individuals_counter][dim_counter], dim_counter)) {
								// prevent out bound
								ss_->OutBoundHandler(Basic::lv_->paehn_archive_, archive_np_counter, &Basic::lv_->population_[np_counter + new_individuals_counter][0], dim_counter);
							}
						}
						Basic::lv_->fitness_[np_counter + new_individuals_counter] = Basic::lv_->fp_->compute(Basic::lv_->population_[np_counter + new_individuals_counter].GetDimensions());
						++Basic::lv_->now_ffe_;
					}
					if (reinitializing_crowding_promising_archive) {
						if (PAEHN::IsArchiveRepeat(archive_np_counter)) {
							if (std::find(reinitial_archive_index.begin(), reinitial_archive_index.end(), archive_np_counter) == reinitial_archive_index.end())
								reinitial_archive_index.push_back(archive_np_counter);
						}
					}

					break; // after creating individuals 
				}
			}
			for (auto index = reinitial_archive_index.begin(); index != reinitial_archive_index.end(); ++index) {
				for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter)
					Basic::lv_->paehn_archive_[*index][dim_counter] = PublicFunction::RandomDouble(Basic::lv_->ev_.GetBoundMin(dim_counter), Basic::lv_->ev_.GetBoundMax(dim_counter));
				Basic::lv_->paehn_archive_fitness_[*index] = Basic::lv_->fp_->compute(Basic::lv_->paehn_archive_[*index].GetDimensions());
				++Basic::lv_->now_ffe_;
			}
			PAEHN::OutputPopulation();
		}

		void PAEHN::SuccessSelectionSHADE_SubpopMemoryCR() {
			if (Basic::lv_->paehn_subpopulations_.size() == 0) return;
			auto now_subpop = Basic::lv_->paehn_subpopulations_.begin();
			auto child = (*now_subpop).children_.begin();
			size_t now_H = (*now_subpop).insert_H_;
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			const size_t H = Basic::lv_->ev_.GetPAEHN_H();

			PAEHN::CalculateSuccessDeltaFitness();
			while (now_subpop != Basic::lv_->paehn_subpopulations_.end()) {
				now_H = (*now_subpop).insert_H_;
				for (child = (*now_subpop).children_.begin(); child != (*now_subpop).children_.end(); ++child) {
					if ((*child).is_less_) {
						(*now_subpop).success_CR_.push_back((*child).GetCR());
					}
				}

				if (!(*now_subpop).success_CR_.empty()) {
					(*now_subpop).memory_CR_[now_H % H] = PublicFunction::WeightedArithmeticMean((*now_subpop).delta_fitness_, (*now_subpop).success_CR_);
					(*now_subpop).insert_H_plus_1_ = true;
				}
				++now_subpop;
			}
		}

		void PAEHN::SuccessSelectionSHADE_SubpopMemoryF() {
			if (Basic::lv_->paehn_subpopulations_.size() == 0) return;
			auto now_subpop = Basic::lv_->paehn_subpopulations_.begin();
			auto child = (*now_subpop).children_.begin();
			size_t now_H = (*now_subpop).insert_H_;
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			const size_t H = Basic::lv_->ev_.GetPAEHN_H();

			PAEHN::CalculateSuccessDeltaFitness();
			while (now_subpop != Basic::lv_->paehn_subpopulations_.end()) {
				now_H = (*now_subpop).insert_H_;
				for (child = (*now_subpop).children_.begin(); child != (*now_subpop).children_.end(); ++child) {
					if ((*child).is_less_) {
						(*now_subpop).success_F_.push_back((*child).GetF());
					}
				}

				if (!(*now_subpop).success_F_.empty()) {
					(*now_subpop).memory_F_[now_H % H] = PublicFunction::WeightedLenmerMean((*now_subpop).delta_fitness_, (*now_subpop).success_F_);
					(*now_subpop).insert_H_plus_1_ = true;
				}
				++now_subpop;
			}
		}

		void PAEHN::CalculateSuccessDeltaFitness() {
			if (Basic::lv_->paehn_is_subpop_delta_fitness_calculated_) return;
			size_t parent = 0;
			const size_t next_gen_subpop_np = Basic::lv_->ev_.GetPAEHN_SubpopulationSize();
			auto now_subpop = Basic::lv_->paehn_subpopulations_.begin();

			while (now_subpop != Basic::lv_->paehn_subpopulations_.end()) {
				for (size_t next_gen_np_counter = 0; next_gen_np_counter < next_gen_subpop_np; ++next_gen_np_counter) {
					if ((*now_subpop).children_[next_gen_np_counter].is_less_) {
						parent = (*now_subpop).children_[next_gen_np_counter].GetComparedParent();
						(*now_subpop).delta_fitness_.push_back((*now_subpop).copied_fitness_[parent] - (*now_subpop).fitness_[next_gen_np_counter]);
					}
				}
				++now_subpop;
			}
			Basic::lv_->paehn_is_subpop_delta_fitness_calculated_ = true;

		}

		void PAEHN::RandomRi() {
			if (Basic::lv_->paehn_is_ri_produced_) return;
			const size_t next_gen_subpop_np = Basic::lv_->ev_.GetPAEHN_SubpopulationSize();
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			const size_t H = Basic::lv_->ev_.GetPAEHN_H();
			auto now_subpop = Basic::lv_->paehn_subpopulations_.begin();
			while (now_subpop != Basic::lv_->paehn_subpopulations_.end()) {
				for (size_t subpop_np_counter = 0; subpop_np_counter < next_gen_subpop_np; ++subpop_np_counter) {
					(*now_subpop).random_ri_[subpop_np_counter] = PublicFunction::RandomInt(0, H - 1);
				}
				++now_subpop;
			}
			Basic::lv_->paehn_is_ri_produced_ = true;

		}

		void PAEHN::SubpopulationSelection(bool random_promising_archive) {
			if (Basic::lv_->paehn_subpopulations_.size() == 0) return;
			const size_t next_gen_subpop_np = Basic::lv_->ev_.GetPAEHN_SubpopulationSize();
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			const size_t stagnation_upper = Basic::lv_->ev_.GetPAEHN_StagnationUpper();
			auto now_subpop = Basic::lv_->paehn_subpopulations_.begin();
			double temp_fitness = 0.0;
			bool is_improving = false;
			double fitness_difference = 0.0;

			while (now_subpop != Basic::lv_->paehn_subpopulations_.end()) {
				if (now_subpop->is_stagnation) {
					++now_subpop;
					continue;
				}
				size_t stagnation_threshold = ceil(((Basic::lv_->worst_fitness_ever_ - now_subpop->GetBestFitness()) / (Basic::lv_->worst_fitness_ever_ - Basic::lv_->best_fitness_ever_)) * stagnation_upper);
				is_improving = false;
				for (size_t subpop_np_counter = 0; subpop_np_counter < next_gen_subpop_np; ++subpop_np_counter) {
					now_subpop->children_[subpop_np_counter].InitialComparisonResult(); // individual
					temp_fitness = Basic::lv_->fp_->compute(now_subpop->children_[subpop_np_counter].GetDimensions()); // individual
					++Basic::lv_->now_ffe_;
					(*now_subpop)[subpop_np_counter].SetFitnessDifference(temp_fitness - (*now_subpop).fitness_[subpop_np_counter]); // individual fitness
					now_subpop->children_[subpop_np_counter].SetComparedParent(subpop_np_counter);

					if (temp_fitness < Basic::lv_->best_fitness_ever_) Basic::lv_->best_fitness_ever_ = temp_fitness;
					else if (temp_fitness > Basic::lv_->worst_fitness_ever_) Basic::lv_->worst_fitness_ever_ = temp_fitness;

					if (temp_fitness < (*now_subpop).fitness_[subpop_np_counter]) {
						// if the new individual is better than its parent
						now_subpop->children_[subpop_np_counter].is_less_ = true;

						if (now_subpop->archive_.size() < 4) {
							now_subpop->archive_.emplace_back(Basic::lv_->ev_.GetDE_Dim());
							now_subpop->archive_.back() = (*now_subpop)[subpop_np_counter].GetDimensions();
						}
						else {
							size_t index = PublicFunction::RandomInt(0, 4 - 1);
							now_subpop->archive_[index] = (*now_subpop)[subpop_np_counter].GetDimensions();
						}

						PAEHN::SuccessSelectionCrowdingArchive((*now_subpop)[subpop_np_counter].GetDimensions(), (*now_subpop).fitness_[subpop_np_counter], random_promising_archive);
						fitness_difference = (*now_subpop).fitness_[subpop_np_counter] - temp_fitness;
						(*now_subpop)[subpop_np_counter] = now_subpop->children_[subpop_np_counter].GetDimensions();
						(*now_subpop).fitness_[subpop_np_counter] = temp_fitness;
						if (fitness_difference > kAccuracyLevel)
							is_improving = true;
						if (temp_fitness < (*now_subpop).GetBestFitness()) {
							(*now_subpop).SetBestFitness(temp_fitness);
							(*now_subpop).best_individual = subpop_np_counter;
						}
						if (temp_fitness < Basic::lv_->best_fitness_ever_)
							Basic::lv_->best_fitness_ever_ = temp_fitness;
					}
					else if (temp_fitness == (*now_subpop).fitness_[subpop_np_counter])
						now_subpop->children_[subpop_np_counter].is_equal_ = true;
					else {
						now_subpop->children_[subpop_np_counter].is_larger_ = true;
						if (temp_fitness > Basic::lv_->worst_fitness_ever_)
							Basic::lv_->worst_fitness_ever_ = temp_fitness;
					}
				}

				if (is_improving) {
					(*now_subpop).stagnation_counter = 0;
				}
				else {
					++((*now_subpop).stagnation_counter);
					if (now_subpop->stagnation_counter == stagnation_threshold) {
						PAEHN::StagnationSubpop(now_subpop);
					}
				}
				++now_subpop;
			}
		}

		void PAEHN::NextGen() {
			auto now_subpop = Basic::lv_->paehn_subpopulations_.begin();
			for (; now_subpop != Basic::lv_->paehn_subpopulations_.end(); ++now_subpop) {
				if (now_subpop->insert_H_plus_1_) {
					++now_subpop->insert_H_;
					now_subpop->insert_H_plus_1_ = false;
				}
				now_subpop->delta_fitness_.clear();
				now_subpop->success_CR_.clear();
				now_subpop->success_F_.clear();
			}
			Basic::lv_->paehn_is_ri_produced_ = false;
			Basic::lv_->paehn_is_subpop_fitness_copied_ = false;
			Basic::lv_->paehn_is_subpop_delta_fitness_calculated_ = false;
		}

		void PAEHN::CleanInTheRadius(size_t base_individual, double distance) {
			const size_t np = Basic::lv_->ev_.GetDE_NP();
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			for (size_t np_counter = 0; np_counter < np; ++np_counter) {
				if (Basic::lv_->population_[np_counter].is_choosen2_) continue;
				if (PublicFunction::CalculateEuclideanDistance(Basic::lv_->population_[np_counter].GetDimensions(), Basic::lv_->population_[base_individual].GetDimensions(), dim) <= distance) {
					if (Basic::lv_->fitness_[np_counter] > Basic::lv_->fitness_[base_individual]) {
						Basic::lv_->population_[np_counter].is_choosen_ = true;
						Basic::lv_->population_[np_counter].is_choosen2_ = true;
						PAEHN::CleanInTheRadius(np_counter, distance);
					}
				}
			}
		}

		void PAEHN::StagnationSubpop(const std::vector<paehn::Subpopulation>::iterator& now_subpop) {
			now_subpop->is_stagnation = true;
		}

		bool PAEHN::AreAllSubpopStagnation() {
			for (auto subpop = Basic::lv_->paehn_subpopulations_.begin(); subpop != Basic::lv_->paehn_subpopulations_.end(); ++subpop) {
				if (subpop->is_stagnation == false)
					return false;
			}
			return true;
		}

		bool PAEHN::IsArchiveRepeat(const size_t base_index) {
			if (Basic::lv_->paehn_subpopulations_.size() == 0) {
				return false;
			}
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			size_t repeat_counter = 0;
			double radius = 0.0;
			for (auto subpop = Basic::lv_->paehn_subpopulations_.begin(); subpop != Basic::lv_->paehn_subpopulations_.end(); ++subpop) {
				radius += subpop->min_internal_distance;
			}
			radius /= Basic::lv_->paehn_subpopulations_.size();

			for (size_t archive_np_counter = 0; archive_np_counter < Basic::lv_->paehn_archive_.GetNowSize(); ++archive_np_counter) {
				if (archive_np_counter == base_index) continue;
				if (Basic::lv_->paehn_archive_fitness_[archive_np_counter] > Basic::lv_->paehn_archive_fitness_[base_index]) continue;
				if (PublicFunction::CalculateEuclideanDistance(Basic::lv_->paehn_archive_[archive_np_counter].GetDimensions(), Basic::lv_->paehn_archive_[base_index].GetDimensions(), dim) < radius) {
					return true;
				}
			}
			return false;
		}

		void PAEHN::ShowAllSubpopBest(std::vector<std::vector<double>>& result) {
			std::cout << "  SUBPOP size: " << Basic::lv_->paehn_subpopulations_.size();
			auto subpop = Basic::lv_->paehn_subpopulations_.begin();
			const size_t subpop_np = Basic::lv_->ev_.GetPAEHN_SubpopulationSize();
			PAEHN::OutputSubpop();
			while (subpop != Basic::lv_->paehn_subpopulations_.end()) {
				for (size_t subpop_np_counter = 0; subpop_np_counter < subpop_np; ++subpop_np_counter) {
					if ((*subpop).fitness_[subpop_np_counter] - Basic::lv_->best_fitness_ever_ < 0.11) {
						result.emplace_back((*subpop)[subpop_np_counter].begin(), (*subpop)[subpop_np_counter].end());
					}
				}
				++subpop;
			}
			for (size_t np_counter = 0; np_counter < Basic::lv_->ev_.GetDE_NP(); ++np_counter) {
				result.emplace_back(Basic::lv_->population_[np_counter].begin(), Basic::lv_->population_[np_counter].end());
			}
			for (size_t archive_counter = 0; archive_counter < Basic::lv_->paehn_archive_.GetNowSize(); ++archive_counter) {
				result.emplace_back(Basic::lv_->paehn_archive_[archive_counter].begin(), Basic::lv_->paehn_archive_[archive_counter].end());
			}
		}

		void PAEHN::InitializeOfstream() {
			if (!kOutput) return;
			string archive_output = Basic::lv_->file_path_ + "/" + Basic::lv_->filename_ + "Archive.csv";
			PAEHN::file_streams_.emplace_back(archive_output, std::ios_base::app);
			PAEHN::file_streams_.emplace_back(Basic::lv_->file_path_ + "/" + Basic::lv_->filename_ + "BFPopulation.csv");
			PAEHN::file_streams_.emplace_back(Basic::lv_->file_path_ + "/" + Basic::lv_->filename_ + "BFArchive.csv");
		}

		void PAEHN::OutputPopulation() {
			if (!kOutput) return;
			const size_t np = Basic::lv_->population_.GetNP();
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			for (size_t np_counter = 0; np_counter < np; ++np_counter) {
				for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
					Basic::lv_->file_stream_ << Basic::lv_->population_[np_counter][dim_counter] << ",";
				}
				Basic::lv_->file_stream_ << Basic::lv_->fitness_[np_counter];
				if (isnan(Basic::lv_->fitness_[np_counter]))
					std::cout << "DIE";
				if (Basic::lv_->paehn_is_this_gen_create_subpopulations_)
					Basic::lv_->file_stream_ << "," << Basic::lv_->paehn_form_subpopulations_gen_counter_ << std::endl;
				else
					Basic::lv_->file_stream_ << ",0" << std::endl;
			}
			Basic::lv_->file_stream_.flush();
		}

		void PAEHN::OutputBeforeAfterPopulation() {
			if (!kOutput) return;
			const size_t np = Basic::lv_->population_.GetNP();
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			for (size_t np_counter = 0; np_counter < np; ++np_counter) {
				for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
					PAEHN::file_streams_[1] << Basic::lv_->population_[np_counter][dim_counter] << ",";
				}
				PAEHN::file_streams_[1] << Basic::lv_->fitness_[np_counter] << std::endl;
			}
			PAEHN::file_streams_[1].flush();
		}

		void PAEHN::OutputBeforAfterArchive() {
			if (!kOutput) return;
			const size_t np = Basic::lv_->paehn_archive_.GetNowSize();
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			for (size_t archive_np_counter = 0; archive_np_counter < np; ++archive_np_counter) {
				for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
					PAEHN::file_streams_[2] << Basic::lv_->paehn_archive_[archive_np_counter][dim_counter] << ",";
				}
				PAEHN::file_streams_[2] << Basic::lv_->paehn_archive_fitness_[archive_np_counter] << std::endl;
			}
			PAEHN::file_streams_[2].flush();
		}

		void PAEHN::OutputArchive() {
			if (!kOutput) return;
			const size_t np = Basic::lv_->paehn_archive_.GetNowSize();
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			for (size_t archive_counter = 0; archive_counter < np; ++archive_counter) {
				for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
					PAEHN::file_streams_[0] << Basic::lv_->paehn_archive_[archive_counter][dim_counter] << ",";
				}
				PAEHN::file_streams_[0] << Basic::lv_->paehn_archive_fitness_[archive_counter] << std::endl;
			}
			PAEHN::file_streams_[0].flush();
		}

		void PAEHN::OutputSubpop() {
			if (!kOutput) return;
			string subpop_path = Basic::lv_->file_path_ + "/" + Basic::lv_->filename_ + "subpop" + std::to_string(++PAEHN::subpop_output_counter) + ".csv";
			string centroid_path = Basic::lv_->file_path_ + "/" + Basic::lv_->filename_ + "centroid" + std::to_string(PAEHN::subpop_output_counter) + ".csv";
			ofstream subpop_file_stream(subpop_path, std::ios_base::app);
			ofstream centroid_file_stream(centroid_path, std::ios_base::app);
			const size_t next_gen_subpop_np = Basic::lv_->ev_.GetPAEHN_SubpopulationSize();
			const size_t dim = Basic::lv_->ev_.GetDE_Dim();
			for (auto subpop = Basic::lv_->paehn_subpopulations_.begin(); subpop != Basic::lv_->paehn_subpopulations_.end(); ++subpop) {
				for (size_t subpop_np_counter = 0; subpop_np_counter < next_gen_subpop_np; ++subpop_np_counter) {
					for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
						subpop_file_stream << (*subpop)[subpop_np_counter][dim_counter] << ",";
					}
					subpop_file_stream << subpop->fitness_[subpop_np_counter];
					if (subpop->is_new_)
						subpop_file_stream << ",1" << std::endl;
					else
						subpop_file_stream << ",0" << std::endl;
				}
				for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
					centroid_file_stream << subpop->centroid[dim_counter] << ",";
				}
				if (subpop->is_new_)
					centroid_file_stream << "1" << std::endl;
				else
					centroid_file_stream << "0" << std::endl;
				subpop->is_new_ = false;
			}
			subpop_file_stream.flush();
			subpop_file_stream.close();

			centroid_file_stream.flush();
			centroid_file_stream.close();
		}
	} // namespace differentialevolution
} // namespace evolutionalgorithm