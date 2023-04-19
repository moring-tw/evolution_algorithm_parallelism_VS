#include "../../../headers/differential_evolution/variants/basic.h"
namespace evolutionalgorithm {
	namespace differentialevolution {

		Basic::Basic()
		{
		}

		Basic::~Basic()
		{
		}

		void Basic::SetLocalVariables(LocalVariables* lv) {
			if (lv == nullptr) {
				std::cout << "Assignment of lv_ failed(passed nullptr)." << std::endl;
				exit(0);
			}
			Basic::lv_ = lv;
			if (Basic::ss_ != nullptr) {
				Basic::ss_->SetLocalVariables((*Basic::lv_));
				Basic::ss_->SetExecutiveVariables(Basic::lv_->ev_);
			} else {
				printf("Before set Locak Variables, Strategies must be set first.");
				exit(0);
			}
		}

		void Basic::SetStrategies(Strategies* ss) {
			Basic::ss_ = ss;
		}

		void Basic::NextRun() {

		}

		void Basic::NextGen() {

		}

		void Basic::CheckLV() {
			if (Basic::lv_ == nullptr) {
				std::cout << "Before using DE object, lv_ must be assigned by DE::SetLocalVariables." << std::endl;
				exit(0);
			}
		}

		void Basic::CheckFP() {
			Basic::CheckLV();
			if (Basic::lv_->fp_ == nullptr) {
				std::cout << "Before using DE object, fp_ must be assigned by DE::SetBenchmarks." << std::endl;
				exit(0);
			}
		}

		void Basic::StrategiesExecuter(Population& population, Individual& next_individual, const size_t base_individual, corematerial::MutantAndCrossover using_strategy, size_t number_of_vector_differences) {
			Basic::CheckLV();
			if ((using_strategy & corematerial::MutantAndCrossover::kRand) == corematerial::MutantAndCrossover::kRand) {
				//std::cout << individual << " using:Rand\n";
				Basic::ss_->DE_Rand(
					population,
					base_individual,
					next_individual,
					number_of_vector_differences
				);
			} else if ((using_strategy & corematerial::MutantAndCrossover::kBest) == corematerial::MutantAndCrossover::kBest) {
				//std::cout << individual << " using:Best\n";
				Basic::ss_->DE_Best(
					population,
					base_individual,
					next_individual,
					number_of_vector_differences
				);
			} else if ((using_strategy & corematerial::MutantAndCrossover::kCurToRand) == corematerial::MutantAndCrossover::kCurToRand) {
				//std::cout << individual << " using:Cur-to-Rand\n";
				Basic::ss_->DE_Cur_To_Rand(
					population,
					base_individual,
					next_individual,
					number_of_vector_differences
				);
			} else if ((using_strategy & corematerial::MutantAndCrossover::kCurToBest) == corematerial::MutantAndCrossover::kCurToBest) {
				//std::cout << individual << " using:Cur-to-Best\n";
				Basic::ss_->DE_Cur_To_Best(
					population,
					base_individual,
					next_individual,
					number_of_vector_differences
				);
			} 
			//else if ((using_strategy & corematerial::MutantAndCrossover::kCurToPBest) == corematerial::MutantAndCrossover::kCurToPBest) {
			//	//std::cout << individual << " using:Cur-to-PBest\n";
			//	next_individual = Basic::ss_->DE_Cur_To_PBest(
			//		population,
			//		base_individual,
			//		Basic::lv_->archive_, 
			//		Basic::lv_->shade_archive_now_size_, 
			//		next_individual,
			//		number_of_vector_differences
			//	);
			//}
			else if ((using_strategy & corematerial::MutantAndCrossover::UNSET) == corematerial::MutantAndCrossover::UNSET) {
				std::cout << "Before executing strategy, using_strategy must had been set." << std::endl;
				exit(0);
			} else {
				std::cout << "Unknown strategy." << std::endl;
				exit(0);
			}
		}

		void Basic::StrategiesExecuter(Population& population, Individual& next_individual, Population& archive, const size_t base_individual, corematerial::MutantAndCrossover using_strategy, size_t number_of_vector_differences) {
			if ((using_strategy & corematerial::MutantAndCrossover::kCurToPBest) == corematerial::MutantAndCrossover::kCurToPBest) {
				Basic::ss_->DE_Cur_To_PBest(
					population,
					base_individual,
					archive,
					next_individual,
					number_of_vector_differences
				);
			}
			else if ((using_strategy & corematerial::MutantAndCrossover::UNSET) == corematerial::MutantAndCrossover::UNSET) {
				std::cout << "Before executing strategy, using_strategy must had been set." << std::endl;
				exit(0);
			}
			else {
				std::cout << "Unknown strategy." << std::endl;
				exit(0);
			}
		}

		size_t Basic::FindStrategyIndex(const corematerial::MutantAndCrossover strategy) {
			std::map<corematerial::MutantAndCrossover, size_t>* strategies_index_map = Basic::lv_->ev_.GetDE_StrategiesIndexMap();
			std::map<corematerial::MutantAndCrossover, size_t>::iterator it = strategies_index_map->find(strategy);
			if (it != strategies_index_map->end()) {
				return it->second;
			} else Basic::UnknownStrategyError(strategy);
			return 0;
		}
	}
}