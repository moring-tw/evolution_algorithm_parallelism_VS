#include "../../headers/differential_evolution/process_controller.h"

namespace evolutionalgorithm {
	namespace differentialevolution {
		ProcessController::ProcessController() {

		}
		ProcessController::~ProcessController() {

		}

		void ProcessController::SetBehaviors(Behavior& b) {
			b.SetDE(&(ProcessController::de_));
			b.SetSaDE(&(ProcessController::sade_));
			b.SetSHADE(&(ProcessController::shade_));
			b.SetRandom(&(ProcessController::random_));
			b.SetPAEHN(&(ProcessController::paehn_));
			b.SetLocalVariables(&(ProcessController::lv_));
		}

		void ProcessController::SetVariants(Basic& b) {
			b.SetStrategies(&(ProcessController::ss_));
			b.SetLocalVariables(&(ProcessController::lv_));
		}

		void ProcessController::NextRun() {
			if (!ProcessController::is_initialed_) {
				ProcessController::SetVariants(ProcessController::de_); // execution instance
				ProcessController::SetVariants(ProcessController::sade_); // execution instance
				ProcessController::SetVariants(ProcessController::shade_); // execution instance
				ProcessController::SetVariants(ProcessController::random_); // execution instance
				ProcessController::SetVariants(ProcessController::paehn_); // execution instance
				ProcessController::SetBehaviors(ProcessController::memory_); // manipulate memories (population, fitness and archive etc.)
				ProcessController::SetBehaviors(ProcessController::parameter_control_); // manipulate algorithms parameters (CR, F etc).
				ProcessController::SetBehaviors(ProcessController::evolve_); // evolve population
				ProcessController::SetBehaviors(ProcessController::select_); // selection population
				ProcessController::SetBehaviors(ProcessController::sample_); // sample population
				ProcessController::is_initialed_ = true;
			}
			ProcessController::de_.NextRun();
			ProcessController::sade_.NextRun();
			ProcessController::shade_.NextRun();
			ProcessController::paehn_.NextRun();
			ProcessController::lv_.NextRun();
		}

		void ProcessController::Execute() {
			ProcessController::BeforeInitial();
			ProcessController::Initial();
			ProcessController::AfterInitial();
			while (!ProcessController::ReachStopCriterion()) {
				//printf("now gen: %d\n", lv_.now_gen_);
				ProcessController::BeforeEvolution();
				ProcessController::Evolution();
				ProcessController::AfterEvolution();

				ProcessController::BeforeSelection();
				ProcessController::Selection();
				ProcessController::AfterSelection();
				ProcessController::NextGen();
				if (ProcessController::debug_)
					std::cout << "temporarily stop here." << std::endl;
			}
			std::vector<std::vector<double>> result;
			paehn_.ShowAllSubpopBest(result);//--------------------------------------------------------------------
			ProcessController::lv_.total_found_optima_.clear();

			for (double accuracy = 0.1; accuracy >= 0.00001; accuracy /= 10) {
				std::vector<std::vector<double>> seeds;
				ProcessController::lv_.total_found_optima_.push_back(how_many_goptima(result, seeds, (CEC2013*)ProcessController::lv_.fp_, accuracy, ((CEC2013*)ProcessController::lv_.fp_)->get_rho()));
				/*cout << accuracy << ": "
					<< ProcessController::lv_.total_found_optima_.back()
					<< "\t";*/
			}
			ProcessController::RecycleMemory();
			//system("pause");
		}

		void ProcessController::SetBenchmark(Benchmarks* fp) {
			ProcessController::lv_.fp_ = fp;
		}

		void ProcessController::SetExecutiveVariables(ExecutiveVariables& ev) {
			ProcessController::lv_.ev_ = ev;
		}

		double ProcessController::GetBest() {
			return ProcessController::lv_.population_.GetBestFitness();
		}

		std::vector<size_t>& ProcessController::GetHowManyOptima() {
			return ProcessController::lv_.total_found_optima_;
		}

		void ProcessController::SetOutputFileStream(const string& filepath, const string& filename, const string& file_ext)
		{
			ProcessController::lv_.file_stream_.close();
			ProcessController::lv_.file_stream_.open(filepath + "/" + filename + "." + file_ext, std::ios_base::app);
			ProcessController::lv_.file_path_ = filepath;
			ProcessController::lv_.filename_ = filename;
		}

		void ProcessController::BeforeInitial() { // malloc

			// Strategies Selection
			if (ProcessController::lv_.ev_.GetDE_NumOfStrategies() >= 2) {
				switch (ProcessController::lv_.ev_.GetDE_StrategiesSelectionPattern()) {
				case corematerial::StrategiesSelection::kSaDE_ss:
					if (ProcessController::debug_)
						std::cout << "MallocSaDE_StrategiesSelection" << std::endl;
					memory_.MallocSaDE_StrategiesSelection();
					break;
				default:
					std::cout << "Unknown Strategies selection pattern." << std::endl;
					exit(0);
				}
			}
			// parameter control CR
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kCR) == corematerial::ParametersControl::kCR) {
				// memory CR
				if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kMemoryCR) == corematerial::ParametersControlMethod::kMemoryCR) {
					switch (ProcessController::lv_.ev_.GetDE_MemoryCR_Pattern()) {
					case(corematerial::MemoryCR::kSaDE_CR):
						if ((ProcessController::lv_.ev_.GetDE_MemoryCR_SettingPattern() & corematerial::MemoryCR_Setting::kMemoryByEachStrategy) == corematerial::MemoryCR_Setting::kMemoryByEachStrategy) {
							if (ProcessController::debug_)
								std::cout << "MallocSaDE_MemoryCR(N)" << std::endl;
							parameter_control_.MallocSaDE_MemoryCR(ProcessController::lv_.ev_.GetDE_NumOfStrategies());
						}
						else {
							if (ProcessController::debug_)
								std::cout << "MallocSaDE_MemoryCR(1)" << std::endl;
							parameter_control_.MallocSaDE_MemoryCR(1);
						}
						break;
					case(corematerial::MemoryCR::kSHADE_CR):
						parameter_control_.MallocSHADE_MemoryCR();
						break;
					case corematerial::MemoryCR::kPAEHN_SHADE_CR:
						break;
					default:
						std::cout << "Unknown memory CR pattern." << std::endl;
						exit(0);
					}
				}
			}
			//parameter control F
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kF) == corematerial::ParametersControl::kF) {
				// memory F
				if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kMemoryF) == corematerial::ParametersControlMethod::kMemoryF) {
					switch (ProcessController::lv_.ev_.GetDE_MemoryF_Pattern()) {
					case(corematerial::MemoryF::kSHADE_F):
						parameter_control_.MallocSHADE_MemoryF();
						break;
					case corematerial::MemoryF::kPAEHN_SHADE_F:
							break;
					default:
						printf("Unknown memory F parrern.");
						exit(0);
					}
				}
			}
		}

		void ProcessController::Initial() { // initial value
			switch (ProcessController::lv_.ev_.GetDE_InitialPopulationPattern()) {
			case corematerial::InitialPopulation::kRandom:
				ProcessController::memory_.InitialDE_Population();
				ProcessController::memory_.InitialDE_NextGeneration();
				ProcessController::memory_.InitialDE_Fitness();
				break;
			case corematerial::InitialPopulation::kLHD:
				break;
			default:
				std::cout << "Unknown initial pattern of DE." << std::endl;
				exit(0);
			}
			// np rank
			if ((ProcessController::lv_.ev_.GetDE_Memory_Pattern() & corematerial::Memory::kNPRank) == corematerial::Memory::kNPRank) {
				memory_.InitialNP_Rank();
			}
		}

		void ProcessController::AfterInitial() {
			// Strategies Selection
			if (ProcessController::lv_.ev_.GetDE_NumOfStrategies() == 1) {
				memory_.InitialDE_Strategy();
			}
			else if (ProcessController::lv_.ev_.GetDE_NumOfStrategies() >= 2) {
				switch (ProcessController::lv_.ev_.GetDE_StrategiesSelectionPattern())
				{
				case corematerial::StrategiesSelection::kSaDE_ss:
					if (ProcessController::debug_)
						std::cout << "InitialSaDE_StrategiesSelection" << std::endl;
					memory_.InitialSaDE_StrategiesSelection();
					break;
				default:
					std::cout << "Unknown Strategies selection pattern." << std::endl;
					exit(0);
				}
			}

			// parameter control CR
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kCR) != corematerial::ParametersControl::kCR) {
				// fixed CR
				parameter_control_.InitialDE_CR();
			}
			else {
				// Memory CR
				if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kMemoryCR) == corematerial::ParametersControlMethod::kMemoryCR) {
					switch (ProcessController::lv_.ev_.GetDE_MemoryCR_Pattern()) {
					case(corematerial::MemoryCR::kSaDE_CR):
						// SaDE 
						if ((ProcessController::lv_.ev_.GetDE_MemoryCR_SettingPattern() & corematerial::MemoryCR_Setting::kMemoryByEachStrategy) == corematerial::MemoryCR_Setting::kMemoryByEachStrategy) {
							if (ProcessController::debug_)
								std::cout << "InitialSaDE_MemoryCR(N)" << std::endl;
							parameter_control_.InitialSaDE_MemoryCR(ProcessController::lv_.ev_.GetDE_NumOfStrategies());
						}
						else {
							if (ProcessController::debug_)
								std::cout << "InitialSaDE_MemoryCR(1)" << std::endl;
							parameter_control_.InitialSaDE_MemoryCR(1);
						}
						break;
					case(corematerial::MemoryCR::kSHADE_CR):
						// SHADE
						parameter_control_.InitialSHADE_MemoryCR();
						break;
					case corematerial::MemoryCR::kPAEHN_SHADE_CR:
						break;
					default:
						std::cout << "Unknown memory CR pattern." << std::endl;
						exit(0);
					}
				}
			}

			// parameter control F
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kF) != corematerial::ParametersControl::kF) {
				// fixed F
				parameter_control_.InitialDE_F();
			}
			else {
				if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kMemoryF) == corematerial::ParametersControlMethod::kMemoryF) {
					switch (ProcessController::lv_.ev_.GetDE_MemoryF_Pattern()) {
					case(corematerial::MemoryF::kSHADE_F):
						// SHADE
						parameter_control_.InitialSHADE_MemoryF();
						break;
					case corematerial::MemoryF::kPAEHN_SHADE_F:
						break;
					default:
						std::cout << "Unknown memory F pattern." << std::endl;
						exit(0);
					}
				}
				else if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kNormalF) == corematerial::ParametersControlMethod::kNormalF) {
					// F with Normal distribution 
				}
			}
			// parameter control K
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kK) != corematerial::ParametersControl::kK) {
				// fixed K
				parameter_control_.InitialDE_K();
			}
			// parameter control P
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kP) != corematerial::ParametersControl::kP) {
				// fixed P
				parameter_control_.InitialDE_P();
			}

			// archive
			if ((ProcessController::lv_.ev_.GetDE_Memory_Pattern() & corematerial::Memory::kArchive) == corematerial::Memory::kArchive) {

				switch (ProcessController::lv_.ev_.GetDE_ArchivePattern()) {
				case(corematerial::Archive::kSHADE_Archive):
					ProcessController::memory_.InitialSHADE_Archive();
					break;
				case(corematerial::Archive::kPAEHN_Archive):
					ProcessController::memory_.InitialPAEHN_Archive();
					break;
				default:
					printf("Unknown archive pattern.");
					std::cout << (unsigned int)ProcessController::lv_.ev_.GetDE_ArchivePattern();
					exit(0);
				}
			}

			// NP Rank
			if ((ProcessController::lv_.ev_.GetDE_Memory_Pattern() & corematerial::Memory::kNPRank) == corematerial::Memory::kNPRank) {
				memory_.SortNP_Rank();
			}
		}

		void ProcessController::BeforeEvolution() {
			// Strategies Selection
			if (ProcessController::lv_.ev_.GetDE_NumOfStrategies() >= 2) {
				switch (ProcessController::lv_.ev_.GetDE_StrategiesSelectionPattern()) {
				case corematerial::StrategiesSelection::kSaDE_ss:
					if (ProcessController::debug_)
						std::cout << "BeforeEvolveSaDE_StrategiesSelection" << std::endl;
					evolve_.BeforeEvolveSaDE_StrategiesSelection();
					break;
				default:
					std::cout << "Unknown Strategies selection pattern when evolving." << std::endl;
					exit(0);
				}
			}
			else if (ProcessController::lv_.ev_.GetDE_NumOfStrategies() != 1) {
				std::cout << "num_of_strategies_ is less than 1, please check." << std::endl;
				exit(0);
			}
			// parameter control CR
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kCR) == corematerial::ParametersControl::kCR) {
				if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kMemoryCR) == corematerial::ParametersControlMethod::kMemoryCR) {
					// Memory CR
					switch (ProcessController::lv_.ev_.GetDE_MemoryCR_Pattern()) {
					case corematerial::MemoryCR::kSaDE_CR:
						// SaDE
						if ((ProcessController::lv_.ev_.GetDE_MemoryCR_SettingPattern() & corematerial::MemoryCR_Setting::kMemoryByEachStrategy) == corematerial::MemoryCR_Setting::kMemoryByEachStrategy) {
							if (ProcessController::debug_)
								std::cout << "BeforeEvolveSaDE_CR_Selection(N)" << std::endl;
							evolve_.BeforeEvolveSaDE_CR_Selection(ProcessController::lv_.ev_.GetDE_NumOfStrategies());
						}
						else {
							if (ProcessController::debug_)
								std::cout << "BeforeEvolveSaDE_CR_Selection(1)" << std::endl;
							evolve_.BeforeEvolveSaDE_CR_Selection(1);
						}
						break;
					case corematerial::MemoryCR::kSHADE_CR:
						// SHADE
						evolve_.BeforeEvolveSHADE_CR_Selection();
						break;
					case corematerial::MemoryCR::kPAEHN_SHADE_CR:
						// paehn subpopulation control CR
						paehn_.SelectSHADE_SubpopCR();//--------------------------------------------------------------------
						break;
					default:
						std::cout << "Unavailable de_memory_CR_pattern setting at Evaluation(). Please check." << std::endl;
						exit(0);
					}
				}
				else parameter_control_error("CR");
			}
			// parameter control F
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kF) == corematerial::ParametersControl::kF) {
				if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kMemoryF) == corematerial::ParametersControlMethod::kMemoryF) {
					// Memory F
					switch (ProcessController::lv_.ev_.GetDE_MemoryF_Pattern()) {
					case corematerial::MemoryF::kSHADE_F:
						// SHADE
						evolve_.BeforeEvolveSHADE_F_Selection();
						break;
					case corematerial::MemoryF::kPAEHN_SHADE_F:
						//subpopulation control F
						paehn_.SelectSHADE_SubpopF();//--------------------------------------------------------------------
						break;
					default:
						std::cout << "Unavailable de_memory_F_pattern setting at Evaluation(). Please check." << std::endl;
						exit(0);
					}
				}
				else if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kNormalF) == corematerial::ParametersControlMethod::kNormalF) {
					// generate F with Normal distribution
				}
				else parameter_control_error("F");
			}
			//parameter control K
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kK) == corematerial::ParametersControl::kK) {
				parameter_control_error("K");
			}
			// parameter control P
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kP) == corematerial::ParametersControl::kP) {
				if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kRandomP) == corematerial::ParametersControlMethod::kRandomP) {
					evolve_.BeforeEvolveDE_RandomP();
				}
				else parameter_control_error("P");
			}
		}

		void ProcessController::Evolution() {
			if (ProcessController::lv_.ev_.GetDE_NumOfStrategies() == 1) {
				if (ProcessController::debug_)
					std::cout << "EvolveDE" << std::endl;
				switch (ProcessController::lv_.ev_.GetSubpopulationPattern()) {
				case corematerial::Subpopulation::kPAEHN_Subpopulation:// subpopulation evolution
					paehn_.SubpopulationMutationAndCrossover(); //--------------------------------------------------------------------
					break;
				case corematerial::Subpopulation::UNSET:
					evolve_.EvolveDE();
					break;
				default:
					std::cout << "Unknown [subpopulation] setting, please check." << std::endl;
					exit(0);
				}
			}
			else if (ProcessController::lv_.ev_.GetDE_NumOfStrategies() >= 2) {
				switch (ProcessController::lv_.ev_.GetDE_StrategiesSelectionPattern()) {
				case corematerial::StrategiesSelection::kSaDE_ss:
					evolve_.EvolveDE();
					break;
				default:
					std::cout << "Unavailable de_strategies_selection_pattern_ setting. Please check." << std::endl;
					exit(0);
				}
			}
			else {
				std::cout << "num_of_strategies_ is less than 1, please check." << std::endl;
				exit(0);
			}

			
		}

		void ProcessController::AfterEvolution() {

		}

		void ProcessController::BeforeSelection() {
			if (ProcessController::lv_.ev_.GetDE_NumOfStrategies() >= 2) {
				switch (ProcessController::lv_.ev_.GetDE_StrategiesSelectionPattern()) {
				case corematerial::StrategiesSelection::kSaDE_ss:
					if (ProcessController::debug_)
						std::cout << "BeforeSelectSaDE_StrategiesSelection" << std::endl;
					ProcessController::select_.BeforeSelectSaDE_StrategiesSelection();
					break;
				}
			}
		}

		void ProcessController::Selection() {
			const corematerial::Selection s = ProcessController::lv_.ev_.GetDE_SelectionPattern();
			//
			// No metter which selection method is using here, individuals must record who is its parent.
			// Not allow any free or delete during selecton, always deallocates in AfterSelection().
			//
			if (s == corematerial::Selection::kCompetition) {
				if (ProcessController::debug_)
					std::cout << "SelectDE" << std::endl;
				switch (ProcessController::lv_.ev_.GetSubpopulationPattern()) {
				case corematerial::Subpopulation::kPAEHN_Subpopulation: {
					// subpopulation selection
					bool RandomPromisingArchive = false;
					if ((ProcessController::lv_.ev_.GetPAEHN_Pattern() & corematerial::PAEHN::kPAEHN_RandomPromisingArchive) == corematerial::PAEHN::kPAEHN_RandomPromisingArchive)
						RandomPromisingArchive = true;
					if ((ProcessController::lv_.ev_.GetSubpopulationPattern() & corematerial::Subpopulation::kPAEHN_Subpopulation) == corematerial::Subpopulation::kPAEHN_Subpopulation)
						paehn_.SubpopulationSelection(RandomPromisingArchive); //--------------------------------------------------------------------
					break;
				}
				case corematerial::Subpopulation::UNSET:
					select_.SelectDE();
					break;
				default:
					std::cout << "Unknown [subpopulation] setting, please check." << std::endl;
					exit(0);
				}
			}
		}

		void ProcessController::AfterSelection() {
			// np_rank
			if ((ProcessController::lv_.ev_.GetDE_Memory_Pattern() & corematerial::Memory::kNPRank) == corematerial::Memory::kNPRank) {
				memory_.SortNP_Rank();
			}

			// strategies selection
			if (ProcessController::lv_.ev_.GetDE_NumOfStrategies() >= 2) {
				switch (ProcessController::lv_.ev_.GetDE_StrategiesSelectionPattern()) {
				case corematerial::StrategiesSelection::kSaDE_ss:
					if (ProcessController::debug_)
						std::cout << "AfterSelectSaDE_SuccessStrategiesSelection" << std::endl;
					if (ProcessController::debug_)
						std::cout << "AfterSelectSaDE_FailStrategiesSelection" << std::endl;
					ProcessController::select_.AfterSelectSaDE_SuccessStrategiesSelection();
					ProcessController::select_.AfterSelectSaDE_FailStrategiesSelection();
					break;
				}
			}

			// parameter control CR
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kCR) == corematerial::ParametersControl::kCR) {
				// memory CR
				if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kMemoryCR) == corematerial::ParametersControlMethod::kMemoryCR) {
					switch (ProcessController::lv_.ev_.GetDE_MemoryCR_Pattern()) {
					case corematerial::MemoryCR::kSaDE_CR:
						if (ProcessController::debug_)
							std::cout << "AfterSelectSaDE_SuccessMemoryCR" << std::endl;
						ProcessController::select_.AfterSelectSaDE_SuccessMemoryCR();
						break;
					case corematerial::MemoryCR::kSHADE_CR:
						ProcessController::select_.AfterSelectSHADE_SuccessMemoryCR();
						break;
					case corematerial::MemoryCR::kPAEHN_SHADE_CR:
						paehn_.SuccessSelectionSHADE_SubpopMemoryCR(); //--------------------------------------------------------------------
						break;
					default:
						printf("Unknown memory CR pattern.");
						exit(0);
					}
				}
			}

			// parameter control F
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kF) == corematerial::ParametersControl::kF) {
				if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kMemoryF) == corematerial::ParametersControlMethod::kMemoryF) {
					switch (ProcessController::lv_.ev_.GetDE_MemoryF_Pattern()) {
					case corematerial::MemoryF::kSHADE_F:
						ProcessController::select_.AfterSelectSHADE_SuccessMemoryF();
						break;
					case corematerial::MemoryF::kPAEHN_SHADE_F:
						paehn_.SuccessSelectionSHADE_SubpopMemoryF();
						break;
					default:
						printf("Unknown memory F pattern.");
						exit(0);
					}
				}
			}

			// archive
			if ((ProcessController::lv_.ev_.GetDE_Memory_Pattern() & corematerial::Memory::kArchive) == corematerial::Memory::kArchive) {
				switch (ProcessController::lv_.ev_.GetDE_ArchivePattern()) {
				case(corematerial::Archive::kSHADE_Archive):
					ProcessController::select_.AfterSelectSHADE_SuccessArchive();
					break;
				case(corematerial::Archive::kPAEHN_Archive):
					//abandoned
					//ProcessController::select_.AfterSelectPAEHN_SuccessArchive();
					break;
				default:
					printf("Unknown archive pattern.");
					exit(0);
				}
			}

			if ((ProcessController::lv_.ev_.GetSubpopulationPattern() & corematerial::Subpopulation::kPAEHN_Subpopulation) == corematerial::Subpopulation::kPAEHN_Subpopulation) {
				//PAEHN create subpopulation
				bool enable_recurrsive_clearing = true, random_new_population = false, reinitializing_crowding_promising_archive = true; // default pattern
				if ((ProcessController::lv_.ev_.GetPAEHN_Pattern() & corematerial::PAEHN::kPAEHN_RecurrsiveClearing) != corematerial::PAEHN::kPAEHN_RecurrsiveClearing) {
					enable_recurrsive_clearing = false;//--------------------------------------------------------------------
				}
				if ((ProcessController::lv_.ev_.GetPAEHN_Pattern() & corematerial::PAEHN::kPAEHN_RandomNewPopulation) == corematerial::PAEHN::kPAEHN_RandomNewPopulation) {
					random_new_population = true;//--------------------------------------------------------------------
				}
				if ((ProcessController::lv_.ev_.GetPAEHN_Pattern() & corematerial::PAEHN::kPAEHN_ReinitializingCrowdingPromisingArchive) != corematerial::PAEHN::kPAEHN_ReinitializingCrowdingPromisingArchive) {
					reinitializing_crowding_promising_archive = false;
				}
				ProcessController::select_.AfterSelectPAEHN_FormSubpopulation(enable_recurrsive_clearing, random_new_population, reinitializing_crowding_promising_archive);
			}

			// 2022-08-13 This might be a little weird here.
			// Is it okay using the same corematerial to clear memories?
			corematerial::Selection s = ProcessController::lv_.ev_.GetDE_SelectionPattern();
			if (s == corematerial::Selection::kCompetition) {
				if (ProcessController::debug_)
					std::cout << "AfterSelectDE" << std::endl;
				// Clear all the individuals in next_generation_.
				// Clear all the individuals which is replaced by some offsprings.
				select_.AfterSelectDE();
			}
		}

		void ProcessController::NextGen() {
			ProcessController::de_.NextGen();
			ProcessController::sade_.NextGen();
			ProcessController::shade_.NextGen();
			ProcessController::paehn_.NextGen();
		}

		bool ProcessController::ReachStopCriterion() {
			const corematerial::StopCriterion sc = ProcessController::lv_.ev_.GetStopCriterionPattern();
			if ((sc & corematerial::StopCriterion::kMaxFE) == corematerial::StopCriterion::kMaxFE) {
				if (ProcessController::lv_.now_ffe_ >= ProcessController::lv_.ev_.GetDE_MaxFE()) {
					return true;
				}
			}
			if ((sc & corematerial::StopCriterion::kMaxGeneration) == corematerial::StopCriterion::kMaxGeneration) {
				if (ProcessController::lv_.now_gen_ >= ProcessController::lv_.ev_.GetDE_MaxGen()) {
					return true;
				}
			}
			if ((sc & corematerial::StopCriterion::kReachOptimal) == corematerial::StopCriterion::kReachOptimal) {
				if (ProcessController::lv_.population_.GetBestFitness() <= ProcessController::lv_.ev_.GetDE_Optimal_()) {
					return true;
				}
			}
			if (sc == corematerial::StopCriterion::UNSET) {
				std::cout << "Unset Stop criterion." << std::endl;
				exit(0);
			}
			return false;
		}

		void ProcessController::RecycleMemory() {
			// initial
			switch (ProcessController::lv_.ev_.GetDE_InitialPopulationPattern()) {
			case corematerial::InitialPopulation::kRandom:
				ProcessController::memory_.RecycleDE_Fitness();
				break;
			case corematerial::InitialPopulation::kLHD:
				break;
			}
			// strategies selection
			if (ProcessController::lv_.ev_.GetDE_NumOfStrategies() >= 2) {
				switch (ProcessController::lv_.ev_.GetDE_StrategiesSelectionPattern()) {
				case corematerial::StrategiesSelection::kSaDE_ss:
					ProcessController::memory_.RecycleSaDE_StrategiesSelection();
					break;
				}
			}
			// parameter control CR
			if ((ProcessController::lv_.ev_.GetParametersControlPattern() & corematerial::ParametersControl::kCR) == corematerial::ParametersControl::kCR) {
				if ((ProcessController::lv_.ev_.GetParametersControlMethodPattern() & corematerial::ParametersControlMethod::kMemoryCR) == corematerial::ParametersControlMethod::kMemoryCR) {
					// Memory CR
					switch (ProcessController::lv_.ev_.GetDE_MemoryCR_Pattern()) {
					case corematerial::MemoryCR::kSaDE_CR:
						if ((ProcessController::lv_.ev_.GetDE_MemoryCR_SettingPattern() & corematerial::MemoryCR_Setting::kMemoryByEachStrategy) == corematerial::MemoryCR_Setting::kMemoryByEachStrategy) {
							ProcessController::parameter_control_.DemallocSaDE_MemoryCR(ProcessController::lv_.ev_.GetDE_NumOfStrategies());
						}
						else {
							ProcessController::parameter_control_.DemallocSaDE_MemoryCR(1);
						}
						break;
					case corematerial::MemoryCR::kSHADE_CR:

						break;
					}
				}
			}

			// shade
			ProcessController::shade_.RecycleMemory();
			ProcessController::lv_.file_stream_.flush();
			ProcessController::lv_.file_stream_.close();
		}
	}
}