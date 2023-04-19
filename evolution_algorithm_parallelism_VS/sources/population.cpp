#include"../headers/population.h"

namespace evolutionalgorithm {
	Population::Population() {
	}

	//Population::Population(Population&& old_pop) noexcept
	//	: furthest_distance(std::move(old_pop.furthest_distance))
	//	, centroid(std::move(old_pop.centroid))
	//{
	//	population_ = old_pop.population_;
	//	old_pop.population_ = nullptr;
	//}

	Population::~Population() {
		//
		//for each allocation, there has to be a deallocation 
		//what was allocated using new should be freed using delete, using new[] should be deallocated using delete[] and using malloc should be deallocated using free
		//
		//if (population_ != nullptr) {
			//for (unsigned i = 0; i < Population::np_; ++i) {
				//delete& (*population_[i]);
				//population_[i] = nullptr;
			//}
			//delete[] population_;
			//population_ = nullptr;
		//}

		//if (population_.size() != 0) {
		//	for (size_t i = 0; i < population_.size(); ++i) {
		//		//delete& (*population_[i]);
		//		//population_[i] = nullptr;
		//	}
		//}
	}

	void Population::Initial(size_t np, size_t dim) {
		if (Population::is_initialed_) {
			for (size_t i = 0; i < Population::np_; ++i) {
				//
				// this loop should be executed before new np being assigned to Population::np_
				// because the new np may different from old np. this will cause memory leak.
				//
				//delete population_[i];
				//population_[i] = nullptr;
			}
			//delete[] population_;
			//population_ = nullptr;
		}
		Population::size = 0;
		Population::np_ = np;
		Population::population_.clear();
		Population::population_.reserve(np);
		for (size_t np_counter = 0; np_counter < np; ++np_counter) {
			population_.emplace_back(dim);
		}
		Population::InitialNpRank();
		if (!Population::is_initialed_)
			Population::is_initialed_ = true;
	}

	//
	//this function won't malloc individual's dimension
	//using this function should assign "double*" to each individual as dimension
	//or will appear undefined behavior cause program crash.
	//
	void Population::Initial(size_t np) {
		//if (Population::is_initialed_) {
		//	for (size_t i = 0; i < Population::np_; ++i) {
		//		//delete population_[i];
		//		//population_[i] = nullptr;
		//	}
		//	//delete[] population_;
		//	//population_ = nullptr;
		//}
		Population::size = 0;
		Population::np_ = np;
		Population::population_.clear();
		Population::population_.reserve(np);
		for (size_t np_counter = 0; np_counter < np; ++np_counter) {
			population_.emplace_back();
		}
		Population::InitialNpRank();
		if (!Population::is_initialed_)
			Population::is_initialed_ = true;
	}

	void Population::SetBestFitness(double best) {
		Population::best_fitness_ = best;
	}

	double Population::GetBestFitness() {
		return Population::best_fitness_;
	}

	std::vector<Individual>::iterator Population::begin() {
		return Population::population_.begin();
	}

	std::vector<Individual>::iterator Population::end() {
		return Population::population_.end();
	}

	Individual Population::emplace_back() {
		Population::population_.emplace_back();
		return Population::population_.back();
	}

	void Population::clear() {
		Population::population_.clear();
	}

	size_t Population::GetNP() {
		return Population::population_.size();
		return Population::np_;
	}

	size_t Population::GetNowSize() {
		return Population::size;
	}

	void Population::AddNowSize()
	{
		++Population::size;
	}

	void Population::InitialNpRank() {
		Population::np_rank_.clear();
		Population::np_rank_.reserve(np_);
		for (size_t np_counter = 0; np_counter < np_; ++np_counter) {
			np_rank_.push_back(np_counter);
		}
	}
}