#include "../headers/individual.h"

namespace evolutionalgorithm {
	Individual::Individual() {
		//Individual::dimensions_ = nullptr;
	}

	Individual::Individual(size_t dim) {
		Individual::dim_ = dim;
		Individual::dimensions_.reserve(dim);
		for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter)
			dimensions_.push_back(0.0);
		//Individual::dimensions_ = (double*)malloc(sizeof(double) * dim);
		//std::fill(dimensions_, dimensions_ + dim, 0);
	}

	Individual::~Individual() {
		/*if (Individual::dimensions_ != nullptr) {
			free(Individual::dimensions_);
			Individual::dimensions_ = nullptr;
		}*/
	}


	double* Individual::GetDimensions() {
		return &Individual::dimensions_[0];
	}

	void Individual::InitialComparisonResult() {
		Individual::is_equal_ = false;
		Individual::is_less_ = false;
		Individual::is_larger_ = false;
	}

	void Individual::SetCR(const double CR) {
		Individual::use_CR_ = CR;
	}

	void Individual::SetF(const double F) {
		Individual::use_F_ = F;
	}

	void Individual::SetK(const double K) {
		Individual::use_K_ = K;
	}

	void Individual::SetP(const double P) {
		Individual::use_P_ = P;
	}
	
	void Individual::SetComparedParent(const size_t compared_parent) {
		Individual::compared_parent_ = compared_parent;
	}

	void Individual::SetCoordinate(const size_t dim_index, const size_t coordinate) {
		Individual::coordinate_[dim_index] = coordinate;
	}

	void Individual::SetFitnessDifference(double difference) {
		Individual::fitness_difference_ = difference;
	}

	void Individual::CopyDimension(double* dimension) {
		for (size_t dim_counter = 0; dim_counter < dim_; ++dim_counter) {
			dimensions_[dim_counter] = dimension[dim_counter];
		}
	}

	double Individual::GetCR() {
		return Individual::use_CR_;
	}

	double Individual::GetF() {
		return Individual::use_F_;
	}

	double Individual::GetK() {
		return Individual::use_K_;
	}

	double Individual::GetP() {
		return Individual::use_P_;
	}

	size_t Individual::GetComparedParent() {
		return Individual::compared_parent_;
	}

	size_t Individual::GetDim()
	{
		return Individual::dim_;
	}

	double Individual::GetFitnessDifference() {
		return Individual::fitness_difference_;
	}

	void Individual::ShowPath() {
		auto path = Individual::record_path_.begin();
		while (path != Individual::record_path_.end()) {
			auto interval = (*path).begin();
			while (interval != (*path).end()) {
				std::cout << (*interval) << ">";
				++interval;
			}
			std::cout << std::endl;
			++path;
		}
	}

	std::vector<double>::iterator Individual::begin() {
		return Individual::dimensions_.begin();
	}

	std::vector<double>::iterator Individual::end() {
		return Individual::dimensions_.end();
	}
} // namespace evolutionalgorithm