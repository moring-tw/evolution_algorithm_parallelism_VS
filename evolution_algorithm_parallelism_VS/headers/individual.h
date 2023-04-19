#pragma once
#ifndef EA_HEADERS_INDIVIDUAL_H_
#define EA_HEADERS_INDIVIDUAL_H_

#include <vector>
#include <iostream>

#include "corematerial.h"

namespace evolutionalgorithm {
	class Individual {
	public:
		Individual();
		Individual(size_t);
		virtual ~Individual();
		double& operator[](int d_index) {
			return dimensions_[d_index];
		}


		// beware the length of pointer arrays are not equal.
		// this function won't change any pointer inside the vector<double> variable, dimensions_.
		void operator=(double* dimensions) {
			if (dimensions != nullptr)
				for (size_t dim_counter = 0; dim_counter < dim_; ++dim_counter)
					dimensions_[dim_counter] = dimensions[dim_counter];
			else
				for (size_t dim_counter = 0; dim_counter < dim_; ++dim_counter)
					dimensions_[dim_counter] = 0.0;
		}

		void operator=(int value) {
			for (size_t dim_counter = 0; dim_counter < dim_; ++dim_counter)
				dimensions_[dim_counter] = value;
		}

		void operator=(double value) {
			for (size_t dim_counter = 0; dim_counter < dim_; ++dim_counter)
				dimensions_[dim_counter] = value;
		}

		double* GetDimensions();
		void InitialComparisonResult();
		void SetCR(const double);
		void SetF(const double);
		void SetK(const double);
		void SetP(const double);
		void SetComparedParent(const size_t);
		void SetCoordinate(const size_t, const size_t);
		void SetFitnessDifference(double);
		void CopyDimension(double*);
		double GetCR();
		double GetF();
		double GetK();
		double GetP();
		size_t GetComparedParent();
		size_t GetDim();
		double GetFitnessDifference();
		void ShowPath();

		std::vector<double>::iterator begin();
		std::vector<double>::iterator end();

		corematerial::MutantAndCrossover use_strategy_ = corematerial::MutantAndCrossover::UNSET;
		size_t use_strategy_index_ = 0;
		int number_of_unimproved_times_ = 0;
		std::vector<std::vector<unsigned>> record_path_;
		std::vector<double> orientation_;
		bool is_out_ = false;
		bool is_equal_ = false;
		bool is_less_ = false;
		bool is_larger_ = false;
		size_t stagnation_counter = 0;
		bool is_choosen_ = false;
		bool is_choosen2_ = false;
		bool form_new_pop_ = false;
		double form_new_pop_fitness_ = 0.0;
		int archive_index = -1;
		double fitness_ = DBL_MAX;
		bool high_priority = false;

		
	private:
		unsigned dim_ = 0;
		std::vector<double> dimensions_;
		//double* dimensions_ = 0;
		double use_CR_ = 0, use_F_ = 0, use_K_ = 0, use_P_ = 0;
		double fitness_difference_ = 0;
		size_t compared_parent_ = 0; // Records which one parent is compared to this individual.
		std::vector<size_t> coordinate_;
	};
}

#endif // !EA_HEADERS_INDIVIDUAL_H_
