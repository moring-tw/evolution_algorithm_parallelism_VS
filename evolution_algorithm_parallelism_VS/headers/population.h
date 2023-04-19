#pragma once
#ifndef EA_HEADERS_POPULATION_H_
#define EA_HEADERS_POPULATION_H_

#include <stdlib.h>

#include <vector>

#include "individual.h"
#include "corematerial.h"

namespace evolutionalgorithm {
	//extern unsigned int ExecutiveVariables::de_initial_pattern_;
	class Population
	{
	public:
		Population();
		//Population(Population&& old_pop) noexcept;
		virtual ~Population();
		Individual& operator[](unsigned i_index) {
			return population_[i_index];
			//return (*population_[i_index]);
		}

		friend bool operator!=(const Population& p, const Individual** ptr) {
			if (ptr == nullptr)
				if (p.population_.size() != 0) {
					return true;
				}
				else return false;
			return false;
		}

		/*friend bool operator!=(const Population& p, const Individual **ptr) {
			return !(p.population_ == ptr);
		}

		friend bool operator==(const Population& p, const Individual** ptr) {
			return (p.population_ == ptr);
		}*/

		void Initial(size_t);
		void Initial(size_t, size_t);
		void InitialNpRank();
		void SetBestFitness(double);
		double GetBestFitness();
		std::vector<Individual>::iterator begin();
		std::vector<Individual>::iterator end();
		Individual emplace_back();
		void clear();
		size_t GetNP();
		size_t GetNowSize();
		void AddNowSize();

		size_t best_individual = 0;
		std::vector<size_t> np_rank_;
	protected:
		std::vector<Individual> population_;
		//Individual** population_ = nullptr;
	private:
		bool is_initialed_ = false;
		unsigned np_ = 0;
		double best_fitness_ = DBL_MAX;
		size_t size = 0;
	};
}
#endif // !EA_HEADERS_POPULATION_H_
