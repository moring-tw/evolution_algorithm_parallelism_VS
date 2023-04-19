#pragma once
#ifndef EA_HEADERS_CELL_H_
#define EA_HEADERS_CELL_H_

#include <map>
#include <string>

#include "differential_evolution/localvariables.h"
#include "individual.h"
#include "population.h"

namespace evolutionalgorithm {
	class Cell {
	public:
		Cell();
		~Cell();
		void AddIndividual(Individual);
		bool DelIndividual(Individual);
		void Initial(size_t);
	private:
		Population population_;
		int individuals_counter_;
		bool is_block_;
		std::vector<double> fitness_;
	};
} // namespace evolutionalgorithm
#endif // EA_HEADERS_CELL_H_
