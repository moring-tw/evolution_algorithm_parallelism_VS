#include "../headers/cell.h"
namespace evolutionalgorithm {
	Cell::Cell() {
		is_block_ = false;
		individuals_counter_ = 0;
	}
	Cell::~Cell() {

	}

	void Cell::Initial(size_t population_size) {
		population_.Initial(population_size);
	}

	void Cell::AddIndividual(Individual i) {
		
	}

	bool Cell::DelIndividual(Individual i) {
		return false;
	}
}