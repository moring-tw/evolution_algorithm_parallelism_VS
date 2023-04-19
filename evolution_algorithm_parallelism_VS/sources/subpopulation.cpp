#include "../headers/subpopulation.h"

namespace evolutionalgorithm {
	namespace differentialevolution {
		namespace paehn {
			Subpopulation::~Subpopulation() {}
			Subpopulation::Subpopulation()
			{
			}

			void Subpopulation::Initial(size_t np, size_t dim) {
				Population::Initial(np, dim);
				children_.Initial(np, dim);
				while (Subpopulation::random_ri_.size() < np)
					Subpopulation::random_ri_.push_back(0);
				for (size_t np_counter = 0; np_counter < np; np_counter++) {
					population_[np_counter].orientation_.reserve(dim);
					for (size_t dim_counter = 0; dim_counter < dim; ++dim_counter) {
						population_[np_counter].orientation_.push_back(0.0);
					}
				}
			}

			void Subpopulation::Initial_Memory_CR(const size_t H, const double default_value) {
				Subpopulation::maxH_ = H;
				for (size_t H_counter = 0; H_counter < H; ++H_counter) {
					if (Subpopulation::memory_CR_.size() <= H_counter)
						Subpopulation::memory_CR_.push_back(default_value);
					else
						Subpopulation::memory_CR_[H_counter] = default_value;
				}
			}

			void Subpopulation::Initial_Memory_F(const size_t H, const double default_value) {
				Subpopulation::maxH_ = H;
				for (size_t H_counter = 0; H_counter < H; ++H_counter) {
					if (Subpopulation::memory_F_.size() <= H_counter)
						Subpopulation::memory_F_.push_back(default_value);
					else
						Subpopulation::memory_F_[H_counter] = default_value;
				}
			}


		}
	}
}