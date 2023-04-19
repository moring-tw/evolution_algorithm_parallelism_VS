#pragma once
#include <vector>

#include "population.h"
namespace evolutionalgorithm {
	namespace differentialevolution {
		namespace paehn {
			class Subpopulation : public Population {
			public:
				Subpopulation();
				~Subpopulation();
				Population children_;
				std::vector<Individual> archive_;
				std::vector<double> fitness_;
				std::vector<double> memory_CR_;
				std::vector<double> memory_F_;
				std::vector<double> copied_fitness_;
				std::vector<double> success_CR_;
				std::vector<double> success_F_;
				std::vector<double> delta_fitness_;
				
				size_t maxH_ = 0, insert_H_ = 0;
				bool insert_H_plus_1_ = false;
				std::vector<size_t> random_ri_;

				double magnification = 0.0;
				double furthest_distance = 0.0;
				double min_internal_distance = DBL_MAX;
				double max_internal_distance = 0.0;
				double extended_distance = 0.0;
				std::vector<double> initial_centroid;
				std::vector<double> centroid;
				size_t stagnation_counter = 0;
				bool is_new_ = true;
				bool is_stagnation = false;
				bool is_merged = false;
				bool from_buffer = false;
				size_t buffer_index = SIZE_MAX;
				bool keep = false;

				void Initial(size_t, size_t);
				void Initial_Memory_CR(const size_t H, const double default_value);
				void Initial_Memory_F(const size_t H, const double default_value);

			protected:
			private:
			};
		}
	}
}