#include "../headers/publicfunction.h"

namespace evolutionalgorithm {
	std::mt19937 PublicFunction::generator_(static_cast<int>(time(0)));
	PublicFunction::PublicFunction() {

	}

	PublicFunction::~PublicFunction() {

	}

	double PublicFunction::RandomDouble(double bound_min, double bound_max) {
		std::uniform_real_distribution<double> distribution(bound_min, bound_max);
		return distribution(PublicFunction::generator_);
	}

	int PublicFunction::RandomInt(int bound_min, int bound_max) { // possible return value includes bound_max
		std::uniform_int_distribution<int> distribution(bound_min, bound_max);
		return distribution(PublicFunction::generator_);
	}

	size_t* PublicFunction::ShufflePopulation(size_t np) {
		if (np <= 0) {
			std::cout << "Can't shuffle population which size is less than or equal to 0." << std::endl;
			exit(0);
		}
		size_t* pool = (size_t*)malloc(sizeof(size_t) * np);
		if (pool == nullptr) {
			std::cout << "Can't malloc pool, this may cause by np equal to ZERO." << std::endl;
			exit(0);
		}
		for (size_t ini = 0; ini < np; ++ini) {
			pool[ini] = ini;
		}
		for (size_t index = 0; index < np; ++index) {
			std::swap(pool[index], pool[index + PublicFunction::RandomInt(0, np - index - 1)]);
		}
		return pool;
	}

	unsigned* PublicFunction::ShufflePopulation(size_t np, size_t base_individual, unsigned how_much) {
		if (how_much > np - 1) {
			std::cout << "Can't shuffle more than np-1" << std::endl;
			exit(0);
		}
		if (base_individual >= np) {
			std::cout << "base_individual out of population range." << std::endl;
			exit(0);
		}
		if (np <= 1) {
			std::cout << "Can't shuffle population when there np is less than 2." << std::endl;
			exit(0);
		}
		unsigned* pool = (unsigned*)malloc(sizeof(unsigned) * np - 1);
		if (pool == nullptr) {
			std::cout << "Can't malloc pool, this may cause by np less than TWO." << std::endl;
			exit(0);
		}
		bool skiped = false;
		for (size_t ini = 0; ini < (np - 1) ; ++ini) {
			if (ini == base_individual) skiped = true;
			if (skiped)
				pool[ini] = ini + 1;
			else
				pool[ini] = ini;
			if (np == 2) break;
		}
		if (np > 2) {
			for (unsigned index = 0; index < how_much; ++index) {
				std::swap(pool[index], pool[index + PublicFunction::RandomInt(0, np - index - 2)]);
			}
		}
		return pool;
	}

	//
	// Quick Sort, time complexity O(n*logn)
	//
	
	/// <summary>
	///		If the population hasn't been sorted, use this function to sort the population.
	///		This function won't change the order of fitness array.
	///		Implementing with Quick sort, time complexity is O(n*logn)
	///		Sort from smaller to larger
	/// </summary>
	/// <param name="fitness"> fitness array of population.</param>
	/// <param name="np_rank"> array that stores population index after sorting.</param>
	/// <param name="np"> population size</param>
	/// <param name="left"></param>
	void PublicFunction::SortPopulation(double* fitness, size_t* np_rank, size_t np, size_t left) {
		if (left < np) {
			size_t i = left, j = np;
			while (true) {
				while (i + 1 < np && fitness[np_rank[++i]] < fitness[np_rank[left]]);
				while (j > 0 && fitness[np_rank[--j]] > fitness[np_rank[left]]);
				if (i >= j) break;
				std::swap(np_rank[i], np_rank[j]);
			}
			std::swap(np_rank[left], np_rank[j]);
			SortPopulation(fitness, np_rank, j, left);
			SortPopulation(fitness, np_rank, np, j + 1);
		}
	}

	/// <summary>
	///		Same as the previous Sortpopulation.
	///		Just change the double array type to "vector&lt;double&gt;."
	/// </summary>
	void PublicFunction::SortPopulation(const std::vector<double>& fitness, std::vector<size_t>& np_rank, size_t np, size_t left) {
		if (left < np) {
			size_t i = left, j = np;
			while (true) {
				while (i + 1 < np && fitness[np_rank[++i]] < fitness[np_rank[left]]);
				while (j > 0 && fitness[np_rank[--j]] > fitness[np_rank[left]]);
				if (i >= j) break;
				std::swap(np_rank[i], np_rank[j]);
			}
			std::swap(np_rank[left], np_rank[j]);
			SortPopulation(fitness, np_rank, j, left);
			SortPopulation(fitness, np_rank, np, j + 1);
		}
	}

	/// <summary>
	///		If the population has been sorted and only one individual is changed,
	///		invoking this function to maintain sorting correctly.
	///		time complexity O(n). 
	///		Sort from smaller to larger
	/// </summary>
	/// <param name="fitness"> fitness array of population.</param>
	/// <param name="np_rank"> array that stores population index after sorting.</param>
	/// <param name="np"> population size</param>
	/// <param name="base_individual"> the individual index that needs to be sorted.</param>
	void PublicFunction::SortIndividual(double* fitness, unsigned* np_rank, size_t np, size_t base_individual) {
		unsigned rank_position = 0;
		for (size_t np_counter = 0; np_counter < np; ++np_counter) {
			if (np_rank[np_counter] == base_individual) {
				rank_position = np_counter;
				break;
			}
		}
		// If fitness of individual is getting better, move its rank to upper place.
		while (rank_position > 0) {
			if (fitness[np_rank[rank_position - 1]] > fitness[np_rank[rank_position]]) {//ASC
				std::swap(np_rank[rank_position - 1], np_rank[rank_position]);
				--rank_position;
			} else {
				break;
			}
		}
		// If fitness of individual is getting worse, move its rank to lower place.
		while (rank_position < np - 1) {
			if (fitness[np_rank[rank_position + 1]] < fitness[np_rank[rank_position]]) {
				std::swap(np_rank[rank_position + 1], np_rank[rank_position]);
				++rank_position;
			} else {
				break;
			}
		}
	}

	double PublicFunction::NormalDistributionRandom(double mean, double stddev) {
		std::normal_distribution<double> distribution(mean, stddev);
		return distribution(generator_);
	}

	double PublicFunction::CauchyDistributionRandom(double mean, double stddev) {
		std::cauchy_distribution<double> distribution(mean, stddev);
		return distribution(generator_);
	}

	double PublicFunction::WeightedArithmeticMean(std::vector<double>& delta_fitness, std::vector<double>& success) {
		double amount_delta = 0, result = 0;
		for (unsigned delta_fitness_counter = 0; delta_fitness_counter < delta_fitness.size(); ++delta_fitness_counter) {
			amount_delta += delta_fitness[delta_fitness_counter];
		}
		for (unsigned success_counter = 0; success_counter < success.size(); ++success_counter) {
			result += (delta_fitness[success_counter] / amount_delta) * success[success_counter];
		}
		return result;
	}

	double PublicFunction::WeightedLenmerMean(std::vector<double>& delta_fitness, std::vector<double>& success) {
		double amount_delta = 0, result = 0, temp = 0;
		for (unsigned delta_fitness_counter = 0; delta_fitness_counter < delta_fitness.size(); ++delta_fitness_counter) {
			amount_delta += delta_fitness[delta_fitness_counter];
		}
		for (unsigned success_counter = 0; success_counter < success.size(); ++success_counter) {
			result += (delta_fitness[success_counter] / amount_delta) * pow(success[success_counter], 2);
			temp += (delta_fitness[success_counter] / amount_delta) * success[success_counter];
		}
		result /= temp;
		return result;
	}

	double PublicFunction::CalculateStandardDeviation(double* values, double mean, unsigned values_size) {
		double stddev = 0;
		for (unsigned size_counter = 0; size_counter < values_size; ++size_counter) {
			stddev += pow(values[size_counter] - mean, 2);
		}
		stddev /= values_size;
		stddev = pow(stddev, 0.5);
		return stddev;
	}

	double PublicFunction::RoundToDecimalPlace(double value, unsigned decimal_place) {
		value = floor(value * pow(10, decimal_place) + 0.5) / pow(10, decimal_place);
		return value;
	}

	double PublicFunction::CalculateEuclideanDistance(const std::vector<double>& individual_a, const std::vector<double>& individual_b) {
		double distance = 0;
		if (individual_a.size() != individual_b.size()) {
			throw std::invalid_argument("It's required two vectors have same sizes when calculation Euclidean Distance.");
		} else {
			for (size_t index = 0; index < individual_a.size(); ++index) {
				distance += pow(individual_a[index] - individual_b[index], 2);
			}
			return sqrt(distance);
		}
	}

	double PublicFunction::CalculateEuclideanDistance(double* individual_a, double* individual_b, size_t length) {
		double distance = 0;
		for (size_t index = 0; index < length; ++index) {
			distance += pow(individual_a[index] - individual_b[index], 2);
		}
		return sqrt(distance);
	}

	double PublicFunction::CalculateEuclideanDistance(double* individual_a, const std::vector<double>& individual_b, size_t length) {
		double distance = 0;
		for (size_t index = 0; index < length; ++index) {
			distance += pow(individual_a[index] - individual_b[index], 2);
		}
		return sqrt(distance);
	}

	double PublicFunction::CalculateManhattanDistance(const std::vector<double>& individual_a, const std::vector<double>& individual_b) {
		double distance = 0;
		if (individual_a.size() != individual_b.size()) {
			throw std::invalid_argument("It's required two vectors have same sizes when calculating Manhattan Distance.");
		} else {
			for (size_t index = 0; index < individual_a.size(); ++index) {
				distance += fabs(individual_a[index] - individual_b[index]);
			}
			return distance;
		}
	}

	double PublicFunction::CalculateManhattanDistance(double* individual_a, double* individual_b, size_t length) {
		double distance = 0;
		for (size_t index = 0; index < length; ++index) {
			distance += fabs(individual_a[index] - individual_b[index]);
		}
		return distance;
	}

	double PublicFunction::FindQuartile(std::vector<double>& v, size_t Q) {
		if (Q > 3 || Q < 1) {
			throw std::invalid_argument("invalid second argument. It should be 1, 2 or 3");
		}
		if (v.size() == 0) return 0;
		else if (v.size() == 1) return v[0];
		else if (v.size() == 2)
			if (Q == 1) return v[0];
			else if (Q == 2) return (v[0] + v[1]) / 2;
			else if (Q == 3) return (v[1]);
		std::sort(v.begin(), v.end());
		std::vector<double> new_v(v.begin(), v.begin() + 2);
		for (double value : new_v) {
			//std::cout << value << std::endl;
		}
		if (Q == 2) return PublicFunction::FindMedian(v);
		else if(Q==1)
			if (v.size() % 2 == 0) {
				return PublicFunction::FindMedian(std::vector<double>(v.begin(), v.begin() + v.size() / 2));
			}
			else {
				return PublicFunction::FindMedian(std::vector<double>(v.begin(), v.begin() + (v.size() - 1) / 2));
			}
		else {
			if (v.size() % 2 == 0) {
				return PublicFunction::FindMedian(std::vector<double>(v.begin() + v.size() / 2, v.end()));
			}
			else {
				return PublicFunction::FindMedian(std::vector<double>(v.begin() + (v.size() + 1) / 2, v.end()));
			}
		}
	}

	double PublicFunction::FindMedian(std::vector<double> v) {
		double res;
		std::sort(v.begin(), v.end());
		if (v.size() % 2 == 0) {
			return res = (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2;
		} else{
			return res = v[(v.size() - 1) / 2];
		}
	}
} // namespace evolutionalgorithm