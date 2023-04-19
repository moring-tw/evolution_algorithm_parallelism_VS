#pragma once
#ifndef EA_HEADERS_CELLS_H_
#define EA_HEADERS_CELLS_H_

#include <unordered_map>
#include <string>

#include "cell.h"
#include "individual.h"

namespace evolutionalgorithm {
	class Cells {
	public:
		Cells();
		~Cells();
		/*bool IsCellExist(std::string);
		Cell GetCell(std::string);
		void AddNewCell(std::string);*/
	private:
		//Cell c;
		std::unordered_map<std::string, Individual> cells_;
	};
} // evolutionalgorithm
#endif // EA_HEADERS_CELLS_H_
