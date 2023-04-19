#include "../headers/cells.h"
namespace evolutionalgorithm {
	Cells::Cells() {

	}

	Cells::~Cells() {
		cells_.clear();
	}

	/*void Cells::AddNewCell(std::string name) {
		if (IsCellExist(name)) {
			Cell c;
			cells_[name] = c;
		}
	}

	bool Cells::IsCellExist(std::string name) {
		return cells_.find(name) != cells_.end();
	}

	Cell Cells::GetPCell(std::string name) {
		if (!IsCellExist(name)) {
			AddNewCell(name);
		}
		return cells_[name];
	}*/
} // namespace evolutionalgorithm