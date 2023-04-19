#pragma once
#ifndef EA_HEADERS_DIFFERENTIAL_EVOLUTION_SAMPLE_H_
#define EA_HEADERS_DIFFERENTIAL_EVOLUTION_SAMPLE_H_
#include "behavior.h"

namespace evolutionalgorithm {
	namespace differentialevolution {
		/// <summary>
		///		checks some variables are valid before invoking corresponding algorithms.
		///		This class is about sampling population.
		/// </summary>
		class Sample : public Behavior {
		public:
			Sample();
			~Sample();
		protected:
		private:
		};
	}
}
#endif // EA_HEADERS_DIFFERENTIAL_EVOLUTION_SAMEPLE_H_
