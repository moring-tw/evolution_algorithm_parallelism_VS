#pragma once
#ifndef EA_HEADERS_DIFFERENTIAL_EVOLUTION_BEHAVIOR_H_
#define EA_HEADERS_DIFFERENTIAL_EVOLUTION_BEHAVIOR_H_

#include "variants/de.h"
#include "variants/sade.h"
#include "variants/shade.h"
#include "variants/random.h"
#include "variants/paehn.h"
#include "localvariables.h"

namespace evolutionalgorithm {
	namespace differentialevolution {
		/// <summary>
		///		For those classes who extend from this class should be considered as a transfer station/hub,
		///		and always be invoked by the class, process_controller.
		///		They should check some relative variables are valid to keep the program from crash. 
		///		Then invoking corresponding algorithm. 
		/// </summary>
		class Behavior {
		public:
			Behavior();
			virtual ~Behavior();
			void SetDE(DE*);
			void SetSaDE(SaDE*);
			void SetSHADE(SHADE*);
			void SetRandom(Random*);
			void SetPAEHN(PAEHN*);
			void SetLocalVariables(LocalVariables*);
			void CheckLV();
		protected:
			DE* de_ = nullptr;
			SaDE* sade_ = nullptr;
			SHADE* shade_ = nullptr;
			Random* random_ = nullptr;
			PAEHN* paehn_ = nullptr;
			LocalVariables* lv_ = nullptr;
		};
	}
}
#endif // !EA_HEADERS_DIFFERENTIAL_EVOLUTION_BEHAVIOR_H_
