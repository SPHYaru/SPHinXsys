#include "composite_elastic_dynamics.h"

#include <numeric>

namespace SPH
{
	//=========================================================================================================//
	namespace solid_dynamics
	{
		//=================================================================================================//
		CompositeIntegration1stHalf::
			CompositeIntegration1stHalf(BaseInnerRelation &inner_relation)
			: Integration1stHalf(inner_relation){}
		//=================================================================================================//
		void CompositeIntegration1stHalf::initialization(size_t index_i, Real dt)
		{
			CompositeMaterial* composite_material = dynamic_cast<CompositeMaterial*>(&particles_->elastic_solid_);
			if (composite_material == nullptr) {
				throw std::runtime_error("Particles_ is not a compositeMaterial object.");
			}
			rho0_ = composite_material->CompositeDensity(index_i);
			inv_rho0_ = 1.0 / rho0_;
			Integration1stHalf::initialization(index_i,dt);
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//
