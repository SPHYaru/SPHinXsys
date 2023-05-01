#include "composite_material.h"
#include "base_particles.hpp"

#include <numeric>

namespace SPH
{
	//=================================================================================================//
	void CompositeMaterial::initializeLocalParameters(BaseParticles* base_particles)
	{
		ElasticSolid::initializeLocalParameters(base_particles);
		for (int i = 0; i < CompositeMaterails_.size(); ++i)
			CompositeMaterails_[i]->initializeLocalParameters(base_particles);

		for (int j = 0; j < CompositeMaterails_.size(); ++j)
			sound_speed_.push_back(CompositeMaterails_[j]->ReferenceSoundSpeed());

		c0_ = *std::max_element(sound_speed_.begin(), sound_speed_.end());
		setContactStiffness(c0_);
		base_particles->registerVariable(materail_id_, "MaterailId");
	}

}
