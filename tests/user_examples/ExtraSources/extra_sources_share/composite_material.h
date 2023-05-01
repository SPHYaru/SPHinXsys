
#ifndef COMPOSITE_MATERIAL_H
#define COMPOSITE_MATERIAL_H

#include "elastic_solid.h"
#include <fstream>

namespace SPH
{
	/**@class CompositeMaterial*/
	class CompositeMaterial : public ElasticSolid
	{
	protected:
		UniquePtrKeepers<ElasticSolid> composite_ptr_keeper_;
		StdVec<ElasticSolid*> CompositeMaterails_;
		/** initialize the local properties, fiber and sheet direction. */
		std::vector<Real> sound_speed_;
	public:
		StdLargeVec<int> materail_id_;

		explicit CompositeMaterial(Real rho0) : ElasticSolid(rho0)
		{
			material_type_name_ = "CompositeMaterial";
		};
		virtual ~CompositeMaterial() {};

		virtual void initializeLocalParameters(BaseParticles* base_particles) override;

		virtual Matd StressPK2(Matd& deformation, size_t particle_index_i) override
		{
			return CompositeMaterails_[materail_id_[particle_index_i]]->StressPK2(deformation, particle_index_i);
		};

		Real CompositeDensity(size_t particle_index_i)
		{
			return CompositeMaterails_[materail_id_[particle_index_i]]->ReferenceDensity();
		};

		virtual Matd StressCauchy(Matd& almansi_strain, Matd& F, size_t particle_index_i) override
		{
			return Matd::Identity();
		};

		virtual Real VolumetricKirchhoff(Real J) override
		{
			return 0.0;
		};

		virtual std::string getRelevantStressMeasureName() override { return "PK2"; };

		template <class AddMaterial, typename... ConstructorArgs>
		void add(ConstructorArgs &&...args)
		{
			CompositeMaterails_.push_back(composite_ptr_keeper_.createPtr<AddMaterial>(std::forward<ConstructorArgs>(args)...));
		};
	};
}
#endif //COMPOSITE_MATERAIL_H