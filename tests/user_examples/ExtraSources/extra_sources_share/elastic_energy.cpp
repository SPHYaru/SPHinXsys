#include "elastic_energy.h"

#include <numeric>

namespace SPH
{
	namespace solid_dynamics
	{
	//=================================================================================================//
	ElasticEnergy::ElasticEnergy(SPHBody& sph_body)
		: LocalDynamicsReduce<Real, ReduceSum<Real>>(sph_body, Real(0)),
		ElasticSolidDataSimple(sph_body), Vol_(particles_->Vol_), F_(particles_->F_),
		stress_PK1_B_(*particles_->getVariableByName<Matd>("CorrectedStressPK1"))
	{
		quantity_name_ = "ElasticEnergy";
	}
	//=================================================================================================//
	Real ElasticEnergy::reduce(size_t index_i, Real dt)
	{
		Matd strain = 0.5 * (F_[index_i].transpose() + F_[index_i]) - Matd::Identity();

		return 0.5 * CalculateDoubleDotProduct(strain, stress_PK1_B_[index_i]) * Vol_[index_i];
	}
	//=================================================================================================//
	SolidKinecticEnergy::SolidKinecticEnergy(SPHBody& sph_body)
		: LocalDynamicsReduce<Real, ReduceSum<Real>>(sph_body, Real(0)),
		ElasticSolidDataSimple(sph_body), vel_(particles_->vel_), mass_(particles_->mass_)
	{
		quantity_name_ = "SolidKinecticEnergy";
	}
	//=================================================================================================//
	Real SolidKinecticEnergy::reduce(size_t index_i, Real dt)
	{

		return 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
	}

	}
}
