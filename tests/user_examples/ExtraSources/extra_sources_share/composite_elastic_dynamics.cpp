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
			Integration1stHalf::initialization(index_i,dt);
			rho_[index_i] = composite_material->CompositeDensity(index_i) / F_[index_i].determinant();
		}
		//=================================================================================================//
		void CompositeIntegration1stHalf::interaction(size_t index_i, Real dt)
		{
			CompositeMaterial* composite_material = dynamic_cast<CompositeMaterial*>(&particles_->elastic_solid_);
			
			Vecd acceleration = Vecd::Zero();
			const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd e_ij = inner_neighborhood.e_ij_[n];
				Real r_ij = inner_neighborhood.r_ij_[n];
				Real dim_r_ij_1 = Dimensions / r_ij;
				Vecd pos_jump = pos_[index_i] - pos_[index_j];
				Vecd vel_jump = vel_[index_i] - vel_[index_j];
				Real strain_rate = dim_r_ij_1 * dim_r_ij_1 * pos_jump.dot(vel_jump);
				Real weight = inner_neighborhood.W_ij_[n] * inv_W0_;
				Matd numerical_stress_ij =
					0.5 * (F_[index_i] + F_[index_j]) * elastic_solid_.PairNumericalDamping(strain_rate, smoothing_length_);
				acceleration += 1/ composite_material->CompositeDensity(index_i) * inner_neighborhood.dW_ijV_j_[n] *
					(stress_PK1_B_[index_i] + stress_PK1_B_[index_j] +
						numerical_dissipation_factor_ * weight * numerical_stress_ij) * e_ij;
			}

			acc_[index_i] = acceleration;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//
