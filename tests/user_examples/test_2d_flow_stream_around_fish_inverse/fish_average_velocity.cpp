#include "fish_average_velocity.h"

namespace SPH
{
	//=========================================================================================================//
	namespace solid_dynamics
	{
		//=================================================================================================//
		FluidVelocityAroundSolid::FluidVelocityAroundSolid(BaseContactRelation& contact_relation, Real frequency)
			: LocalDynamics(contact_relation.getSPHBody()), FSIContactData(contact_relation), dynamicArray(), frequency_(frequency),
			compute_total_momentum_(contact_relation.getSPHBody(), "Velocity"),
			vel_(particles_->vel_), pos_(particles_->pos_),
			relative_velocity_(*particles_->getVariableByName<Vecd>("RelativeAverageVelocity")),
			solid_ave_vel_(*particles_->getVariableByName<Vecd>("SolidAverageVelocity")),
			ave_acc_(*particles_->getVariableByName<Vecd>("RelativeAverageAcceleration")),
			mass_(particles_->mass_)
		{
			for (size_t k = 0; k != contact_particles_.size(); ++k)
			{
				contact_vel_n_.push_back(&(contact_particles_[k]->vel_));
			}

			sum_ = Vecd::Zero();
			new_relative_velocity_ = Vecd::Zero();
			previous_relative_velocity_ = Vecd::Zero();
			acceleration_ = Vecd::Zero();

			ReduceDynamics<QuantitySummation<Real>> compute_total_mass_(contact_relation.getSPHBody(), "MassiveMeasure");
			total_mass_ = compute_total_mass_.exec();
		}
		//=================================================================================================//
		void FluidVelocityAroundSolid::setupDynamics(Real dt)
		{
			Vecd fluid_ave_vel = Vecd::Zero();
			Vecd solid_ave_vel = Vecd::Zero();
			std::set<size_t> index_set;

			/** select the fluid particle around the solid. */
			for (size_t i = 0; i < particles_->total_real_particles_; ++i)
			{
				for (size_t k = 0; k < contact_configuration_.size(); ++k)
				{
					Neighborhood& contact_neighborhood = (*contact_configuration_[k])[i];
					for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
					{
						size_t index_j = contact_neighborhood.j_[n];
						index_set.insert(index_j);
					}
				}
			}

			/** average fluid velocity around the solid. */
			for (size_t k = 0; k < contact_configuration_.size(); ++k)
			{
				StdLargeVec<Vecd>& vel_n_k = *(contact_vel_n_[k]);
				for (std::set<size_t>::iterator n = index_set.begin(); n != index_set.end(); ++n)
				{
					fluid_ave_vel += vel_n_k[*n] / index_set.size();
				}
			}

			/** average solid velocity. */
			solid_ave_vel = compute_total_momentum_.exec(dt) / total_mass_;

			/** relative velocity. */
			//Vecd rel_vel = -(solid_ave_vel - fluid_ave_vel);
			Vecd rel_vel = -solid_ave_vel;

			/** average within a contraction frequency */
			dynamicArray.emplace_back(rel_vel, GlobalStaticVariables::physical_time_);
			Real delta_t = dynamicArray.back().second - dynamicArray.front().second;

			if (delta_t >= (1 / frequency_))
			{
				sum_ = Vecd::Zero();
				for (size_t jj = 0; jj < dynamicArray.size(); jj++)
				{
					sum_ += dynamicArray[jj].first;
				}
				sum_ /= dynamicArray.size();

				dynamicArray.clear();
			}

			/** new relative_velocity_ */
			new_relative_velocity_ = Vec2d(sum_.norm(), 0.0);

			/** calculate acceleration_ */
			acceleration_ = (new_relative_velocity_ - previous_relative_velocity_) / dt;

			previous_relative_velocity_ = new_relative_velocity_;

			solid_ave_vel_[0] = solid_ave_vel;
		}
		//=================================================================================================//
		void FluidVelocityAroundSolid::update(size_t index_i, Real dt)
		{
			relative_velocity_[index_i] = new_relative_velocity_;
			solid_ave_vel_[index_i] = solid_ave_vel_[0];
			ave_acc_[index_i] = acceleration_;
		}
	}

//=================================================================================================//
} 
  //=================================================================================================//