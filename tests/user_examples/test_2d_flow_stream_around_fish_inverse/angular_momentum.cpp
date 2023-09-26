#include "angular_momentum.h"

namespace SPH
{
	namespace solid_dynamics
	{
		//=================================================================================================//
		SolidTorque::SolidTorque(SPHBody& sph_body, Real frequency, Real covariance)
			: SolidDataSimple(sph_body), LocalDynamics(sph_body), frequency_(frequency), covariance_(covariance),
			omega(0.0), omega_filter(0.0), theta(0.0), dynamicArray(), sum(0.0), avg_direction(0.0), theta_filter(0.0),
			vel_(particles_->vel_), pos_(particles_->pos_), pos0_(particles_->pos0_), mass_(particles_->mass_),
			theta_(*particles_->getVariableByName<Real>("Theta")),
			omega_(*particles_->getVariableByName<Real>("AngularVelocity")),
			theta_filter_(*particles_->getVariableByName<Real>("Theta_filter")) 
			{}
		//=================================================================================================//
		void SolidTorque::setupDynamics(Real dt)
		{
			Real rotational_inertia_(0.0);
			Vec3d angular_momentum_ = Vec3d::Zero();
			Vecd mass_center_ = Vecd::Zero();
			Vecd vel_center_ = Vecd::Zero();

			for (size_t i = 0; i != particles_->total_real_particles_; ++i)
			{
				mass_center_ += pos_[i] / particles_->total_real_particles_;
				vel_center_ += vel_[i] / particles_->total_real_particles_;
			}

			for (size_t j = 0; j != particles_->total_real_particles_; ++j)
			{
				Vecd rij = pos_[j] - mass_center_;
				Vecd vij = vel_[j] - vel_center_;
				rotational_inertia_ += mass_[j] * rij.dot(rij);
				angular_momentum_ += Vec3d(rij[0], rij[1], 0.0).cross(mass_[j] * Vec3d(vij[0], vij[1], 0.0));  //cross product of vectors
			}

			omega = angular_momentum_[2] / rotational_inertia_;

            dynamicArray.emplace_back(omega, GlobalStaticVariables::physical_time_);
			Real delta_t = dynamicArray.back().second - dynamicArray.front().second;

			if (delta_t >= (1 / frequency_))
			{
				avg_direction = 0.0;

				for (size_t jj = 0; jj < dynamicArray.size(); jj++)
				{
					avg_direction += dynamicArray[jj].first;
				}
				avg_direction /= dynamicArray.size();

				dynamicArray.clear();
			}
			omega_filter = avg_direction;
		}
		//=================================================================================================//
		void SolidTorque::update(size_t index_i, Real dt)
		{
			omega_[index_i] = omega;
			theta_[index_i] += omega * dt;
			theta_filter_[index_i] += omega_filter * dt;
		}
		//=================================================================================================//
	}
	//=================================================================================================//
}
//=================================================================================================//