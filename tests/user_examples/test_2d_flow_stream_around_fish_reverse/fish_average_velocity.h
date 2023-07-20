
#ifndef FISH_AVERAGE_VELOCITY_H
#define FISH_AVERAGE_VELOCITY_H

#include "fluid_dynamics_complex.h"
#include "elastic_dynamics.h"

namespace SPH
{
	namespace solid_dynamics
	{

		//----------------------------------------------------------------------
		//		for general solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<SolidParticles> SolidDataSimple;
		typedef DataDelegateContact<SolidParticles, BaseParticles> FSIContactData;

		/**
			 * @class ViscousForceFromFluid
			 * @brief Computing the viscous force from the fluid
			 */
		class FluidVelocityAroundSolid : public LocalDynamics, public FSIContactData
		{
		private:
			Real total_mass_;
			Vecd sum_;
			std::vector<std::pair<Vecd, Real>> dynamicArray;
			Real frequency_;
			ReduceDynamics<QuantityMoment<Vecd>> compute_total_momentum_;
			Vecd new_relative_velocity_, previous_relative_velocity_, acceleration_;

		public:
			explicit FluidVelocityAroundSolid(BaseContactRelation& contact_relation, Real frequency);
			virtual ~FluidVelocityAroundSolid() {};


			Vecd getNewRelativeVelocity() const {return new_relative_velocity_;}

			Vecd getAcceleration() const {return acceleration_;}

			virtual void setupDynamics(Real dt) override;

			void update(size_t index_i, Real dt = 0.0);

		protected:
			StdLargeVec<Vecd>& vel_, & pos_;
			StdVec<StdLargeVec<Vecd>*> contact_vel_n_;
			StdLargeVec<Vecd>& relative_velocity_, & solid_ave_vel_, & ave_acc_;
			StdLargeVec<Real>& mass_;
		};
	}

} 
#endif // FISH_AVERAGE_VELOCITY_H
