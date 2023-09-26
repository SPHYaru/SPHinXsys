
#ifndef ANGULAR_MOMENTUM_H
#define ANGULAR_MOMENTUM_H

#include "elastic_dynamics.h"

namespace SPH
{
	namespace solid_dynamics
	{
		typedef DataDelegateSimple<SolidParticles> SolidDataSimple;

		/**
		* @class SolidTorque
		*/
		class SolidTorque : public  SolidDataSimple, public LocalDynamics
		{

		private:
			Real  frequency_, covariance_, omega, omega_filter, theta;
			std::vector<std::pair<Real, Real>> dynamicArray;
			Real sum, avg_direction, theta_filter;
			int mm;
			
		protected:
			StdLargeVec<Vecd>& vel_, & pos_, & pos0_;
			StdLargeVec<Real>& mass_;
			StdLargeVec<Real>& theta_, & omega_, & theta_filter_;

			virtual void setupDynamics(Real dt) override;

		public:
			explicit SolidTorque(SPHBody& sph_body, Real frequency, Real covariance);
			virtual ~SolidTorque() {};

			void update(size_t index_i, Real dt = 0.0);
		
		};
	}
} // 
#endif // ANGULAR_MOMENTUM_H
