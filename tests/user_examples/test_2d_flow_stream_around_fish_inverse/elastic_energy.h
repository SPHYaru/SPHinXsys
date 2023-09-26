
#ifndef ELASTIC_ENERGY_H
#define ELASTIC_ENERGY_H

#include "solid_body.h"
#include "solid_particles.h"
#include "elastic_solid.h"

namespace SPH
{
	namespace solid_dynamics
	{
		//----------------------------------------------------------------------
		//		for general solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<ElasticSolidParticles> ElasticSolidDataSimple;

		/**
        * @class ElasticEnergy
        */
		class ElasticEnergy
			: public LocalDynamicsReduce<Real, ReduceSum<Real>>,
			public ElasticSolidDataSimple
		{
		protected:
			StdLargeVec<Real>& Vol_;
			StdLargeVec<Matd>& F_, &stress_PK1_B_, & active_strain_;
			StdLargeVec<int>& materail_id_;

		public:
			ElasticEnergy(SPHBody& sph_body);
			virtual ~ElasticEnergy() {};

			Real reduce(size_t index_i, Real dt = 0.0);
		};

		/**
	    * @class KineticEnergy
	    */
		class SolidKinecticEnergy
			: public LocalDynamicsReduce<Real, ReduceSum<Real>>,
			public ElasticSolidDataSimple
		{
		protected:
			StdLargeVec<Vecd>&vel_;
			StdLargeVec<Real>& mass_;

		public:
			SolidKinecticEnergy(SPHBody& sph_body);
			virtual ~SolidKinecticEnergy() {};

			Real reduce(size_t index_i, Real dt = 0.0);
		};
	}
}
#endif //ELASTIC_ENERGY_H