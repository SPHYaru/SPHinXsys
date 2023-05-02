/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	elastic_dynamics.h
 * @brief 	Here, we define the algorithm classes for elastic solid dynamics.
 * @details 	We consider here a weakly compressible solids.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef COMPOSITE_ELASTIC_DYNAMICS_H
#define COMPOSITE_ELASTIC_DYNAMICS_H

#include "elastic_dynamics.h"
#include "composite_material.h"

namespace SPH
{
	namespace solid_dynamics
	{
		/**
		 * @class CompositeIntegration1stHalf
		 */
		class CompositeIntegration1stHalf : public Integration1stHalf
		{
		public:
			explicit CompositeIntegration1stHalf(BaseInnerRelation &inner_relation);
			virtual ~CompositeIntegration1stHalf(){};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction (size_t index_i, Real dt = 0.0);
		};

	}
}
#endif // COMPOSITE_ELASTIC_DYNAMICS_H