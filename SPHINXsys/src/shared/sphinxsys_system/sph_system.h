/**
 * @file sph_system.h
 * @brief The SPH_System managing objects in the system level.
 * @details Note that the system operation prefer these are application independent.
 * @author  Xiangyu Hu, Luhui Han and Chi Zhang
 */
#pragma once

#include "base_data_package.h"
#include "sph_data_conainers.h"

namespace SPH 
{
	/**
	 * @brief Preclaimed classes.
	 */
	class SPHBody;

	/**
	 * @class SPHSystem
	 * @brief The SPHsystem managing objects in the system level.
	 */
	class SPHSystem
	{

	public:
		/**
		 * @brief Default constructor.
		 * Note that the lower and upper domian bounds are used to build mesh cell linked list, 
		 * therefore, should be bigger than the lower and upper bounds of 
		 * all bodies plus the boundary width. Otherwise, it may come to the situation that 
		 * a particle may can not found a valid cell to build the cell linked list.
		 * @param[in] lower_bound Lower bound of the system computational domain.
		 * @param[in] upper_bound Upper bound of the system computational domain.
		 * @param[in] particle_spacing_ref Reference particle spacing.
		 * @param[in] smoothinglength_ratio The Referncen ratio of smoothing length to particle spacing.
		 */
 		SPHSystem(Vecd lower_bound, Vecd upper_bound, Real particle_spacing_ref, 
			int number_of_threads = tbb::task_scheduler_init::automatic);
		virtual ~SPHSystem();

		Vecd lower_bound_, upper_bound_;	/**< Lower and Upper domain bound. */
		Real particle_spacing_ref_;			/**< Refernce initial particle spacing. */
		/** start the simulation with relaxed particles*/
		bool reload_particles_;
		/** restart step*/
		int restart_step_;
		/** computing from roeload particles from files. */
		bool run_particle_relaxation_;

		task_scheduler_init tbb_init_;		/**< TBB library. */

		SPHBodyVector bodies_;			/**< All sph bodies. */
		SPHBodyVector fictitious_bodies_;/**< The bodies without inner particle configuration. */
		SPHBodyVector real_bodies_;		/**< The bodies with inner particle configuration. */
		SPHBodyTopology* body_topology_;	/**< SPH body topology. */

		/** Add a new body to the SPH system. */
		void AddBody(SPHBody* body);
		/** Add a new body to the SPH real bodies. */
		void AddRealBody(SPHBody* body);
		/** Add a new body to the SPH fictitious bodies. */
		void AddFictitiousBody(SPHBody* body);
		/** Set up the body topology and allocate memeries 
		  * for mesh cell linked lists and particle configurations. */
		void SetBodyTopology(SPHBodyTopology* body_topology);
		/** Initialize cell linked lists. */
		void InitializeSystemCellLinkedLists();
		/** Initialize particle interacting configurations. */
		void InitializeSystemConfigurations();
		/** Set up cell-linked list and configuration for simulation. */
		void SetupSPHSimulation();
	};
}