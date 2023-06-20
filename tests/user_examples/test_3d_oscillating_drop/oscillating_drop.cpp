/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D dambreak example                        *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

/** Set the file path to the stl file. */
std::string full_path_to_stl_file = "./input/qiu.stl";

// general parameters for geometry
Real resolution_ref = 0.02;	  // particle spacing
Real BW = resolution_ref * 4; // boundary width
Real DL = 5.366;			  // tank length
Real DH = 2.0;				  // tank height
Real DW = 2.0;				  // tank width
Real LL = 2.0;				  // liquid length
Real LH = 2.0;				  // liquid height
Real LW = 2.0;				  // liquid width

// for material properties of the fluid
Real rho0_f = 1000.0;
//Real gravity_g = 1.0;
Real U_f = 1.0;
Real c_f = 10.0 * U_f;

//	define the water block shape

class WaterBlock : public ComplexShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(full_path_to_stl_file, Vecd::Zero(), 0.001);
	}
};

//	define an observer particle generator
class WaterObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit WaterObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
	{
		// add observation points
		positions_.push_back(Vecd(DL, 0.01, 0.5 * DW));
	}
};

/**
 * application dependent initial velocity
 */
class InitialVelocity
	: public fluid_dynamics::FluidInitialCondition
{
public:
	InitialVelocity(SPHBody& sph_body)
		: fluid_dynamics::FluidInitialCondition(sph_body) {};

	void update(size_t index_particle_i, Real dt)
	{
		/** initial velocity profile */

		vel_[index_particle_i][0] = -1.0 * pos_[index_particle_i][0];
		vel_[index_particle_i][1] = 1.0 * pos_[index_particle_i][1];		
	}
};

/**
 * application dependent
 */
class ExternalField
	: public fluid_dynamics::FluidInitialCondition
{
public:
	ExternalField(SPHBody& sph_body)
		: fluid_dynamics::FluidInitialCondition(sph_body),
		 fluid_particles_(dynamic_cast<FluidParticles*>(&sph_body.getBaseParticles())),
		 acc_prior_(fluid_particles_->acc_prior_){};

	void update(size_t index_particle_i, Real dt)
	{
		acc_prior_[index_particle_i][0] = -1.0 * 1.0 * pos_[index_particle_i][0];
		acc_prior_[index_particle_i][1] = -1.0 * 1.0 * pos_[index_particle_i][1];
		acc_prior_[index_particle_i][2] = -1.0 * 1.0 * pos_[index_particle_i][2];
	}

protected:
	FluidParticles* fluid_particles_;
	StdLargeVec<Vecd>& acc_prior_;
};

class ExternalEnergy : public TotalMechanicalEnergy
{
public:
	ExternalEnergy(SPHBody& sph_body) :
		TotalMechanicalEnergy(sph_body)
	{
		quantity_name_ = "ExternalEnergy";
	}

	Real reduce(size_t index_particle_i, Real dt)
	{
		return 1.0 * 1.0 * 0.5 * mass_[index_particle_i] * pos_[index_particle_i].squaredNorm();
	}
};

class KineticEnergy : public TotalMechanicalEnergy
{
public:
	KineticEnergy(SPHBody& sph_body) :
		TotalMechanicalEnergy(sph_body)
	{
		quantity_name_ = "KineticEnergy";
	}

	Real reduce(size_t index_particle_i, Real dt)
	{
		return 0.5 * mass_[index_particle_i] * vel_[index_particle_i].squaredNorm();
	}
};

class MaximumCoordinate : public MaximumSpeed
{
public:
	MaximumCoordinate(SPHBody& sph_body) :
		MaximumSpeed(sph_body), 
		fluid_particles_(dynamic_cast<FluidParticles*>(&sph_body.getBaseParticles())), 
		pos_(fluid_particles_->pos_)
	{
		quantity_name_ = "MaximumCoordinate";
	}

	Real reduce(size_t index_particle_i, Real dt)
	{
		return pos_[index_particle_i][1];
	}

protected:
	FluidParticles* fluid_particles_;
	StdLargeVec<Vecd>&pos_;


};


// the main program with commandline options
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up an SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vecd(-DL, - DH, - DW), Vecd(DL + BW, DH + BW, DW + BW));
	SPHSystem system(system_domain_bounds, resolution_ref);

	system.setRunParticleRelaxation(true);
	// Tag for reload initially relaxed particles.
	system.setReloadParticles(false);

	system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
	water_block.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	(!system.RunParticleRelaxation() && system.ReloadParticles())
		? water_block.generateParticles<ParticleGeneratorReload>(io_environment, water_block.getName())
		: water_block.generateParticles<ParticleGeneratorLattice>();

	ObserverBody fluid_observer(system, "FluidObserver");
	fluid_observer.generateParticles<WaterObserverParticleGenerator>();

	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation water_block_inner(water_block);
	ContactRelation fluid_observer_contact(fluid_observer, {&water_block});

	//----------------------------------------------------------------------
	//	check whether run particle relaxation for body fitted particle distribution.
	//----------------------------------------------------------------------
	if (system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		// Random reset the insert body particle position.
		SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(water_block);
		// Write the particle reload files.
		ReloadParticleIO write_particle_reload_files(io_environment, water_block);
		// A  Physics relaxation step.
		relax_dynamics::RelaxationStepInner relaxation_step_inner(water_block_inner);
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_inserted_body_particles.exec(0.25);
		relaxation_step_inner.SurfaceBounding().exec();
		write_states.writeToFile(0);
		//----------------------------------------------------------------------
		//	Particle relaxation loop.
		//----------------------------------------------------------------------
		int ite_p = 0;
		while (ite_p < 1000)
		{
			relaxation_step_inner.exec();
			ite_p += 1;
			if (ite_p % 200 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
				write_states.writeToFile(ite_p);
			}
		}
		std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
		// Output particles position for reload.
		write_particle_reload_files.writeToFile(0);
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define the numerical methods used in the simulation.
	//	Note that there may be data dependence on the sequence of constructions.
	//----------------------------------------------------------------------
	/** external force */
	SimpleDynamics<ExternalField> drop_external_field(water_block);
	/** Initial velocity field */
	SimpleDynamics<InitialVelocity> drop_initial_velocity(water_block);

	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
	Dynamics1Level<fluid_dynamics::Integration1stHalfRiemann> pressure_relaxation(water_block_inner);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemann> density_relaxation(water_block_inner);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> update_density_by_summation(water_block_inner);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

	water_block.addBodyStateForRecording<Real>("Pressure");
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_water_block_states(io_environment, system.real_bodies_);
	ReducedQuantityRecording<ReduceDynamics<KineticEnergy>>
		write_water_kinetic_energy(io_environment, water_block);
	ReducedQuantityRecording<ReduceDynamics<ExternalEnergy>>
		write_water_external_energy(io_environment, water_block);
	ReducedQuantityRecording<ReduceDynamics<MaximumCoordinate>>
		write_coordinate(io_environment, water_block);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
		write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	drop_initial_velocity.exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = system.RestartStep();
	int screen_output_interval = 100;
	Real end_time = 20.0;
	Real output_interval = 0.05;
	Real dt = 0.0; // default acoustic time step sizes
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_water_block_states.writeToFile(0);
	write_water_kinetic_energy.writeToFile(0);
	write_water_external_energy.writeToFile(0);
	write_coordinate.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_interval)
		{
			initialize_a_fluid_step.exec();
			Real Dt = get_fluid_advection_time_step_size.exec();
			update_density_by_summation.exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				drop_external_field.exec();
				pressure_relaxation.exec(dt);
				density_relaxation.exec(dt);
				dt = get_fluid_time_step_size.exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;

			water_block.updateCellLinkedListWithParticleSort(100);
			water_block_inner.updateConfiguration();
			fluid_observer_contact.updateConfiguration();
			write_recorded_water_pressure.writeToFile(number_of_iterations);
		}

		write_water_kinetic_energy.writeToFile(number_of_iterations);
		write_water_external_energy.writeToFile(number_of_iterations);
		write_coordinate.writeToFile(number_of_iterations);

		TickCount t2 = TickCount::now();
		write_water_block_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	//if (system.generate_regression_data_)
	//{
	//	write_water_mechanical_energy.generateDataBase(1.0e-3);
	//	write_recorded_water_pressure.generateDataBase(1.0e-3);
	//}
	//else
	//{
	//	write_water_mechanical_energy.testResult();
	//	write_recorded_water_pressure.testResult();
	//}

	return 0;
}
