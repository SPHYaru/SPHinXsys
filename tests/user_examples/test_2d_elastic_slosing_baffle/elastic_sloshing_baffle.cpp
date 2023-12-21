
/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h"
#define PI 3.1415926
 /**
* @brief Namespace cite here.
*/
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 1.0; 								           /**< Tank length. */
Real DH = 0.7; 							              /**< Tank height. */
Real L_W = 1.0;                                      /**< water width. */
Real L_H = 0.15;                                    /**< water depth. */
Real Gate_x = 0.5 * L_W;						   /**< Width of the gate. */
Real Gate_width = 0.006;						  /**< Width of the gate. */
Real Gate_height = 0.174;						 /**< Height of the gate. */
Real particle_spacing_ref = Gate_width / 2; 	    /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4.0; 		   /**< Extending width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));

// parameters for liquid sloshing Case:Xue Mian
Real a_0 = 0.01; 					     //amplitude of the sloshing
Real w_0 = 1.6 * 4.142814038; 			//frequency of the sloshing  0.583*8.9556
Real k_0 = a_0 * w_0;            // parameter of sloshing x = k_0 * sin(2*pi*f_0*t)
Real k_1 = a_0 * w_0 * w_0;            // parameter of sloshing x = k_0 * sin(2*pi*f_0*t)
Real fre_ = w_0 /(2 * PI);            // parameter of sloshing x = k_0 * sin(2*pi*f_0*t)

/** create a water block shape */
std::vector<Vecd> CreatWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.0, 0.0));
	water_block_shape.push_back(Vecd(0.0, L_H));
	water_block_shape.push_back(Vecd(L_W, L_H));
	water_block_shape.push_back(Vecd(L_W, 0.0));
	water_block_shape.push_back(Vecd(0.0, 0.0));

	return water_block_shape;
}

/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1000.0;							/**< Reference density of fluid. */
Real gravity_g = 9.81; 					        /**< Value of gravity. */
Real U_max = 2.0*sqrt(gravity_g*L_H);		/**< Characteristic velocity. */
Real c_f = 10.0* U_max;					        /**< Reference sound speed. */
/**
 * @brief Material properties of the elastic gate.
 */
Real rho0_s = 1250.0;
Real Youngs_modulus1 = 1.0e5;
Real Youngs_modulus2 = 1.0e5;
Real Youngs_modulus3 = 1.0e5;
Real poisson = 0.47;

// Observer location
StdVec<Vecd> observation_location = { Vecd(DL / 2, 0.1 * Gate_height),Vecd(DL / 2, 0.2 * Gate_height),
									  Vecd(DL / 2, 0.3 * Gate_height),Vecd(DL / 2, 0.4 * Gate_height),
									  Vecd(DL / 2, 0.5 * Gate_height),Vecd(DL / 2, 0.6 * Gate_height),
									  Vecd(DL / 2, 0.6 * Gate_height),Vecd(DL / 2, 0.7 * Gate_height),
									  Vecd(DL / 2, 0.8 * Gate_height),Vecd(DL / 2, 0.9 * Gate_height),
									  Vecd(DL / 2, Gate_height) };

/**
 * @brief 	wall body definition.
 */
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> outer_wall_shape;
		outer_wall_shape.push_back(Vecd(-BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
		outer_wall_shape.push_back(Vecd(DL + BW, -BW));
		outer_wall_shape.push_back(Vecd(-BW, -BW));

		std::vector<Vecd> inner_wall_shape;
		inner_wall_shape.push_back(Vecd(0.0, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, DH));
		inner_wall_shape.push_back(Vecd(DL, DH));
		inner_wall_shape.push_back(Vecd(DL, 0.0));
		inner_wall_shape.push_back(Vecd(0.0, 0.0));

		multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
	}
};

/**
* @brief create a gate shape
*/
std::vector<Vecd> CreatGateShape()
{
	std::vector<Vecd> gate_shape;
	gate_shape.push_back(Vecd(Gate_x - 0.0262, 0.0));
	gate_shape.push_back(Vecd(Gate_x - 0.0262, 0.006));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075, 0.006));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075, 0.0125));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008, 0.0125));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008, 0.026));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002, 0.026));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002, 0.026 + Gate_height));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006, 0.026 + Gate_height));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006, 0.026));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006 + 0.002 + 0.006, 0.026));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006 + 0.002 + 0.006, 0.0125));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006 + 0.002 + 0.006 + 0.008, 0.0125));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006 + 0.002 + 0.006 + 0.008, 0.006));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006 + 0.002 + 0.006 + 0.008 + 0.0075, 0.006));
	gate_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.006 + 0.002 + 0.006 + 0.002 + 0.006 + 0.008 + 0.0075, 0.0));
	gate_shape.push_back(Vecd(Gate_x - 0.0262, 0.0));

	return gate_shape;
}

/**
* @brief create a Gate constrain shape
*/
MultiPolygon CreateGateConstrainShape()
{
	std::vector<Vecd> gate_constrain_shape;
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262, 0.0));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262, 0.006));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075, 0.006));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075, 0.0125));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008, 0.0125));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008, 0.026));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022, 0.026));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022, 0.0125));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022 + 0.008, 0.0125));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022 + 0.008, 0.006));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022 + 0.008 + 0.0075, 0.006));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262 + 0.0075 + 0.008 + 0.022 + 0.008 + 0.0075, 0.0));
	gate_constrain_shape.push_back(Vecd(Gate_x - 0.0262, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(gate_constrain_shape, ShapeBooleanOps::add);
	return multi_polygon;
}
/**
 * @brief  Define the elastic gate body.
 */
class Gate : public MultiPolygonShape
{
public:
	explicit Gate(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(CreatGateShape(), ShapeBooleanOps::add);
	}
};

/**
 * @brief 	Fluid body definition.
 */
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(CreatWaterBlockShape(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(CreatGateShape(), ShapeBooleanOps::sub);
	}
};

//----------------------------------------------------------------------
//	Case dependent composite material
//----------------------------------------------------------------------
class BaffleComposite : public CompositeSolid
{
public:
	BaffleComposite() : CompositeSolid(rho0_s)
	{
		add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus1, poisson);
		add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus2, poisson);
		add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus3, poisson);
	};
};

class BaffelMaterialInitialization
	: public MaterialIdInitialization
{
public:
	explicit BaffelMaterialInitialization(SolidBody& solid_body)
		: MaterialIdInitialization(solid_body) {};

	void update(size_t index_i, Real dt = 0.0)
	{
		Real y = pos0_[index_i][1];

		if (y >= 0.10)
		{
			material_id_[index_i] = 0; // region for muscle
		}
		else if (y <= 0.05)
		{
			material_id_[index_i] = 2;
		}
		else
		{
			material_id_[index_i] = 1;
		}
	};
};

class TimeDependentAcceleration : public Gravity
{
	Real t_ref_, u_ref_, du_ave_dt_;

public:
	explicit TimeDependentAcceleration(Vecd gravity_vector)
		: Gravity(gravity_vector){}

	virtual Vecd InducedAcceleration(const Vecd& position) override
	{
		
		global_acceleration_[0]= -k_1 * sin(w_0 * GlobalStaticVariables::physical_time_);
		global_acceleration_[1] = -gravity_g;

		return  global_acceleration_;
	}
};

Real h = 1.3 * particle_spacing_ref;
Real probe_x1 = DL - 0.032;
MultiPolygon CreateWaveProbeShape1()
{
	std::vector<Vecd> pnts;
	pnts.push_back(Vecd(probe_x1 - h, 0.0));
	pnts.push_back(Vecd(probe_x1 - h, DH));
	pnts.push_back(Vecd(probe_x1 + h, DH));
	pnts.push_back(Vecd(probe_x1 + h, 0.0));
	pnts.push_back(Vecd(probe_x1 - h, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
	return multi_polygon;
}

Real probe_x2 = 0.032;
MultiPolygon  CreateWaveProbeShape2()
{
	std::vector<Vecd> pnts;
	pnts.push_back(Vecd(probe_x2 - h, 0.0));
	pnts.push_back(Vecd(probe_x2 - h, DH));
	pnts.push_back(Vecd(probe_x2 + h, DH));
	pnts.push_back(Vecd(probe_x2 + h, 0.0));
	pnts.push_back(Vecd(probe_x2 - h, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
	return multi_polygon;
}

Real probe_x3 = DL / 2.0;
MultiPolygon  CreateWaveProbeShape3()
{
	std::vector<Vecd> pnts;
	pnts.push_back(Vecd(probe_x3 - h, 0.0));
	pnts.push_back(Vecd(probe_x3 - h, DH));
	pnts.push_back(Vecd(probe_x3 + h, DH));
	pnts.push_back(Vecd(probe_x3 + h, 0.0));
	pnts.push_back(Vecd(probe_x3 - h, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
	return multi_polygon;
}

/**
 * @brief 	Main program starts here.
 */
int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
   //	Build up the environment of a SPHSystem.
   //----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
	water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
	water_block.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();

	SolidBody gate(sph_system, makeShared<Gate>("Gate"));
	gate.defineAdaptationRatios(1.15, 2.0);
	gate.defineParticlesAndMaterial<ElasticSolidParticles, BaffleComposite>();
	gate.generateParticles<ParticleGeneratorLattice>();

	ObserverBody gate_observer(sph_system, "Observer");
	gate_observer.defineAdaptationRatios(1.15, 2.0);
	gate_observer.generateParticles<ObserverParticleGenerator>(observation_location);

	/** topology */
	InnerRelation water_block_inner(water_block);
	ContactRelation water_block_contact(water_block, RealBodyVector{ &wall_boundary, &gate });
	InnerRelation gate_inner(gate);
	ContactRelation gate_water_contact(gate, { &water_block });
	ContactRelation gate_observer_contact(gate_observer, { &gate });
	//----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
	ComplexRelation water_block_complex(water_block_inner, water_block_contact);
	//----------------------------------------------------------------------
   //	Define the main numerical methods used in the simulation.
   //	Note that there may be data dependence on the constructors of these methods.
   //----------------------------------------------------------------------
	InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
		free_surface_indicator(water_block_inner, water_block_contact);
	//----------------------------------------------------------------------
   //	Algorithms of fluid dynamics.
   //----------------------------------------------------------------------
	InteractionWithUpdate<KernelCorrectionMatrixComplex> corrected_configuration_fluid(ConstructorArgs(water_block_inner, 0.3), water_block_contact);
	Dynamics1Level<fluid_dynamics::Integration1stHalfCorrectionWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
	InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(water_block_inner, water_block_contact);
	SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, makeShared<TimeDependentAcceleration>(Vec2d::Zero()));
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
	/** Computing vorticity in the flow. */
	InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(water_block_inner);
	//----------------------------------------------------------------------
	//	Algorithms of FSI.
	//----------------------------------------------------------------------
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	SimpleDynamics<NormalDirectionFromBodyShape> gate_normal_direction(gate);
	InteractionWithUpdate<KernelCorrectionMatrixInner> gate_corrected_configuration(gate_inner);
	InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluidRiemann> fluid_pressure_force_on_gate(gate_water_contact);
	solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(gate);

	//----------------------------------------------------------------------
   //	Algorithms of Elastic dynamics.
   //----------------------------------------------------------------------
	SimpleDynamics<BaffelMaterialInitialization> composite_material_id(gate);
	Dynamics1Level<solid_dynamics::Integration1stHalfPK2> gate_stress_relaxation_first_half(gate_inner);
	Dynamics1Level<solid_dynamics::Integration2ndHalf> gate_stress_relaxation_second_half(gate_inner);
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> gate_computing_time_step_size(gate);

	BodyRegionByParticle gate_constraint_part(gate, makeShared<MultiPolygonShape>(CreateGateConstrainShape()));
	SimpleDynamics<solid_dynamics::FixBodyPartConstraint> gate_constraint(gate_constraint_part);
	SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> gate_update_normal(gate);

	water_block.addBodyStateForRecording<Real>("Pressure");
	water_block.addBodyStateForRecording<int>("Indicator");
	gate.addBodyStateForRecording<int>("MaterialID");
	//----------------------------------------------------------------------
	 //	Define the methods for I/O operations and observations of the simulation.
	 //----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_real_body_states_to_vtp(io_environment, sph_system.real_bodies_);

	/** WaveProbes. */
	BodyRegionByCell wave_probe_buffer_no_1(water_block, makeShared<MultiPolygonShape>(CreateWaveProbeShape1(), "WaveProbe_01"));
	ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
		wave_probe_1(io_environment, wave_probe_buffer_no_1, "FreeSurfaceHeight");

	BodyRegionByCell wave_probe_buffer_no_2(water_block, makeShared<MultiPolygonShape>(CreateWaveProbeShape2(), "WaveProbe_02"));
	ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
		wave_probe_2(io_environment, wave_probe_buffer_no_2, "FreeSurfaceHeight");

	BodyRegionByCell wave_probe_buffer_no_3(water_block, makeShared<MultiPolygonShape>(CreateWaveProbeShape3(), "WaveProbe_03"));
	ReducedQuantityRecording<UpperFrontInAxisDirection<BodyPartByCell>>
		wave_probe_3(io_environment, wave_probe_buffer_no_3, "FreeSurfaceHeight");

	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_gate_displacement("Position", io_environment, gate_observer_contact);
	/**
	 * @brief The time stepping starts here.
	 */
	 /**
	  * @brief Prepare quantities will be used once only and initial condition.
	  */
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_boundary_normal_direction.exec();
	gate_normal_direction.exec();
	gate_corrected_configuration.exec();
	/** initialize material ids for the baffle. */
	composite_material_id.exec();

	//write_real_body_states_to_plt.writeToFile(0);
	write_real_body_states_to_vtp.writeToFile(0);
	wave_probe_1.writeToFile(0);
	wave_probe_2.writeToFile(0);
	wave_probe_3.writeToFile(0);
	write_gate_displacement.writeToFile(0);

	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	Real end_time = 10.0;
	Real output_interval = 0.01;
	Real dt = 0.0;   /**< Default acoustic time step sizes. */
	Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
	/** statistics for computing CPU time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	TimeInterval interval_computing_time_step;
	TimeInterval interval_computing_pressure_relaxation;
	TimeInterval interval_updating_configuration;
	TickCount time_instance;

	/**
	 * @brief Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < output_interval)
		{
			/** Acceleration due to viscous force and gravity. */
			initialize_a_fluid_step.exec();
			Real Dt = get_fluid_advection_time_step_size.exec();
			free_surface_indicator.exec();
			update_density_by_summation.exec();
			corrected_configuration_fluid.exec();
			/** Update normal direction at elastic body surface. */
			gate_update_normal.exec();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				/** Fluid relaxation and force computation. */
				pressure_relaxation.exec(dt);
				fluid_pressure_force_on_gate.exec();
				density_relaxation.exec(dt);
				/** Solid dynamics time stepping. */
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.exec();
				while (dt_s_sum < dt)
				{
					dt_s = gate_computing_time_step_size.exec();
					if (dt - dt_s_sum < dt_s)
						dt_s = dt - dt_s_sum;
					gate_stress_relaxation_first_half.exec(dt_s);
					gate_constraint.exec();
					gate_stress_relaxation_second_half.exec(dt_s);
					dt_s_sum += dt_s;
				}
				average_velocity_and_acceleration.update_averages_.exec(dt);
				dt = get_fluid_time_step_size.exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			water_block.updateCellLinkedListWithParticleSort(100);
			gate.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			gate_water_contact.updateConfiguration();
			
		}
		TickCount t2 = TickCount::now();
		/** write run-time observation into file */
		compute_vorticity.exec();
		write_real_body_states_to_vtp.writeToFile();
		wave_probe_1.writeToFile(number_of_iterations);
		wave_probe_2.writeToFile(number_of_iterations);
		wave_probe_3.writeToFile(number_of_iterations);
		/** write run-time observation into file */
		write_gate_displacement.writeToFile(number_of_iterations);
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();
	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
