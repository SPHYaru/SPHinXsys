/**
 * @file	standingwave.cpp
 * @brief	2D standingwave example.
 * @author	Yaru Ren, Chi Zhang and Xiangyu Hu
 */
#include "fluid_dynamics_complex_wkgc.hpp"
#include "fluid_dynamics_inner_wkgc.hpp"
#include "general_dynamics_wkgc.h"
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real LL = 1.0;                      /**< Liquid column length. */
Real LH = 1.0;                      /**< Liquid column height. */
Real particle_spacing_ref = 0.005;   /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Extending width for boundary conditions. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(LL + BW, LH + BW));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                    /**< Reference density of fluid. */
//Real gravity_g = 9.81;                   /**< Gravity force of fluid. */
Real U_max = 1.0; /**< Characteristic velocity. */
Real c_f = 5.0 * sqrt(2.0);                 /**< Reference sound speed. */

Vec2d DamP_lb(-LL / 2, -LH / 2);	 /**< Left bottom. */
Vec2d DamP_lt(-LL / 2,  LH / 2);	 /**< Left top. */
Vec2d DamP_rt(LL / 2, LH / 2); /**< Right top. */
Vec2d DamP_rb(LL / 2, -LH / 2);	 /**< Right bottom. */

//----------------------------------------------------------------------
//	Geometry definition.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(DamP_lb);
    water_block_shape.push_back(DamP_lt);
    water_block_shape.push_back(DamP_rt);
    water_block_shape.push_back(DamP_rb);
    water_block_shape.push_back(DamP_lb);

    return water_block_shape;
}
class WaterBlock : public MultiPolygonShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
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
        : fluid_dynamics::FluidInitialCondition(sph_body),
        fluid_particles_(dynamic_cast<FluidParticles*>(&sph_body.getBaseParticles())),
        p_(fluid_particles_->p_) {};

    void update(size_t index_i, Real dt)
    {
        Real omega = 1.0;
        /** initial velocity profile */
        vel_[index_i][0] =  omega * pos_[index_i][1];
        vel_[index_i][1] = -omega * pos_[index_i][0];

        for (size_t m = 1; m!= 10; ++m)
        {
            for (size_t n = 1; n!= 10; ++n)
            {
                if (m % 2 == 1 && n % 2 == 1)
                {
                    Real x_star = pos_[index_i][0] + LL / 2;
                    Real y_star = pos_[index_i][1] + LL / 2;
                    Real coefficient1 = m * n * PI * PI * (pow((m * PI / LL), 2) + pow((n * PI / LL), 2));
                    p_[index_i] += rho0_f * (-32 * omega * omega) / coefficient1 * sin(m * PI * x_star / LL)
                        * sin(n * PI * y_star / LL);
                }
                
            }
        }
    }

protected:
    FluidParticles* fluid_particles_;
    StdLargeVec<Real>&p_;
};
//----------------------------------------------------------------------
//	wave gauge
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
MultiPolygon createWaveProbeShape()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(0.0, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-LL, -LL), Vec2d(LL + 10 * BW, LH + 10 * BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(false);
    sph_system.generate_regression_data_ = true;
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineAdaptation<SPHAdaptation>(1.3, 1.0);
    water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    // Using relaxed particle distribution if needed
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<ParticleGeneratorReload>(io_environment, water_block.getName())
        : water_block.generateParticles<ParticleGeneratorLattice>();

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = { Vecd(0.0, 0.0) };
    fluid_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    //-------------------------------------------------------------------------------------------------------------------------------------------  
    InnerRelation water_body_inner(water_block);
    ComplexRelation water_block_complex(water_body_inner, {});
    ContactRelation fluid_observer_contact(fluid_observer, { &water_block });
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    SimpleDynamics<InitialVelocity> initial_condition(water_block);
    InteractionWithUpdate<CorrectionMatrixInner> corrected_configuration_fluid(water_body_inner, 2, 0.3);

    /** time-space method to detect surface particles. */
    InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationInner>
        free_surface_indicator(water_body_inner);
    /** modify the velocity of boundary particles with free-stream velocity. */
    /** Apply transport velocity formulation. */
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionInner> transport_velocity_correction(water_body_inner);
    
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemann> fluid_pressure_relaxation_correct(water_body_inner);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfRiemann> fluid_density_relaxation(water_body_inner);
    //InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> fluid_density_by_summation(water_block_complex);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> fluid_density_by_summation(water_block_complex);
    SimpleDynamics<TimeStepInitialization> fluid_step_initialization(water_block);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
    /** We can output a method-specific particle data for debug */
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<int>("SurfaceIndicator");
    water_block.addBodyStateForRecording<Real>("PositionDivergence");
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RestartIO restart_io(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
        write_water_mechanical_energy(io_environment, water_block);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>>
        write_recorded_water_pressure("Pressure", io_environment, fluid_observer_contact);
  
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    initial_condition.exec();
    //----------------------------------------------------------------------
    //	Load restart file if necessary.
    //----------------------------------------------------------------------
    if (sph_system.RestartStep() != 0)
    {
        GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
        water_block.updateCellLinkedList();
        water_body_inner.updateConfiguration();
    }
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 10.0;
    Real output_interval = 0.1;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_water_mechanical_energy.writeToFile(number_of_iterations);
    write_recorded_water_pressure.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            fluid_step_initialization.exec();
            Real advection_dt = fluid_advection_time_step.exec();
            free_surface_indicator.exec();
            fluid_density_by_summation.exec();
            //corrected_configuration_fluid.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
            while (relaxation_time < advection_dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                acoustic_dt = fluid_acoustic_time_step.exec();
                fluid_pressure_relaxation_correct.exec(acoustic_dt);
                fluid_density_relaxation.exec(acoustic_dt);
                relaxation_time += acoustic_dt;
                integration_time += acoustic_dt;
                GlobalStaticVariables::physical_time_ += acoustic_dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	advection_dt = " << advection_dt << "	acoustic_dt = " << acoustic_dt << "\n";
                
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            write_water_mechanical_energy.writeToFile(number_of_iterations);
            write_recorded_water_pressure.writeToFile(number_of_iterations);

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            water_block.updateCellLinkedListWithParticleSort(100);
            water_body_inner.updateConfiguration();
            interval_updating_configuration += TickCount::now() - time_instance;

            
        }
        body_states_recording.writeToFile();
        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    //if (sph_system.generate_regression_data_)
    //{
        //write_water_mechanical_energy.generateDataBase(1.0e-3);
        //wave_probe.generateDataBase(1.0e-3);
    //}
    //else if (sph_system.RestartStep() == 0)
    //{
        //write_water_mechanical_energy.testResult();
        //wave_probe.testResult();
    //}

    return 0;
};
