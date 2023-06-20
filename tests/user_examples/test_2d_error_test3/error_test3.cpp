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
Real particle_spacing_ref = 0.05;   /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Extending width for boundary conditions. */
BoundingBox system_domain_bounds(Vec2d(-LL / 2, -LH / 2), Vec2d(LL / 2, LH / 2));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                    /**< Reference density of fluid. */
//Real gravity_g = 9.81;                   /**< Gravity force of fluid. */
Real U_max = 1.0; /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                 /**< Reference sound speed. */

/** Water block shape definition */
class WaterBlock : public MultiPolygonShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(Vecd(0.0, 0.0));
        water_block_shape.push_back(Vecd(0.0, LH));
        water_block_shape.push_back(Vecd(LL, LH));
        water_block_shape.push_back(Vecd(LL, 0.0));
        water_block_shape.push_back(Vecd(0.0, 0.0));
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};

Vec2d outer_wall_halfsize = Vec2d(0.5 * LL + BW, 0.5 * LH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * LL, 0.5 * LH);
Vec2d inner_wall_translation = inner_wall_halfsize;

//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform2d(outer_wall_translation), outer_wall_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform2d(inner_wall_translation), inner_wall_halfsize);
    }
};
/**
 * application dependent initial velocity
 */
class Gradient
    : public CorrectionMatrixInner
{
public:
    Gradient(BaseInnerRelation& inner_relation, int beta, Real alpha)
        : CorrectionMatrixInner(inner_relation,beta,alpha),pos_(particles_->pos_),
        B_(*particles_->registerSharedVariable<Matd>("WeightedCorrectionMatrix"))
    {
        particles_->registerVariable(gradient_value_, "GradientValue");
        particles_->registerVariable(func_, "GiveFunction");
        particles_->registerVariable(gradient_func_, "AnalyticalGradient");
    };

    void interaction(size_t index_i, Real dt)
    {
        Vecd gradient_ = Vecd::Zero();
        const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];

            //gradient_ += (func_[inner_neighborhood.j_[n]]-func_[index_i]) * B_[index_i] * gradW_ij;
            gradient_ += (func_[inner_neighborhood.j_[n]] - func_[index_i]) *  gradW_ij;
        }

        gradient_value_[index_i] = gradient_;        

     
    }

    void update(size_t index_i, Real dt)
    {
        func_[index_i] = exp(-0.2 * pos_[index_i][1])*cos(4 * pos_[index_i][1]);
        gradient_func_[index_i] = -4 * exp(-0.2 * pos_[index_i][1]) * sin(4 * pos_[index_i][1]) + (-0.2) * exp(-0.2 * pos_[index_i][1]) * cos(4 * pos_[index_i][1]);
    }

protected:
    StdLargeVec<Vecd> gradient_value_;
    StdLargeVec<Real> func_,gradient_func_;
    StdLargeVec<Vecd> &pos_;
    StdLargeVec<Matd> &B_;
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-LL / 2, -LH / 2), Vec2d(LL/2, LH/2));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    /** Tag for run particle relaxation for the initial body fitted distribution. */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(false);
    //sph_system.generate_regression_data_ = true;
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<FluidParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    // Using relaxed particle distribution if needed
    sph_system.ReloadParticles()
        ? water_block.generateParticles<ParticleGeneratorReload>(io_environment, water_block.getName())
        : water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");

    ObserverBody fluid_observer(sph_system, "FluidObserver");
    StdVec<Vecd> observation_location = { Vecd(0.0, 0.0) };
    fluid_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    //-------------------------------------------------------------------------------------------------------------------------------------------  
    InnerRelation water_body_inner(water_block);
    ComplexRelation water_block_complex(water_block, { &wall_boundary });
    //----------------------------------------------------------------------
/** check whether run particle relaxation for body fitted particle distribution. */
    if (sph_system.RunParticleRelaxation())
    {
        /**
          * @brief 	Methods used for particle relaxation.
          */
          /** Random reset the insert body particle position. */
        SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_block);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
        /** Write the particle reload files. */
        ReloadParticleIO write_real_body_particle_reload_files(io_environment, sph_system.real_bodies_);

        /** A  Physics relaxation step. */
        relax_dynamics::RelaxationStepInner relaxation_step_inner(water_body_inner, true);
        /**
         * @brief 	Particle relaxation starts here.
         */
        random_water_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_real_body_states.writeToFile(0);

        /** relax particles of the insert body. */
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_real_body_states.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;

        /** Output results. */
        write_real_body_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
     /** time-space method to detect surface particles. */
    InteractionWithUpdate<fluid_dynamics::SpatialTemporalFreeSurfaceIdentificationInner>
        free_surface_indicator(water_body_inner);
    InteractionWithUpdate<Gradient> initial_condition(water_body_inner, 2,0.3);
    InteractionWithUpdate<CorrectionMatrixInner> corrected_configuration_fluid(water_body_inner, 2, 0.3);
    /** modify the velocity of boundary particles with free-stream velocity. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
    PeriodicConditionUsingCellLinkedList periodic_condition_x(water_block, water_block.getBodyShapeBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_y(water_block, water_block.getBodyShapeBounds(), yAxis);
    /** We can output a method-specific particle data for debug */
    water_block.addBodyStateForRecording<Vecd>("GradientValue");
    water_block.addBodyStateForRecording<Real>("GiveFunction");
    water_block.addBodyStateForRecording<Real>("AnalyticalGradient");
    water_block.addBodyStateForRecording<int>("SurfaceIndicator");
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_x.update_cell_linked_list_.exec();
    //periodic_condition_y.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 10.0;
    Real output_interval = 0.01;
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
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */

        /** outer loop for dual-time criteria time-stepping. */
        time_instance = TickCount::now();       
        interval_computing_time_step += TickCount::now() - time_instance;

        time_instance = TickCount::now();
        Real relaxation_time = 0.0;
        Real acoustic_dt = 0.0;
           
        free_surface_indicator.exec();
        //corrected_configuration_fluid.exec();
        initial_condition.exec();

        acoustic_dt = fluid_acoustic_time_step.exec();
        
        relaxation_time += acoustic_dt;
        integration_time += acoustic_dt;
        GlobalStaticVariables::physical_time_ += acoustic_dt;
            
        interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

        /** screen output, write body reduced values and restart files  */
        if (number_of_iterations % screen_output_interval == 0)
        {
           std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
           << GlobalStaticVariables::physical_time_<< "	acoustic_dt = " << acoustic_dt << "\n";      
        }
        number_of_iterations++;

        /** Update cell linked list and configuration. */
        periodic_condition_x.bounding_.exec();
        //periodic_condition_y.bounding_.exec();
        water_block.updateCellLinkedList();
        periodic_condition_x.update_cell_linked_list_.exec();
        //periodic_condition_y.update_cell_linked_list_.exec();
        water_body_inner.updateConfiguration();
        interval_updating_configuration += TickCount::now() - time_instance;
       
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

    return 0;
};
