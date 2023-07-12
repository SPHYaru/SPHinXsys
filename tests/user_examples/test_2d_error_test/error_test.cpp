/**
 * @file	standingwave.cpp
 * @brief	2D standingwave example.
 * @author	Yaru Ren, Chi Zhang and Xiangyu Hu
 */
#include "fluid_dynamics_complex_wkgc.hpp"
#include "fluid_dynamics_inner_wkgc.hpp"
#include "general_dynamics_wkgc.h"
#include <algorithm>
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real LL = 1.0;                      /**< Liquid column length. */
Real LH = 1.0;                      /**< Liquid column height. */
Real particle_spacing_ref = 0.025;   /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Extending width for boundary conditions. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(LL + BW, LH + BW));
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                    /**< Reference density of fluid. */
//Real gravity_g = 9.81;                   /**< Gravity force of fluid. */
Real U_max = 1.0; /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                 /**< Reference sound speed. */

Vec2d insert_circle_center(0.0, 0.0);	/**< Location of the cylinder center. */
Real insert_circle_radius = 1.0;		/**< Radius of the cylinder. */

/** Water block shape definition */
class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        MultiPolygon multi_polygon;
        multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
        add<MultiPolygonShape>(multi_polygon);
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
        : CorrectionMatrixInner(inner_relation, beta, alpha),pos_(particles_->pos_)
    {
        particles_->registerVariable(gradient_value_, "GradientValue");
        particles_->registerVariable(func_, "GiveFunction");
        particles_->registerVariable(ana_gradient_, "AnalyticalGradient");
    };

    void interaction(size_t index_i, Real dt)
    {
        Vecd gradient_ = Vecd::Zero();
        const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];

            gradient_ += (func_[inner_neighborhood.j_[n]]- func_[index_i]) * B_[index_i] * gradW_ij;
            //gradient_ += (func_[inner_neighborhood.j_[n]] - func_[index_i]) *  gradW_ij;
        }

        gradient_value_[index_i] = gradient_;        
    }

    void update(size_t index_i, Real dt)
    {
        func_[index_i] = exp(-pow(pos_[index_i][0],2)/0.1);
        ana_gradient_[index_i][0] = - 2.0 * (pos_[index_i][0])/0.1 * exp(-pow(pos_[index_i][0], 2) / 0.1);
        ana_gradient_[index_i][1] = 0.0;
        //func_[index_i] = exp(-0.2 * pos_[index_i][0]) * cos(4.0 * pos_[index_i][0]);
        //ana_gradient_[index_i][0] = -exp(-0.2 * pos_[index_i][0]) * (0.2 * cos(4.0 * pos_[index_i][0]) + 4.0 * sin(4.0 * pos_[index_i][0]));
        //ana_gradient_[index_i][1] = 0.0;
    }

protected:
    StdLargeVec<Vecd> gradient_value_, ana_gradient_;
    StdLargeVec<Real> func_;
    StdLargeVec<Vecd> &pos_;
};

class OriginalB
    : public CorrectionMatrixInner
{
public:
    OriginalB(BaseInnerRelation& inner_relation, int beta, Real alpha)
        : CorrectionMatrixInner(inner_relation, beta, alpha), pos_(particles_->pos_)
    {
        particles_->registerVariable(original_B_, "OriginalB");
        particles_->registerVariable(matrix_error_, "MatrixError");
        particles_->registerVariable(matrix_difference_, "MatrixDifference");
        particles_->registerVariable(matrix_determinant_, "MatrixDeterminant");
    };

    void interaction(size_t index_i, Real dt)
    {
        Matd local_configuration = Eps * Matd::Identity();

        const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            Vecd gradW_ij = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
            local_configuration -= r_ji * gradW_ij.transpose();
        }

        original_B_[index_i] = local_configuration.inverse();
        matrix_determinant_[index_i] = original_B_[index_i].determinant();
    }

    void update(size_t index_i, Real dt)
    {
        matrix_difference_[index_i] = B_[index_i] - original_B_[index_i];
        matrix_error_[index_i] = sqrt(pow(matrix_difference_[index_i](0, 0),2) + pow(matrix_difference_[index_i](0, 1), 2) + pow(matrix_difference_[index_i](1, 0), 2) + pow(matrix_difference_[index_i](1, 1), 2));
    }

protected:
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd> original_B_, matrix_difference_;
    StdLargeVec<Real> matrix_error_, matrix_determinant_;
};

class Error
    : public fluid_dynamics::FluidInitialCondition
{
public:
    Error(SPHBody& sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body),
        fluid_particles_(dynamic_cast<FluidParticles*>(&sph_body.getBaseParticles())),
        surface_indicator_(fluid_particles_->surface_indicator_),
        gradient_value_(*fluid_particles_->registerSharedVariable<Vecd>("GradientValue")),
        ana_gradient_(*fluid_particles_->registerSharedVariable<Vecd>("AnalyticalGradient")),
        matrix_error_(*fluid_particles_->registerSharedVariable<Real>("MatrixError"))
    {
        fluid_particles_->registerVariable(L_2_gradient_m, "L2normGradientM");
        fluid_particles_->registerVariable(L_2_gradient_ana, "L2normGradientAna");
        fluid_particles_->registerVariable(L_2_gradient_inf, "L2normGradientInfinite");
        fluid_particles_->registerVariable(relative_difference, "RelativeDifference");
        fluid_particles_->registerVariable(L_2_matrix_, "L2normMatrix");
    };

    void setupDynamics(Real dt)
    {
        Real L_2_1(0.0),L_2_2(0.0),L_2_ana(0.0);
        std::vector<Real> L_2_inf;

        int m = 0;
        for (size_t i = 0; i != surface_indicator_.size(); ++i)
        {
            relative_difference[i] = abs(gradient_value_[i][0] - ana_gradient_[i][0]);

            if (surface_indicator_[i]==0)
            {
                L_2_1 += pow((gradient_value_[i][0] - ana_gradient_[i][0]), 2);
                L_2_ana += pow(ana_gradient_[i][0], 2);
                L_2_inf.push_back(relative_difference[i]);
                L_2_2 += matrix_error_[i];
                m++;
            }          
        }

        L_2_gradient_m[0]   = sqrt(L_2_1) / m;  //L_2 norm divide the particle number
        L_2_gradient_ana[0] = sqrt(L_2_1) / sqrt(L_2_ana); //L_2 norm divide the L_2 anlaytical results
        L_2_gradient_inf[0] = *std::max_element(L_2_inf.begin(), L_2_inf.end());
        L_2_matrix_[0] = L_2_2 / m;
    }

    void update(size_t index_i, Real dt)
    {
        L_2_gradient_m[index_i] = L_2_gradient_m[0];
        L_2_gradient_ana[index_i] = L_2_gradient_ana[0];
        L_2_gradient_inf[index_i] = L_2_gradient_inf[0];
        L_2_matrix_[index_i] = L_2_matrix_[0];
    }

protected:
    FluidParticles* fluid_particles_;
    StdLargeVec<int> & surface_indicator_;
    StdLargeVec<Vecd> &gradient_value_, &ana_gradient_;
    StdLargeVec<Real> L_2_gradient_m, L_2_gradient_ana, L_2_gradient_inf, relative_difference, L_2_matrix_;
    StdLargeVec<Real> & matrix_error_;
};

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
    //sph_system.generate_regression_data_ = true;
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineBodyLevelSetShape();
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
    InteractionWithUpdate<Gradient> calculate_gradient(water_body_inner, 2, 0.0);
    InteractionWithUpdate<OriginalB> calculate_matrix_difference(water_body_inner, 2, 0.0);
    InteractionWithUpdate<CorrectionMatrixInner> corrected_configuration_fluid(water_body_inner, 2, 0.0);
    SimpleDynamics<Error> calculate_error(water_block);
    /** modify the velocity of boundary particles with free-stream velocity. */

    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
    /** We can output a method-specific particle data for debug */
    water_block.addBodyStateForRecording<Vecd>("GradientValue");
    water_block.addBodyStateForRecording<Real>("GiveFunction");
    water_block.addBodyStateForRecording<int>("SurfaceIndicator");
    water_block.addBodyStateForRecording<Vecd>("AnalyticalGradient");
    water_block.addBodyStateForRecording<Matd>("WeightedCorrectionMatrix");
    water_block.addBodyStateForRecording<Matd>("OriginalB");
    water_block.addBodyStateForRecording<Real>("MatrixDeterminant");
    water_block.addBodyStateForRecording<Matd>("MatrixDifference");
    water_block.addBodyStateForRecording<Real>("L2normGradientM");
    water_block.addBodyStateForRecording<Real>("L2normGradientAna");
    water_block.addBodyStateForRecording<Real>("L2normGradientInfinite"); 
    water_block.addBodyStateForRecording<Real>("RelativeDifference");
    water_block.addBodyStateForRecording<Real>("L2normMatrix");
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
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();       
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real acoustic_dt = 0.0;
           
            free_surface_indicator.exec();
            corrected_configuration_fluid.exec();
            calculate_gradient.exec();
            calculate_matrix_difference.exec();
            calculate_error.exec();

            acoustic_dt = fluid_acoustic_time_step.exec();
        
            relaxation_time += acoustic_dt;
            integration_time += acoustic_dt;
            GlobalStaticVariables::physical_time_ += acoustic_dt;
            
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	acoustic_dt = " << acoustic_dt << "\n";      
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            //water_block.updateCellLinkedListWithParticleSort(100);
            //water_body_inner.updateConfiguration();
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

    return 0;
};
