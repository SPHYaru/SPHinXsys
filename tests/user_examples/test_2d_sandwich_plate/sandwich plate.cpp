/* ---------------------------------------------------------------------------*
 *            SPHinXsys: 2D oscillation beam example-one body version           *
 * ----------------------------------------------------------------------------*
 * This is the one of the basic test cases, also the first case for            *
 * understanding SPH method for solid simulation.                              *
 * In this case, the constraint of the beam is implemented with                *
 * internal constrained subregion.                                             *
 * ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
#include "composite_material.h"
#include "elastic_energy.h"
#include "composite_elastic_dynamics.h"
using namespace SPH;
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.2;	// beam length
Real PH = 0.02; // for thick plate; =0.01 for thin plate
Real SL = 0.06; // depth of the insert
// reference particle spacing
Real resolution_ref = PH / 40.0;
Real BW = resolution_ref * 4; // boundary width, at least three particles
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0),
								 Vec2d(PL + 3.0 * BW, PL / 2.0));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_s1 = 1.0e2;		 //reference density
Real Youngs_modulus1 = 1.0e6; //reference Youngs modulus
Real poisson1 = 0.3975;		 //Poisson ratio

Real rho0_s2 = 1.0e3;		 //reference density
Real Youngs_modulus2 = 1.0e6; //reference Youngs modulus
Real poisson2 = 0.3975;		 //Poisson ratio

Real equivalent_Youngs = (7.0 * Youngs_modulus1 + Youngs_modulus2) / 8;
Real k_0 = equivalent_Youngs / (3 - 6 * poisson1);
Real c_s = sqrt(k_0 / rho0_s1);
//----------------------------------------------------------------------
//	Parameters for initial condition on velocity
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.01;
Real R = PL / (0.5 * Pi);
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
// a beam base shape
std::vector<Vecd> beam_base_shape{
	Vecd(-SL - BW, -PH / 2 - BW), Vecd(-SL - BW, PH / 2 + BW), Vecd(0.0, PH / 2 + BW),
	Vecd(0.0, -PH / 2 - BW), Vecd(-SL - BW, -PH / 2 - BW)};
// a beam shape
std::vector<Vecd> beam_shape{
	Vecd(-SL, -PH / 2), Vecd(-SL, PH / 2), Vecd(PL, PH / 2), Vecd(PL, -PH / 2), Vecd(-SL, -PH / 2)};
// Beam observer location
StdVec<Vecd> observation_location = {Vecd(PL, 0.0)};
//----------------------------------------------------------------------
//	Define the beam body
//----------------------------------------------------------------------
class Beam : public MultiPolygonShape
{
public:
	explicit Beam(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(beam_shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class BeamInitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	explicit BeamInitialCondition(SPHBody &sph_body)
		: solid_dynamics::ElasticDynamicsInitialCondition(sph_body){};

	void update(size_t index_i, Real dt)
	{
		/** initial velocity profile */
		Real x = pos_[index_i][0] / PL;
		if (x > 0.0)
		{
			vel_[index_i][1] = vf * c_s *
							   (M * (cos(kl * x) - cosh(kl * x)) - N * (sin(kl * x) - sinh(kl * x))) / Q;
		}
	};
};

//----------------------------------------------------------------------
class SolidBodyMaterial : public CompositeMaterial
{
public:
	SolidBodyMaterial() : CompositeMaterial(rho0_s1)
	{
		add<SaintVenantKirchhoffSolid>(rho0_s1, Youngs_modulus1, poisson1);
		add<SaintVenantKirchhoffSolid>(rho0_s2, Youngs_modulus2, poisson2);
	};
};

//	Setup material ID
//----------------------------------------------------------------------
class MaterialId
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	explicit  MaterialId(SolidBody& solid_body)
		: solid_dynamics::ElasticDynamicsInitialCondition(solid_body) ,
		solid_particles_(dynamic_cast<SolidParticles*>(&solid_body.getBaseParticles())),
		materail_id_(*solid_particles_->getVariableByName<int>("MaterailId"))
	{};
	virtual void update(size_t index_i, Real dt = 0.0)
	{
		if (abs(pos_[index_i][1]) < 0.005 )
		{
			materail_id_[index_i] = 1;
		}
		else
		{
			materail_id_[index_i] = 0;
		}

	};

protected:
	SolidParticles* solid_particles_;
	StdLargeVec<int>& materail_id_;
};

//----------------------------------------------------------------------
//	define the beam base which will be constrained.
//----------------------------------------------------------------------
MultiPolygon createBeamConstrainShape()
{
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(beam_shape, ShapeBooleanOps::sub);
	return multi_polygon;
};
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
#ifdef BOOST_AVAILABLE
	// handle command line arguments
	system.handleCommandlineOptions(ac, av);
#endif //----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody beam_body(system, makeShared<Beam>("BeamBody"));
	beam_body.defineParticlesAndMaterial<ElasticSolidParticles, SolidBodyMaterial>();
	beam_body.generateParticles<ParticleGeneratorLattice>();

	ObserverBody beam_observer(system, "BeamObserver");
	beam_observer.defineAdaptationRatios(1.15, 2.0);
	beam_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation beam_body_inner(beam_body);
	ContactRelation beam_observer_contact(beam_observer, {&beam_body});
	//-----------------------------------------------------------------------------
	// this section define all numerical methods will be used in this case
	//-----------------------------------------------------------------------------
	SimpleDynamics<MaterialId> CompositematerialID(beam_body);
	SimpleDynamics<BeamInitialCondition> beam_initial_velocity(beam_body);
	// corrected strong configuration
	InteractionDynamics<solid_dynamics::CorrectConfiguration> beam_corrected_configuration(beam_body_inner);
	// time step size calculation
	ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(beam_body);
	// stress relaxation for the beam
	Dynamics1Level<solid_dynamics::Integration1stHalf> stress_relaxation_first_half(beam_body_inner);
	Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(beam_body_inner);
	// clamping a solid body part. This is softer than a direct constraint
	BodyRegionByParticle beam_base(beam_body, makeShared<MultiPolygonShape>(createBeamConstrainShape()));
	SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_beam_base(beam_base);

	beam_body.addBodyStateForRecording<int>("MaterailId");
	//-----------------------------------------------------------------------------
	// outputs
	//-----------------------------------------------------------------------------
	IOEnvironment io_environment(system);
	BodyStatesRecordingToVtp write_beam_states(io_environment, system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<Vecd>>
		write_beam_tip_displacement("Position", io_environment, beam_observer_contact);
	/** Elastic Energy of beam. */
	RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::ElasticEnergy>>>
		write_beam_elastic_energy(io_environment, beam_body);
	RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::SolidKinecticEnergy>>>
		write_beam_kinetic_energy(io_environment, beam_body);
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	beam_initial_velocity.exec();
	beam_corrected_configuration.exec();
	CompositematerialID.exec();
	//----------------------------------------------------------------------
	//	Setup computing time-step controls.
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 1.0;
	Real end_time = T0;
	// time step size for output file
	Real output_interval = 0.0001 * T0;
	Real Dt = 0.1 * output_interval; /**< Time period for data observing */
	Real dt = 0.0;					 // default acoustic time step sizes

	// statistics for computing time
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	//-----------------------------------------------------------------------------
	// from here the time stepping begins
	//-----------------------------------------------------------------------------
	write_beam_states.writeToFile(0);
	write_beam_tip_displacement.writeToFile(0);

	// computation loop starts
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		// integrate time (loop) until the next output time
		while (integration_time < output_interval)
		{

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				stress_relaxation_first_half.exec(dt);
				constraint_beam_base.exec();
				stress_relaxation_second_half.exec(dt);

				ite++;
				dt = computing_time_step_size.exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: "
							  << dt << "\n";
				}
			}
		}

		write_beam_tip_displacement.writeToFile(ite);
		write_beam_elastic_energy.writeToFile(ite);
		write_beam_kinetic_energy.writeToFile(ite);

		TickCount t2 = TickCount::now();
		write_beam_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	if (system.generate_regression_data_)
	{
		// The lift force at the cylinder is very small and not important in this case.
		write_beam_tip_displacement.generateDataBase(Vec2d(1.0e-2, 1.0e-2), Vec2d(1.0e-2, 1.0e-2));
	}
	else
	{
		write_beam_tip_displacement.newResultTest();
	}

	return 0;
}
