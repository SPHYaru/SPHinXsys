
#include "2d_fish_and_bones.h"
#include "composite_material.h"
#include "sphinxsys.h"
#define PI 3.1415926
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.8;                                /**< Channel length. */
Real DH = 0.4;                                /**< Channel height. */
Real particle_spacing_ref = 0.0025;           /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = particle_spacing_ref * 4.0;         /**< Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));

Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;

Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
Vec2d emitter_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d emitter_buffer_translation = Vec2d(-DL_sponge, 0.0) + emitter_buffer_halfsize;
Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                /**< Density. */
Real U_f = 1.0;                      /**< freestream velocity. */
Real c_f = 10.0 * U_f;               /**< Speed of sound. */
Real Re = 30000.0;                   /**< Reynolds number. */
Real mu_f = rho0_f * U_f * 0.3 / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
//----------------------------------------------------------------------
Real cx = 0.3 * DL;           /**< Center of fish in x direction. */
Real cy = DH / 2;             /**< Center of fish in y direction. */
Real fish_length = 0.2;       /**< Length of fish. */
Real fish_thickness = 0.03;   /**< The maximum fish thickness. */
Real muscle_thickness = 0.02; /**< The maximum fish thickness. */
Real head_length = 0.03;      /**< Length of fish bone. */
Real bone_thickness = 0.003;  /**< Length of fish bone. */
Real fish_shape_resolution = particle_spacing_ref * 0.5;

Real rho0_s = 1050.0;
Real Youngs_modulus1 = 0.8e6;
Real Youngs_modulus2 = 0.5e6;
Real Youngs_modulus3 = 1.1e6;
Real poisson = 0.49;

Real a1 = 1.22 * fish_thickness / fish_length;
Real a2 = 3.19 * fish_thickness / fish_length / fish_length;
Real a3 = -15.73 * fish_thickness / pow(fish_length, 3);
Real a4 = 21.87 * fish_thickness / pow(fish_length, 4);
Real a5 = -10.55 * fish_thickness / pow(fish_length, 5);

Real b1 = 1.22 * muscle_thickness / fish_length;
Real b2 = 3.19 * muscle_thickness / fish_length / fish_length;
Real b3 = -15.73 * muscle_thickness / pow(fish_length, 3);
Real b4 = 21.87 * muscle_thickness / pow(fish_length, 4);
Real b5 = -10.55 * muscle_thickness / pow(fish_length, 5);

Vec2d rotation_center_(cx + fish_length / 2, DH / 2);
//----------------------------------------------------------------------
//	SPH bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

    return water_block_shape;
}

/**
 * Fish body with tethering constraint.
 */
class FishBody : public MultiPolygonShape
{

  public:
    explicit FishBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> fish_shape = CreatFishShape(cx, cy, fish_length, fish_shape_resolution);
        multi_polygon_.addAPolygon(fish_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Define case dependent bodies material, constraint and boundary conditions.
//----------------------------------------------------------------------
/** Fluid body definition */
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
        MultiPolygon fish(CreatFishShape(cx, cy, fish_length, fish_shape_resolution));
        subtract<MultiPolygonShape>(fish);
    }
};

//----------------------------------------------------------------------
//	Free-stream velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType& boundary_condition)
        : u_ref_(0.0), t_ref_(2.0) {}

    Vecd operator()(Vecd& position, Vecd& velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real run_time = GlobalStaticVariables::physical_time_;
        target_velocity[0] = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        return target_velocity;
    }
};

// ----------------------------------------------------------------------
//	Free-stream velocity
//----------------------------------------------------------------------
class FreeStreamVelocity : public fluid_dynamics::EmitterInflowCondition
{

    Real u_ave_, u_ref_, t_ref_;

public:
    FreeStreamVelocity(BodyAlignedBoxByParticle& aligned_box_part, SolidBody& solid_body)
        : EmitterInflowCondition(aligned_box_part),
        solid_particles_(dynamic_cast<SolidParticles*>(&solid_body.getBaseParticles())),
        relative_velocity_(*solid_particles_->getVariableByName<Vecd>("RelativeAverageVelocity")),
        theta_(*solid_particles_->getVariableByName<Real>("Theta")),
        theta_filter_(*solid_particles_->getVariableByName<Real>("Theta_filter")),
        u_ave_(0), u_ref_(1.0), t_ref_(2.0)
    {
        rotation_angle_ = 0.0;
    }

protected:
    Real rotation_angle_;
    SolidParticles* solid_particles_;
    StdLargeVec<Vecd>& relative_velocity_;
    StdLargeVec<Real>& theta_, & theta_filter_;

    virtual void updateTransform()override
    {
        rotation_angle_ = 0.0;//theta_filter_[0];

        old_transform_ = updated_transform_;

        Vec2d center = rotation_center_;
        Vec2d emitter_center_offset = emitter_buffer_translation - center;

        Rotation new_rotation(rotation_angle_);
        Vec2d new_translation(new_rotation.toRotationMatrix()*(emitter_center_offset) + center);

        Transform new_transform(new_rotation, new_translation);

        updated_transform_ = new_transform;

    }

    virtual Vecd getTargetVelocity(Vecd& position, Vecd& velocity) override
    {
        Vecd target_velocity = Vecd::Zero();

        target_velocity[0] = relative_velocity_[0][0];

        return target_velocity;
    }

};

//----------------------------------------------------------------------
//freestream correction
//----------------------------------------------------------------------
class FreeStreamCorrection : public fluid_dynamics::FreeStreamVelocityCorrection<InflowVelocity>
{
    Real u_ave_, u_ref_, t_ref_;

public:
    explicit FreeStreamCorrection(SPHBody& sph_body, SolidBody& solid_body, const Transform& transform = Transform())
        : FreeStreamVelocityCorrection<InflowVelocity>(sph_body, transform),
        solid_particles_(dynamic_cast<SolidParticles*>(&solid_body.getBaseParticles())),
        relative_velocity_(*solid_particles_->getVariableByName<Vecd>("RelativeAverageVelocity")),
        theta_(*solid_particles_->getVariableByName<Real>("Theta")),
        theta_filter_(*solid_particles_->getVariableByName<Real>("Theta_filter")),
        u_ave_(0), u_ref_(1.0), t_ref_(2.0)
    {
        rotation_angle_ = 0.0;
    }

    void setupDynamics(Real dt = 0.0)
    {
        rotation_angle_ = 0.0;//theta_filter_[0];

        Vec2d center = rotation_center_;
        Vec2d emitter_center_offset = Vecd(DL / 2, DH / 2) - rotation_center_;

        Rotation new_rotation(rotation_angle_);
        Vec2d new_translation(new_rotation.toRotationMatrix() * (emitter_center_offset) + center);

        transform_ = Transform(new_rotation, new_translation);
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        if (surface_indicator_[index_i] == 1)
        {
            Vecd frame_position = transform_.shiftBaseStationToFrame(pos_[index_i]);
            Vecd frame_velocity = transform_.xformBaseVecToFrame(vel_[index_i]);
            Real frame_u_stream_direction = frame_velocity[0];
            Real u_freestream = relative_velocity_[0][0];
            frame_velocity[0] = u_freestream + (frame_u_stream_direction - u_freestream) *
                SMIN(rho_sum_[index_i], rho0_) / rho0_;
            vel_[index_i] = transform_.xformFrameVecToBase(frame_velocity);
        }
    };

protected:
    Real rotation_angle_;
    SolidParticles* solid_particles_;
    StdLargeVec<Vecd>& relative_velocity_;
    StdLargeVec<Real>& theta_, & theta_filter_;
};

//outflow
class ModifiedDisposerOutflowDeletion2 : public LocalDynamics
{
public:
    ModifiedDisposerOutflowDeletion2(FluidBody& fluid_body, SolidBody& solid_body, int axis, const Transform& transform = Transform())
        : LocalDynamics(fluid_body), transform_(transform),
        fluid_particles_(dynamic_cast<BaseParticles*>(&fluid_body.getBaseParticles())),
        solid_particles_(dynamic_cast<SolidParticles*>(&solid_body.getBaseParticles())),
        pos_(fluid_particles_->pos_),
        theta_(*solid_particles_->getVariableByName<Real>("Theta")),
        theta_filter_(*solid_particles_->getVariableByName<Real>("Theta_filter"))
    {
        rotation_angle_ = 0.0;
    }

    void setupDynamics(Real dt = 0.0)
    {
        rotation_angle_ = 0.0;//theta_filter_[0];

        Vec2d center_ = rotation_center_;
        Vec2d disposer_center_offset = disposer_translation - center_;

        Rotation new_rotation(rotation_angle_);
        Vec2d new_translation(new_rotation.toRotationMatrix() * (disposer_center_offset) + center_);

        transform_ = Transform(new_rotation, new_translation);
    }
   
    void update(size_t index_i, Real dt)
    {
        mutex_switch_to_buffer_.lock();
        while ((transform_.shiftBaseStationToFrame(pos_[index_i])[0] > disposer_halfsize[0])
            && index_i < fluid_particles_->total_real_particles_)
        {
            fluid_particles_->switchToBufferParticle(index_i);
        }
        mutex_switch_to_buffer_.unlock();
    }

protected:
    std::mutex mutex_switch_to_buffer_;
    Transform transform_;
    Real rotation_angle_;
    BaseParticles* fluid_particles_;
    SolidParticles* solid_particles_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Real>& theta_, & theta_filter_;
};


//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class TimeDependentAcceleration : public Gravity
{
    Real t_ref_, u_ref_, du_ave_dt_;

  public:
    explicit TimeDependentAcceleration(Vecd gravity_vector, SolidBody& solid_body)
        : Gravity(gravity_vector), solid_particles_(dynamic_cast<SolidParticles*>(&solid_body.getBaseParticles())), 
        t_ref_(2.0), u_ref_(0.0), du_ave_dt_(0)
    {
        solid_particles_->registerVariable(ave_acc_, "RelativeAverageAcceleration");
        solid_particles_->registerVariable(relative_velocity_, "RelativeAverageVelocity");
        solid_particles_->registerVariable(solid_ave_vel_, "SolidAverageVelocity");
        solid_particles_->registerVariable(theta_, "Theta");
        solid_particles_->registerVariable(omega_, "AngularVelocity");
        solid_particles_->registerVariable(theta_filter_, "Theta_filter");
    }

    virtual Vecd InducedAcceleration(Vecd &position) override
    {
        //Real run_time_ = GlobalStaticVariables::physical_time_;
        //du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);
        //return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
        return  solid_ave_vel_[0];
    }

protected:
    //Real rotation_angle_;
    SolidParticles* solid_particles_;
    StdLargeVec<Vecd> ave_acc_, relative_velocity_, solid_ave_vel_;
    StdLargeVec<Real> theta_, omega_, theta_filter_;
};

// Material ID
class SolidBodyMaterial : public CompositeMaterial
{
  public:
    SolidBodyMaterial() : CompositeMaterial(rho0_s)
    {
        add<ActiveModelSolid>(rho0_s, Youngs_modulus1, poisson);
        add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus2, poisson);
        add<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus3, poisson);
    };
};
//----------------------------------------------------------------------
//	Case dependent initialization material ids
//----------------------------------------------------------------------
class FishMaterialInitialization
    : public MaterialIdInitialization
{
  public:
    explicit FishMaterialInitialization(SolidBody &solid_body)
        : MaterialIdInitialization(solid_body){};

    void update(size_t index_i, Real dt = 0.0)
    {
        Real x = pos0_[index_i][0] - cx;
        Real y = pos0_[index_i][1];

        Real x1 = abs(pos0_[index_i][0] - (cx + fish_length));

        Real y1 = a1 * pow(x1, 0 + 1) + a2 * pow(x1, 1 + 1) + a3 * pow(x1, 2 + 1) + a4 * pow(x1, 3 + 1) + a5 * pow(x1, 4 + 1);
        if (x >=  head_length && y > (y1 - 0.004 + cy) && y > (cy + bone_thickness / 2))
        {
            material_id_[index_i] = 0; // region for muscle
        }
        else if (x >= head_length && y < (-y1 + 0.004 + cy) && y < (cy - bone_thickness / 2))
        {
            material_id_[index_i] = 0; // region for muscle
        }
        else if ((x < head_length) || ((y < (cy + bone_thickness / 2)) && (y > (cy - bone_thickness / 2))))
        {
            material_id_[index_i] = 2;
        }
        else
        {
            material_id_[index_i] = 1;
        }
    };
};

// imposing active strain to fish muscle
class ImposingActiveStrain
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit ImposingActiveStrain(SolidBody &solid_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(solid_body),
          material_id_(*particles_->getVariableByName<int>("MaterialID")),
          pos0_(*particles_->getVariableByName<Vecd>("InitialPosition"))
    {
        particles_->registerVariable(active_strain_, "ActiveStrain");
    };
    virtual void update(size_t index_i, Real dt = 0.0)
    {
        if (material_id_[index_i] == 0)
        {
            Real x = abs(pos0_[index_i][0] - (cx + fish_length));
            Real y = pos0_[index_i][1];

            Real Am = 0.12;
            Real frequency = 4.0;
            Real w = 2 * Pi * frequency;
            Real lambda = 3.0 * fish_length;
            Real wave_number = 2 * Pi / lambda;
            Real hx = -(pow(x, 2) - pow(fish_length, 2)) / pow(fish_length, 2);
            Real start_time = 0.2;
            Real current_time = GlobalStaticVariables::physical_time_;
            Real strength = 1 - exp(-current_time / start_time);

            Real phase_shift = y > (cy + bone_thickness / 2) ? 0 : Pi / 2;
            active_strain_[index_i](0, 0) =
                -Am * hx * strength * pow(sin(w * current_time / 2 + wave_number * x / 2 + phase_shift), 2);
        }
    };

  protected:
    StdLargeVec<int> &material_id_;
    StdLargeVec<Vecd> &pos0_;
    StdLargeVec<Matd> active_strain_;
};
