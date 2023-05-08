/**
* @author 	Yaru Ren, Chi Zhang and Xiangyu Hu
*/

#include "sphinxsys.h"
#include "2d_fish_and_bones.h"
#define PI 3.1415926
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;								     /**< Channel length. */
Real DH = 0.5;								    /**< Channel height. */
Real particle_spacing_ref = 0.005;			   /**< Initial reference particle spacing. */
Real DL_sponge = particle_spacing_ref * 20.0; /**< Sponge region to impose inflow condition. */
Real BW = particle_spacing_ref * 4.0;		 /**< Boundary width, determined by specific layer of boundary particles. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -0.25 * DH), Vec2d(DL, 1.25 * DH));

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
Real rho0_f = 1.0;										     /**< Density. */
Real U_f = 1.0;												    /**< freestream velocity. */
Real c_f = 10.0 * U_f;										   /**< Speed of sound. */
Real Re = 30000.0;											  /**< Reynolds number. */
Real mu_f = rho0_f * U_f * 0.03/ Re;                  /**< Dynamics viscosity. */

//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
Real cx = 0.5 * DL;				   /**< Center of fish in x direction. */
Real cy = DH / 2;					  /**< Center of fish in y direction. */
Real fish_length = 0.2;			     /**< Length of fish. */
Real fish_thickness = 0.03;		    /**< The maximum fish thickness. */
Real muscel_thickness = 0.02;	   /**< The maximum fish thickness. */
Real head_length = 0.03;		  /**< Length of fish bone. */
Real bone_thickness = 0.003;	 /**< Length of fish bone. */
Real fish_shape_resolution = particle_spacing_ref * 0.5;

Real rho0_s = 1.05;
Real Ae = 1.27e4; /**< Normalized Youngs Modulus. */
Real a01[4] = { 2.68e2, 0.0, 0.0, 0.0 };
Real b0[4] = { 8.023, 0.0, 0.0, 0.0 };
Vecd fiber_direction(1.0, 0.0);
Vecd sheet_direction(0.0, 1.0);
Real bulk_modulus1 = Ae * rho0_s * 1.0 * 1.0;
Real Youngs_modulus = Youngs_modulus = Ae * rho0_s * 1.0 * 1.0 * 0.08; // 1.1e6;
Real poisson = 0.49;

Real a1 = 1.22 * fish_thickness / fish_length;
Real a2 = 3.19 * fish_thickness / fish_length / fish_length;
Real a3 = -15.73 * fish_thickness / pow(fish_length, 3);
Real a4 = 21.87 * fish_thickness / pow(fish_length, 4);
Real a5 = -10.55 * fish_thickness / pow(fish_length, 5);

Real b1 = 1.22 * muscel_thickness / fish_length;
Real b2 = 3.19 * muscel_thickness / fish_length / fish_length;
Real b3 = -15.73 * muscel_thickness / pow(fish_length, 3);
Real b4 = 21.87 * muscel_thickness / pow(fish_length, 4);
Real b5 = -10.55 * muscel_thickness / pow(fish_length, 5);

// Observation locations
Vec2d point_coordinate_1(cx + 0.2 * fish_length, 0.263);
Vec2d point_coordinate_2(cx + 0.4 * fish_length, 0.264);
Vec2d point_coordinate_3(cx + 0.6 * fish_length, 0.262);
Vec2d point_coordinate_4(cx + 0.8 * fish_length, 0.257);

Vec2d point_coordinate_5(cx + 0.2 * fish_length, 0.237);
Vec2d point_coordinate_6(cx + 0.4 * fish_length, 0.236);
Vec2d point_coordinate_7(cx + 0.6 * fish_length, 0.238);
Vec2d point_coordinate_8(cx + 0.8 * fish_length, 0.243);

StdVec<Vecd> observation_locations = { point_coordinate_1, point_coordinate_2, point_coordinate_3,
									   point_coordinate_4, point_coordinate_5, point_coordinate_6,
                                       point_coordinate_7, point_coordinate_8};

//----------------------------------------------------------------------
//	SPH bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
	water_block_shape.push_back(Vecd(-DL_sponge, DH));
	water_block_shape.push_back(Vecd(DL, DH));
	water_block_shape.push_back(Vecd(DL, 0.0));
	water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

	return water_block_shape;
}

/** create outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));

	return outer_wall_shape;
}
/**
* @brief create inner wall shape
*/
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, DH));
	inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, DH));
	inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));
	inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));

	return inner_wall_shape;
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
class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> fish_shape = CreatFishShape(cx, cy, fish_length, fish_shape_resolution);
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(fish_shape, ShapeBooleanOps::sub);
	}
};

/* Definition of the solid body. */
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		std::vector<Vecd> outer_shape = createOuterWallShape();
		std::vector<Vecd> inner_shape = createInnerWallShape();
		multi_polygon_.addAPolygon(outer_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(inner_shape, ShapeBooleanOps::sub);
	}
};

/** create a inflow buffer shape. */
MultiPolygon createInflowBufferShape()
{
	std::vector<Vecd> inflow_buffer_shape;
	inflow_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));
	inflow_buffer_shape.push_back(Vecd(-DL_sponge, DH));
	inflow_buffer_shape.push_back(Vecd(0.0, DH));
	inflow_buffer_shape.push_back(Vecd(0.0, 0.0));
	inflow_buffer_shape.push_back(Vecd(-DL_sponge, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(inflow_buffer_shape, ShapeBooleanOps::add);
	return multi_polygon;
}

/** Case dependent inflow boundary condition. */
struct FreeStreamVelocity
{
	Real u_ref_, t_ref_;

	template <class BoundaryConditionType>
	FreeStreamVelocity(BoundaryConditionType& boundary_condition)
		: u_ref_(0.02), t_ref_(2.0) {}

	Vecd operator()(Vecd& position, Vecd& velocity)
	{
		Vecd target_velocity = Vecd::Zero();
		Real run_time = GlobalStaticVariables::physical_time_;
		target_velocity[0] = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
		return target_velocity;
	}
};

//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class TimeDependentAcceleration : public Gravity
{
	Real t_ref_, u_ref_, du_ave_dt_;

public:
	explicit TimeDependentAcceleration(Vecd gravity_vector)
		: Gravity(gravity_vector), t_ref_(2.0), u_ref_(0.02), du_ave_dt_(0) {}

	virtual Vecd InducedAcceleration(Vecd& position) override
	{
		Real run_time_ = GlobalStaticVariables::physical_time_;
		du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);

		return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
	}
};

/**
 * Assign case dependent muscle activation histroy
 */
class MyocardiumActivation
	: public active_muscle_dynamics::MuscleActivation
{
public:
	explicit MyocardiumActivation(SolidBody &myocardium)
		: active_muscle_dynamics::MuscleActivation(myocardium) {};

	void update(size_t index_i, Real dt)
	{
        size_t n = pos0_.size();
        Real x = abs(pos0_[index_i][0] - pos0_[n - 1][0]);
        Real x1 = pos0_[index_i][0] - pos0_[0][0];
        Real y = pos0_[index_i][1];
		Real y1(0);
		Real y2(0);

		y1 = a1 * pow(x, 0 + 1) + a2 * pow(x, 1 + 1) + a3 * pow(x, 2 + 1) + a4 * pow(x, 3 + 1) + a5 * pow(x, 4 + 1);
		y2 = b1 * pow(x, 0 + 1) + b2 * pow(x, 1 + 1) + b3 * pow(x, 2 + 1) + b4 * pow(x, 3 + 1) + b5 * pow(x, 4 + 1);

		Real xx = abs(pos0_[index_i][0] - (cx+ fish_length));
		Real am = 300;
		Real frequency = 4.0;
		Real w = 2 * PI * frequency;
		Real lamda = fish_length;
		Real wave_number = 2 * PI / lamda;
		Real hx = (pow(xx, 2) - pow(0.2, 2)) / pow(0.2, 2);
		Real ta = 0.2;
		Real st = 1 - exp(-GlobalStaticVariables::physical_time_ / ta);

		if (x1 > (head_length) && y >(y2 + cy))
		{
			active_contraction_stress_[index_i] = am * hx * st * pow(sin(w * GlobalStaticVariables::physical_time_/2+ wave_number * xx/2 ),2);
		}

		if (x1 > (head_length) && y < (-y2 + cy))
		{
			active_contraction_stress_[index_i] = am * hx * st * pow(sin(w * GlobalStaticVariables::physical_time_/2 + wave_number * xx /2 + PI/ 2),2);
		}
	};
};