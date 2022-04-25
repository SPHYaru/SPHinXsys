/**
 * @file total_artificial_heart.cpp
 * @brief This is the example of total artificial heart implantation path simulation
 * @author John Benjamin, Bence Rochlitz - Virtonomy GmbH
 */
#include "bernoulli_beam.h"

#ifdef __EMSCRIPTEN__

#include <emscripten.h>
#include <emscripten/bind.h>

EMSCRIPTEN_BINDINGS(SPHINXSYS)
{
    emscripten::value_array<std::array<double, 3>>("ArrayDouble3")
        .element(emscripten::index<0>())
        .element(emscripten::index<1>())
        .element(emscripten::index<2>());

    emscripten::value_object<StlData>("StlData")   
        .field("name", &StlData::name)
        .field("ptr", &StlData::ptr);
    
    emscripten::register_vector<std::string>("StringVector");
    emscripten::register_vector<double>("DoubleVector");
    emscripten::register_vector<IndexVector::value_type>("UIntVector");
    emscripten::register_vector<StlData>("StlsList");
    emscripten::register_vector<IndexVector>("BodiesList");
    emscripten::register_map<std::string, std::string>("VtuData");

    emscripten::value_object<BernoulliBeamInput>("BernoulliBeamInput")
        .field("material_model_name", &BernoulliBeamInput::material_model_name)
        .field("scale_stl", &BernoulliBeamInput::scale_stl)
        .field("resolution", &BernoulliBeamInput::resolution)
        .field("rho_0", &BernoulliBeamInput::rho_0)
        .field("poisson", &BernoulliBeamInput::poisson)
        .field("Youngs_modulus", &BernoulliBeamInput::Youngs_modulus)
        .field("Youngs_modulus", &BernoulliBeamInput::Youngs_modulus)
        .field("physical_viscosity", &BernoulliBeamInput::physical_viscosity)
        .field("translation", &BernoulliBeamInput::translation)
        .field("stls", &BernoulliBeamInput::stls)
        .field("relative_input_path", &BernoulliBeamInput::relative_input_path)
        .field("contacting_bodies_list", &BernoulliBeamInput::contacting_bodies_list);

    emscripten::class_<BernoulliBeamJS>("SimTotalArtificialHeart")
        .constructor<BernoulliBeamInput>()
        .function("runSimulation", &BernoulliBeamJS::runSimulation)
        .function("onError", &BernoulliBeamJS::onError)
        .property("vtuData", &BernoulliBeamJS::getVtuData);
}

#else

int main()
{
    BernoulliBeamInput input;
    input.scale_stl = 0.001;
    input.resolution = { 1.5 };
    input.rho_0 = 1e3;
    input.poisson = 0.3;
    input.Youngs_modulus = 5e8;
    input.physical_viscosity = 5e6;
    input.translation = { 0, 0, 0 };

    string beam_stl = "bernoulli_beam_20x.stl";

    input.stls = { beam_stl };
    input.relative_input_path = "./input/";

    /* DOWNLOAD STLs files at this point */

    try
    {
        // set up the simulation
        BernoulliBeamJS simTotalArtificialHeart(input);
        int number_of_steps = 700;
        std::cout << "About to run the simulation" << std::endl;
        simTotalArtificialHeart.runSimulation(number_of_steps);
    }
    catch (const std::exception& e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}

#endif