#include <ancse/fvm_rate_of_change.hpp>

#include <Eigen/Dense>
#include <ancse/numerical_flux.hpp>
#include <ancse/reconstruction.hpp>
#include <fmt/format.h>
#include <ancse/limiters.hpp>

#define REGISTER_NUMERICAL_FLUX(token, FluxType, flux)                         \
    if (config["flux"] == (token)) {                                           \
        return std::make_shared<FVMRateOfChange<FluxType, Reconstruction>>(    \
            grid, model, flux, reconstruction);                                       \
    }

template <class Reconstruction>
std::shared_ptr<RateOfChange>
deduce_numerical_flux(const nlohmann::json &config,
                      const Grid &grid,
                      const std::shared_ptr<Model> &model,
                      const std::shared_ptr<SimulationTime> &simulation_time,
                      const Reconstruction &reconstruction)
{
    REGISTER_NUMERICAL_FLUX("central_flux", CentralFlux, CentralFlux(model))
    REGISTER_NUMERICAL_FLUX("laxfriedrichs",LaxFriedrichs, LaxFriedrichs(grid,model,simulation_time))
    REGISTER_NUMERICAL_FLUX("rusanov", Rusanov, Rusanov(model))
    REGISTER_NUMERICAL_FLUX("roe", Roe, Roe(model))
    REGISTER_NUMERICAL_FLUX("HLL", HLL, HLL(model))

    // Register the other numerical fluxes.

    if (model->get_name().compare("euler") == 0) {
        auto model_euler = std::dynamic_pointer_cast<Euler>(model);

        REGISTER_NUMERICAL_FLUX("hllc", HLLCEuler, HLLCEuler(model_euler))
    }


    throw std::runtime_error(
        fmt::format("Unknown numerical flux. {}", std::string(config["flux"])));
}
#undef REGISTER_NUMERICAL_FLUX

#define REGISTER_RECONSTRUCTION(token, reconstruction)                         \
    if (config["reconstruction"] == token) {                                   \
        return deduce_numerical_flux(                                          \
            config, grid, model, simulation_time, reconstruction);                     \
    }

std::shared_ptr<RateOfChange> make_fvm_rate_of_change(
    const nlohmann::json &config,
    const Grid &grid,
    const std::shared_ptr<Model> &model,
    const std::shared_ptr<SimulationTime> &simulation_time)
{
    REGISTER_RECONSTRUCTION("o1", PWConstantReconstruction{})
    REGISTER_RECONSTRUCTION("minmod-c", (PWLinearReconstruction<MinMod,Conserved>(MinMod{})))
    REGISTER_RECONSTRUCTION("minmod-p", (PWLinearReconstruction<MinMod,Primitive>(model, MinMod{})))
    REGISTER_RECONSTRUCTION("superbee-c", (PWLinearReconstruction<SuperBee,Conserved>( SuperBee{})))
    REGISTER_RECONSTRUCTION("superbee-p", (PWLinearReconstruction<SuperBee,Primitive>(model,SuperBee{})))
    REGISTER_RECONSTRUCTION("MC-c", (PWLinearReconstruction<MonotonizedCentral,Conserved>( MonotonizedCentral{})))
    REGISTER_RECONSTRUCTION("MC-p", (PWLinearReconstruction<MonotonizedCentral,Primitive>(model,MonotonizedCentral{})))



    // Register piecewise linear reconstructions.

    throw std::runtime_error(fmt::format(
        "Unknown reconstruction. [{}]", std::string(config["reconstruction"])));
}

#undef REGISTER_RECONSTRUCTION
