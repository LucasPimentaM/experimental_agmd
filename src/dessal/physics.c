#include "../properties/properties.h"
#include "../entrydata/entrydata.h"
#include "physics.h"

/*
Maxwell's model for the thermal conductivity of the membrane and empirical correlation for the Nusselt number in spacer-filled channels

Reference: I. Hitsov, K. De Sitter, C. Dotremont, P. Cauwenberg, I. Nopens, Full-scale validated Air Gap Membrane Distillation (AGMD) model
           without calibration parameters. J. Membrane Sci. 533 (2017) 309-320. https://doi.org/10.1016/j.memsci.2017.04.002
*/

PetscReal MembraneConductivity(MoistAirProperties *pore_air_prop,
                               PetscReal polymer_conductivity,
                               PetscReal membrane_porosity)
{
    PetscReal air_conductivity = pore_air_prop->thermal_conductivity, beta;

    beta = (polymer_conductivity - air_conductivity) / (polymer_conductivity + 2.0 * air_conductivity);

    return 0.93 * air_conductivity * (1.0 + 2.0 * beta * (1.0 - membrane_porosity)) / (1.0 - beta * (1.0 - membrane_porosity));
}

PetscReal ChannelMassTransfCoef(SaltWaterProperties *bulk_water_prop,
                                SaltWaterProperties *wall_water_prop,
                                PetscReal mass_flow_rate,
                                PetscReal channel_height,
                                PetscReal channel_width,
                                PetscReal membrane_area,
                                PetscReal number_channels,
                                PetscReal spacer_porosity)
{
    // Properties
    PetscReal dyn_viscosity = bulk_water_prop->dyn_viscosity,
              density = bulk_water_prop->density,
              schmidt = dyn_viscosity / (density * salt_diffusivity);

    PetscReal channel_length, hydraulic_diameter, mass_velocity, reynolds, inverse_graetz, sherwood;

    channel_length = membrane_area / (number_channels * channel_width);

    hydraulic_diameter = 2.0 * channel_height;

    mass_velocity = mass_flow_rate / (number_channels * channel_height * channel_width);

    reynolds = mass_velocity * channel_height / dyn_viscosity;

    inverse_graetz = channel_length / (reynolds * schmidt * hydraulic_diameter);

    if (inverse_graetz <= 0.001)
    {
        sherwood = 2.236 * PetscPowReal(inverse_graetz, -1.0 / 3.0);
    }
    else if (inverse_graetz > 0.001 && inverse_graetz <= 0.01)
    {
        sherwood = 2.236 * PetscPowReal(inverse_graetz, -1.0 / 3.0) + 0.9;
    }
    else
    {
        sherwood = 8.235 + 0.0364 / inverse_graetz;
    }

    return salt_diffusivity * sherwood / channel_height;
}

PetscReal ChannelHeatTransfCoef(SaltWaterProperties *bulk_water_prop,
                                PetscReal mass_flow_rate,
                                PetscReal channel_height,
                                PetscReal channel_width,
                                PetscReal membrane_area,
                                PetscInt number_channels,
                                PetscReal spacer_porosity)
{
    // Properties
    PetscReal dyn_viscosity = bulk_water_prop->dyn_viscosity,
              thermal_conductivity = bulk_water_prop->thermal_conductivity,
              prandtl = bulk_water_prop->prandtl;

    PetscReal channel_length, hydraulic_diameter, mass_velocity, reynolds, inverse_graetz, nusselt;

    channel_length = membrane_area / (number_channels * channel_width);

    hydraulic_diameter = 2.0 * channel_height;

    mass_velocity = mass_flow_rate / (number_channels * channel_height * channel_width);

    reynolds = mass_velocity * hydraulic_diameter / dyn_viscosity;

    inverse_graetz = channel_length / (reynolds * prandtl * hydraulic_diameter);

    if (inverse_graetz <= 0.001)
    {
        nusselt = 2.236 * PetscPowReal(inverse_graetz, -1.0 / 3.0);
    }
    else if (inverse_graetz > 0.001 && inverse_graetz <= 0.01)
    {
        nusselt = 2.236 * PetscPowReal(inverse_graetz, -1.0 / 3.0) + 0.9;
    }
    else
    {
        nusselt = 8.235 + 0.0364 / inverse_graetz;
    }

    return thermal_conductivity * nusselt / channel_height;
}

PetscReal MembraneSalinity(SaltWaterProperties *bulk_water_prop,
                           PetscReal bulk_salinity,
                           PetscReal mass_transfer_coef,
                           PetscReal mass_flux)
{
    // Properties
    PetscReal density = bulk_water_prop->density;

    return bulk_salinity * PetscExpReal(mass_flux / (density * mass_transfer_coef));
}

/*
Water mass flux across the membrane using approximate membrane and air gap permeabilities

Reference: K.M. Lisboa, D.B. Moraes, C.P. Naveira-Cotta, R.M. Cotta, Analysis of the membrane effects on the energy efficiency of water
           desalination in a direct contact membrane distillation (DCMD) system with heat recovery. Appl. Thermal Eng. 182 (2021) 116063
           https://doi.org/10.1016/j.applthermaleng.2020.116063
*/

PetscReal MolecularDiffusion(PetscReal membrane_porosity, PetscReal membrane_tortuosity, PetscReal temperature)
{
    PetscReal molecular_diffusivity;

    molecular_diffusivity = 4.46e-6 * membrane_porosity / membrane_tortuosity;
    molecular_diffusivity *= PetscPowReal(temperature, 2.334);

    return molecular_diffusivity;
}

PetscReal KnudsenDiffusion(PetscReal membrane_porosity,
                           PetscReal membrane_tortuosity,
                           PetscReal pore_diameter,
                           PetscReal temperature)
{
    PetscReal knudsen_diffusivity;

    knudsen_diffusivity = pore_diameter / 3.0;
    knudsen_diffusivity *= membrane_porosity / membrane_tortuosity;
    knudsen_diffusivity *= PetscSqrtReal(8.0 * gas_constant * temperature / (M_PI * water_molar_mass));

    return knudsen_diffusivity;
}

WaterProduction MassFlux(PetscReal membrane_porosity,
                         PetscReal membrane_tortuosity,
                         PetscReal membrane_thickness,
                         PetscReal pore_diameter,
                         PetscReal gap_spacer_porosity,
                         PetscReal air_gap_thickness,
                         PetscReal temperature_membrane,
                         PetscReal temperature_gap,
                         PetscReal feed_membrane_pressure,
                         PetscReal film_boundary_pressure,
                         PetscReal membrane_gap_pressure,
                         PetscReal vacuum_pressure,
                         PetscReal mass_flux,
                         PetscReal film_thickness)
{
    WaterProduction water_production;
    PetscReal molecular_diffusivity, knudsen_diffusivity, effective_diffusivity,
              total_pressure = atm_pressure + vacuum_pressure;
    
    temperature_membrane += 273.15;
    temperature_gap += 273.15;

    molecular_diffusivity = MolecularDiffusion(1.0, 1.0, temperature_gap);

    membrane_gap_pressure = -mass_flux * gas_constant * temperature_gap * (air_gap_thickness - film_thickness) / (water_molar_mass * molecular_diffusivity);
    membrane_gap_pressure = total_pressure - (total_pressure - film_boundary_pressure) * PetscExpReal(membrane_gap_pressure);

    molecular_diffusivity = MolecularDiffusion(membrane_porosity, membrane_tortuosity, temperature_membrane);
    knudsen_diffusivity = KnudsenDiffusion(membrane_porosity, membrane_tortuosity, pore_diameter, temperature_membrane);

    effective_diffusivity = molecular_diffusivity * knudsen_diffusivity / (molecular_diffusivity + (atm_pressure + vacuum_pressure) * knudsen_diffusivity);

    mass_flux = (molecular_diffusivity - effective_diffusivity * membrane_gap_pressure) / (molecular_diffusivity - effective_diffusivity * feed_membrane_pressure);
    mass_flux = water_molar_mass * molecular_diffusivity * PetscLogReal(mass_flux) / (gas_constant * temperature_membrane * membrane_thickness);

    water_production.mass_flux = mass_flux;
    water_production.membrane_gap_pressure = membrane_gap_pressure;

    return water_production;
}