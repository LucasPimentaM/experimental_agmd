#include "../entrydata/entrydata.h"
#include "physics.h"

/*
Mass and energy balance in the desalination module
*/

PetscErrorCode DessalBalance(DessalData *dessal_data)
{
    PetscFunctionBeginUser;

    // Operational, properties, and geometrical data
    PetscReal feed_mass_flow_rate = dessal_data->feed_mass_flow_rate,
              cool_mass_flow_rate = dessal_data->cool_mass_flow_rate,
              entry_temperature_feed = dessal_data->entry_temperature_feed,
              entry_temperature_cool = dessal_data->entry_temperature_cool,
              entry_salinity_feed = dessal_data->entry_salinity_feed,
              entry_salinity_cool = dessal_data->entry_salinity_cool,
              vacuum_pressure = dessal_data->vacuum_pressure,
              membrane_area = dessal_data->membrane_area,
              membrane_thickness = dessal_data->membrane_thickness,
              membrane_porosity = dessal_data->membrane_porosity,
              pore_diameter = dessal_data->pore_diameter,
              feed_channel_height = dessal_data->feed_channel_height,
              cool_channel_height = dessal_data->cool_channel_height,
              channel_width = dessal_data->channel_width,
              spacer_porosity = dessal_data->spacer_porosity,
              gap_spacer_porosity = dessal_data->gap_spacer_porosity,
              air_gap_thickness = dessal_data->air_gap_thickness,
              wall_thickness = dessal_data->wall_thickness,
              polymer_conductivity = dessal_data->polymer_conductivity,
              spacer_conductivity = dessal_data->spacer_conductivity,
              wall_conductivity = dessal_data->wall_conductivity;
    PetscInt number_channels = dessal_data->number_channels;

    // Iterative data
    PetscReal out_temperature_feed = dessal_data->out_temperature_feed,
              out_temperature_cool = dessal_data->out_temperature_cool,
              feed_membrane_temperature = dessal_data->feed_membrane_temperature,
              gap_membrane_temperature = dessal_data->gap_membrane_temperature,
              film_boundary_temperature = dessal_data->film_boundary_temperature,
              film_wall_temperature = dessal_data->film_wall_temperature,
              cool_wall_temperature = dessal_data->cool_wall_temperature,
              membrane_salinity = dessal_data->membrane_salinity,
              out_salinity_feed = dessal_data->out_salinity_feed,
              membrane_gap_pressure = dessal_data->membrane_gap_pressure,
              mass_flux = dessal_data->mass_flux,
              film_thickness = dessal_data->film_thickness,
              heat_flux = dessal_data->heat_flux,
              vapor_heat_flux = dessal_data->vapor_heat_flux,
              feed_outflow_rate = dessal_data->feed_outflow_rate;

    // Feed
    SaltWaterProperties feed_prop, feed_memb_prop;
    PetscReal feed_resistance, feed_heat_transf_coef, feed_mass_trasf_coef;
    PetscReal avg_feed_temperature = 0.5 * (entry_temperature_feed + out_temperature_feed),
              avg_feed_salinity = PetscMax(0.5 * (entry_salinity_feed + out_salinity_feed), 0.0);

    SaltWaterPropBuild(&feed_prop, avg_feed_temperature, avg_feed_salinity);
    SaltWaterPropBuild(&feed_memb_prop, feed_membrane_temperature, membrane_salinity);

    feed_heat_transf_coef = ChannelHeatTransfCoef(&feed_prop,
                                                  feed_mass_flow_rate,
                                                  feed_channel_height,
                                                  channel_width,
                                                  membrane_area,
                                                  number_channels,
                                                  spacer_porosity);
    feed_resistance = 1.0 / feed_heat_transf_coef;

    feed_mass_trasf_coef = ChannelMassTransfCoef(&feed_prop,
                                                 &feed_memb_prop,
                                                 feed_mass_flow_rate,
                                                 feed_channel_height,
                                                 channel_width,
                                                 membrane_area,
                                                 number_channels,
                                                 spacer_porosity);

    // Updating the salinity at the feed-membrane interface
    membrane_salinity = MembraneSalinity(&feed_prop,
                                         avg_feed_salinity,
                                         feed_mass_trasf_coef,
                                         mass_flux);

    // Heat conduction in the membrane
    MoistAirProperties pore_air_prop;
    PetscReal membrane_resistance, membrane_conductivity;

    MoistAirPropBuild(&pore_air_prop, 0.5 * (feed_membrane_temperature + gap_membrane_temperature));

    membrane_conductivity = MembraneConductivity(&pore_air_prop, polymer_conductivity, membrane_porosity);
    membrane_resistance = membrane_thickness / membrane_conductivity;

    // Mass flux in the air gap
    SaltWaterProperties film_prop;
    WaterProduction water_production;
    PetscReal latent_resistance;

    SaltWaterPropBuild(&film_prop, film_boundary_temperature, 0.0);

    water_production = MassFlux(membrane_porosity,
                                membrane_tortuosity,
                                membrane_thickness,
                                pore_diameter,
                                gap_spacer_porosity,
                                air_gap_thickness,
                                0.5 * (feed_membrane_temperature + gap_membrane_temperature),
                                0.5 * (gap_membrane_temperature + film_boundary_temperature),
                                feed_memb_prop.vapor_pressure,
                                film_prop.vapor_pressure,
                                membrane_gap_pressure,
                                vacuum_pressure,
                                mass_flux,
                                film_thickness);

    mass_flux = water_production.mass_flux;
    membrane_gap_pressure = water_production.membrane_gap_pressure;

    mass_flux = PetscMax(mass_flux, 0.0);
    mass_flux = PetscMin(mass_flux, 0.2 * feed_mass_flow_rate / membrane_area);

    vapor_heat_flux = mass_flux * feed_memb_prop.latent_heat_vaporization;

    latent_resistance = (feed_membrane_temperature - film_boundary_temperature) / vapor_heat_flux;

    // Heat conduction in the distillate film
    MoistAirProperties gap_air_prop;
    PetscReal film_resistance;
    PetscReal effective_conductivity;

    MoistAirPropBuild(&gap_air_prop, 0.5 * (gap_membrane_temperature + film_boundary_temperature));

    film_thickness = film_prop.density * (film_prop.density - gap_air_prop.density) * acceleration_gravity;
    film_thickness = 3.0 * mass_flux * film_prop.dyn_viscosity * channel_width / film_thickness;
    film_thickness = 0.8 * PetscPowReal(film_thickness, 1.0 / 3.0);

    effective_conductivity = gap_spacer_porosity * film_prop.thermal_conductivity + (1.0 - gap_spacer_porosity) * spacer_conductivity;

    film_resistance = film_thickness / effective_conductivity;

    // Heat conduction in the air gap
    PetscReal gap_resistance;

    effective_conductivity = gap_spacer_porosity * gap_air_prop.thermal_conductivity + (1.0 - gap_spacer_porosity) * spacer_conductivity;

    gap_resistance = PetscMax(0.0, air_gap_thickness - film_thickness) / effective_conductivity;

    // Heat conduction in the wall
    PetscReal wall_resistance;

    wall_resistance = wall_thickness / wall_conductivity;

    // Coolant
    SaltWaterProperties cool_prop, cool_wall_prop;
    PetscReal cool_resistance, cool_heat_transf_coef;
    PetscReal avg_cool_temperature = 0.5 * (entry_temperature_cool + out_temperature_cool);

    SaltWaterPropBuild(&cool_prop, avg_cool_temperature, entry_salinity_cool);
    SaltWaterPropBuild(&cool_wall_prop, cool_wall_temperature, entry_salinity_cool);

    cool_heat_transf_coef = ChannelHeatTransfCoef(&cool_prop,
                                                  cool_mass_flow_rate,
                                                  cool_channel_height,
                                                  channel_width,
                                                  membrane_area,
                                                  number_channels,
                                                  spacer_porosity);
    cool_resistance = 1.0 / cool_heat_transf_coef;

    // Calculating the total heat flux
    PetscReal equiv_resistance;

    equiv_resistance = latent_resistance * (membrane_resistance + gap_resistance) / (latent_resistance + gap_resistance + membrane_resistance);
    equiv_resistance += feed_resistance + film_resistance + wall_resistance + cool_resistance;

    heat_flux = (avg_feed_temperature - avg_cool_temperature) / equiv_resistance;

    // Calculating the temperature at the interface between the feed and the membrane
    feed_membrane_temperature = avg_feed_temperature - heat_flux * feed_resistance;

    // Calculating the temperature at the interface between the feed and the air gap
    gap_membrane_temperature = feed_membrane_temperature - membrane_resistance * (heat_flux - vapor_heat_flux);

    // Calculating the temperature at the interface between the wall and the coolant
    cool_wall_temperature = avg_cool_temperature + heat_flux * cool_resistance;

    // Calculating the temperature at the interface between the wall and the distillate film
    film_wall_temperature = cool_wall_temperature + heat_flux * wall_resistance;

    // Calculating the temperature at the interface between the gap and the distillate film
    film_boundary_temperature = film_wall_temperature + heat_flux * film_resistance;

    // Calculating the mass outflow rate in the feed side
    feed_outflow_rate = feed_mass_flow_rate - mass_flux * membrane_area;

    // Calculating the salinity of the feed at the outlet
    out_salinity_feed = entry_salinity_feed * feed_mass_flow_rate / feed_outflow_rate;

    // Calculating the temperature of the feed at the outlet
    out_temperature_feed = entry_temperature_feed - heat_flux * membrane_area / (feed_mass_flow_rate * feed_prop.specific_heat);

    // Calculating the temperature of the coolant at the outlet
    out_temperature_cool = entry_temperature_cool + heat_flux * membrane_area / (cool_mass_flow_rate * cool_prop.specific_heat);

    // Updating the iterative data
    dessal_data->out_temperature_feed = out_temperature_feed;
    dessal_data->out_temperature_cool = out_temperature_cool;
    dessal_data->feed_membrane_temperature = feed_membrane_temperature;
    dessal_data->gap_membrane_temperature = gap_membrane_temperature;
    dessal_data->film_boundary_temperature = film_boundary_temperature;
    dessal_data->film_wall_temperature = film_wall_temperature;
    dessal_data->cool_wall_temperature = cool_wall_temperature;
    dessal_data->membrane_salinity = membrane_salinity;
    dessal_data->out_salinity_feed = out_salinity_feed;
    dessal_data->membrane_gap_pressure = membrane_gap_pressure;
    dessal_data->mass_flux = mass_flux;
    dessal_data->film_thickness = film_thickness;
    dessal_data->heat_flux = heat_flux;
    dessal_data->vapor_heat_flux = vapor_heat_flux;
    dessal_data->feed_outflow_rate = feed_outflow_rate;

    return 0;
}