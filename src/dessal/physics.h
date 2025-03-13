#ifndef PHYSICS

#define PHYSICS

#include "../properties/properties.h"

typedef struct
{
    PetscReal membrane_gap_pressure, mass_flux;
} WaterProduction;

// Function to calculate the mass transfer coefficient in the feed water channel
PetscReal ChannelMassTransfCoef(SaltWaterProperties *bulk_water_prop,
                                SaltWaterProperties *wall_water_prop,
                                PetscReal mass_flow_rate,
                                PetscReal channel_height,
                                PetscReal channel_width,
                                PetscReal membrane_area,
                                PetscReal number_channels,
                                PetscReal spacer_porosity);

// Function to calculate the heat transfer coefficients in the water channels
PetscReal ChannelHeatTransfCoef(SaltWaterProperties *bulk_water_prop,
                                PetscReal mass_flow_rate,
                                PetscReal channel_height,
                                PetscReal channel_width,
                                PetscReal membrane_area,
                                PetscInt number_channels,
                                PetscReal spacer_porosity);

// Function to calculate the salinity at the feed-membrane interface
PetscReal MembraneSalinity(SaltWaterProperties *bulk_water_prop,
                           PetscReal bulk_salinity,
                           PetscReal mass_transfer_coef,
                           PetscReal mass_flux);

// Function to calculate the effective thermal conductivity of the membrane
PetscReal MembraneConductivity(MoistAirProperties *pore_air_prop,
                               PetscReal polymer_conductivity,
                               PetscReal membrane_porosity);

// Function to calculate the distillate mass flux across the membrane
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
                         PetscReal film_thickness);

#endif