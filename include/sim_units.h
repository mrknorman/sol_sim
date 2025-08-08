#include "constants.h"

#ifndef SIM_UNITS_H
#define SIM_UNITS_H

// Astronomical constants
#define AU_METERS              1.495978707e11        // 1 AU      [m]
#define SOLAR_MASS_KG          1.98847e30            // 1 M☉      [kg]
#define SIDEREAL_YEAR_SECONDS  31558149.7635456      // 1 yr      [s]

// Conversion factors (direct constants for simplicity)
#define METERS_TO_AU           (1.0 / AU_METERS)
#define KG_TO_SOLAR_MASS       (1.0 / SOLAR_MASS_KG)
#define SECONDS_TO_YEAR        (1.0 / SIDEREAL_YEAR_SECONDS)

// Derived conversion factors
#define METERS_PER_SECOND_TO_AU_PER_YEAR     (METERS_TO_AU / SECONDS_TO_YEAR)
#define METERS_PER_SECOND_SQUARED_TO_AU_PER_YEAR_SQUARED  (METERS_TO_AU / (SECONDS_TO_YEAR * SECONDS_TO_YEAR))
#define KG_METERS_PER_SECOND_TO_SOLAR_MASS_AU_PER_YEAR   (KG_TO_SOLAR_MASS * METERS_TO_AU / SECONDS_TO_YEAR)

// Gravitational constant in code units
#define GRAVITATIONAL_CONSTANT_CODE   (GRAVITATIONAL_CONSTANT * SOLAR_MASS_KG * \
                                       SIDEREAL_YEAR_SECONDS * SIDEREAL_YEAR_SECONDS / \
                                       (AU_METERS * AU_METERS * AU_METERS))             // ≈ 4π²

// Optional sanity check at compile-time (C11)
_Static_assert(GRAVITATIONAL_CONSTANT_CODE > 39.47 && GRAVITATIONAL_CONSTANT_CODE < 39.49, "Gravitational constant code out of range");

static inline void convert_bodies_to_astronomical_units(
    bodies_t       *destination,
    const bodies_t *source
) {
    for (size_t i = 0; i < NUM_ATTRACTORS; ++i) {
        // positions ----------------------------------------------------
        destination->position.x[i].value = source->position.x[i].value * METERS_TO_AU;
        destination->position.y[i].value = source->position.y[i].value * METERS_TO_AU;
        destination->position.z[i].value = source->position.z[i].value * METERS_TO_AU;

        // velocities ---------------------------------------------------
        destination->velocity.x[i].value = source->velocity.x[i].value * METERS_PER_SECOND_TO_AU_PER_YEAR;
        destination->velocity.y[i].value = source->velocity.y[i].value * METERS_PER_SECOND_TO_AU_PER_YEAR;
        destination->velocity.z[i].value = source->velocity.z[i].value * METERS_PER_SECOND_TO_AU_PER_YEAR;

        // accelerations ------------------------------------------------
        destination->acceleration.x[i].value = source->acceleration.x[i].value * METERS_PER_SECOND_SQUARED_TO_AU_PER_YEAR_SQUARED;
        destination->acceleration.y[i].value = source->acceleration.y[i].value * METERS_PER_SECOND_SQUARED_TO_AU_PER_YEAR_SQUARED;
        destination->acceleration.z[i].value = source->acceleration.z[i].value * METERS_PER_SECOND_SQUARED_TO_AU_PER_YEAR_SQUARED;

        // momenta ------------------------------------------------------
        destination->momentum.x[i].value = source->momentum.x[i].value * KG_METERS_PER_SECOND_TO_SOLAR_MASS_AU_PER_YEAR;
        destination->momentum.y[i].value = source->momentum.y[i].value * KG_METERS_PER_SECOND_TO_SOLAR_MASS_AU_PER_YEAR;
        destination->momentum.z[i].value = source->momentum.z[i].value * KG_METERS_PER_SECOND_TO_SOLAR_MASS_AU_PER_YEAR;

        // radii & masses ----------------------------------------------
        destination->radius[i].value = source->radius[i].value * METERS_TO_AU;
        destination->mass[i].value = source->mass[i].value * KG_TO_SOLAR_MASS;
    }
}

#endif
