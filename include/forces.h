#include <stddef.h>
#include <math.h>


#include "units.h"
#include "planets.h"
#include "vectors.h"
#include "constants.h"

three_acceleration gravitational_acceleration(
    const three_position position,
    const mass_t *attractor_masses,
    const attractor_positions_t attractor_positions
) {
    three_acceleration acceleration = THREE_ZERO;

    for (size_t i = 0; i < NUM_ATTRACTORS; i++) {
        // Δx = attractor.x - position.x
        const double dx = attractor_positions.x[i].meters - position.x[0].meters;
        const double dy = attractor_positions.y[i].meters - position.y[0].meters;
        const double dz = attractor_positions.z[i].meters - position.z[0].meters;

        const double r_squared = dx * dx + dy * dy + dz * dz;
        if (r_squared == 0.0) continue;

		const double r = sqrt(r_squared);

        // Acceleration magnitude: a = G * m / r^2
        const double a_mag = GRAVITAIONAL_CONSTANT * attractor_masses[i].kilograms / r_squared;

        // Normalize direction and scale
        acceleration.x[0].meters_per_second_squared += a_mag * dx / r;
        acceleration.y[0].meters_per_second_squared += a_mag * dy / r;
        acceleration.z[0].meters_per_second_squared += a_mag * dz / r;
    }

    return acceleration;
}