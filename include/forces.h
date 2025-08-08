#include <stddef.h>

#include "units.h"
#include "planets.h"
#include "vectors.h"
#include "constants.h"

static inline __attribute__((always_inline, const))
three_acceleration calculate_gravitational_acceleration(
    const three_position                  position,
    const mass_t                *restrict attractor_masses,
    const attractor_positions_t           attractor_positions)   /* ← pointer + restrict */
{
    three_acceleration acceleration = THREE_ZERO;

    #pragma clang loop vectorize(enable)
    #pragma clang loop unroll(disable)
    for (size_t attractor_index = 0; attractor_index < NUM_ATTRACTORS; ++attractor_index) {

        three_position displacement_vector = three_position_subtract(
            (three_position) EXTRACT_VEC(attractor_positions, attractor_index),
            position
        );

        length_t distance = three_vector_length(displacement_vector);
        if (distance.value == 0) continue;

        inv_length_t inv_distance = inv_length(distance); // 1/r
        inv_length_cubed_t inv_distance_cubed = inv_length_cubed(inv_distance); // 1/r^3

        double scale = convert_units.G_code  * attractor_masses[attractor_index].value * inv_distance_cubed.value;
        
        acceleration.x[0].value += scale * displacement_vector.x[0].value;
        acceleration.y[0].value += scale * displacement_vector.y[0].value;
        acceleration.z[0].value += scale * displacement_vector.z[0].value;
    }

    return acceleration;
}
