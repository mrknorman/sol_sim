#include <stddef.h>
#include <stdio.h>

#include "units.h"
#include "planets.h"
#include "vectors.h"
#include "forces.h"

#define NUM_FRAGMENTS 1000

typedef struct {
    kinematics_t(NUM_FRAGMENTS) kinematics;
} fragments_t;

int main(int argc, char** argv) {

    for (size_t i = 0; i < NUM_ATTRACTORS; i++) {
        three_acceleration test = gravitational_acceleration(
            (three_position) {
                .x = solar_system.position.x[i], 
                .y = solar_system.position.y[i], 
                .z = solar_system.position.z[i]
            },
            solar_system.mass,
            solar_system.position
        );
        printf("Gravitational Acceleration: %e, %e, %e \n", test.x[0].meters_per_second_squared, test.y[0].meters_per_second_squared, test.z[0].meters_per_second_squared);
    }

    return 0; 
}

