#include <stddef.h>
#include <stdio.h>
#include <omp.h> 
#include <math.h>

#include "units.h"
#include "planets.h"
#include "sim_units.h"
#include "vectors.h"
#include "forces.h"
#include "integrator.h"


#define SIMULATION_DURATION_YEARS 1
#define SIMULATION_DURATION_SECONDS (SIMULATION_DURATION_YEARS * SECONDS_PER_SIDEREAL_YEAR)
#define TIME_STEP_DURATION_SECONDS SECONDS_PER_DAY

#define SECONDS_PER_DAY 86400.002
#define DAYS_PER_SIDEREAL_YEAR 365.256363004
#define SECONDS_PER_SIDEREAL_YEAR SECONDS_PER_DAY*DAYS_PER_SIDEREAL_YEAR

#define NUM_FRAGMENTS 1000
#define NUM_TIME_STEPS (size_t) round(SIMULATION_DURATION_SECONDS/TIME_STEP_DURATION_SECONDS)

time_t TIME_STEP_DURATION;

int main(int argc, char** argv) {

    TIME_STEP_DURATION.value = SECONDS_PER_DAY * convert_units.T_from_SI;

    bodies_t solar_system;
    convert_bodies_to_sim_units(&solar_system, &initial_solar_system);

    attractor_accelerations_t current_acceleration;

    #pragma omp parallel default(none)            \
                    shared(solar_system,          \
                           TIME_STEP_DURATION,    \
                           current_acceleration)
    {

    FILE *csv_fp = NULL;
    #pragma omp master
    {
        csv_fp = fopen("positions.csv", "w");
        if (!csv_fp) { perror("positions.csv"); exit(1); }
        fprintf(csv_fp, "step,body,x,y,z\n");
    }

    for (size_t time_index = 0; time_index < NUM_TIME_STEPS; time_index++) {
        
            #pragma omp for simd schedule(static)
            for (size_t attractor_index = 0; attractor_index < NUM_ATTRACTORS; attractor_index++) {
                three_position additional_distance = calculate_distance_travelled_verlet(
                    (three_velocity) EXTRACT_VEC(solar_system.velocity, attractor_index),
                    (three_acceleration) EXTRACT_VEC(solar_system.acceleration, attractor_index),
                    TIME_STEP_DURATION
                );
                ADD_LENGTH_VEC(solar_system.position, attractor_index, additional_distance);
            }

            #pragma omp for simd schedule(static)
            for (size_t attractor_index = 0; attractor_index < NUM_ATTRACTORS; attractor_index++) {
                const three_acceleration acceleration = calculate_gravitational_acceleration(
                    (three_position) EXTRACT_VEC(solar_system.position, attractor_index),
                    solar_system.mass,
                    solar_system.position
                );
                INSERT_VEC(
                    current_acceleration, 
                    attractor_index,             
                    acceleration
                );
            }

            #pragma omp for simd schedule(static)
            for (size_t attractor_index = 0; attractor_index < NUM_ATTRACTORS; attractor_index++) {
                three_acceleration previous_acceleration = 
                    (three_acceleration) EXTRACT_VEC(solar_system.acceleration, attractor_index);

                three_velocity additional_velocity = calculate_velocity_verlet(
                    (three_acceleration) EXTRACT_VEC(solar_system.acceleration, attractor_index),
                    (three_acceleration) EXTRACT_VEC(current_acceleration, attractor_index),
                    TIME_STEP_DURATION
                );
                ADD_VELOCITY_VEC(solar_system.velocity, attractor_index, additional_velocity);

                COPY_VEC(solar_system.acceleration, attractor_index, current_acceleration);
            }

            #pragma omp barrier
            #pragma omp master
            {
                for (size_t i = 0; i < NUM_ATTRACTORS; ++i) {
                    fprintf(
                        csv_fp, "%zu,%zu,%.15e,%.15e,%.15e\n",
                        time_index, i,
                        solar_system.position.x[i].value,
                        solar_system.position.y[i].value,
                        solar_system.position.z[i].value
                    );
                }
            }
        }
    }

    return 0; 
}

