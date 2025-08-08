#ifndef INTEGRATOR
#define INTEGRATOR

static inline length_t calculate_acceleration_correction_verlet(acceleration_t acceleration, time_squared_t travel_time_squared) {
    return (length_t) {travel_time_squared.value * acceleration.value};
}

static inline three_position calculate_distance_travelled_verlet(
    three_velocity velocity, 
    three_acceleration acceleration,
    time_t duration
) {
    const time_squared_t half_duration_squared = (time_squared_t) {(duration.value * duration.value) / 2.0};

    return (three_position) {
        .x = {add_length(
            calculate_distance_travelled(velocity.x[0], duration),
            calculate_acceleration_correction_verlet(acceleration.x[0], half_duration_squared)
        )},
        .y = {add_length(
            calculate_distance_travelled(velocity.y[0], duration),
            calculate_acceleration_correction_verlet(acceleration.y[0], half_duration_squared)
        )},
        .z = {add_length(
            calculate_distance_travelled(velocity.z[0], duration),
            calculate_acceleration_correction_verlet(acceleration.z[0], half_duration_squared)
        )}
    };
}

static inline three_velocity calculate_velocity_verlet(
    three_acceleration previous_acceleration,
    three_acceleration current_acceleration,
    time_t duration
) {
    const time_t half_duration = (time_t) {duration.value / 2.0};

    // Update velocity
    return (three_velocity) {
        .x = calculate_velocity_increase(
            add_acceleration(previous_acceleration.x[0], current_acceleration.x[0]), half_duration),
        .y = calculate_velocity_increase(
            add_acceleration(previous_acceleration.y[0], current_acceleration.y[0]), half_duration),
        .z = calculate_velocity_increase(
            add_acceleration(previous_acceleration.z[0], current_acceleration.z[0]), half_duration)
    };
}

#endif
