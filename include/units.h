
#ifndef UNITS
#define UNITS

typedef struct {
    double value;
} time_t;

typedef struct {
    double value;
} time_squared_t;

typedef struct {
    double value;
} velocity_t;

typedef struct {
    double value;
} length_t;

typedef struct {
    double value;
} inv_length_t;

typedef struct {
    double value;
} inv_length_cubed_t;

typedef struct {
    double value;
} length_squared_t;

typedef struct {
    double value;
} mass_t;

typedef struct {
    double value;
} momentum_t;

typedef struct {
    double value;
} acceleration_t;

static inline time_squared_t square_time(time_t time_a, time_t time_b) {
    return (time_squared_t) {time_a.value * time_b.value};
}

static inline velocity_t calculate_velocity_increase(acceleration_t acceleration, time_t acceleration_time) {
    return (velocity_t) {acceleration.value * acceleration_time.value};
}

static inline length_t
length_sum_impl(size_t count, const length_t *values)
{
    length_t total = { 0.0 };
    for (size_t i = 0; i < count; ++i)
        total.value += values[i].value;
    return total;
}

#define LENGTH_SUM(...)                                                     \
    length_sum_impl(                                                        \
        sizeof((length_t[]){ __VA_ARGS__ }) / sizeof(length_t),             \
        (length_t[]){ __VA_ARGS__ } )
    
#define add_length(...)  LENGTH_SUM(__VA_ARGS__)

static inline length_t subtract_length(length_t length_a, length_t length_b) {
    return (length_t) {length_a.value - length_b.value};
}

static inline length_t mult_length(length_t length_a, length_t length_b) {
    return (length_t) {length_a.value * length_b.value};
}

static inline inv_length_t inv_length(length_t length) {
    return (inv_length_t) {1.0 / length.value};
}

static inline inv_length_cubed_t inv_length_cubed(inv_length_t inv_length) {
    return (inv_length_cubed_t) {inv_length.value*inv_length.value*inv_length.value};
}

static inline length_squared_t square_length(length_t length) {
    return (length_squared_t) {length.value * length.value};
}

static inline length_t calculate_distance_travelled(velocity_t velocity, time_t travel_time) {
    return (length_t) {travel_time.value * velocity.value};
}

static inline acceleration_t add_acceleration(acceleration_t acceleration_a, acceleration_t acceleration_b) {
    return (acceleration_t) {acceleration_a.value + acceleration_b.value};
}


static inline length_squared_t
length_squared_sum_impl(size_t count, const length_squared_t *values)
{
    length_squared_t total = { 0.0 };
    for (size_t i = 0; i < count; ++i)
        total.value += values[i].value;
    return total;
}

#define LENGTH_SQUARED_SUM(...)                                                     \
    length_squared_sum_impl(                                                        \
        sizeof((length_squared_t[]){ __VA_ARGS__ }) / sizeof(length_squared_t),     \
        (length_squared_t[]){ __VA_ARGS__ } )
    
#define add_length_squared(...)  LENGTH_SQUARED_SUM(__VA_ARGS__)

#endif
