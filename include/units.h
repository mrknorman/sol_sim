
#ifndef UNITS
#define UNITS

typedef struct {
    double meters;
} length_t;

typedef struct {
    double meters_per_second;
} velocity_t;

typedef struct {
    double kilograms;
} mass_t;

typedef struct {
    double newton_seconds;
} momentum_t;

typedef struct {
    double meters_per_second_squared;
} acceleration_t;

#endif