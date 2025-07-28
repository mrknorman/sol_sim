#include <stddef.h>
#define NUM_FRAGMENTS 1000
#define NUM_PLANETS 10
#define MAX_NAME_LENGTH 32

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
    length_t x_position[NUM_FRAGMENTS];
    length_t y_position[NUM_FRAGMENTS];
    length_t z_position[NUM_FRAGMENTS];
    velocity_t x_velocity[NUM_FRAGMENTS];
    velocity_t y_velocity[NUM_FRAGMENTS];
} fragments_t;

typedef struct {
    char name[NUM_PLANETS][MAX_NAME_LENGTH];
    length_t x_position[NUM_PLANETS];
    length_t y_position[NUM_PLANETS];
    length_t z_position[NUM_PLANETS];
    velocity_t x_velocity[NUM_PLANETS];
    velocity_t y_velocity[NUM_PLANETS];
    velocity_t z_velocity[NUM_PLANETS];
    momentum_t x_momentum[NUM_PLANETS];
    momentum_t y_momentum[NUM_PLANETS];
    momentum_t z_momentum[NUM_PLANETS];
    length_t radius[NUM_PLANETS];
    mass_t mass[NUM_PLANETS];
} planets_t;

planets_t solar_system = {
    /* names (≤32 chars each) */
    .name = {
        "Sol", "Mercury", "Venus", "Earth", "Moon",
        "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"
    },

    /* --- Positions (m) --------------------------------------------------- */
    .x_position = {
        {0.0}, {4.600e10}, {1.0748e11}, {1.4709e11}, {1.47366e11},
        {2.0662e11}, {7.4052e11},
        {1.357554e12}, {2.732696e12}, {4.47105e12}
    },
    .y_position = { {0},{0},{0},{0},{0},{0},{0},{0},{0},{0} },
    .z_position = { {0},{0},{0},{0},{0},{0},{0},{0},{0},{0} },

    /* --- Velocities (m s‑¹) --------------------------------------------- */
    .x_velocity = { {0},{0},{0},{0},{0},{0},{0},{0},{0},{0} },
    .y_velocity = {
        {-14.245104}, {58540.0}, {35200.0}, {30290.0}, {30362.0},
        {26386.0}, {14717.0},
        {10140.0}, {7130.0}, {5470.0}
    },
    .z_velocity = {
        {-0.304442}, {7188.0}, {2085.0}, {0.0}, {1076.0},
        {852.0}, {312.0},
        {0.0}, {0.0}, {0.0}
    },

    /* --- Momenta (kg m s‑¹) --------------------------------------------- */
    .x_momentum = { {0},{0},{0},{0},{0},{0},{0},{0},{0},{0} },
    .y_momentum = {
        {-2.83263894e31}, {1.93246394e28}, {1.71336e29},
        {1.80903996e29}, {2.23039252e27}, {1.693216006e28},
        {2.793566223e31},
        {5.7627648e30}, {6.1896243e29}, {5.6017723e29}
    },
    .z_momentum = {
        {-6.05382628e29}, {2.37283068e27}, {1.01487375e28},
        {0.0}, {7.904296e25}, {5.4673692e26}, {5.9223528e29},
        {0.0}, {0.0}, {0.0}
    },

    /* --- Radii (m) ------------------------------------------------------- */
    .radius = {
        {6.957e8},  {2.4397e6}, {6.0518e6}, {6.3781e6}, {1.7381e6},
        {3.3962e6}, {6.6854e7}, {5.4364e7}, {2.4973e7}, {2.4341e7}
    },

    /* --- Masses (kg) ----------------------------------------------------- */
    .mass = {
        {1.9885e30}, {3.3011e23}, {4.8675e24}, {5.9724e24}, {7.346e22},
        {6.4171e23}, {1.89819e27}, {5.6832e26}, {8.6811e25}, {1.02409e26}
    }
};

int main(int argc, char** argv) {

    return 0; 
}

