#define NUM_PLANETS 10
#define MAX_NAME_LEN 32

typedef struct {
    length_t x_position[NUM_PLANETS];
    length_t y_position[NUM_PLANETS];
    length_t z_position[NUM_PLANETS];
    velocity_t x_velocity[NUM_PLANETS];
    velocity_t y_velocity[NUM_PLANETS];
    velocity_t z_velocity[NUM_PLANETS];
    mass_t mass[NUM_PLANETS];
} planets_t;

planets_t solar_system = {
    .planets = {
        {
            .name = "Sun",
            .mass = {1.9885e30},
            .radius = {695700000},
            .position = {{0}, {0}, {0}},
            .velocity = {{0}, {-14.2451}, {-0.3044}},
            .momentum = {{0}, {-2.832638941798e31}, {-6.0538262806e29}}
        },
        {
            .name = "Mercury",
            .mass = {3.3011e23},
            .radius = {2439700},
            .position = {{4.6e10}, {0}, {0}},
            .velocity = {{0}, {58540}, {7188}},
            .momentum = {{0}, {1.93246394e28}, {2.37283068e27}}
        },
        {
            .name = "Venus",
            .mass = {4.8675e24},
            .radius = {6051800},
            .position = {{1.0748e11}, {0}, {0}},
            .velocity = {{0}, {35200}, {2085}},
            .momentum = {{0}, {1.71336e29}, {1.01487375e28}}
        },
        {
            .name = "Earth",
            .mass = {5.9724e24},
            .radius = {6378100},
            .position = {{1.4709e11}, {0}, {0}},
            .velocity = {{0}, {30290}, {0}},
            .momentum = {{0}, {1.80903996e29}, {0}}
        },
        {
            .name = "Moon",
            .mass = {7.346e22},
            .radius = {1738100},
            .position = {{1.47366e11}, {0}, {0}},
            .velocity = {{0}, {30362}, {1076}},
            .momentum = {{0}, {2.23039252e27}, {7.904296e25}}
        },
        {
            .name = "Mars",
            .mass = {6.4171e23},
            .radius = {3396200},
            .position = {{2.0662e11}, {0}, {0}},
            .velocity = {{0}, {26386}, {852}},
            .momentum = {{0}, {1.693216006e28}, {5.4673692e26}}
        },
        {
            .name = "Jupiter",
            .mass = {1.89819e27},
            .radius = {66854000},
            .position = {{7.4052e11}, {0}, {0}},
            .velocity = {{0}, {14717}, {312}},
            .momentum = {{0}, {2.793566223e31}, {5.9223528e29}}
        }
    }
};