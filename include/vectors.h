#include "units.h"
#include <math.h>

#ifndef VECTORS
#define VECTORS

#define THREE_ZERO { .x = { {0} }, .y = { {0} }, .z = { {0} } }


#define length_v(N) struct { length_t x[N]; length_t y[N]; length_t z[N]; }
typedef length_v(1) three_position;


#define velocity_v(N) struct { velocity_t x[N]; velocity_t y[N]; velocity_t z[N]; }
typedef velocity_v(1) three_velocity;

#define momentum_v(N) struct { momentum_t x[N]; momentum_t y[N]; momentum_t z[N]; }
#define acceleration_v(N) struct { acceleration_t x[N]; acceleration_t y[N]; acceleration_t z[N]; }

typedef acceleration_v(1) three_acceleration;

#define INSERT_VEC(DESTINATION, INDEX, SOURCE)        \
    do {                                                          \
        (DESTINATION).x[(INDEX)] = (SOURCE).x[0];                             \
        (DESTINATION).y[(INDEX)] = (SOURCE).y[0];                             \
        (DESTINATION).z[(INDEX)] = (SOURCE).z[0];                             \
    } while (0)

#define COPY_VEC(DESTINATION, INDEX, SOURCE)        \
    do {                                                          \
        (DESTINATION).x[(INDEX)] = (SOURCE).x[(INDEX)];                             \
        (DESTINATION).y[(INDEX)] = (SOURCE).y[(INDEX)];                             \
        (DESTINATION).z[(INDEX)] = (SOURCE).z[(INDEX)];                             \
    } while (0)

#define EXTRACT_VEC(SOURCE, INDEX)        \
    {                                                          \
        .x = (SOURCE).x[(INDEX)],                             \
        .y = (SOURCE).y[(INDEX)],                            \
        .z = (SOURCE).z[(INDEX)]                             \
    }

/* Add the (x,y,z)[0] components of SRC to DEST[IDX] ------------------- */
#define ADD_LENGTH_VEC(DEST, IDX, SRC)             \
    do {                                              \
        (DEST).x[(IDX)].value += (SRC).x[0].value;  \
        (DEST).y[(IDX)].value += (SRC).y[0].value;  \
        (DEST).z[(IDX)].value += (SRC).z[0].value;  \
    } while (0)

/* Add the (x,y,z)[0] components of SRC to DEST[IDX] ------------------- */
#define ADD_VELOCITY_VEC(DEST, IDX, SRC)             \
    do {                                              \
        (DEST).x[(IDX)].value += (SRC).x[0].value;  \
        (DEST).y[(IDX)].value += (SRC).y[0].value;  \
        (DEST).z[(IDX)].value += (SRC).z[0].value;  \
    } while (0)

static inline length_t three_vector_length(three_position p)
{
    const double dx = p.x[0].value;
    const double dy = p.y[0].value;
    const double dz = p.z[0].value;
    return (length_t){ sqrt(dx*dx + dy*dy + dz*dz) };
}

static inline three_position three_position_subtract(three_position position_a, three_position position_b) {
        return (three_position) {
        .x = { subtract_length(position_a.x[0], position_b.x[0]) },
        .y = { subtract_length(position_a.y[0], position_b.y[0]) },
        .z = { subtract_length(position_a.z[0], position_b.z[0]) }
    };
}

#endif
