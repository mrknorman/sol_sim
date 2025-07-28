#include "units.h"

#ifndef VECTORS
#define VECTORS

#define length_v(N) struct { length_t x[N]; length_t y[N]; length_t z[N]; }

typedef length_v(1) three_position;

#define THREE_ZERO {.x=0.0,.y=0.0,.z=0.0}

#define velocity_v(N) struct { velocity_t x[N]; velocity_t y[N]; velocity_t z[N]; }
#define momentum_v(N) struct { momentum_t x[N]; momentum_t y[N]; momentum_t z[N]; }
#define accelleration_v(N) struct { acceleration_t x[N]; acceleration_t y[N]; acceleration_t z[N]; }

typedef accelleration_v(1) three_acceleration;

#define kinematics_t(X) \
    struct { \
        length_v(X) position; \
        velocity_v(X) velocity; \
    }

#endif