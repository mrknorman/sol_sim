#include "constants.h"

#ifndef SIM_UNITS_H
#define SIM_UNITS_H

#define AU_METERS              1.495978707e11        /* 1 AU      [m] */
#define MSUN_KG                1.98847e30            /* 1 M☉      [kg] */
#define SIDEREAL_YEAR_SECONDS  31558149.7635456      /* 1 yr      [s] */
#define G_CODE   (GRAVITATIONAL_CONSTANT * MSUN_KG * \
                  SIDEREAL_YEAR_SECONDS * SIDEREAL_YEAR_SECONDS / \
                  (AU_METERS * AU_METERS * AU_METERS))             /* ≈ 4π² */

typedef struct {
    double L_to_SI,  L_from_SI;
    double M_to_SI,  M_from_SI;
    double T_to_SI,  T_from_SI;
    double G_code;
} unit_system_t;

const unit_system_t convert_units = {
    .L_to_SI = AU_METERS,
    .L_from_SI = 1.0 / AU_METERS,

    .M_to_SI = MSUN_KG,
    .M_from_SI = 1.0 / MSUN_KG,

    .T_to_SI = SIDEREAL_YEAR_SECONDS,
    .T_from_SI = 1.0 / SIDEREAL_YEAR_SECONDS,

    .G_code  = G_CODE
};

/* Optional sanity check at compile-time (C11) */
_Static_assert(G_CODE > 39.47 && G_CODE < 39.49, "G_code out of range");

static inline void convert_bodies_to_sim_units(
	bodies_t       *dst,
	const bodies_t *src
) {
    const double L = convert_units.L_from_SI;                       /* 1 / AU */
    const double T = convert_units.T_from_SI;                       /* 1 / yr */
    const double M = convert_units.M_from_SI;                       /* 1 / M☉ */

	const double Vfac = L / T;            /* (m/s) → AU/yr   */
	const double Afac = L / (T * T);      /* m/s² → AU/yr²   */
	const double Pfac = M * L / T;        /* kg·m/s → M☉·AU/yr */
    const double Rfac   = L;                /* m → AU          */

    for (size_t i = 0; i < NUM_ATTRACTORS; ++i) {
        /* positions ---------------------------------------------------- */
        dst->position.x[i].value = src->position.x[i].value * L;
        dst->position.y[i].value = src->position.y[i].value * L;
        dst->position.z[i].value = src->position.z[i].value * L;

        /* velocities --------------------------------------------------- */
        dst->velocity.x[i].value = src->velocity.x[i].value * Vfac;
        dst->velocity.y[i].value = src->velocity.y[i].value * Vfac;
        dst->velocity.z[i].value = src->velocity.z[i].value * Vfac;

        /* accelerations ------------------------------------------------ */
        dst->acceleration.x[i].value = src->acceleration.x[i].value * Afac;
        dst->acceleration.y[i].value = src->acceleration.y[i].value * Afac;
        dst->acceleration.z[i].value = src->acceleration.z[i].value * Afac;

        /* momenta ------------------------------------------------------ */
        dst->momentum.x[i].value = src->momentum.x[i].value * Pfac;
        dst->momentum.y[i].value = src->momentum.y[i].value * Pfac;
        dst->momentum.z[i].value = src->momentum.z[i].value * Pfac;

        /* radii & masses ---------------------------------------------- */
        dst->radius[i].value = src->radius[i].value * Rfac;
        dst->mass[i].value   = src->mass[i].value * M;
    }
}



#endif
