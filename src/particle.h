#ifndef DEMOTRACK_PARTICLE_H__
#define DEMOTRACK_PARTICLE_H__

#if !defined( PARTICLE_ARGPTR_DEC )
    #define PARTICLE_ARGPTR_DEC __global
#endif /* !defined( PARTICLE_ARGPTR_DEC ) */

typedef struct Particle {
  double x;
  double px;
  double y;
  double py;
  double zeta;
  double delta;
  double rpp;
  double rvv;
  double psigma;
  double chi;
  double charge_ratio;
  double q0;
  double mass0;
  double beta0;
  double gamma0;
  double p0c;
  long state;
  long at_element;
  long at_turn;
  long id;
} Particle;

DEMOTRACK_STATIC DEMOTRACK_FN double Particle_energy0(
    PARTICLE_ARGPTR_DEC const Particle *const p );

DEMOTRACK_STATIC DEMOTRACK_FN void Particle_add_to_energy(
    PARTICLE_ARGPTR_DEC Particle* p, double const delta_energy );

DEMOTRACK_STATIC DEMOTRACK_FN void ParticleSet_setup(
    PARTICLE_ARGPTR_DEC Particle* pset,
    long int const num_particles, double const p0c, double const mass0 );

/* ************************************************************************* */

DEMOTRACK_INLINE double Particle_energy0(
    PARTICLE_ARGPTR_DEC const Particle *const p )
{
    return sqrt( p->p0c * p->p0c + p->mass0 * p->mass0 );
}

DEMOTRACK_INLINE void Particle_add_to_energy( PARTICLE_ARGPTR_DEC Particle* p,
    double const delta_energy )
{
    double const delta_beta0 = p->delta * p->beta0;
    double const ptau_beta0 = delta_energy / Particle_energy0( p ) +
        sqrt( delta_beta0 * delta_beta0 + ( double )2 * delta_beta0 * p->beta0
                + ( double )1 ) - ( double )1;

    double const ptau = ptau_beta0 / p->beta0;
    double const psigma = ptau / p->beta0;
    double const delta = sqrt( ptau * ptau + ( double )2 * psigma +
        ( double )1 ) - ( double )1;

    double const one_plus_delta = delta + ( double )1;
    double const rvv = one_plus_delta / ( ( double )1 + ptau_beta0 );

    p->delta = delta;
    p->psigma = psigma;
    p->zeta *= rvv / p->rvv;
    p->rvv = rvv;
    p->rpp = ( double )1 / one_plus_delta;
}

DEMOTRACK_INLINE void ParticleSet_setup( PARTICLE_ARGPTR_DEC Particle* pset,
    long int const num_particles, double const p0c, double const mass0 )
{
    #if defined( __cplusplus )
    using std::sqrt;
    #endif /* ADL for cmath methods */

    double const E0 = sqrt( p0c * p0c + mass0 * mass0 );
    double const beta0 = p0c / E0;
    double const gamma0 = E0 / mass0;

    for( long int idx = 0 ; idx < num_particles ; ++idx )
    {
        PARTICLE_ARGPTR_DEC Particle* p = &pset[ idx ];

        double const delta = 0.0;
        double ptau = sqrt( delta * delta + 2. * delta + 1./( beta0 * beta0 ) );
        ptau -= 1.0 / beta0;
        double const psigma = ptau / beta0;

        p->x = 0.0; p->y = 0.0; p->px = 0.0; p->py = 0.0;
        p->zeta = 0.0; p->delta = delta;
        p->psigma = psigma;

        p->chi = 1.0;
        p->charge_ratio = 1.0;
        p->rvv = 1.0;
        p->rpp = 1.0;

        p->q0 = 1.0;
        p->mass0 = mass0;
        p->beta0 = beta0;
        p->gamma0 = gamma0;
        p->p0c = p0c;

        p->state = 1;
        p->at_element = 0;
        p->at_turn = 0;
        p->id = idx;
    }
}

#endif /* DEMOTRACK_PARTICLE_H__ */
