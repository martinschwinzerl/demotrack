#ifndef DEMOTRACK_PARTICLE_H__
#define DEMOTRACK_PARTICLE_H__

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

#if !defined( PARTICLE_ARGPTR_DEC )
    #define PARTICLE_ARGPTR_DEC __global 
#endif /* !defined( PARTICLE_ARGPTR_DEC ) */

double Particle_energy0( PARTICLE_ARGPTR_DEC const Particle *const restrict p );

void Particle_add_to_energy( PARTICLE_ARGPTR_DEC Particle* restrict p,
    double const delta_energy );


/* ************************************************************************* */

double Particle_energy0( PARTICLE_ARGPTR_DEC const Particle *const restrict p )
{
    return sqrt( p->p0c * p->p0c + p->mass0 * p->mass0 );
}

void Particle_add_to_energy( PARTICLE_ARGPTR_DEC Particle* restrict p,
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

#endif /* DEMOTRACK_PARTICLE_H__ */
