#ifndef DEMOTRACK_TRACK_H__
#define DEMOTRACK_TRACK_H__

#include "beam_elements.h"
#include "space_charge.h"
#include "particle.h"

#if !defined( BE_ARGPTR_DEC )
    #define BE_ARGPTR_DEC __global 
#endif /* !defined( BE_ARGPTR_DEC ) */

#if !defined( PARTICLE_ARGPTR_DEC )
    #define PARTICLE_ARGPTR_DEC __global 
#endif /* !defined( PARTICLE_ARGPTR_DEC ) */

void Drift_track( BE_ARGPTR_DEC const Drift *const restrict drift,
    PARTICLE_ARGPTR_DEC Particle* restrict p );

void LimitGlobal_track( PARTICLE_ARGPTR_DEC Particle* restrict p );

void Multipole_track(
    BE_ARGPTR_DEC const Multipole *const restrict quad,
    PARTICLE_ARGPTR_DEC Particle* restrict p );

void Cavity_track(
    BE_ARGPTR_DEC const Cavity *const restrict cavity,
    PARTICLE_ARGPTR_DEC Particle* restrict p );

void SpaceChargeCoasting_track(
    BE_ARGPTR_DEC const SpaceChargeCoasting *const restrict cavity,
    PARTICLE_ARGPTR_DEC Particle* restrict p );

/* ------------------------------------------------------------------------- */

unsigned int Track_beam_element_dispatcher(
    PARTICLE_ARGPTR_DEC Particle* restrict p, 
    BE_ARGPTR_DEC double const* restrict beam_elements_buffer,
    unsigned int const current_slot_idx, 
    unsigned int const num_slots_in_beam_elements_buffer );

void Track_particle_until_turn(
    PARTICLE_ARGPTR_DEC Particle* restrict p,
    BE_ARGPTR_DEC double const* restrict beam_elements_buffer, 
    unsigned int const num_slots_in_beam_elements_buffer,
    long int const until_turn );

/* ************************************************************************* */

void Drift_track( BE_ARGPTR_DEC const Drift *const restrict drift,
    PARTICLE_ARGPTR_DEC Particle* restrict p )
{
    double const one_plus_delta = p->delta + ( double )1;
    double const lpzi = ( double )1 / sqrt( one_plus_delta * one_plus_delta 
        - ( p->px * p->px + p->py * p->py ) );

    p->x += p->px * drift->length * lpzi;
    p->y += p->py * drift->length * lpzi;    
    p->zeta += p->rvv * drift->length - one_plus_delta * lpzi;
}

void LimitGlobal_track( PARTICLE_ARGPTR_DEC Particle* restrict p ) 
{
    double const sign_x = ( double )( ( ( double )0 < p->x ) - 
                          ( double )( p->x < ( double )0 ) );
    
    double const sign_y = ( double )( ( ( double )0 < p->y ) - 
                          ( double )( p->y < ( double )0 ) );

    p->state &= ( long int )( ( ( sign_x * p->x ) < ( double )1 ) &
                              ( ( sign_y * p->y ) < ( double )1 ) );
}

void Multipole_track( BE_ARGPTR_DEC const Multipole *const restrict mp,
    PARTICLE_ARGPTR_DEC Particle* restrict p ) 
{
    long int index_x = 2 * ( long int )mp->order;
    long int index_y = index_x + 1;
    
    double d_px = mp->bal[ index_x ];
    double d_py = mp->bal[ index_y ];
    
    while( index_x > 0 )
    {
        index_x -= 2;
        index_y -= 2;

        d_px = d_px * p->x - d_py * p->y + mp->bal[ index_x ];
        d_py = d_px * p->y + d_py * p->x + mp->bal[ index_y ];
    }

    d_px = -p->chi * d_px;
    d_py =  p->chi * d_py;

    p->px += d_px;
    p->py += d_py;
}

void Cavity_track( BE_ARGPTR_DEC const Cavity *const restrict cavity,
    PARTICLE_ARGPTR_DEC Particle* restrict p )
{
    double const PI = ( double )3.141592653589793;    
    double const phase = ( PI / ( real_t )180 ) * cavity->lag -
        ( ( ( real_t )2 * PI ) / ( double )299792458.0 ) * 
        ( p->zeta / ( p->beta0 * p->rvv ) ) * cavity->frequency;
           
    Particle_add_to_energy( p, 
        cavity->voltage * p->q0 * p->charge_ratio * sin( phase ) );
}

void SpaceChargeCoasting_track(
    BE_ARGPTR_DEC const SpaceChargeCoasting *const restrict sc_elem,
    PARTICLE_ARGPTR_DEC Particle* restrict p )
{
    double Ex;
    double Ey;
    double Gx;
    double Gy;

    double const charge = ( double )1.602176634e-19 * p->q0;
    double const p0c = ( double )1.602176634e-19 * p->p0c;
    
    double fact_kick = sc_elem->number_of_particles * 
        sc_elem->length * p->chi * p->charge_ratio * charge * charge *
        ( ( double )1 - p->beta0 * beta0 );

    fact_kick /= sc_elem->circumference * p0c * p->beta0 * p->rvv;

    NS(get_Ex_Ey_Gx_Gy_gauss)( p->x - sc_elem->x_co, p->y - sc_elem->y_co,
        sc_elem->sigma_x, sc_elem->sigma_y, sc_elem->min_sigma_diff, 
        ( int )1, &Ex, &Ey, &Gx, &Gy );

    p->px += fact_kick * Ex;
    p->py += fact_kick * Ey;
}    

/* ------------------------------------------------------------------------- */

unsigned int Track_beam_element_dispatcher(
    PARTICLE_ARGPTR_DEC Particle* restrict p, 
    BE_ARGPTR_DEC double const* restrict beam_elements_buffer,
    unsigned int const current_slot_idx, 
    unsigned int const num_slots_in_beam_elements_buffer )
{
    unsigned int next_slot_idx = current_slot_idx;
    
    beam_element_type const type_id = ( beam_element_type 
        )beam_elements_buffer[ current_slot_idx ];
        
    switch( type_id )
    {
        case BEAM_ELEMENT_DRIFT:
        {
            BE_ARGPTR_DEC Drift const* belem = ( BE_ARGPTR_DEC Drift const* 
                )( uintptr_t )&beam_elements_buffer[ current_slot_idx ];
                
            Drift_track( belem, p );            
            LimitGlobal_track( p );
            
            next_slot_idx = ( p->state == 1 )
                ? next_slot_idx + ( unsigned int )NUM_SLOTS_DRIFT
                : num_slots_in_beam_elements_buffer;
            
            break;
        }
        
        case BEAM_ELEMENT_MULTIPOLE:
        {
            BE_ARGPTR_DEC Multipole const* belem = ( BE_ARGPTR_DEC Multipole 
                const* )( uintptr_t )&beam_elements_buffer[ current_slot_idx ];
                
            Multipole_track( belem, p );
            next_slot_idx += ( unsigned int )NUM_SLOTS_MULTIPOLE;
            break;
        }
        
        case BEAM_ELEMENT_CAVITY:
        {
            BE_ARGPTR_DEC Cavity const* belem = ( BE_ARGPTR_DEC Cavity 
                const* )( uintptr_t )&beam_elements_buffer[ current_slot_idx ];
                
            Cavity_track( belem, p );
            next_slot_idx += ( unsigned int )NUM_SLOTS_CAVITY;
            break;
        }
        
        case BEAM_ELEMENT_SC_COASTING:
        {
            BE_ARGPTR_DEC SpaceChargeCoasting const* belem = ( BE_ARGPTR_DEC 
                SpaceChargeCoasting const* )( uintptr_t 
                    )&beam_elements_buffer[ current_slot_idx ];
                
            SpaceChargeCoasting_track( belem, p );
            next_slot_idx += ( unsigned int )NUM_SLOTS_SC_COASTING;
            break;
        }
        
        default:
        {
            next_slot_idx = num_slots_in_beam_elements_buffer;
        }
    };
    
    if( next_slot_idx < num_slots_in_beam_elements_buffer )
    {
        p->at_element_id++;
    }
    else
    {
        next_slot_idx = num_slots_in_beam_elements_buffer;
    }
    
    return next_slot_idx;
}

void Track_particle_until_turn( PARTICLE_ARGPTR_DEC Particle* restrict p,
    BE_ARGPTR_DEC double const* restrict beam_elements_buffer, 
    unsigned int const num_slots_in_beam_elements_buffer,
    long int const until_turn )
{
    typedef long int nelements_t;
    typedef long int index_t;
    
    while( ( p->at_turn < until_turn ) && ( p->state == 1 ) )
    {
        unsigned int idx = ( unsigned int )0;
        p->at_element_id = 0;
        
        while( ( p->state == 1 ) && 
               ( idx < num_slots_in_beam_elements_buffer ) )
        {
            idx = Track_beam_element_dispatcher( p, beam_elements_buffer, 
                    idx, num_slots_in_beam_elements_buffer );
        }

        if( p->state == 1 )
        {
            p->at_element_id = 0;
            p->at_turn++;
        }
    }
}

#endif /* DEMOTRACK_TRACK_H__ */
