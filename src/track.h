#ifndef DEMOTRACK_TRACK_H__
#define DEMOTRACK_TRACK_H__

#include "beam_elements.h"
#include "space_charge.h"
#include "particle.h"

DEMOTRACK_STATIC DEMOTRACK_FN void Drift_track(
    BE_ARGPTR_DEC const Drift *const  drift,
    PARTICLE_ARGPTR_DEC Particle*  p );

DEMOTRACK_STATIC DEMOTRACK_FN void LimitGlobal_track(
    PARTICLE_ARGPTR_DEC Particle*  p );

DEMOTRACK_STATIC DEMOTRACK_FN void Multipole_track(
    BE_ARGPTR_DEC const Multipole *const  quad,
    PARTICLE_ARGPTR_DEC Particle*  p );

DEMOTRACK_STATIC DEMOTRACK_FN void Cavity_track(
    BE_ARGPTR_DEC const Cavity *const  cavity,
    PARTICLE_ARGPTR_DEC Particle*  p );

/* ------------------------------------------------------------------------- */

DEMOTRACK_STATIC DEMOTRACK_FN unsigned int Track_beam_element_dispatcher(
    PARTICLE_ARGPTR_DEC Particle*  p,
    BE_ARGPTR_DEC double const*  beam_elements_buffer,
    unsigned int const current_slot_idx,
    unsigned int const num_slots_in_beam_elements_buffer );

DEMOTRACK_STATIC DEMOTRACK_FN void Track_particle_until_turn(
    PARTICLE_ARGPTR_DEC Particle*  p,
    BE_ARGPTR_DEC double const*  beam_elements_buffer,
    unsigned int const num_slots_in_beam_elements_buffer,
    long int const until_turn );

/* ************************************************************************* */

DEMOTRACK_INLINE void Drift_track(
    BE_ARGPTR_DEC const Drift *const  drift,
    PARTICLE_ARGPTR_DEC Particle*  p )
{
    double const one_plus_delta = p->delta + ( double )1;
    double const lpzi = ( double )1 / sqrt( one_plus_delta * one_plus_delta
        - ( p->px * p->px + p->py * p->py ) );

    p->x += p->px * drift->length * lpzi;
    p->y += p->py * drift->length * lpzi;
    p->zeta += p->rvv * drift->length - one_plus_delta * lpzi;
}

DEMOTRACK_INLINE void LimitGlobal_track(
    PARTICLE_ARGPTR_DEC Particle*  p )
{
    double const sign_x = ( double )( ( ( double )0 < p->x ) -
                          ( double )( p->x < ( double )0 ) );

    double const sign_y = ( double )( ( ( double )0 < p->y ) -
                          ( double )( p->y < ( double )0 ) );

    p->state &= ( long int )( ( ( sign_x * p->x ) < ( double )1 ) &
                              ( ( sign_y * p->y ) < ( double )1 ) );
}

DEMOTRACK_INLINE void Multipole_track(
    BE_ARGPTR_DEC const Multipole *const  mp,
    PARTICLE_ARGPTR_DEC Particle*  p )
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

DEMOTRACK_INLINE void Cavity_track(
    BE_ARGPTR_DEC const Cavity *const  cavity,
    PARTICLE_ARGPTR_DEC Particle*  p )
{
    double const PI = ( double )3.141592653589793;
    double const phase = ( PI / ( double )180 ) * cavity->lag -
        ( ( ( double )2 * PI ) / ( double )299792458.0 ) *
        ( p->zeta / ( p->beta0 * p->rvv ) ) * cavity->frequency;

    Particle_add_to_energy( p,
        cavity->voltage * p->q0 * p->charge_ratio * sin( phase ) );
}

/* ------------------------------------------------------------------------- */

DEMOTRACK_INLINE unsigned int Track_beam_element_dispatcher(
    PARTICLE_ARGPTR_DEC Particle*  p,
    BE_ARGPTR_DEC double const*  beam_elements_buffer,
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

        #if defined( DISABLE_SPACECHARGE ) && ( DISABLE_SPACECHARGE == 0 )
        case BEAM_ELEMENT_SC_COASTING:
        {
            BE_ARGPTR_DEC SpaceChargeCoasting const* belem = ( BE_ARGPTR_DEC
                SpaceChargeCoasting const* )( uintptr_t
                    )&beam_elements_buffer[ current_slot_idx ];

            SpaceChargeCoasting_track( belem, p );
            next_slot_idx += ( unsigned int )NUM_SLOTS_SC_COASTING;
            break;
        }
        #endif /* defined( DISABLE_SPACECHARGE ) && ( DISABLE_SPACECHARGE == 0 ) */

        default:
        {
            next_slot_idx = num_slots_in_beam_elements_buffer;
        }
    };

    if( next_slot_idx < num_slots_in_beam_elements_buffer )
    {
        p->at_element++;
    }
    else
    {
        next_slot_idx = num_slots_in_beam_elements_buffer;
    }

    return next_slot_idx;
}

DEMOTRACK_INLINE void Track_particle_until_turn( PARTICLE_ARGPTR_DEC Particle*  p,
    BE_ARGPTR_DEC double const*  beam_elements_buffer,
    unsigned int const num_slots_in_beam_elements_buffer,
    long int const until_turn )
{
    typedef long int nelements_t;
    typedef long int index_t;

    while( ( p->at_turn < until_turn ) && ( p->state == 1 ) )
    {
        unsigned int idx = ( unsigned int )0;
        p->at_element = 0;

        while( ( p->state == 1 ) &&
               ( idx < num_slots_in_beam_elements_buffer ) )
        {
            unsigned int const next_idx = Track_beam_element_dispatcher(
                p, beam_elements_buffer, idx,
                num_slots_in_beam_elements_buffer );

            idx = next_idx;
        }

        if( p->state == 1 )
        {
            p->at_element = 0;
            p->at_turn++;
        }
    }
}

#endif /* DEMOTRACK_TRACK_H__ */
