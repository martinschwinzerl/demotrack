import numpy as np
import pyopencl
import pyopencl.tools
import os

from .particle import Particle

def setup_pyopencl_ctx(platform_id=None, device_id=0, interactive=True):
    ctx = None
    cl_Particle = None
    cl_Particle_cdecl = None

    if platform_id is not None:
        platforms = pyopencl.get_platforms()
        ctx = pyopencl.Context(
            dev_type=pyopencl.device_type.ALL,
            properties=[
                (pyopencl.context_properties.PLATFORM,
                 platforms[ platform_id ] )])
    else:
        ctx = pyopencl.Context(interactive=interactive)

    if ctx is not None:
        cl_Particle, cl_Particle_cdecl = pyopencl.tools.match_dtype_to_c_struct(
            ctx.devices[device_id], "Particle", Particle )

        cl_Particle = pyopencl.tools.get_or_register_dtype( "Particle", cl_Particle )

    return ( ctx, cl_Particle, cl_Particle_cdecl )


def get_track_particle_program(ctx, particles="global", has_spacecharge=False):
    import demotrack as dt
    base_path = os.path.dirname( dt.__file__ )
    path_to_src_dir = os.path.join( base_path, os.pardir, "src" )

    src = ""
    if particles == "global":
        src += f"""
        #define PARTICLE_ARGPTR_DEC __global
        """
    else:
        src += """
        #define PARTICLE_ARGPTR_DEC __private
        """

    spacecharge_flag = 1 if has_spacecharge else 0

    src += f"""
    #define ENABLE_SPACECHARGE {spacecharge_flag}
    #define BE_ARGPTR_DEC __global

    #include "particle.h"
    #include "beam_elements.h"
    #include "space_charge.h"
    #include "track.h"

    __kernel void Track_particle_global(
        PARTICLE_ARGPTR_DEC Particle* restrict pset,
        long int const num_particles,
        BE_ARGPTR_DEC double const* restrict beam_elements_buffer,
        long int const num_slots_in_buffer,
        long int const until_turn )"""

    if particles == "global":
        src += """
        {
            long int part_idx = ( long int )get_global_id( 0 );

            for( ; part_idx < num_particles ; part_idx += get_global_size( 0 ) )
            {
                __global Particle* p = &pset[ part_idx ];
                Track_particle_until_turn( p, beam_elements_buffer,
                    num_slots_in_buffer, until_turn );
            }
        }
        """
    else:
        src += """
        {
            long int part_idx = ( long int )get_global_id( 0 );

            for( ; part_idx < num_particles ; part_idx += get_global_size( 0 ) )
            {
                __private Particle p = pset[ part_idx ];
                Track_particle_until_turn( &p, beam_elements_buffer,
                    num_slots_in_buffer, until_turn );

                pset[ part_idx ] = p;
            }
        }
        """

    prg = pyopencl.Program( ctx, src )
    prg.build( options=[ f"-I {path_to_src_dir}" ] )

    return prg
