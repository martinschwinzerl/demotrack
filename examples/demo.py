import demotrack as dt
import pyopencl
import pyopencl.tools
import pyopencl.array
import numpy as np
import os

def get_track_particle_program_global(ctx, cl_Particle_cdecl):
    base_path = os.path.dirname( dt.__file__ )
    path_to_src_dir = os.path.join( base_path, os.pardir, "src" )
    
    src = f"""
    #define BE_ARGPTR_DEC __global 
    #define PARTICLE_ARGPTR_DEC __global
    
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
    
    prg = pyopencl.Program( ctx, src )
    prg.build( options=[ f"-I {path_to_src_dir}" ] )
    
    return prg

if __name__ == '__main__':
    ctx, cl_Particle, cl_Particle_cdecl = dt.setup_pyopencl_ctx(
        platform_id=1, device_id=0)
    
    queue = pyopencl.CommandQueue(ctx)
    
    # Setup lattice
    lattice = dt.Lattice()    
    num_cells = 100
    
    for jj in range( 0, num_cells):
        lattice.add_drift(length=0.0)
        lattice.add_drift(length=0.012499999999999956)
        lattice.add_multipole(length=0.025, order=3, knl=[0., 0.70763461, 0.,  0., ])
        lattice.add_drift(length=0.025000000000000022)
        lattice.add_multipole(length=0.025, order=3, knl=[0., 0.70763461, 0.,  0., ])
        lattice.add_drift(length=0.21250000000000002)
        lattice.add_cavity(voltage=0.0, frequency=0.0, lag=0.0 )
        lattice.add_drift(length=0.21250000000000002)
        lattice.add_multipole(length=0.025, order=3, knl=[0., -0.70763461, 0.,  0., ])
        lattice.add_drift(length=0.012499999999999956)
        lattice.add_drift(length=0.0)
        lattice.add_drift(length=0.012499999999999956)
        lattice.add_multipole(length=0.025, order=3, knl=[0., -0.70763461, 0.,  0., ])
        lattice.add_drift(length=0.025000000000000022)
        lattice.add_multipole(length=0.025, order=3, knl=[0., -0.70763461, 0.,  0., ])
        lattice.add_drift(length=0.012499999999999956)
        lattice.add_drift(length=0.0)
        lattice.add_drift(length=0.012499999999999956)
        lattice.add_multipole(length=0.025, order=3, knl=[0., -0.70763461, 0.,  0., ])
        lattice.add_drift(length=0.21250000000000013)
        lattice.add_cavity(voltage=0.0, frequency=0.0, lag=0.0 )
        lattice.add_drift(length=0.2124999999999999)
        lattice.add_multipole(length=0.025, order=3, knl=[0., 0.70763461, 0.,  0., ])
        lattice.add_drift(length=0.025000000000000133)
        lattice.add_multipole(length=0.025, order=3, knl=[0., 0.70763461, 0.,  0., ])
        lattice.add_drift(length=0.012499999999999956)
        
    lattice_buffer = lattice.pack()
    lattice_arg = pyopencl.array.to_device(queue,lattice_buffer)
    num_slots_in_buffer = np.int64(len(lattice_buffer))
        
    
    # Setup particle beam -> argument for OpenCL kernels:
    num_particles = np.int64(50 * 1024)
    beam_buffer = dt.ParticleSet( num_particles=num_particles, p0c=1e9 )
    beam_arg = pyopencl.array.to_device(queue, beam_buffer.particle_set_buffer)
    
    # Prepare compiling the programs
    track_global_prg = get_track_particle_program_global( ctx, cl_Particle_cdecl)
           
    until_turn=np.int64(100)
    track_ev = track_global_prg.Track_particle_global(
        queue, (50*1024,), None, 
        beam_arg.data, num_particles, 
        lattice_arg.data, num_slots_in_buffer, until_turn, )
    
    
                                 
    
    
