import demotrack as dt
import pyopencl
import pyopencl.tools
import pyopencl.array
import numpy as np
import time


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
    particles_memory = "global"
    has_spacecharge = False
    track_prg = dt.get_track_particle_program(
        ctx, particles=particles_memory, has_spacecharge=has_spacecharge)

    until_turn=np.int64(100)

    start_time = time.perf_counter()
    track_ev = track_prg.Track_particle(
        queue, (50*1024,), None,
        beam_arg.data, num_particles,
        lattice_arg.data, num_slots_in_buffer, until_turn, )

    track_ev.wait()
    stop_time = time.perf_counter()

    print( f"Elapsed time: {stop_time - start_time} sec" )
    beam_arg.get(queue,beam_buffer.particle_set_buffer)


    for ii in range( 0, num_particles):
        if beam_buffer.particle_set_buffer[ ii ][ "state" ] == 1:
            assert beam_buffer.particle_set_buffer[ ii ][ "at_turn" ] == until_turn
            assert beam_buffer.particle_set_buffer[ ii ][ "at_element" ] == 0








