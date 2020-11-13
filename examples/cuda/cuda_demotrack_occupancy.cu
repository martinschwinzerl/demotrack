#include <cuda.h>
// #include <cudart.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>
#include <sys/time.h>

#if !defined( PARTICLE_ARGPTR_DEC )
    #define PARTICLE_ARGPTR_DEC
#endif /* !defined( PARTICLE_ARGPTR_DEC ) */

#if !defined( BE_ARGPTR_DEC )
    #define BE_ARGPTR_DEC
#endif /* !defined( BE_ARGPTR_DEC ) */

#if !defined( DISABLE_SPACECHARGE )
    #define DISABLE_SPACECHARGE 1
#endif /* !defined( DISABLE_SPACECHARGE ) */

#if !defined( DEMOTRACK_STATIC )
    #define DEMOTRACK_STATIC static
#endif /* !defined( DEMOTRACK_STATIC ) */

#if !defined( DEMOTRACK_FN )
    #define DEMOTRACK_FN __device__ __host__
#endif /* !defined( DEMOTRACK_FN ) */

#if !defined( DEMOTRACK_INLINE )
    #define DEMOTRACK_INLINE inline
#endif /* !defined( DEMOTRACK_INLINE ) */

#include "src/particle.h"
#include "src/beam_elements.h"
#include "src/space_charge.h"
#include "src/track.h"

__global__ void track_particles_until_turn_nonoptimised(
    PARTICLE_ARGPTR_DEC Particle* pset, long int const num_particles,
    BE_ARGPTR_DEC double const* beam_elements_buffer,
    long int const num_slots_in_buffer, long int const until_turn )
{
    /* Stride = ( threads/block ) * ( number of blocks ) */
    long int const STRIDE = blockDim.x * gridDim.x;
    long int idx = threadIdx.x + blockIdx.x * blockDim.x;

    for( ; idx < num_particles ; idx += STRIDE )
    {
        /* Copy particle to thread-local memory  */
        Particle p = pset[ idx ];
        Track_particle_until_turn( &p,
            beam_elements_buffer, num_slots_in_buffer, until_turn );

        /* Copy results back to global mmemory */
        pset[ idx ] = p;
    }
}

int main()
{
    /* Source: https://developer.nvidia.com/blog/cuda-pro-tip-occupancy-api-simplifies-launch-configuration/ */

    using std::size_t;

    /* ********************************************************************* */
    /* Prepare particle set to track */
    long int const NUM_PARTICLES = 50 * 1024;
    double const P0_C = 1.0e9; /* Kinetic energy, [eV]  */
    double const MASS0 = 938.272081e6; /* Proton rest mass, [eV] */

    std::vector< ::Particle > pset( NUM_PARTICLES, ::Particle{} );
    ::ParticleSet_setup( pset.data(), pset.size(), P0_C, MASS0 );

    /* ********************************************************************* */
    /* Read lattice from prepared and saved dump */

    size_t const LATTICE_NUM_SLOTS = 153600 / sizeof( double );
    std::vector< double > lattice( LATTICE_NUM_SLOTS, 0.0 );

    char const PATH_TO_LATTICE[] = "./demo_lattice.bin";
    FILE* fp = std::fopen( PATH_TO_LATTICE, "rb" );

    if( fp != nullptr )
    {
        size_t const ret = std::fread( lattice.data(), sizeof( double ),
            lattice.size(), fp );

        if( ret != LATTICE_NUM_SLOTS )
        {
            std::cerr << "Unable to read lattice" << std::endl;
            return 0;
        }

        std::fclose( fp );
        fp = nullptr;
    }

    /* ******************************************************************** */
    /* Prepare device memory */

    ::Particle* particles_dev = nullptr;
    double* lattice_dev = nullptr;

    auto status = ::cudaMalloc( ( void** )&particles_dev,
        sizeof( ::Particle ) * NUM_PARTICLES );
    assert( status == CUDA_SUCCESS );

    status = ::cudaMalloc( ( void** )&lattice_dev,
        sizeof( double ) * LATTICE_NUM_SLOTS );
    assert( status == CUDA_SUCCESS );

    /* Copy particle and lattice data to device */

    status = ::cudaMemcpy( lattice_dev, lattice.data(),
        LATTICE_NUM_SLOTS * sizeof( double ), ::cudaMemcpyHostToDevice );

    assert( status == CUDA_SUCCESS );

    status = ::cudaMemcpy( particles_dev, pset.data(),
        pset.size() * sizeof( ::Particle ), ::cudaMemcpyHostToDevice );

    assert( status == CUDA_SUCCESS );

    /* ******************************************************************** */
    /* Estimate block size */

    int BLOCK_SIZE = 0;
    int MIN_GRID_SIZE = 0;

    status = ::cudaOccupancyMaxPotentialBlockSize(
        &MIN_GRID_SIZE, /* -> minimum grid size needed for max occupancy */
        &BLOCK_SIZE, /* -> estimated optimal block size */
        track_particles_until_turn_nonoptimised, /* the kernel */
        0u, /* -> dynamic shared memory per block required [bytes] */
        0u /* -> max block size limit for the kernel; 0 == no limit */ );

    assert( status == CUDA_SUCCESS );

    assert( BLOCK_SIZE > 0 );
    int const GRID_SIZE = ( NUM_PARTICLES + BLOCK_SIZE - 1 ) / BLOCK_SIZE;

    std::cout << "NUM_OF_BLOCKS     : " << GRID_SIZE << "\r\n"
              << "THREADS_PER_BLOCK : " << BLOCK_SIZE << std::endl;

    /* ******************************************************************** */
    /* Perform calculation */

    long int const TRACK_UNTIL_TURN = 10;

    struct timeval  start;
    struct timeval  stop;

    ::gettimeofday( &start, NULL );

    track_particles_until_turn_nonoptimised<<< GRID_SIZE, BLOCK_SIZE >>>(
        particles_dev, NUM_PARTICLES, lattice_dev,
            LATTICE_NUM_SLOTS, TRACK_UNTIL_TURN );

    status = ::cudaDeviceSynchronize();
    assert( status == CUDA_SUCCESS );

    ::gettimeofday( &stop, NULL );

    std::cout << "Elapsed time: "
              << static_cast< double >( stop.tv_usec - start.tv_usec ) / 1000000.0 +
                 static_cast< double >( stop.tv_sec - start.tv_sec )
              << " seconds" << std::endl;

    /* Calculate theoretical Occupancy */
    int MAX_ACTIVE_BLOCKS;
    status = ::cudaOccupancyMaxActiveBlocksPerMultiprocessor(
        &MAX_ACTIVE_BLOCKS, /* -> returned number of max active blocks */
        track_particles_until_turn_nonoptimised, /* -> Kernel */
        BLOCK_SIZE, /* -> threads / block for kernel */
        0u /* -> dynamic shared memory per block required [bytes] */ );

    assert( status == CUDA_SUCCESS );

    int device;
    ::cudaDeviceProp props;
    status = ::cudaGetDevice( &device );
    assert( status == CUDA_SUCCESS );

    status = ::cudaGetDeviceProperties( &props, device );
    assert( status == CUDA_SUCCESS );

    double const occupancy = ( MAX_ACTIVE_BLOCKS * BLOCK_SIZE / props.warpSize) /
                    (double)(props.maxThreadsPerMultiProcessor /
                            props.warpSize);

    std::cout << "Theoretical occupancy: " << occupancy << "\r\n"
              << "max num active blocks: " << MAX_ACTIVE_BLOCKS << "\r\n"
              << std::endl;


    /* ******************************************************************** */
    /* Fetch particle data */

    status = ::cudaMemcpy( pset.data(), particles_dev,
        pset.size() * sizeof( ::Particle ), ::cudaMemcpyDeviceToHost );

    assert( status == CUDA_SUCCESS );
    ( void )status;

    ::cudaFree( lattice_dev );
    ::cudaFree( particles_dev );

    return 0;
}
