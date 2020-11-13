configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/demo_lattice.bin"
                "${CMAKE_CURRENT_BINARY_DIR}/demo_lattice.bin" COPYONLY )

# ------------------------------------------------------------------------------

add_executable( cuda_demotrack_nonopt cuda_demotrack_nonopt.cu )

target_include_directories( cuda_demotrack_nonopt PRIVATE ${CMAKE_SOURCE_DIR} )
target_compile_options( cuda_demotrack_nonopt PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-cudart shared -lineinfo>)

set_target_properties( cuda_demotrack_nonopt PROPERTIES
    CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 11 CUDA_STANDARD_REQUIRED 1 )

target_link_libraries( cuda_demotrack_nonopt PRIVATE ${DEMOTRACK_CUDA_LIBRARIES} )

# ------------------------------------------------------------------------------

add_executable( cuda_demotrack_local_particles cuda_demotrack_local_particles.cu )

target_include_directories( cuda_demotrack_local_particles PRIVATE ${CMAKE_SOURCE_DIR} )
target_compile_options( cuda_demotrack_local_particles PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-cudart shared -lineinfo>)

set_target_properties( cuda_demotrack_local_particles PROPERTIES
    CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 11 CUDA_STANDARD_REQUIRED 1 )

target_link_libraries( cuda_demotrack_local_particles PRIVATE ${DEMOTRACK_CUDA_LIBRARIES} )

# ------------------------------------------------------------------------------

add_executable( cuda_demotrack_occupancy cuda_demotrack_occupancy.cu )

target_include_directories( cuda_demotrack_occupancy PRIVATE ${CMAKE_SOURCE_DIR} )
target_compile_options( cuda_demotrack_occupancy PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-cudart shared -lineinfo>)

set_target_properties( cuda_demotrack_occupancy PROPERTIES
    CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON
    CUDA_STANDARD 11 CUDA_STANDARD_REQUIRED 1 )

target_link_libraries( cuda_demotrack_occupancy PRIVATE ${DEMOTRACK_CUDA_LIBRARIES} )
