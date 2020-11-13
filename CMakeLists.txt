cmake_minimum_required( VERSION 3.8 FATAL_ERROR )
project( demotrack LANGUAGES C CXX )

# Enable CUDA support

include( CheckLanguage )
check_language( CUDA )

if( CMAKE_CUDA_COMPILER )
    message( STATUS "-- enabling CUDA language" )
    enable_language( CUDA )
    message( STATUS "-- CUDA compiler version ${CMAKE_CUDA_COMPILER_VERSION}" )
endif()

# Find CUDA runtime library -> CMake 3.17 or higher has a much more elegant
# way of doing this

set( DEMOTRACK_CUDA_INCLUDE_DIRS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES} )

find_library( CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES} )

if( CUDART_LIBRARY )
    set( DEMOTRACK_CUDA_LIBRARIES ${CUDART_LIBRARY} )
else()
    message( WARNING "------ Unable to find libcudart in
        ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES}" )
endif()

# Build examples

add_subdirectory( examples )
