import numpy as np

Particle = np.dtype(
    [
        ( "x", np.float64 ), ( "px", np.float64 ),
        ( "y", np.float64 ), ( "py", np.float64 ),
        ( "zeta", np.float64 ), ( "delta", np.float64 ),
        ( "rpp", np.float64 ), ( "rvv", np.float64 ),
        ( "psigma", np.float64 ), ( "chi", np.float64 ),
        ( "charge_ratio", np.float64 ),
        ( "q0", np.float64 ), ( "mass0", np.float64 ),
        ( "beta0", np.float64 ), ( "gamma0", np.float64 ),
        ( "p0c", np.float64 ),
        ( "state", np.int64 ), ( "at_element", np.int64 ),
        ( "at_turn", np.int64 ), ( "id", np.int64 ),
    ] )


class ParticleSet(object):
    def __init__(self, num_particles=1, p0c=1e9, q0=1, mass0=938.272081e6):
        E0 = np.sqrt( p0c * p0c + mass0 * mass0 )
        beta0 = p0c / E0
        gamma0 = E0 / mass0
        _delta = float( 0.0 )

        self.particle_set_buffer = np.empty( num_particles, dtype=Particle )
        self.particle_set_buffer[ "x" ].fill( 0.0 )
        self.particle_set_buffer[ "px" ].fill( 0.0 )
        self.particle_set_buffer[ "y" ].fill( 0.0 )
        self.particle_set_buffer[ "py" ].fill( 0.0 )
        self.particle_set_buffer[ "zeta" ].fill( 0.0 )
        self.particle_set_buffer[ "delta" ].fill( _delta )
        self.particle_set_buffer[ "chi" ].fill( 1.0 )
        self.particle_set_buffer[ "charge_ratio" ].fill( 1.0 )
        self.particle_set_buffer[ "rvv" ].fill( 1.0 )
        self.particle_set_buffer[ "rpp" ].fill( 1.0 )
        self.particle_set_buffer[ "q0" ].fill( q0 )
        self.particle_set_buffer[ "mass0" ].fill( mass0 )
        self.particle_set_buffer[ "beta0" ].fill( beta0 )
        self.particle_set_buffer[ "gamma0" ].fill( gamma0 )
        self.particle_set_buffer[ "p0c" ].fill( p0c )
        self.particle_set_buffer[ "state" ].fill( 1 )
        self.particle_set_buffer[ "at_element" ].fill( 0 )
        self.particle_set_buffer[ "at_turn" ].fill( 0 )

        for ii in range( 0, num_particles ):
            self.particle_set_buffer[ ii ][ "id" ] = ii
            self.set_delta( _delta, ii )

    def set_delta(self, _delta, ii ):
        beta0 = self.particle_set_buffer[ "beta0" ][ ii ]
        ptau = np.sqrt( _delta * _delta + 2 * _delta + 1 / ( beta0 * beta0 ) )
        ptau -= 1 / beta0
        self.particle_set_buffer[ "delta" ][ ii ] = _delta
        self.particle_set_buffer[ "psigma" ][ ii ] = ptau / beta0
