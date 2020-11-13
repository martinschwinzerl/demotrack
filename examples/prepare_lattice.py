import demotrack as dt
import numpy as np

if __name__ == '__main__':
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
    lattice_buffer.astype( 'float64' ).tofile( './demo_lattice.bin' )

