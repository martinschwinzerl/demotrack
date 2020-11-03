import numpy as np
from scipy.special import factorial

class Lattice(object):
    def __init__(self):
        self._slots = []
        self.packed_buffer = None
        
    def add_drift(self, length):
        self._slots.append( float( 0 ) )
        self._slots.append( float( length ) )
    
    def add_multipole( self, order, length=0.0, knl=None, ksl=None ):
        if knl is None:
            knl = []
        if ksl is None:
            ksl = []
        if order is None:
            order = 0

        n = max((order + 1), max(len(knl), len(ksl)))
        assert n > 0

        _knl = np.array(knl)
        nknl = np.zeros(n, dtype=_knl.dtype)
        nknl[: len(knl)] = knl
        knl = nknl
        del _knl
        assert len(knl) == n

        _ksl = np.array(ksl)
        nksl = np.zeros(n, dtype=_ksl.dtype)
        nksl[: len(ksl)] = ksl
        ksl = nksl
        del _ksl
        assert len(ksl) == n

        order = n - 1
        bal = np.zeros(2 * order + 2)
        assert len( bal ) <= 16

        idx = np.array([ii for ii in range(0, len(knl))])
        inv_factorial = 1.0 / factorial(idx, exact=True)
        bal[0::2] = knl * inv_factorial
        bal[1::2] = ksl * inv_factorial
        
        stored_bal = np.zeros( 16 )
        stored_bal[0:len(bal)] = bal

        self._slots.append( float( 1 ) )
        self._slots.append( float( order ) )
        self._slots.append( float( length ) )
        self._slots.append( b for b in stored_bal )
    
    def add_cavity(self, voltage, frequency, lag ):
        self._slots.append( float( 2 ) )
        self._slots.append( float( voltage ) )
        self._slots.append( float( frequency ) )
        self._slots.append( float( lag ) )
        
    def pack(self):
        return np.array( self._slots )
    
    def num_slots(self):
        return len( self._slots ) if self._slots is not None else 0
    
    def clear(self):
        self._slots = []

