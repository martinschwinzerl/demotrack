import numpy as np
import pyopencl
import pyopencl.tools

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
    
