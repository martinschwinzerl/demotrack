# demotrack

## About

`demotrack` is a stripped-down, stand-alone, proof-of-concept GPU accelerated beam-dynamics particle tracking tool similar to [SixTrackLib](https://github.com/SixTrack/sixtracklib) or [SixTrack](https://github.com/SixTrack/sixtrack). It implements a particle model similar to `SixTrackLib` and a small list of beam-elements
- drifts
- multipole
- cavity
- coasting SpaceCharge

These elements are sufficent to track a simple FODO lattice. `demotrack` relies on [PyOpenCL](https://documen.tician.de/pyopencl/) for accelerated parallel problems. Cf the `examples/demo.py` file for an usage example.

## Installation

It should be sufficient to run 
```
pip install -e .
``` 

from the main directory of your working copy. 
