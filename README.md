# OpticalRayTrace

Smallish program to trace rays through the optical setup from this [paper](https://doi.org/10.1039/D0AY01101K).

Currently models a point source from centre of bottle, and a ring of point sources on bottle's surface.


## Usage

To run one a single core/thread
```
make && ./raytrace
```

To run on all available threads
```
make mp && ./raytrace
```

## Requirments
Gfortran 7.5+ or 9+ if you want to supress openMP warnings whilst using default(none)
GNU Make 4.1+
