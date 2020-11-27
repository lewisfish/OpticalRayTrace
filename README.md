# OpticalRayTrace

Smallish program to trace rays through the optical setup from this [paper](https://doi.org/10.1039/D0AY01101K).

Currently models a point source from centre of bottle, and a ring of point sources on bottle's surface.


![Image of Rays from ring of sources](https://github.com/lewisfish/OpticalRayTrace/raw/main/images/ring-rays.png)
![Image of rays from point source](https://github.com/lewisfish/OpticalRayTrace/raw/main/images/point-rays.png)



## Usage

To run one a single core/thread
```
./install.sh -n 1
```

To run on all available threads
```
./install.sh
```

## Requirments
Gfortran 7.5+ or 9+ if you want to supress openMP warnings whilst using default(none)
GNU Make 4.1+
