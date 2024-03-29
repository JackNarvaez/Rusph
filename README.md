Rusph
=====
> Implementation of the SPH algorithm in the Rust programming language.

Rusph is a Smoothed Particle Hydrodynamics code implemented in Rust, a low-level language, that uses an ownership system for all data in memory that intrinsically guarantees memory safety and concurrency safety. It is parallelized using the Rayon library and operates in 3D.

Test systems
------------

- Sedov blast wave
- Sod shock tube
- Kelvin-Helmholtz instability
- Toy star
- Turbulent Gas

Future work
-----------

So far, Rusph includes only hydrodynamical equations for non-viscous fluids. However, efforts are underway to incorporate additional factors in the future, such as Gravity, Real Viscosity, and Magnetohydrodynamics.

License
-----------

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

