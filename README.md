Rusph
=====
> Implementation of the SPH algorithm in the Rust programming language.

Rusph is a Smoothed Particle Hydrodynamics code implemented in Rust, a low-level programming language, that uses an ownership system for all data in memory that intrinsically guarantees memory safety and concurrency safety. It is parallelized using the Rayon library and operates in 3D.

Test systems
------------

- Sedov blast wave
- Sod shock tube
- Kelvin-Helmholtz instability
- Toy star
- Turbulent Gas
- Accretion Disc

Future work
-----------

So far, Rusph includes only hydrodynamical equations for non-viscous fluids. However, efforts are underway to incorporate more physics and additional factors in the future, such as self-gravity, magnetohydrodynamics, and dust particle dynamics.

License
-----------

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

