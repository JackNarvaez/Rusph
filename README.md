Rusph
=====
> Implementation of the SPH algorithm in the Rust programming language.

Rusph is a Smoothed Particle Hydrodynamics code implemented in Rust, a low-level language, that uses an ownership system for all data in memory that intrinsically guarantees memory safety and concurrency safety. It is parallelized using the Rayon library and operates in 3D.

Test systems
------------

- The Sedov blast wave problem
- The Kelvin-Helmholtz problem
- Toy star
- Turbulent Gas
