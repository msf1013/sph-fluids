# SPH-based Fluid Simulation

Final project for CS 395T Physical Simulation.

This simulation implements the Smoothed Particle Hydrodynamics method for fluids described in http://matthias-mueller-fischer.ch/publications/sca03.pdf.

## Team members

Abheek Ghosh and Mario Fuentes

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `fluids_bin` binary.

## Run

From within the `build` directory just issue:

    ./fluids_bin

A glfw app should launch displaying a GUI.

## Dependencies

The only dependencies are [libigl](libigl.github.io/libigl/) and the dependencies of its GUI (glfw and opengl).

We recommend you to install libigl using git by cloning the repository at https://github.com/libigl/libigl.

Note from Abheek: I submitted the survey. Thanks. 