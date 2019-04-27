# Furious Birds Milestone II Start Code

A simple simulation framework using libigl and cmake. Based on Alec Jacobson's libigl example project. This project contains some boilerplate that sets up a physical simulation to run in its own thread, with rendering provided by libigl.

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `birds_bin` binary.

## Run

From within the `build` directory just issue:

    ./birds_bin

A glfw app should launch displaying a GUI. The code will try to load scene files from either
the ./scenes or ../scenes folder, so you need to run the binary from either the project root
folder, or a build subdirectory.

## Dependencies

The only dependencies are [libigl](libigl.github.io/libigl/) and the dependencies of its GUI (glfw and opengl).

We recommend you to install libigl using git by cloning the repository at https://github.com/libigl/libigl.