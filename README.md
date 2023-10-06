# Ptycho.jl

[![Build Status](https://github.com/xpwu/Ptycho.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/xpwu/Ptycho.jl/actions/workflows/CI.yml?query=branch%3Amain)

Ptycho.jl is a Julia library for ptychographic reconstruction.

Currently the package has only implemented the ePIE (extended Ptychographic Iterative Engine) algorithm, and is still under development.

Ptycho.jl relies on KernelAbstractions.jl to provide an unified interface for different hardware architectures, right now it should work on any KA compatible backends other than Metal, due to the lack of FFT library for the Apple GPUs.

To try out the package, please see the ePIE example script in the examples folder.