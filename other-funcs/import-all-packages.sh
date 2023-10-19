#!/bin/bash

allpackages=("Statistics","NBInclude","JLD2","ITensors","LinearAlgebra","Test","KrylovKit","Printf","Random")
remainingpackages=("NBInclude","JLD2")

julia import-julia-packages.jl "$allpackages"









