#!/bin/bash

allpackages=("Statistics","NBInclude","JLD2","ITensors","LinearAlgebra","Test","Revise","KrylovKit","Printf","BenchmarkTools","Random")

julia import-julia-packages.jl "$allpackages"









