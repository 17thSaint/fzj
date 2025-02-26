#####################################################
#=

This file contains little overwriting functions over the existing TTN
or TTN-dependency connected functions.

Depends on:
    TTN

=#
######################################################


TTN.linear_ind(lat::TTN.SimpleLattice{D}, here::Vector{Int}) where D = TTN.linear_ind(lat, (here[1],here[2]))








































"fin"