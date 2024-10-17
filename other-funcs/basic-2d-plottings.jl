#####################################################
#=

This file contains useful plotting functions for 2D physics

=#
######################################################

using PyPlot,LaTeXStrings

function range_cdwsf_angles(points_count::Int64,dens_corr_mat::Array{Float64},radius::Int64=1; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)
    plot_title::String = get(kwargs,:plot_title,"")
    if_multiple_lines::Bool = get(kwargs,:if_multiple_lines,false)
    line_label::String = get(kwargs,:line_label,"")

    angles = range(0.0,2*pi,length=points_count)
    cdwsfs::Vector{ComplexF64} = zeros(ComplexF64,points_count)
    for (idx,angle) in enumerate(angles)
        qvec::Vector{Float64} = [radius*cos(angle),radius*sin(angle)]
        cdwsfs[idx] = get_cdwsf(qvec,denscorrs)
    end

    if if_plot
        if if_multiple_lines
            plot(angles ./ (2*pi),abs.(cdwsfs),"-p",label=line_label)
            legend()
        else
            fig = figure()
            plot(angles ./ (2*pi),abs.(cdwsfs),"-p")
        end
        xlabel("Angle")
        ylabel("CDW Structure Factor")
        title(plot_title)
    end

    return angles,cdwsfs
end

function range_cdwsf_radii(points_count::Int64,dens_corr_mat::Array{ComplexF64}; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)
    plot_title::String = get(kwargs,:plot_title,"")

    Lphys::Int = size(dens_corr_mat)[1]
    
    min_radius::Float64 = get(kwargs,:min_radius,0.1)
    max_radius::Float64 = get(kwargs,:max_radius,Lphys/2)
    angles::Vector{Float64} = get(kwargs,:angles,[0.0,pi/2,pi,3*pi/2])

    radii = range(min_radius,max_radius,length=points_count)
    cdwsfs::Vector{Vector{ComplexF64}} = [zeros(ComplexF64,points_count) for i in 1:length(angles)]
    for (idx2,angle) in enumerate(angles)
        for (idx,radius) in enumerate(radii)
            cdwsfs[idx2][idx] = get_cdwsf([radius*cos(angle),radius*sin(angle)],denscorrs)
        end

        if if_plot
            plot(radii,abs.(cdwsfs[idx2]),"-p",label="$(round(angle/(2*pi),digits=3))π")
            xlabel("Radius")
            ylabel("CDW Structure Factor")
            title(plot_title)
            legend()
        end
    end

    return radii,cdwsfs
end
































"fin"