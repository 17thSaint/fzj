#####################################################
#=

This file contains useful plotting functions for 2D physics

=#
######################################################

using LaTeXStrings,PyPlot

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

function plot_spectrum(xxs::Vector,nrgs::Vector,idx::Int,nev::Int,xstring::AbstractString="x",if_diff::Bool=true; kwargs...)
    plot_title = get(kwargs,:plot_title,"")

    cols = ["b","g","r","m","c"]
    if nev > length(cols)
        cols = repeat(cols,ceil(Int,nev/length(cols)))
    end

    for i in 1:nev
        change = abs(xxs[1] - xxs[2])
        xval = xxs[idx]
        shift = (i - nev/2) * ((0.1*change)/(nev/2))
        display(xval)
        display(nrgs[i] - if_diff*nrgs[1])
        scatter(xval + shift,nrgs[i] - if_diff*nrgs[1],c=cols[i])
    end
    xlabel(xstring)
    ystring = if_diff ? "NRG - E0" : "NRG"
    ylabel(ystring)
    title("Energy Spectrum"*plot_title)

    return
end
plot_spectrum(xxs::StepRangeLen,nrgs::Vector,idx::Int,nev::Int,xstring::AbstractString,if_diff::Bool; kwargs...) = plot_spectrum(collect(xxs),nrgs,idx,nev,xstring,if_diff; kwargs...)

function plot_fullspectrum(xs::Vector,nrgs::Dict,xstring::AbstractString="x",if_diff::Bool=true; kwargs...)
    fig = figure()
    for (idx,x) in enumerate(xs)
        local_nrgs = filter(x->x<1000.0,[nrgs[string(i)][idx] for i in 1:length(nrgs)])
        plot_spectrum(xs,local_nrgs,idx,length(local_nrgs),xstring,if_diff; kwargs...)
    end
    return
end

function plot_gamma(theta_xs::Vector{Float64},theta_ys::Vector{Float64},gammas::Matrix{ComplexF64},which_gamma::Int; kwargs...)
    plot_title = get(kwargs,:plot_title,"")

    plotting_gamma = reverse(transpose(abs.(gammas)),dims=1)

    fig = figure()
    imshow(plotting_gamma; cmap="viridis", extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)], vmin=0.0)
    colorbar()
    which_gamma == 1 ? title("Magnitude of "*L"\Lambda_1 "*plot_title) : title("Magnitude of "*L"\Lambda_2 "*plot_title)
    ylabel(L"\theta_y / 2\pi")
    xlabel(L"\theta_x / 2\pi")
end
plot_gamma(theta_xs::StepRangeLen,theta_ys::StepRangeLen,gammas::Matrix{ComplexF64},which_gamma::Int; kwargs...) = plot_gamma(collect(theta_xs),collect(theta_ys),gammas,which_gamma; kwargs...)

function plot_omega(theta_xs::Vector{Float64},theta_ys::Vector{Float64},omegas::Matrix{ComplexF64}; kwargs...)
    plot_title::String = get(kwargs,:plot_title,"")
    if_mag::Bool = get(kwargs,:if_mag,false)
    if_perfect_grid::Bool = get(kwargs,:if_perfect_grid,true)
    if_count_chern::Bool = get(kwargs,:if_count_chern,true)

    omegas_phase::Matrix{Float64} = zeros(Float64,length(theta_xs),length(theta_ys))
    for i in 1:length(theta_xs)
        for j in 1:length(theta_ys)
            omegas_phase[i,j] = angle(omegas[i,j]) + pi
        end
    end
    #omegas_phase[1,1] = 0.0

    plotting_omega = reverse(transpose(omegas_phase),dims=1)
    fig = figure()
    if if_perfect_grid
        imshow(plotting_omega; cmap="hsv", extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)], vmax=2*pi, vmin=0)
        colorbar()
    end
    xs = transpose(repeat(theta_xs,1,length(theta_ys)))
    ys = reverse(repeat(theta_ys,1,length(theta_xs)),dims=1)
    us = cos.(plotting_omega)
    vs = sin.(plotting_omega)
    quiver(xs, ys, us, vs)
    title("Phase of "*L"\Omega"*" "*plot_title)
    xlabel(L"\theta_x / 2\pi")
    ylabel(L"\theta_y / 2\pi")
    #xlim([minimum(theta_xs),maximum(theta_xs)])
    #ylim([minimum(theta_ys),maximum(theta_ys)])

    if if_mag
        fig = figure()
        imshow(abs.(omegas); cmap="viridis", extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)])
        colorbar()
        title("Magnitude of Omega"*plot_title)
        xlabel("Theta_x / 2pi")
        ylabel("Theta_y / 2pi")
        xlim([minimum(theta_xs),maximum(theta_xs)])
        ylim([minimum(theta_ys),maximum(theta_ys)])
    end

    if_count_chern ? count_chern_number(theta_xs,theta_ys,plotting_omega; kwargs...) : nothing
end
plot_omega(theta_xs::StepRangeLen,theta_ys::StepRangeLen,omegas::Matrix{ComplexF64}; kwargs...) = plot_omega(collect(theta_xs),collect(theta_ys),omegas; kwargs...)

function count_chern_number(theta_xs::Vector{Float64},theta_ys::Vector{Float64},plotting_omega::Matrix{Float64}; kwargs...)
    plot_title::String = get(kwargs,:plot_title,"")
    if_plot::Bool = get(kwargs,:if_plot,true)

    reverse!(theta_ys)

    vortex_counting::Matrix{Float64} = zeros(Float64,length(theta_xs),length(theta_ys))
    for i in 2:length(theta_xs)-1
        topsite = i + 1
        downsite = i - 1
        for j in 2:length(theta_ys)-1
            rightsite = j + 1
            leftsite = j - 1
            alldiffs = Matrix{Float64}(undef,2,2)
            alldiffs[1,1] = plotting_omega[i,leftsite] - plotting_omega[topsite,j]
            alldiffs[2,1] = plotting_omega[topsite,j] - plotting_omega[i,rightsite]
            alldiffs[2,2] = plotting_omega[i,rightsite] - plotting_omega[downsite,j]
            alldiffs[1,2] = plotting_omega[downsite,j] - plotting_omega[i,leftsite]
            for ii in 1:2
                for jj in 1:2
                    if abs(alldiffs[ii,jj]) > pi
                        alldiffs[ii,jj] -= sign(alldiffs[ii,jj]) * 2 * pi
                    end
                end
            end
            #display(round.(alldiffs ./ pi,digits=3))
            vortex_counting[i,j] = sum(alldiffs)
        end
    end
    

    if if_plot
        fig = figure()
        imshow(vortex_counting ./ (1*pi); cmap="viridis", extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)], vmax=1.0, vmin=-1.0)
        colorbar()
        #title(L"\sum \Omega"*plot_title)
        title("Chern Number "*L"U_{ir}"*"=1000")
        xlabel(L"\theta_x / 2\pi")
        ylabel(L"\theta_y / 2\pi")
    end

    return sum(vortex_counting) / (2*pi)
end

# plotting the minimum amount the GSs go into the excited states as a function of interaction strength
function plot_twistflatness_vs_intstren(intstrens,flatnesses; kwargs...)
    plot_title = get(kwargs,:plot_title,"")

    fig = figure()
    scatter(intstrens,flatnesses,c="b")
    xlabel("Interaction Strength")
    ylabel("Maximum Twist Flatness")
    ylim(-0.02,1.05)
    title(plot_title)

end

# plot spectrum under twisting
function plot_twisting_spectrum(tw1s::Vector{Float64},tw2s::Vector{Float64},nrgs::Dict{String,Vector{Float64}}; kwargs...)
    plot_title::String = get(kwargs,:plot_title,"")
    if_gap::Bool = get(kwargs,:if_gap,true)

    fig = figure()
    scatter3D(tw1s,tw2s,nrgs["1"] .- nrgs["1"].*if_gap,c="b")
    scatter3D(tw1s,tw2s,nrgs["2"].- nrgs["1"].*if_gap,c="g")
    scatter3D(tw1s,tw2s,nrgs["3"].- nrgs["1"].*if_gap,c="r")
    xlabel(L"\theta_x / 2\pi")
    ylabel(L"\theta_y / 2\pi")
    zlabel("Energy")
    title("Spectrum "*plot_title)
end

# plot distance density-density correlation function
function plot_distdensdenscorrs(distances::Vector{Int64},denscorrs::Vector{Float64},which_direction::String; kwargs...)
    plot_title::String = get(kwargs,:plot_title,"")

    fig = figure()
    plot(distances,denscorrs,"-p")
    xlabel("Distance")
    ylabel("Density-Density Correlation")
    title("$which_direction Density-Density Correlation "*plot_title)
end

function plot_occupancy(exp_occ; kwargs...)
    vmax = get(kwargs,:vmax,nothing)
    #exp_occ = transpose(exp_occ)
	fig = figure()
	isnothing(vmax) ? imshow(exp_occ,origin="lower",vmin=0) : imshow(exp_occ,origin="lower",vmin=0,vmax=vmax)
	colorbar()
	plot_title = get(kwargs, :plot_title, "")
	title_string = "Occupancy, " * plot_title
	title(title_string)
	ylabel("Synthetic")
	xlabel("Physical")

    return nothing
end

function plot_fourpointcorrelator(fourpointcorrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    vmax = get(kwargs,:vmax,nothing)
    fig = figure()
    isnothing(vmax) ? imshow(fourpointcorrs,vmin=0.0,vmax=0.02,origin="lower") : imshow(fourpointcorrs,vmin=0.0,vmax=vmax)
    colorbar()
    title("Four Point Correlator " * plot_title)
    xlabel("Physical")
    ylabel("Synthetic")
    return nothing
end

function plot_pairdistribution(pairdist::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    vmax = get(kwargs,:vmax,nothing)
    fig = figure()
    isnothing(vmax) ? imshow(pairdist,vmin=0.0) : imshow(pairdist,vmin=0.0,vmax=vmax)
    colorbar()
    title("Pair Distribution " * plot_title)
    xlabel("Physical")
    ylabel("Synthetic")
    return nothing
end

function plot_four_point(results::Matrix{Float64}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    if_2pt::Bool = get(kwargs,:if_2pt,false)
    if_2pt ? pt_title = "2pt " : pt_title = "4pt "
    fig = figure()
    imshow(results,vmin=0.0,origin="lower")
    colorbar()
    title("$pt_title Momentum "*plot_title)
    ylabel("m")
    xlabel("m'")
    return nothing
end

function plot_four_point(results::Vector{Float64},mp::Int64; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    plot(collect(0:length(results)-1),results,"-p")
    title("4pt Momentum "*plot_title)
    ylabel("4pt Momentum")
    xlabel("Momentum k = m / Ly, m' = $mp")
    return nothing
end

function plot_adiabatic_condition(xs::Vector{Float64},ys::Vector{Float64},f1s::Matrix{Float64},f2s::Matrix{Float64}; kwargs...)
    plot_title::String = get(kwargs,:plot_title,"")

    fs_angle::Matrix{Float64} = atan.(f1s,f2s)
    fs_mag::Matrix{Float64} = sqrt.(f1s.^2 .+ f2s.^2)
    
    fig = figure()
    
    pcolormesh(xs,ys,log10.(fs_mag); shading="auto")
    #imshow(log10.(fs_mag); cmap="viridis", origin="lower", extent=[minimum(xs),maximum(xs),minimum(ys),maximum(ys)])
    colorbar()

    #us = cos.(fs_angle)
    #vs = sin.(fs_angle)
    #quiver(xs[1:end-1], ys[1:end-1], us, vs)
    title("Adiabatic Condition "*plot_title)
    xlabel("Physical Hopping, "*L"t_x")
    ylabel("IR Interaction Strength, "*L"U_{ir}")

end

































"fin"