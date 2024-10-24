#####################################################
#=

This file contains plotting functions for 1D effective MPS results

Depends on:
    other-funcs/basic-2d-stuff.jl
    two-dimensions.jl
    observables.jl
    hatsugai-mbcn.jl

=#
######################################################

using PyPlot,LaTeXStrings


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
        scatter(xval + shift,nrgs[i] - if_diff*nrgs[1],c=cols[i])
    end
    xlabel(xstring)
    ystring = if_diff ? "NRG - E0" : "NRG"
    ylabel(ystring)
    title("Energy Spectrum"*plot_title)

    return
end
plot_spectrum(xxs::StepRangeLen,nrgs::Vector,idx::Int,nev::Int,xstring::AbstractString,if_diff::Bool; kwargs...) = plot_spectrum(collect(xxs),nrgs,idx,nev,xstring,if_diff; kwargs...)

function plot_omega(theta_xs::Vector{Float64},theta_ys::Vector{Float64},omegas::Matrix{ComplexF64}; kwargs...)
    plot_title::String = get(kwargs,:plot_title,"")
    if_mag::Bool = get(kwargs,:if_mag,false)
    if_perfect_grid::Bool = get(kwargs,:if_perfect_grid,true)

    omegas_phase::Matrix{Float64} = zeros(Float64,length(theta_xs),length(theta_ys))
    for i in 1:length(theta_xs)
        for j in 1:length(theta_ys)
            omegas_phase[i,j] = angle(omegas[i,j]) + pi
        end
    end
    omegas_phase[1,1] = 0.0

    plotting_omega = reverse(transpose(omegas_phase),dims=1)

    fig = figure()
    if if_perfect_grid
        imshow(plotting_omega; cmap="hsv", extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)])
        colorbar()
    end
    xs = transpose(repeat(theta_xs,1,length(theta_ys)))
    ys = reverse(repeat(theta_ys,1,length(theta_xs)),dims=1)
    us = cos.(plotting_omega)
    vs = sin.(plotting_omega)
    quiver(xs, ys, us, vs)
    title("Phase of "*L"\Omega"*plot_title)
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
end
plot_omega(theta_xs::StepRangeLen,theta_ys::StepRangeLen,omegas::Matrix{ComplexF64}; kwargs...) = plot_omega(collect(theta_xs),collect(theta_ys),omegas; kwargs...)

function plot_gamma(theta_xs::Vector{Float64},theta_ys::Vector{Float64},gammas::Matrix{ComplexF64},which_gamma::Int; kwargs...)
    plot_title = get(kwargs,:plot_title,"")

    plotting_gamma = reverse(transpose(abs.(gammas)),dims=1)

    fig = figure()
    imshow(plotting_gamma; cmap="viridis", extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)])
    colorbar()
    which_gamma == 1 ? title("Magnitude of "*L"\Lambda_1 "*plot_title) : title("Magnitude of "*L"\Lambda_2 "*plot_title)
    ylabel(L"\theta_y / 2\pi")
    xlabel(L"\theta_x / 2\pi")
end
plot_gamma(theta_xs::StepRangeLen,theta_ys::StepRangeLen,gammas::Matrix{ComplexF64},which_gamma::Int; kwargs...) = plot_gamma(collect(theta_xs),collect(theta_ys),gammas,which_gamma; kwargs...)

function plot_greenfunc(all_greens,hopping_direction; kwargs...)
	virt_edge_length,phys_edge_length = size(all_greens)
	title_string = "$hopping_direction Greens Function " * get(kwargs, :plot_title, "")
	if_lines = get(kwargs, :if_lines, false)
	fig = figure()
	if if_lines
		if hopping_direction == "virt"
			for i in 1:phys_edge_length
				plot(1:virt_edge_length,all_greens[:,i],"-p",label="$i")
			end
			xlabel("Synthetic")
		else
			for i in 1:virt_edge_length
				plot(1:phys_edge_length,all_greens[i,:],"-p",label="$i")
			end
			xlabel("Physical")
		end
	else
		imshow(abs.(all_greens), vmin=0.0, vmax=1.0)
		colorbar()
		xlabel("Physical")
		ylabel("Synthetic")
	end
	title(title_string)
	
	#=if_save_fig = get(kwargs, :if_save_fig, false)
	if if_save_fig
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "$hopping_direction-dir-GF")
		filename = check_plot_label(filename,"$hopping_direction-dir-GF")
	end
	if_save_fig ? save_figure(filename; location=location) : nothing=#
end



































"fin"