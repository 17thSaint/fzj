using Test
include("long-range-ttn.jl")


do_all = true


if do_all | false
@testset "Flat with fixed edge" begin;

	edge_length = 8
	end_length = 4
	strength = 1.0
	expected = [strength,strength,strength,strength,strength,0,0,0]
	expected_onsite = [strength,0,0,0,0,0,0,0]
	calc_onsite = long_range_scaling(0,edge_length,strength; scaling="flat",cliff=false)
	calc_norm = long_range_scaling(end_length,edge_length,strength; scaling="flat",cliff=false)
	calc_cliff = long_range_scaling(end_length,edge_length,strength; scaling="flat",cliff=true)
	
	@test expected == calc_norm
	@test expected == calc_cliff
	@test expected_onsite == calc_onsite
end;
end

if do_all | false
@testset "Exponential decay" begin;

	edge_length = 8
	end_length = 4
	strength = 1.0
	trunc_error = 0.006
	
	expected_open = strength .* [exp(-log(10^3)*(i-1)/end_length) for i in 1:edge_length] 
	expected_cliff = [expected_open[1:end_length+1]; [0 for i in 1:edge_length - end_length - 1]]
	expected_trunc = [expected_open[1:findfirst(x -> x <= trunc_error,expected_open)-1]; [0,0,0,0,0]]
	
	calc_open = long_range_scaling(end_length,edge_length,strength; scaling="exp",rounding=false)
	calc_cliff = long_range_scaling(end_length,edge_length,strength; scaling="exp",rounding=false,cliff=true)
	calc_trunc = long_range_scaling(end_length,edge_length,strength; scaling="exp",rounding=true,trunc_min=trunc_error)
	
	@test isapprox(expected_open,calc_open)
	@test isapprox(expected_cliff,calc_cliff)
	@test isapprox(expected_trunc,calc_trunc)
end;
end
































"fin"
