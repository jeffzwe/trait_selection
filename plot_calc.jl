using DataFrames
using SparseArrays
using HypothesisTests
using Statistics
using StatsBase
using Plots
using CSV
using JLD
using FileIO
using Images
using Distributed
using ProgressMeter
using Combinatorics
using SearchSortedNearest
using KernelDensity
using StatsPlots
using Noise
using SharedArrays
using Distances
using Printf
using LinearAlgebra
using Polynomials
using CurveFit
using TiffImages
using Random
using Plots.PlotMeasures



function fill_grid_mat(n_trait_comb, trait_comb, grid_res_hull, grid_res_hull_norm, grid_res_rao, grid_res_rao_norm, n_scal_spec)
	for t in 1:n_trait_comb
		temp = trait_comb[t]
		results = load(@sprintf("community_qhull_jld/grid_res_hull_%d_%d_%d_%d_%d.jld", temp[1], temp[2], temp[3], temp[4], temp[5]))
		results = collect(values(results))
		results = results[1]
		grid_res_hull[:,:,t] = results
		for s in 1:n_scal_spec+1
			grid_res_hull_norm[s,:,t] = results[s,:]./norm(results[s,:])
		end
	end
	
	for t in 1:n_trait_comb
		temp = trait_comb[t]
		results = load(@sprintf("community_rao_jld/grid_res_rao_%d_%d_%d_%d_%d.jld", temp[1], temp[2], temp[3], temp[4], temp[5]))
		results = collect(values(results))
		results = results[1]
		grid_res_rao[:,:,t] = results
		for s in 1:n_scal_spec+1
			grid_res_rao_norm[s,:,t] = results[s,:]./norm(results[s,:])
		end
	end
end

function fill_grid_var(n_grid, n_scal_spec, n_trait_comb, res_mat, var_mat, st_t_2, st_t_3, st_t_4, st_t_5)

	for g in 1:n_grid
	
		if res_mat[1, g, 1] == 0 continue end
		
		# vector for analysis of specific scaling over all trait combinations: column 1 = diversity metric, column 2 = group affiliation
		results = zeros(n_trait_comb, 2)

		# add group affiliation to each combination
		for t in 1:n_trait_comb
			
			results[t, 1] = res_mat[1,g,t]

			if (1 <= t <= st_t_2-1) results[t,2] = 1 end
			if (st_t_2 <= t <= st_t_3-1) results[t,2] = 2 end
			if (st_t_3 <= t <= st_t_4-1) results[t,2] = 3 end
			if (st_t_4 <= t <= st_t_5-1) results[t,2] = 4 end
			if (st_t_5 <= t <= n_trait_comb) results[t,2] = 5 end
		end

		temp1 = [results[x,1] for x in collect(1:n_trait_comb) if results[x,2] == 1]
		temp2 = [results[x,1] for x in collect(1:n_trait_comb) if results[x,2] == 2]
		temp3 = [results[x,1] for x in collect(1:n_trait_comb) if results[x,2] == 3]
		temp4 = [results[x,1] for x in collect(1:n_trait_comb) if results[x,2] == 4]
		temp5 = [results[x,1] for x in collect(1:n_trait_comb) if results[x,2] == 5]


		# Calc of coefficient variances of total and per group
		var_mat[g, 1] = std(results[:,1])/mean(results[:,1])

		if length(temp1) != 0 var_mat[g, 2] = std(temp1)/mean(temp1) end
		if length(temp2) != 0 var_mat[g, 3] = std(temp2)/mean(temp2) end
		if length(temp3) != 0 var_mat[g, 4] = std(temp3)/mean(temp3) end
		if length(temp4) != 0 var_mat[g, 5] = std(temp4)/mean(temp4) end
		if length(temp5) != 0 var_mat[g, 6] = std(temp5)/mean(temp5) end
	end
end

function map_coeff_var(n_scal_spec, n_grid, scal_spec, grid_raster_sp_indices, grid_raster, var_mat, metric)
	for i in 1:6
		fill!(grid_raster, NaN)
		for g in 1:n_grid
			x = getindex(grid_raster_sp_indices[g], 1)
			y = getindex(grid_raster_sp_indices[g], 2)
			grid_raster[x,y] = var_mat[g, i]
		end

		tiff = reinterpret(Gray{Float64}, grid_raster)
		tiff = TiffImages.DenseTaggedImage(tiff)
		TiffImages.save(@sprintf("plots/%s/map_group_%d.tif", metric,i), tiff)
	end	
end

function plot_rich_cv(n_grid, grid_species, var_mat, group, colors, metric)

	p = Plots.plot(xlabel="Species Richness",ylabel="Coefficient of Variation", size=(1800,1200), left_margin = [10mm 0mm], bottom_margin = [10mm 0mm], grid = false, xtickfontsize=18, ytickfontsize=18, xguidefontsize=18, yguidefontsize=18, legendfontsize=18)

	for i in [4, 6, 2, 3, 5]
		temp = zeros(n_grid, 2)
		temp[:,1] = length.(grid_species)
		temp[:,2] = var_mat[:,i]

		temp = temp[temp[:,2] .!= 0, :]
		temp = temp[sortperm(temp[:, 1]), :]
		sub = sort(sample(collect(1:size(temp, 1)), Int(floor(size(temp, 1)/3)), replace=false))
		Plots.plot!(temp[sub,1], temp[sub,2], label=group[i], lc=colors[i], lw=3, la= 0.6, ylims=(0,2))
	end

	png(p, @sprintf("plots/%s/rich_cv.png", metric))
end

function plot_cv_dens(n_grid, grid_species, var_mat, group, colors, metric)
	p = Plots.plot(xlabel="Coefficient of Variation",ylabel="Density", size=(1800,1200), left_margin = [10mm 0mm], right_margin = [5mm 0mm], bottom_margin = [10mm 0mm], grid = false, xtickfontsize=18, ytickfontsize=18, xguidefontsize=18, yguidefontsize=18, legendfontsize=18)

	for i in 2:6
		temp = filter(x -> x != 0, var_mat[:,i])
		dens = kde(temp)

		StatsPlots.plot!(dens, label=group[i], lw=3, lc=colors[i], xlims=(0, 1.5))
	end

	png(p, @sprintf("plots/%s/dens_cv.png", metric))
end


function plot_scaling_fd(n_grid, n_scal_spec, n_trait_comb, res_mat, scal_spec, metric, st_t_2, st_t_3, st_t_4, st_t_5, colors, group)
	for g in [4030, 1105, 9970, 18015, 20530]
		# if res_mat[end, g, 1] == 0 continue end  #only look at large communities
		l = length([x for x in res_mat[2:end, g, 1] if x != 0])
		a = findfirst(x->x!=0, res_mat[2:end, g, 1])

		if metric == "qhull" || metric == "qhull_norm"
			p = Plots.plot(xlabel="Species Richness Scaling",ylabel="Functional Richness", size=(1800,1200), left_margin = [10mm 0mm], bottom_margin = [10mm 0mm], grid = false, legend=:topleft, xtickfontsize=18, ytickfontsize=18, xguidefontsize=18, yguidefontsize=18, legendfontsize=18)
		else
			p = Plots.plot(xlabel="Species Richness Scaling",ylabel="Functional Divergence", size=(1800,1200), left_margin = [10mm 0mm], bottom_margin = [10mm 0mm], grid = false, legend=:topleft, xtickfontsize=18, ytickfontsize=18, xguidefontsize=18, yguidefontsize=18, legendfontsize=18)
		end
		
		vec_mean = zeros(l)

		
		for t in shuffle(collect(1:n_trait_comb))

			temp = zeros(n_scal_spec,2)
			temp[:,1] = res_mat[2:end, g, t]
			temp[:,2] = scal_spec

			temp = [temp[x,:] for x in collect(1:n_scal_spec) if temp[x,1] != 0]
			temp = reduce(hcat,temp)'

			if (1 <= t <= st_t_2-1) c = 2 end
			if (st_t_2 <= t <= st_t_3-1) c = 3 end
			if (st_t_3 <= t <= st_t_4-1) c = 4 end
			if (st_t_4 <= t <= st_t_5-1) c = 5 end
			if (st_t_5 <= t <= n_trait_comb) c = 6 end

			if (t==1 || t==st_t_2 || t==st_t_3 || t==st_t_4 || t==st_t_5) 
				p = Plots.plot!(temp[:,2], temp[:,1], lc=colors[c], alpha=0.4, label=group[c])
			else
				p = Plots.plot!(temp[:,2], temp[:,1], lc=colors[c], alpha=0.4, label="")
			end

			vec_mean = vec_mean + temp[:,1]./n_trait_comb
		end
		Plots.plot!(scal_spec[a:a+l-1], vec_mean, lc="blue", lw=3, label="mean")


		if isdir(@sprintf("plots/%s/scaling_fd", metric))==false mkdir(@sprintf("plots/%s/scaling_fd", metric)) end
		png(p, @sprintf("plots/%s/scaling_fd/grid_%d.png", metric, g))
	end
end


function plot_spec_fd_fit(n_grid, grid_species, res_mat, st_t_2, st_t_3, st_t_4, st_t_5, group, colors, metric)

	if metric == "qhull" || metric == "qhull_norm"
		p = Plots.plot(xlabel="Species Richness",ylabel="Functional Richness", size=(1800,1200), left_margin = [10mm 0mm], bottom_margin = [10mm 0mm], grid = false, legend=:bottomright, xtickfontsize=18, ytickfontsize=18, xguidefontsize=18, yguidefontsize=18, legendfontsize=18)
	else
		p = Plots.plot(xlabel="Species Richness",ylabel="Functional Divergence", size=(1800,1200), left_margin = [10mm 0mm], bottom_margin = [10mm 0mm], grid = false, legend=:bottomright, xtickfontsize=18, ytickfontsize=18, xguidefontsize=18, yguidefontsize=18, legendfontsize=18)
	end
	

	res_1 = []
	res_2 = []
	res_3 = []
	res_4 = []
	res_5 = []
 
	for g in 1:n_grid

		if res_mat[1, g, 1] == 0 continue end

		temp_1 = [(length(grid_species[g]), res_mat[1, g, x]) for x in collect(1:st_t_2-1)]
		temp_2 = [(length(grid_species[g]), res_mat[1, g, x]) for x in collect(st_t_2:st_t_3-1)]
		temp_3 = [(length(grid_species[g]), res_mat[1, g, x]) for x in collect(st_t_3:st_t_4-1)]
		temp_4 = [(length(grid_species[g]), res_mat[1, g, x]) for x in collect(st_t_4:st_t_5-1)]
		temp_5 = [(length(grid_species[g]), res_mat[1, g, x]) for x in collect(st_t_5:size(res_mat, 3))]

		append!(res_1, temp_1)
		append!(res_2, temp_2)
		append!(res_3, temp_3)
		append!(res_4, temp_4)
		append!(res_5, temp_5)
	end

	f1 = log_fit(getindex.(res_1,1), getindex.(res_1,2))
	f_1(x) = f1[1] + f1[2]*log(x)

	f2 = log_fit(getindex.(res_2,1), getindex.(res_2,2))
	f_2(x) = f2[1] + f2[2]*log(x)

	f3 = log_fit(getindex.(res_3,1), getindex.(res_3,2))
	f_3(x) = f3[1] + f3[2]*log(x)

	f4 = log_fit(getindex.(res_4,1), getindex.(res_4,2))
	f_4(x) = f4[1] + f4[2]*log(x)

	f5 = log_fit(getindex.(res_5,1), getindex.(res_5,2))
	f_5(x) = f5[1] + f5[2]*log(x)


	Plots.plot!(f_1, extrema(getindex.(res_1,1))..., lc=colors[2], lw=3, label=group[2])
	Plots.plot!(f_2, extrema(getindex.(res_2,1))..., lc=colors[3], lw=3, label=group[3])
	Plots.plot!(f_3, extrema(getindex.(res_3,1))..., lc=colors[4], lw=3, label=group[4])
	Plots.plot!(f_4, extrema(getindex.(res_4,1))..., lc=colors[5], lw=3, label=group[5])
	Plots.plot!(f_5, extrema(getindex.(res_5,1))..., lc=colors[6], lw=3, label=group[6])


	png(p, @sprintf("plots/%s/spec_fd_fit.png", metric))

end



function plot_fd_density_cell(n_grid, n_trait_comb, res_mat, st_t_2, st_t_3, st_t_4, st_t_5, group, colors, metric)
	# Plotting of density per grid cell
	for g in [4030, 1105, 9970, 18015, 20530]
		
		# vector for analysis of specific scaling over all trait combinations: column 1 = diversity metric, column 2 = group affiliation
		results = zeros(n_trait_comb, 2)

		# add group affiliation to each combination
		for t in 1:n_trait_comb
			
			results[t, 1] = res_mat[1,g,t]

			if (1 <= t <= st_t_2-1) results[t,2] = 1 end
			if (st_t_2 <= t <= st_t_3-1) results[t,2] = 2 end
			if (st_t_3 <= t <= st_t_4-1) results[t,2] = 3 end
			if (st_t_4 <= t <= st_t_5-1) results[t,2] = 4 end
			if (st_t_5 <= t <= n_trait_comb) results[t,2] = 5 end
		end

		# Calc Standard Deviation and mean for total
		res_std = std(results[:,1])
		res_mean = mean(results[:,1])

		# plot density all results
		dens = kde(results[:,1])

		if metric == "qhull" || metric == "qhull_norm"
			p = StatsPlots.plot(xlabel="Functional Richness",ylabel="Density", size=(1500,1000), left_margin = [10mm 0mm], bottom_margin = [10mm 0mm], grid = false, xtickfontsize=18, ytickfontsize=18, xguidefontsize=18, yguidefontsize=18, legendfontsize=18)
		else
			p = StatsPlots.plot(xlabel="Functional Divergence",ylabel="Density", size=(1500,1000), left_margin = [10mm 0mm], bottom_margin = [10mm 0mm], grid = false, xtickfontsize=18, ytickfontsize=18, xguidefontsize=18, yguidefontsize=18, legendfontsize=18)
		end
		
		p = StatsPlots.plot!(dens, label="density")

		# add total mean and std lines
		vline!([res_mean - res_std, res_mean, res_mean + res_std], label="mean +- std", lc="violet")
		
		# get highest point of density calculation: Needed for plotting height of points in groups
		h = maximum(dens.density)/40

		# create groups
		temp5 = [(results[x,1],h*2) for x in collect(1:n_trait_comb) if results[x,2] == 5]
		temp4 = [(results[x,1],h*6) for x in collect(1:n_trait_comb) if results[x,2] == 4]
		temp3 = [(results[x,1],h*10) for x in collect(1:n_trait_comb) if results[x,2] == 3]
		temp2 = [(results[x,1],h*14) for x in collect(1:n_trait_comb) if results[x,2] == 2]
		temp1 = [(results[x,1],h*18) for x in collect(1:n_trait_comb) if results[x,2] == 1]

		# Calc std and mean per group
		std_temp1 = std(getindex.(temp1,1))
		mean_temp1 = mean(getindex.(temp1,1))
		std_temp2 = std(getindex.(temp2,1))
		mean_temp2 = mean(getindex.(temp2,1))
		std_temp3 = std(getindex.(temp3,1))
		mean_temp3 = mean(getindex.(temp3,1))
		std_temp4 = std(getindex.(temp4,1))
		mean_temp4 = mean(getindex.(temp4,1))
		std_temp5 = std(getindex.(temp5,1))
		mean_temp5 = mean(getindex.(temp5,1))
		
		# plotting position of groups at different heights
		StatsPlots.scatter!(temp1, label=group[2], mc=colors[2], ms=6)	
		StatsPlots.scatter!(temp2, label=group[3], mc=colors[3], ms=6)
		StatsPlots.scatter!(temp3, label=group[4], mc=colors[4], ms=6)
		StatsPlots.scatter!(temp4, label=group[5], mc=colors[5], ms=6)
		StatsPlots.scatter!(temp5, label=group[6], mc=colors[6], ms=6)

		# plotting std and mean per group
		if length(temp1) != 0 StatsPlots.plot!([mean_temp1-std_temp1, mean_temp1-std_temp1], [h*17, h*19], lc=colors[2], lw=3, label=""); StatsPlots.plot!([mean_temp1, mean_temp1], [h*17, h*19], lc=colors[2], lw=3, label=""); StatsPlots.plot!([mean_temp1+std_temp1, mean_temp1+std_temp1], [h*17, h*19], lc=colors[2], lw=3, label="") end
		if length(temp2) != 0 StatsPlots.plot!([mean_temp2-std_temp2, mean_temp2-std_temp2], [h*13, h*15], lc=colors[3], lw=3, label=""); StatsPlots.plot!([mean_temp2, mean_temp2], [h*13, h*15], lc=colors[3], lw=3, label=""); StatsPlots.plot!([mean_temp2+std_temp2, mean_temp2+std_temp2], [h*13, h*15], lc=colors[3], lw=3, label="") end
		if length(temp3) != 0 StatsPlots.plot!([mean_temp3-std_temp3, mean_temp3-std_temp3], [h*9, h*11], lc=colors[4], lw=3, label=""); StatsPlots.plot!([mean_temp3, mean_temp3], [h*9, h*11], lc=colors[4], lw=3, label=""); StatsPlots.plot!([mean_temp3+std_temp3, mean_temp3+std_temp3], [h*9, h*11], lc=colors[4], lw=3, label="") end
		if length(temp4) != 0 StatsPlots.plot!([mean_temp4-std_temp4, mean_temp4-std_temp4], [h*5, h*7], lc=colors[5], lw=3, label=""); StatsPlots.plot!([mean_temp4, mean_temp4], [h*5, h*7], lc=colors[5], lw=3, label=""); StatsPlots.plot!([mean_temp4+std_temp4, mean_temp4+std_temp4], [h*5, h*7], lc=colors[5], lw=3, label="") end
		if length(temp5) != 0 StatsPlots.plot!([mean_temp5-std_temp5, mean_temp5-std_temp5], [h, h*3], lc=colors[6], lw=3, label=""); StatsPlots.plot!([mean_temp5, mean_temp5], [h, h*3], lc=colors[6], lw=3, label=""); StatsPlots.plot!([mean_temp5+std_temp5, mean_temp5+std_temp5], [h, h*3], lc=colors[6], lw=3, label="") end

		# saving plot
		if isdir(@sprintf("plots/%s/grid_idx", metric))==false mkdir(@sprintf("plots/%s/grid_idx", metric)) end
		png(p, @sprintf("plots/%s/grid_idx/grid_%d.png", metric, g))
	end
end

function trait_buckets(trait_names, trait_comb, st_t_2, st_t_3, st_t_4, st_t_5, group, colors)

	dic_t = Dict(zip(collect(1:18), trait_names))
	temp = []
	push!(temp, vcat(collect(1:18), reduce(vcat,trait_comb[1: st_t_2-1])))
	push!(temp, vcat(collect(1:18), reduce(vcat,trait_comb[st_t_2: st_t_3-1])))
	push!(temp, vcat(collect(1:18), reduce(vcat,trait_comb[st_t_3: st_t_4-1])))
	push!(temp, vcat(collect(1:18), reduce(vcat,trait_comb[st_t_4: st_t_5-1])))
	push!(temp, vcat(collect(1:18), reduce(vcat,trait_comb[st_t_5: end])))

	plots = []

	for i in 1:5
		dic = countmap(temp[i])
		res = zeros(18,2)
		res[:,1] = collect(keys(dic))
		res[:,2] = collect(values(dic)).-1
		res = res[sortperm(res[:,2], rev=true),:]
		p = Plots.plot(ylabel="Frequency", size=(1800,1200), left_margin = [10mm 0mm], bottom_margin = [15mm 0mm], grid = false, title=@sprintf("Group %s", group[i+1]), xtickfontsize=11, ytickfontsize=11, xguidefontsize=15, yguidefontsize=15)
		Plots.bar!(get.(Ref(dic_t), res[:,1], missing), res[:,2], xticks=(0.5:1:18.5, get.(Ref(dic_t), res[:,1], missing)), xrotation=60, label="", color=colors[i+1])
		png(p, @sprintf("plots/trait_usage_group%d", i+1))
		push!(plots, p)
	end

	p = Plots.plot(plots[1], plots[2], plots[3], plots[4], plots[5], layout = 5)
	png(p, "plots/trait_usage_groups")
end


#Main!

traits = DataFrame(CSV.File("data/Traits.csv"))
trait_names = names(traits[:, 2:end])
traits = Matrix(traits[:, 2:end])
n_traits = size(traits, 1)

# loading world grid 	sparse() <-> Array()
grid_raster = load("data/Grid_map.tif")
grid_raster = convert(Matrix{Float64}, grid_raster)
grid_raster_sp = sparse(grid_raster)
grid_raster_sp_indices = findall(!iszero, grid_raster_sp)
n_grid = size(grid_raster_sp_indices,1)
# n_grid = 21366

# Trait Combinations for comparison
trait_corr = DataFrame(CSV.File("data/trait_correlation_table.csv"))
combs = collect(combinations(1:18,5))
ncombs = length(combs)

# Calculation of specific specs for a trait combination
combs_spec = zeros(ncombs, 3)
combs_spec[:,1] = collect(1:ncombs)

for i in 1:ncombs
    temp = collect(combinations(combs[i],2))
    temp2 = [trait_corr[x[1],x[2]] for x in temp]
    combs_spec[i, 2] = sum(temp2)
    combs_spec[i, 3] = sum(abs.(temp2))
end

# Listing of specific trait combinations we are interested in
# one from each group
trait_comb1 = sort.(collect.(vec(collect(Iterators.product([1, 11],[2], [3, 4, 8, 9, 10, 12, 13, 17, 18], [7, 14, 15, 16], [5, 6])))))

# sorted by column sum
combs_spec = combs_spec[sortperm(combs_spec[:, 2]), :]
# most negatively correlated
trait_comb2 = [combs[Int(x)] for x in combs_spec[1:150, 1]]
# most positively correlated
trait_comb3 = [combs[Int(x)] for x in combs_spec[end-149:end, 1]]

# sorted by column abs sum
combs_spec = combs_spec[sortperm(combs_spec[:, 3]), :]
# least correlated
trait_comb4 = [combs[Int(x)] for x in combs_spec[1:150, 1]]
# greatest sum of absolute values in correlation
trait_comb5 = [combs[Int(x)] for x in combs_spec[end-149:end, 1]]

# concate them all together
trait_comb = vcat(trait_comb1, trait_comb2)
trait_comb = vcat(trait_comb, trait_comb3)
trait_comb = vcat(trait_comb, trait_comb4)
trait_comb = vcat(trait_comb, trait_comb5)

n_trait_comb = size(trait_comb,1)

# get starting indices for each group of traits
st_t_2 = 145
st_t_3 = 295
st_t_4 = 445
st_t_5 = 595

scal_spec = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 150, 200, 300, 500, 800, 1200, 1800]
n_scal_spec = size(scal_spec,1)

# Init grid_species: load from file

# Weird format -> transformations
grid_species = load("data/grid_species.jld")
grid_species = collect(values(grid_species))
grid_species = grid_species[1]


# Init result matrix to assemble files from metric_calculations.jl
grid_res_hull = zeros(n_scal_spec+1, n_grid, n_trait_comb)
grid_res_rao = zeros(n_scal_spec+1, n_grid, n_trait_comb)

grid_res_hull_norm = zeros(n_scal_spec+1, n_grid, n_trait_comb)
grid_res_rao_norm = zeros(n_scal_spec+1, n_grid, n_trait_comb)

grid_var_hull = zeros(n_grid, 6)
grid_var_rao = zeros(n_grid, 6)

grid_var_hull_norm = zeros(n_grid, 6)
grid_var_rao_norm = zeros(n_grid, 6)

# Definitions of plotting attributes
group = ["total", "one/cluster", "most neg", "most pos", "least", "greatest"]
colors = ["red", "green", "red", "purple", "orange", "turquoise"]
metric = ["qhull", "rao", "qhull_norm", "rao_norm"]

# Filling grid res matrices
fill_grid_mat(n_trait_comb, trait_comb, grid_res_hull, grid_res_hull_norm, grid_res_rao, grid_res_rao_norm, n_scal_spec)

# Filling grid var matrices
fill_grid_var(n_grid, n_scal_spec, n_trait_comb, grid_res_hull, grid_var_hull, st_t_2, st_t_3, st_t_4, st_t_5)
fill_grid_var(n_grid, n_scal_spec, n_trait_comb, grid_res_rao, grid_var_rao, st_t_2, st_t_3, st_t_4, st_t_5)

fill_grid_var(n_grid, n_scal_spec, n_trait_comb, grid_res_hull_norm, grid_var_hull_norm, st_t_2, st_t_3, st_t_4, st_t_5)
fill_grid_var(n_grid, n_scal_spec, n_trait_comb, grid_res_rao_norm, grid_var_rao_norm, st_t_2, st_t_3, st_t_4, st_t_5)


# Plots

# Creating trait bucket plots

trait_buckets(trait_names, trait_comb, st_t_2, st_t_3, st_t_4, st_t_5, group, colors)

# World map plots

# Plot Coeff Var as World Map
map_coeff_var(n_scal_spec, n_grid, scal_spec, grid_raster_sp_indices, grid_raster, grid_var_hull, metric[1])
map_coeff_var(n_scal_spec, n_grid, scal_spec, grid_raster_sp_indices, grid_raster, grid_var_rao, metric[2])

# Plot Normed Coeff Var per 
map_coeff_var(n_scal_spec, n_grid, scal_spec, grid_raster_sp_indices, grid_raster, grid_var_hull_norm, metric[3])
map_coeff_var(n_scal_spec, n_grid, scal_spec, grid_raster_sp_indices, grid_raster, grid_var_rao_norm, metric[4])

# Line Plots

# Plot #species vs Coeff Var in total species per grid
plot_rich_cv(n_grid, grid_species, grid_var_hull, group, colors, metric[1])
plot_rich_cv(n_grid, grid_species, grid_var_rao, group, colors, metric[2])

plot_rich_cv(n_grid, grid_species, grid_var_hull_norm, group, colors, metric[3])
plot_rich_cv(n_grid, grid_species, grid_var_rao_norm, group, colors, metric[4])

# Plot CV density
plot_cv_dens(n_grid, grid_species, grid_var_hull, group, colors, metric[1])
plot_cv_dens(n_grid, grid_species, grid_var_rao, group, colors, metric[2])

plot_cv_dens(n_grid, grid_species, grid_var_hull_norm, group, colors, metric[3])
plot_cv_dens(n_grid, grid_species, grid_var_rao_norm, group, colors, metric[4])

# Plot #species in scaling vs FD
plot_scaling_fd(n_grid, n_scal_spec, n_trait_comb, grid_res_hull, scal_spec, metric[1], st_t_2, st_t_3, st_t_4, st_t_5, colors, group)
plot_scaling_fd(n_grid, n_scal_spec, n_trait_comb, grid_res_rao, scal_spec, metric[2], st_t_2, st_t_3, st_t_4, st_t_5, colors, group)

plot_scaling_fd(n_grid, n_scal_spec, n_trait_comb, grid_res_hull_norm, scal_spec, metric[3], st_t_2, st_t_3, st_t_4, st_t_5, colors, group)
plot_scaling_fd(n_grid, n_scal_spec, n_trait_comb, grid_res_rao_norm, scal_spec, metric[4], st_t_2, st_t_3, st_t_4, st_t_5, colors, group)

# Plot #species vs FD LSQ
plot_spec_fd_fit(n_grid, grid_species, grid_res_hull, st_t_2, st_t_3, st_t_4, st_t_5, group, colors, metric[1])
plot_spec_fd_fit(n_grid, grid_species, grid_res_rao, st_t_2, st_t_3, st_t_4, st_t_5, group, colors, metric[2])

plot_spec_fd_fit(n_grid, grid_species, grid_res_hull_norm, st_t_2, st_t_3, st_t_4, st_t_5, group, colors, metric[3])
plot_spec_fd_fit(n_grid, grid_species, grid_res_rao_norm, st_t_2, st_t_3, st_t_4, st_t_5, group, colors, metric[4])

# Plot FD density in grid cell
plot_fd_density_cell(n_grid, n_trait_comb, grid_res_hull, st_t_2, st_t_3, st_t_4, st_t_5, group, colors, metric[1])
plot_fd_density_cell(n_grid, n_trait_comb, grid_res_rao, st_t_2, st_t_3, st_t_4, st_t_5, group, colors, metric[2])

plot_fd_density_cell(n_grid, n_trait_comb, grid_res_hull_norm, st_t_2, st_t_3, st_t_4, st_t_5, group, colors, metric[3])
plot_fd_density_cell(n_grid, n_trait_comb, grid_res_rao_norm, st_t_2, st_t_3, st_t_4, st_t_5, group, colors, metric[4])