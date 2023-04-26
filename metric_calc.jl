# # Used onlyu once in prepocessing to create Comm_map.csv !!!!
# traits = DataFrame(CSV.File("Traits.csv"))
# n_spec = size(traits, 1)

# # loading species for mapping from data
# spec_map = DataFrame(CSV.File("Comm_map_old.csv"))
# temp = zeros(size(spec_map,1))

# for i in 1:size(spec_map,1)
# 	name = spec_map[i, 3]
# 	temp[i] = findall(x -> x == name, traits.accepted_bin)[1]
# end

# spec_map[!,"spec_id"] = hoi
# CSV.write("/Users/jeffreyzweidler/Desktop/Thesis/trait_selection/Comm_map.csv", spec_map)

# Load in main thread
using DataFrames
using SparseArrays
using Statistics
using StatsBase
using Plots
using CSV
using JLD
using FileIO
using Images
using Distributed
using Combinatorics
using KernelDensity

addprocs(length(Sys.cpu_info())-1)

# Load in all threads
@everywhere begin
	using QHull
	using Noise
	using SharedArrays
	using Distances
	using Printf
	using SearchSortedNearest
	using ProgressMeter
end

@everywhere function calc_metrics(spec_traits, scal_spec, n_scal_spec, grid_idx)

	n_species_total = size(spec_traits, 1)

	if n_species_total > 1

		# filter duplicates and add noise to remaining, this helps Qhull calculations for points which are almost identical
		temp_spec_traits = unique(spec_traits,dims=1)
		add_gauss!(temp_spec_traits, 0.0001)
		n_species = size(temp_spec_traits,1)

		# shorten scaling to #species in grid cell + exact #species in grid cell
		scal_spec_temp = filter(x -> x < n_species, scal_spec)
		push!(scal_spec_temp, n_species)

		# call functions
		ret1 = calc_chull(scal_spec_temp, n_scal_spec, temp_spec_traits, size(spec_traits, 2), n_species_total, scal_spec)
		ret2 = calc_rao(scal_spec_temp, n_scal_spec, temp_spec_traits, n_species_total, scal_spec)
		# create return vector
		return vcat(ret1,ret2)
	end

	return zeros((n_scal_spec+1))#*2)
end


@everywhere function calc_rao(scal_spec_temp, n_scal_spec, temp_spec_traits, n_species_total, scal_spec)
	
	# initialize the results
	outvec = zeros(n_scal_spec+1) 

	# pairwise distance matrix
	n_scal_spec_temp = size(scal_spec_temp,1)
	d1 = pairwise(SqEuclidean(), transpose(temp_spec_traits))   

	for s in 1:n_scal_spec_temp
		if s == n_scal_spec_temp
			outvec[1] = sum(d1[1:scal_spec_temp[s], 1:scal_spec_temp[s]])*(1/scal_spec_temp[s]^2)   
		else 
			outvec[s+1] = sum(d1[1:scal_spec_temp[s], 1:scal_spec_temp[s]])*(1/scal_spec_temp[s]^2)  
		end
	end
	for s in n_scal_spec_temp:n_scal_spec
		if scal_spec[s] <= n_species_total
			outvec[s+1] = outvec[1]
		end
	end

	return outvec
end


@everywhere function calc_chull(scal_spec_temp, n_scal_spec, temp_spec_traits, n_traits, n_species_total, scal_spec)
	# initialize the output
	outvec = zeros(n_scal_spec+1)
	# we need more points then dims
	if size(temp_spec_traits, 1) > n_traits 
		n_scal_spec_temp = size(scal_spec_temp,1)
		# check where to start since QHull sometimes throws error if points are too close
		safe = false
		iter = findfirst(x -> x > n_traits, scal_spec_temp)
		ch = []
		hull = []

		while !safe && iter <= n_scal_spec_temp
			try
				ch = chull(temp_spec_traits[1:scal_spec_temp[iter], :])
				hull = ch.points[sort(ch.vertices),:]
				safe = true
			catch e
				iter+=1
			end
		end

		if !safe return outvec end

		outvec[iter+1] = ch.volume

		# We found a starting point
		for s in (iter+1):n_scal_spec_temp
			hull = vcat(hull, temp_spec_traits[scal_spec_temp[s-1]+1:scal_spec_temp[s], :])
			ch = chull(hull)
			hull = ch.points[sort(ch.vertices),:]
			if s == n_scal_spec_temp
				outvec[1] = ch.volume
			else
				outvec[s+1] = ch.volume
			end
		end
		for s in n_scal_spec_temp:n_scal_spec
			if scal_spec[s] <= n_species_total
				outvec[s+1] = outvec[1]
			end
		end
	end
	return outvec
end



# Main!

# Loading Trait Data
traits = DataFrame(CSV.File("data/Traits.csv"))
traits_norm = Matrix(mapcols(x -> (x.-minimum(x))./(maximum(x) - minimum(x)), traits[:,2:end]))
traits = Matrix(traits[:, 2:end])
n_traits = size(traits, 1)

# loading world grid 	sparse() <-> Array()
grid_raster = load("Grid_map.tif")
grid_raster = convert(Matrix{Float64}, grid_raster)
grid_raster_sp = sparse(grid_raster)
grid_raster_sp_indices = findall(!iszero, grid_raster_sp)
grid_nr = (getindex.(grid_raster_sp_indices,1).-1).*size(grid_raster,2) .+ getindex.(grid_raster_sp_indices, 2)
n_grid = size(grid_raster_sp_indices,1)
# n_grid = 21366


# Mapping of Community species to world grid

# Init grid_species, load from file or create file

if isfile("data/grid_species.jld")
	# Load it from file, weird format -> transformations
	grid_species_temp = JLD.load("data/grid_species.jld")
	grid_species = collect(values(grid_species_temp))
	grid_species = grid_species[1]
else 
	# create file for later re-runs
	# Loading Community
	comm_map = DataFrame(CSV.File("data/Comm_map.csv"))
	n_comm = size(comm_map,1)

	grid_species = [zeros(0) for _ in 1:n_grid]
	for i in 1:n_comm
		idx = indexin(comm_map[i,1], grid_nr)
		push!(grid_species[idx[1]], comm_map[i,4])
	end
	JLD.save("data/grid_species.jld", "grid_species", grid_species)
end


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
# most not correlated at all
trait_comb4 = [combs[Int(x)] for x in combs_spec[1:150, 1]]
# greatest sum of absolute values in correlation
trait_comb5 = [combs[Int(x)] for x in combs_spec[end-149:end, 1]]

# concat them all together
trait_comb = vcat(trait_comb1, trait_comb2)
trait_comb = vcat(trait_comb, trait_comb3)
trait_comb = vcat(trait_comb, trait_comb4)
trait_comb = vcat(trait_comb, trait_comb5)

trait_comb = unique(trait_comb)
n_trait_comb = size(trait_comb,1)

# scaling of species per grid cell
scal_spec = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 150, 200, 300, 500, 800, 1200, 1800]
n_scal_spec = size(scal_spec,1)


# run metric calculations
for t in collect(1:n_trait_comb)
	@printf("Iteration %d out of %d\n", t, n_trait_comb)

	temp = trait_comb[t]
	# if isfile(@sprintf("community_qhull_jld/grid_res_hull_%d_%d_%d_%d_%d.jld", temp[1], temp[2], temp[3], temp[4], temp[5])) continue end
	
	# create results matrix
	grid_res_hull = zeros(n_scal_spec+1, n_grid)
	grid_res_rao = zeros(n_scal_spec+1, n_grid)
	
	# in parallel run the metric calculations
	arrayOfTuples = [(t,y) for y in collect(1:n_grid)]
	out = @showprogress pmap(arrayOfTuples) do i
		calc_metrics(traits_norm[Int.(grid_species[i[2]]),trait_comb[i[1]]], scal_spec, n_scal_spec, i[2])
	end

	out = reduce(hcat, out)
	grid_res_hull = out[1:n_scal_spec+1,:]
	grid_res_rao = out[n_scal_spec+2:end,:]

	JLD.save(@sprintf("community_qhull_jld/grid_res_hull_%d_%d_%d_%d_%d.jld", temp[1], temp[2], temp[3], temp[4], temp[5]), "grid_res_hull", grid_res_hull)
	JLD.save(@sprintf("community_rao_jld/grid_res_rao_%d_%d_%d_%d_%d.jld", temp[1], temp[2], temp[3], temp[4], temp[5]), "grid_res_rao", grid_res_rao)
end