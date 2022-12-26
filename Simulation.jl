const USE_GPU = false
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@static if USE_GPU
	@init_parallel_stencil(CUDA, Float64, 2)
else
	@init_parallel_stencil(Threads, Float64, 2)
end
using DataFrames
using Arrow
using Polyhedra
using QHull
using SparseArrays
using ImplicitGlobalGrid

using HypothesisTest
using Statistics



# PARALLELSTENCIL EXPLICIT INDICES WHEN CALLING @parallel (1:size(A,1), 1:size(A,3)) bc_y!(A)
# SEE IF PARALLELSTENCIL CAN RUN IN 1D !!!!



# orig_gf = convert.(Int64,DataFrame(Arrow.Table( "prep_files/lookup_tab.feather")))

# println(first(orig_gf, 5))

# traits = convert.(Float64, DataFrame(Arrow.Table("prep_files/trait_tab_norm.feather")))

# println(first(traits, 5))

# ll = DataFrame(Arrow.Table("prep_files/lat_long_mat.feather"))

# println(first(ll, 22))

# println(ll[22, [2, 1]])

# plots = sort(unique(orig_gf.plot_id[orig_gf.use_plot.==1]))
# println(first(plots, 5))

# max ammount of species we want to compare per plot
n_samp_species = 

# stepsize of species we want to compare at a time
scaling_species[] = 

# loading species from data
species_in =

# abundancy matrix (in case i even need one)
abund =

# loading world grid 	sparse() <-> Array()
grid_raster = load("prep_files/SR_all_trees_observed.tif")
convert!(Matrix{Float64}, grid_raster)
grid_raster_sp = sparse(grid_raster)
grid_raster_sp_indices = findall(!iszero, grid_raster_sp)
nindices = size(grid_raster_sp_indices,1)


grid_species = zeros(Float32, nindices, n_samp_species+2)
grid_species[:, 1] = getindex.(grid_raster_sp_indices, 1)
grid_species[:, 2] = getindex.(grid_raster_sp_indices, 2)

# this size depends on how many metric calculations we want to do???????
grid_results = zeros(Float32, nindices, size(scaling_species)+2)
grid_species[:, 1] = getindex.(grid_raster_sp_indices, 1)
grid_species[:, 2] = getindex.(grid_raster_sp_indices, 2)


# fill grid_species
@parallel (1:nindices) calc_null_mod(grid_species, species_in, n_samp_species, frec_tab)
# fill grid_results
@parallel (1:nindices) calc_metrics(grid_species, traits, trait_comb, results, abund, scaling_species)

# do some fancy analysis on these results

# map back to map and plot

@parallel_indices (ix) function calc_metrics(grid_species, traits, trait_comb, results, abund, scaling_species)
	results[somewhat undtil somewhere] = calc_rao()
	results[somewhere until somewhat different] = calc_chull()
end


# Goal of this function is to perform the random sampling for the null model and fill the matrix grid_species
@parallel_indices (ix) function calc_null_mod(grid_species, species_in, n_samp_species, frec_tab)
	grid_species[ix, 3:end] =  sample(species_in.id, frec_tab, n_samp_species; replace=false, ordered=true)
end



function calc_rao(tmp_traits::AbstractArray{Float64}, tmp_abund::Vector{Float64}, tmp_abund_mat::Matrix{Float64}, sub_traits::Matrix{Float64}, scaling::Vector{Float64})
	
	# subtraits Dimension 1 = amount of subtrait configs, Dimension 2 = subtrait configs

	# create a mock matrix of ones for non-abundance calc
	one_abund_mat = ones(Float64, size(sub_traits, 2), size(tmp_traits, 1))
	# initialize the outcomes
	outcomes = zeros(Float64, size(sub_traits, 1), 16 * size(scaling)) 
	# if we have at least one species
	if size(tmp_abund, 1) > 1
		# pairwise distance matrix
		d1 = pairwise(SqEuclidean(), transpose(tmp_traits), dims=2)

		for i in 1:size(scaling)
			for j in 1:size(sub_traits, 2)
				t_selec = sub_traits[j]
				index = (i-1)*16 + 1								#INDEX NOMAL CHECKE!!!!
				# get metrics, once with abundance once without
				outcomes[j, index:index+6] = calc_mpd(d1[t_selec,1:scaling[i]], tmp_abund_mat[t_selec, 1:scaling[i]])
				outcomes[j, index+7:index+13] = calc_mpd(d1[t_selec,1:scaling[i]], one_abund_mat)

				# get convex hull volume
				outcomes[j, index+14:index+15] = calc_chull(tmp_traits[1:scaling[i], t_selec], tmp_abund_mat[1:scaling[i], t_selec])   #might need to check dim of tmp_abund_mat!!!!!!!
			end
		end
		
	end
	return outcomes
end



# calculate the pairwise distance metrics
function calc_mpd(tmp_dist::Matrix{Float64}, tmp_abund_mat::Matrix{Float64})
	# initialize the output
	outvec = zeros(Float64, 7)

	# get weighted pairwise distances and summarize
	rq_mat = tmp_dist .* tmp_abund_mat
	outvec[1] = sum(rq_mat)
	outvec[2] = maximum(rq_mat[tril!(trues(size(rq_mat)), -1)])
	outvec[3] = minimum(rq_mat[tril!(trues(size(rq_mat)), -1)])
	outvec[4] = median(rq_mat[tril!(trues(size(rq_mat)), -1)])    

	# get minimum distances
	rq_mat[diagind(rq_mat)] .= Inf
	mv = minimum(rq_mat, dims=1)[1, :]
	outvec[5] = mean(mv)

	# get maximum distances
	rq_mat[diagind(rq_mat)] .= -Inf
	mv = maximum(rq_mat, dims=1)[1, :]
	outvec[6] = mean(mv)

	# now dendrogram
	rq_mat[diagind(rq_mat)] .= 0
	outvec[7] = sum_branch_dendro(sqrt.(tmp_dist) .* tmp_abund_mat)

	return outvec
end



function calc_chull(tmp_traits::Matrix{Float64}, tmp_abund_mat::Matrix{Float64})

	# initialize the output
	outvec = zeros(Float64, 2)
	ab_mat = tmp_traits .* tmp_abund_mat	#MAYBE?? Chamer das überhaupt so mache wenn mer kei Distance Matrix het?
	
	v = vrep(tmp_traits)
	v_ab = vrep(ab_mat)

	p_qhull = polyhedron(v, QHull.Library())
	# removevredundancy!(p_qhull)		#bruchts demfall ned wueki

	p_ab_hull = polyhedron(v_ab, QHull.Library())
	# removevredundancy!(p_ab_hull)	#bruchts demfall ned wueki

	outvec[1] = volume(p_hull)
	outvec[2] = volume(p_ab_hull)

	return outvec

end
