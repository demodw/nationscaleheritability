include("HeritabilityModule.jl")
using .HeritabilityModule
using ArgParse
using DataFrames
using CSV
using MatrixMarket
using BenchmarkTools

function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin
		"--fe"
			default = "static_covars.tsv"
			help = "File containing fixed effects in the Strata column"		
		"--Y"
			default = "Ysim.tsv"
			help = "File containing the phenotype in the Y column"
		"--covariance_matrices"
			default = "A.mtx"
			help = "Covariance matrix"
	end
	return parse_args(s)
end


function main()
	parsed_args = parse_commandline()

	# Read the phenotype
	Y = CSV.read(parsed_args["Y"], DataFrame; delim='\t')[:, "Y"];

	# Read covariance matricx
	mat_list = [MatrixMarket.mmread(parsed_args["covariance_matrices"])];

	# Fixed effects
	fixed_groups = CSV.read(parsed_args["fe"], DataFrame; delim='\t')[:, "Strata"];

	# Estimate h2
	# takes 54s
	@time h2 = calculate_cov(phe=Y, fixed=fixed_groups, mat_list=mat_list);
	print(h2)

	# Estimation can be sped up greatly by pre-processing the covariance matrix
	# Takes 30s after preprocessing
	long_list = preprocessing(mat_list=copy(mat_list));
	@time h2 = calculate_cov(phe=Y, fixed=fixed_groups, long=long_list);

end

main()
