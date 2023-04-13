include("HeritabilityModule.jl")
using .HeritabilityModule
using ArgParse
using DataFrames
using CSV
using MatrixMarket

function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin
		"--phenotype_dir"
			default = "/data/workdata/222190/heritability/data/article_data/extended_all/ICD10"
			help = "Directory containing phenotypes in format <phenotype>/"
		"--phenotype_list"
			help = "List of phenotypes to calculate from folder"
			default = "O20"
			#default = "UTI,Lymphatic_Disorder,Cataract,Lupus_Erythematosus,Metabolic_Syndrome_X,Overweight,Eustachian_Tube_Disorder"
		"--precomputed_dir"
			default = "/data/workdata/222190/heritability/data/article_data/extended_all/pedigree_clusters"
			help = "Directory containing covariates and pre-computed matrices in long format"
		"--model"
			default = "ASDCF1F2"
			help = "Model. Default is to only include A. Choices: A, AS, AC, AF, ASC, ASF, ACF, ASCF"
		"--outname"
			default = "out3.csv"
			help = "Output filename"
	end
	return parse_args(s)
end


function main()
	parsed_args = parse_commandline()

	# If the output exists, gracefully exit
	if isfile(parsed_args["outname"])
		return
	end

	sub_populations = readdir(parsed_args["precomputed_dir"]);

	# Could probably be made smarter but here we are
	df = DataFrame(Phenotype=String[], Subpopulation=String[]);
	mat_list = String[]
	if contains(parsed_args["model"], "A")
		df = hcat(df, DataFrame(A=Float64[]));
		push!(mat_list, "A")
	end
	if contains(parsed_args["model"], "S")
		df = hcat(df, DataFrame(Sib=Float64[]));
		push!(mat_list, "Sibling")
	end
	if contains(parsed_args["model"], "D")
		df = hcat(df, DataFrame(D=Float64[]));
		push!(mat_list, "D")
	end
	if contains(parsed_args["model"], "C")
		df = hcat(df, DataFrame(Sp=Float64[]));
		push!(mat_list, "Spouse")
	end
	if contains(parsed_args["model"], "F1")
		df = hcat(df, DataFrame(F1=Float64[]));
		push!(mat_list, "Fam1")
	end
	if contains(parsed_args["model"], "F2")
		df = hcat(df, DataFrame(F2=Float64[]));		
		push!(mat_list, "Fam2")
	end
	if contains(parsed_args["model"], "Fa")
		df = hcat(df, DataFrame(Fa=Float64[]));		
		push!(mat_list, "NewFam")
	end

	# Add individual 
	df = hcat(df, DataFrame(E=Float64[]));
	push!(mat_list, "long_mat.csv");

	pheno_list = split(parsed_args["phenotype_list"], ",");

	# Loop-de-loop
	for s in sub_populations
		# Matrix to load
		long_mat_filename = joinpath(parsed_args["precomputed_dir"], s, join(mat_list, "_"));
			
		# Load sub-population files
		fixed_groups = CSV.read(joinpath(parsed_args["precomputed_dir"], s, "static_covars.tsv"), DataFrame; delim='\t')[:,"Strata"];
		long_list = CSV.read(long_mat_filename, DataFrame, types=Dict(:row=>Int32, :col=>Int32, :mat_val=>String));
		long_list = transform(long_list, [:mat_val] => ByRow(x->[Tuple(parse.(Float32, split(chop(x, head=1, tail=1), ',')))])=>[:mat_val]);

		# Loop over phenotypes
		for pheno in pheno_list
			# Load phenotype
			Y = CSV.read(joinpath(parsed_args["phenotype_dir"], pheno, s, "Y.tsv"), DataFrame; delim='\t')[:, "Y"];
			h2 = try
				 calculate_cov(phe=Y, fixed=fixed_groups, mat_list=nothing, long=long_list);
			catch e
				fill(-1, length(first(long_list.mat_val)));
			end

			push!(df, vcat([pheno, s], h2));
		end
	end

	# Save it
	CSV.write(parsed_args["outname"], df)
end

main()
