include("HeritabilityModule.jl")
using .HeritabilityModule
using ArgParse
using DataFrames
using CSV
using MatrixMarket

function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin
		"--fe"
			default = "static_covars.tsv"
			help = "File containing fixed effects in the Strata column"		
		"--cv_chol"
			default = "A_chol.mtx"
			help = "List of covariance matrices in Cholesky form, seperated by comma"
        "--fileout"
			default = "Ysim.tsv"
			help = "Simulated phenotype"            
	end
	return parse_args(s)
end


function main()
	parsed_args = parse_commandline()

    c = split(parsed_args["cv_chol"], ",");
	mat_list = [];
	for i in c
		m = MatrixMarket.mmread(i);
		append!(mat_list, m);
	end

    fixed_groups = CSV.read(parsed_args["fe"], DataFrame; delim='\t');
    fixed_groups_u = unique(fixed_groups, "Strata");

    # Simulate a phenotype, assuming there is no sex effect
    fe_groups = sort!(fixed_groups_u, "Strata");
    fe_groups[:, "Age"] = rownumber.(eachrow(fe_groups));
    # Corresponds to a prevalence of approximately 5%
    fe_groups[:, "FE"] = -0.14617*fe_groups[:, "Age"];

    fixed_groups_w_fe = leftjoin(fixed_groups, fe_groups[:, ["Strata", "FE"]], on=:Strata);
    phe = make_pheno(
        var=[0.5], 
        fixed=fixed_groups_w_fe[!, "FE"],
        chol_list=mat_list
    );
    CSV.write(parsed_args["fileout"],DataFrame(Y=phe))
end

main()
