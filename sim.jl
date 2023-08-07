include("/home/xn0/kdrev/Metode/Challenge/xn0/heritability/src/nationscaleheritability-main/HeritabilityModule.jl")
using .HeritabilityModule
using MatrixMarket
using LinearAlgebra
using CSV
using DataFrames
using Distributions
using SparseArrays
using Base.Threads
import Base.Threads.@threads

function load_data(dir)
    # structural data
    Sibling_chol = MatrixMarket.mmread(joinpath(dir, "Sibling_chol.mtx"))
    Sibling = MatrixMarket.mmread(joinpath(dir, "Sibling.mtx"))

    Spouse_chol = MatrixMarket.mmread(joinpath(dir, "Spouse_chol.mtx"))
    Spouse = MatrixMarket.mmread(joinpath(dir, "Spouse.mtx"))

    Fam1_chol = MatrixMarket.mmread(joinpath(dir, "Fam1_chol.mtx"))
    Fam1 = MatrixMarket.mmread(joinpath(dir, "Fam1.mtx"))

    Fam2_chol = MatrixMarket.mmread(joinpath(dir, "Fam2_chol.mtx"))
    Fam2 = MatrixMarket.mmread(joinpath(dir, "Fam2.mtx"))

    A_chol = MatrixMarket.mmread(joinpath(dir, "A_chol.mtx"))
    A = MatrixMarket.mmread(joinpath(dir, "A.mtx"))

    fixed_groups = CSV.read(joinpath(dir, "static_covars.tsv"), DataFrame, delim="\t")

    return([A, A_chol, Sibling, Sibling_chol, Spouse, Spouse_chol, Fam1, Fam1_chol, Fam2, Fam2_chol, fixed_groups])
end

function simulate(chollist, fixed_groups)
    fixed_groups_u = unique(fixed_groups, "Strata")
    
    if lenght(chollist) == 1
        ageF = rand(Uniform(-0.025, 0.01))
        ageM = rand(Uniform(-0.025, 0.01))
        offset = rand(Uniform(-3.5, -0.2))j
    else
        ageF = rand(Uniform(-0.01, 0.01))
        ageM = rand(Uniform(-0.01, 0.01))
        offset = rand(Uniform(-1, -0.05))
    end

    # Get variance components
    var_components = rand(Uniform(0.05,0.5), size(1, 1), length(chollist)+1)
    VAR = var_components ./ sum(var_components, dims=2)

    male_groups = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"1"),:], "Strata")
    male_groups[:, "Age"] = rownumber.(eachrow(male_groups))
    male_groups[:, "FE"] = ageM*male_groups[:, "Age"] .+ offset
    female_groups = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"2"),:], "Strata")
    female_groups[:, "Age"] = rownumber.(eachrow(female_groups))
    female_groups[:, "FE"] = ageF*female_groups[:, "Age"] .+ offset
    new_fe = vcat(male_groups, female_groups)
    fixed_groups_w_fe = leftjoin(fixed_groups, new_fe[:, ["Strata", "FE"]], on=:Strata)

    fe = fixed_groups_w_fe[!, "FE"]
    fg = fixed_groups_w_fe[!, "Strata"]

    phe = make_pheno(var=VAR[1, 1:length(chollist)], 
        fixed=fe,
        chol_list=chollist
    )
    prev = sum(phe)/size(phe)[1]
    return([phe, prev, VAR])
end

function Acholesky(X)
    X = X + Diagonal(X)*0.00001
    F = cholesky(X, perm=1:size(X)[1])
    L = sparse(F.L)
    return(L)
end

function main()

    A, A_chol, Sibling, Sibling_chol, Spouse, Spouse_chol, FamD1, FamD1_chol, FamD2, FamD2_chol, fixed_groups = load_data("/home/xn0/kdrev/c2_exchange/heritability_2023-04-20/")
    FamA1 = MatrixMarket.mmread("/home/xn0/kdrev/Metode/Challenge/xn0/heritability/out/mat/Famawj_1.mtx")
    FamA2 = MatrixMarket.mmread("/home/xn0/kdrev/Metode/Challenge/xn0/heritability/out/mat/Famawj_2.mtx")
    FamA1_chol = Acholesky(FamA1)
    FamA2_chol = Acholesky(FamA2)
    fg = fixed_groups[:, "Strata"]
    
    # Simualtion 1 - only genetic factors
    chollist=[A_chol]
    Threads.@threads for ii in range(1, 101)
        
        phe, prev, VAR = simulate(chollist, fixed_groups)

        h2 = calculate_cov(phe=phe, fixed=fg, mat_list=[A])
        
        out = DataFrame(
            prev = [prev, prev],
            v_A = [h2[1], VAR[1, 1]],
        )

        CSV.write(joinpath("/home/xn0/kdrev/Metode/Challenge/xn0/heritability/out/sim/_1/", string(ii)), out)

    # Simualtion 2 - only genetic factors for siulation but estimation with environmental factors
    Threads.@threads for ii in range(1, 101)
    
        phe, prev, VAR = simulate(chollist, fixed_groups)

        h2 = calculate_cov(phe=phe, fixed=fg, mat_list=[A, Sibling, Spouse, FamD1, FamD2])
        
        out = DataFrame(
            prev = [prev, prev],
            v_A = [h2[1], VAR[1, 1]]
        )

        CSV.write(joinpath("/home/xn0/kdrev/Metode/Challenge/xn0/heritability/out/sim/_2/", string(ii)), out)


    # Simulation 3 - genetic + environment for simualtion and estimation
    chollist=[A_chol, Sibling_chol, Spouse_chol, FamA1_chol, FamA2_chol]

    Threads.@threads for ii in range(81, 120)
        
        phe, prev, VAR = simulate(chollist, fixed_groups)

        h2 = calculate_cov(phe=phe, fixed=fg, mat_list=[A, Sibling, Spouse, FamD1, FamD2])
        
        out = DataFrame(
            prev = [prev, prev],
            v_A = [h2[1], VAR[1, 1]],
            v_sib = [h2[2], VAR[1, 2]],
            v_spo = [h2[3], VAR[1, 3]],
            v_fam1 = [h2[4], VAR[1, 4]],
            v_fam2 = [h2[5], VAR[1, 5]]
        )

        CSV.write(joinpath("/home/xn0/kdrev/Metode/Challenge/xn0/heritability/out/sim/_3/", string(ii)), out)

    end
end

main()
