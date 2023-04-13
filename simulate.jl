include("/data/workdata/222190/heritability/scripts/FHJ/article_scripts/main.jl")

using Base.Threads
import Base.Threads.@threads

in_folder = "/data/workdata/222190/heritability/data/article_data/extended_all/pedigree_clusters/"
out_folder = "/data/workdata/222190/heritability/data/article_data/extended_all/simulated_phenotypes_ASCF1F2"

# I manually found these values zzzz, see below
prevs_to_sim = DataFrame(
    Prevalence = [0.001, 0.01, 0.1], 
    AgeH = [-0.618642, -0.343759, -0.14617],
    SexH = [0.0877188, -0.0352351, -0.085785]);
sim_df = repeat(prevs_to_sim, 300)
sim_df[!, "Sim"] = 1:size(sim_df, 1)
# Get variance components
var_components = rand(Uniform(0.01,0.5), size(sim_df, 1), 6);
var_components_norm = var_components ./ sum(var_components, dims=2);
sim_df[!, "A"] = var_components_norm[:, 1];
sim_df[!, "Sibling"] = var_components_norm[:, 2];
sim_df[!, "Spouse"] = var_components_norm[:, 3];
sim_df[!, "Fam1"] = var_components_norm[:, 4];
sim_df[!, "Fam2"] = var_components_norm[:, 5];
sim_df[!, "res"] = var_components_norm[:, 6];

# Read in existing
sim_df = DataFrame(CSV.File("/data/workdata/222190/heritability/data/article_data/extended_all/simulated_info_ASCF1F2.csv"))

# For each prevalence, simulate phenotypes with varying heritability
CSV.write("/data/workdata/222190/heritability/data/article_data/extended_all/simulated_info_ASCF1F2.csv", sim_df)
for p in 1:10
    p_folder = joinpath(in_folder, string(p))
    A_chol = MatrixMarket.mmread(joinpath(p_folder, "A_chol.mtx"));
    A = MatrixMarket.mmread(joinpath(p_folder, "A.mtx"));
    Sib_chol = MatrixMarket.mmread(joinpath(p_folder, "Sibling_chol.mtx"));
    Sib = MatrixMarket.mmread(joinpath(p_folder, "Sibling.mtx"));
    Spouse_chol = MatrixMarket.mmread(joinpath(p_folder, "Spouse_chol.mtx"));
    Spouse = MatrixMarket.mmread(joinpath(p_folder, "Spouse.mtx"));
    Fam1_chol = MatrixMarket.mmread(joinpath(p_folder, "Fam1_chol.mtx"));
    Fam1 = MatrixMarket.mmread(joinpath(p_folder, "Fam1.mtx"));
    Fam2_chol = MatrixMarket.mmread(joinpath(p_folder, "Fam2_chol.mtx"));
    Fam2 = MatrixMarket.mmread(joinpath(p_folder, "Fam2.mtx"));

    # Fixed effect groups
    fixed_groups = CSV.read(joinpath(p_folder, "static_covars.tsv"), DataFrame; delim='\t');
    fixed_groups_u = unique(fixed_groups, "Strata");

    #for i in reachrow(sim_df)
    Threads.@threads for i in eachrow(sim_df)
        # Simulate a phenotype
        male_groups = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"1"),:], "Strata");
        male_groups[:, "Age"] = rownumber.(eachrow(male_groups));
        male_groups[:, "FE"] = i[2]*male_groups[:, "Age"] .+ i[3];
    
        female_groups = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"2"),:], "Strata");
        female_groups[:, "Age"] = rownumber.(eachrow(female_groups));
        female_groups[:, "FE"] = i[2]*female_groups[:, "Age"];

        new_fe = vcat(male_groups, female_groups);
        fixed_groups_w_fe = leftjoin(fixed_groups, new_fe[:, ["Strata", "FE"]], on=:Strata);
        calc_pe = calculate_pre(fixed_groups_w_fe[!, "FE"])

        # OK, now lets simulate the phenotype
        #A_h2 = h_to_sim[sample(1:size(h_to_sim, 1))]
        #Sib_h2 = h_to_sim[sample(1:size(h_to_sim, 1))]
        #Spouse_h2 = h_to_sim[sample(1:size(h_to_sim, 1))]
        #Fam1_h2 = h_to_sim[sample(1:size(h_to_sim, 1))]
        #Fam2_h2 = h_to_sim[sample(1:size(h_to_sim, 1))]
        #res_h2 =  h_to_sim[sample(1:size(h_to_sim, 1))]
        #var_components = rand(Uniform(0.01,0.5), 6);
        #var_components_norm = var_components / sum(var_components)
        phe = make_pheno(var=i[["A", "Sibling", "Spouse", "Fam1", "Fam2"]], 
            fixed=fixed_groups_w_fe[!, "FE"],
            chol_list=[A_chol, Sib_chol, Spouse_chol, Fam1_chol, Fam2_chol]
        );

        # Lets save it
        exact_out = joinpath(out_folder, "prev_"*string(i["Prevalence"])*"_"*"Sim_"*string(i["Sim"]), string(p))
        mkpath(exact_out)
        CSV.write(joinpath(exact_out, "Y.tsv"),DataFrame(Y=phe))
    end
    #=
     #for h2 in h_to_sim
        for prev in eachrow(prevs_to_sim)
            # Simulate a phenotype
            male_groups = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"1"),:], "Strata");
            male_groups[:, "Age"] = rownumber.(eachrow(male_groups));
            male_groups[:, "FE"] = prev[2]*male_groups[:, "Age"] .+ prev[3];
        
            female_groups = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"2"),:], "Strata");
            female_groups[:, "Age"] = rownumber.(eachrow(female_groups));
            female_groups[:, "FE"] = prev[2]*female_groups[:, "Age"];

            new_fe = vcat(male_groups, female_groups);
            fixed_groups_w_fe = leftjoin(fixed_groups, new_fe[:, ["Strata", "FE"]], on=:Strata);
            calc_pe = calculate_pre(fixed_groups_w_fe[!, "FE"])

            # OK, now lets simulate the phenotype
            phe = make_pheno(var=[h2],fixed=fixed_groups_w_fe[!, "FE"],chol_list=[A_chol]);

            # Lets save it
            exact_out = joinpath(out_folder, "prev_"*string(prev[1])*"_"*"h2_"*string(h2), string(p))
            mkpath(exact_out)
            CSV.write(joinpath(exact_out, "Y.tsv"),DataFrame(Y=phe))
        end
    =#
end



# OK, now lets simulate genetic correlations for diseases with h2 > 0.1. Randomly select ~ 1000 pairs from sim_df
n_gencorr_sim = 1000
sim_gencorr_df_a = sim_df[sample(1:size(sim_df, 1), n_gencorr_sim), :]
rename!(sim_gencorr_df_a, names(sim_gencorr_df_a) .* "_A");
sim_gencorr_df_b = sim_df[sample(1:size(sim_df, 1), n_gencorr_sim), :]
rename!(sim_gencorr_df_b, names(sim_gencorr_df_b) .* "_B");
# Select genetic correlation
gen_cor_sim = DataFrame(rand(Uniform(-1, 1), n_gencorr_sim, 5), [:rgA, :rgSib, :rgSp, :rgFam1, :rgFam2]);
gen_cor_sim[:, "rgRes"] .= 0;
# Combine it
sim_genncor =  hcat(sim_gencorr_df_a, sim_gencorr_df_b, gen_cor_sim);
sim_genncor[:, "SimNum"] = 1:n_gencorr_sim

# Save it
out_folder = "/data/workdata/222190/heritability/data/article_data/extended_all/simulated_phenotypes_gc_ASCF1F2"
sim_genncor = DataFrame(CSV.File("/data/workdata/222190/heritability/data/article_data/extended_all/simulated_info_gc_ASCF1F2.csv"))

# For each prevalence, simulate phenotypes with varying heritability
CSV.write("/data/workdata/222190/heritability/data/article_data/extended_all/simulated_info_gc_ASCF1F2.csv", sim_genncor)

for p in 1:10
    p_folder = joinpath(in_folder, string(p))
    A_chol = MatrixMarket.mmread(joinpath(p_folder, "A_chol.mtx"));
    A = MatrixMarket.mmread(joinpath(p_folder, "A.mtx"));
    Sib_chol = MatrixMarket.mmread(joinpath(p_folder, "Sibling_chol.mtx"));
    Sib = MatrixMarket.mmread(joinpath(p_folder, "Sibling.mtx"));
    Spouse_chol = MatrixMarket.mmread(joinpath(p_folder, "Spouse_chol.mtx"));
    Spouse = MatrixMarket.mmread(joinpath(p_folder, "Spouse.mtx"));
    Fam1_chol = MatrixMarket.mmread(joinpath(p_folder, "Fam1_chol.mtx"));
    Fam1 = MatrixMarket.mmread(joinpath(p_folder, "Fam1.mtx"));
    Fam2_chol = MatrixMarket.mmread(joinpath(p_folder, "Fam2_chol.mtx"));
    Fam2 = MatrixMarket.mmread(joinpath(p_folder, "Fam2.mtx"));

    # Fixed effect groups
    fixed_groups = CSV.read(joinpath(p_folder, "static_covars.tsv"), DataFrame; delim='\t');
    fixed_groups_u = unique(fixed_groups, "Strata");

    p_id = sparse(I,size(A, 1), size(A,1 ));

    #for i in reachrow(sim_df)
    Threads.@threads for i in eachrow(sim_genncor)
        # Phenotype a
        male_groups_a = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"1"),:], "Strata");
        male_groups_a[:, "Age"] = rownumber.(eachrow(male_groups_a));
        male_groups_a[:, "FE"] = i["AgeH_A"]*male_groups_a[:, "Age"] .+ i["SexH_A"];
        female_groups_a = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"2"),:], "Strata");
        female_groups_a[:, "Age"] = rownumber.(eachrow(female_groups_a));
        female_groups_a[:, "FE"] = i["AgeH_A"]*female_groups_a[:, "Age"];
        new_fe_a = vcat(male_groups_a, female_groups_a);
        fixed_groups_w_fe_a = leftjoin(fixed_groups, new_fe_a[:, ["Strata", "FE"]], on=:Strata);
        calc_pe_a = calculate_pre(fixed_groups_w_fe_a[!, "FE"])

        # Phenotype b
        male_groups_b = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"1"),:], "Strata");
        male_groups_b[:, "Age"] = rownumber.(eachrow(male_groups_b));
        male_groups_b[:, "FE"] = i["AgeH_B"]*male_groups_b[:, "Age"] .+ i["SexH_B"];
        female_groups_b = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"2"),:], "Strata");
        female_groups_b[:, "Age"] = rownumber.(eachrow(female_groups_b));
        female_groups_b[:, "FE"] = i["AgeH_B"]*female_groups_b[:, "Age"];
        new_fe_b = vcat(male_groups_b, female_groups_b);
        fixed_groups_w_fe_b = leftjoin(fixed_groups, new_fe_b[:, ["Strata", "FE"]], on=:Strata);
        calc_pe_b = calculate_pre(fixed_groups_w_fe_b[!, "FE"])

        # OK, now lets simulate the phenotype
        phe1, phe2 = make_gc_phe(
            var1 = i[["A_A", "Sibling_A", "Spouse_A", "Fam1_A", "Fam2_A", "res_A"]],
            var2 = i[["A_B", "Sibling_B", "Spouse_B", "Fam1_B", "Fam2_B", "res_B"]],
            chol = [A_chol, Sib_chol, Spouse_chol, Fam1_chol, Fam2_chol, p_id],
            cor_list = i[["rgA", "rgSib", "rgSp", "rgFam1", "rgFam2", "rgRes"]],
            fixed1 = fixed_groups_w_fe_a[!, "FE"],
            fixed2 = fixed_groups_w_fe_b[!, "FE"]
        );

        # Lets save it
        o1 = joinpath(out_folder, "Sim_"*string(i["SimNum"])*"_A", string(p))
        mkpath(o1)
        CSV.write(joinpath(o1, "Y.tsv"),DataFrame(Y=phe1))
        o2 = joinpath(out_folder, "Sim_"*string(i["SimNum"])*"_B", string(p))
        mkpath(o2)
        CSV.write(joinpath(o2, "Y.tsv"),DataFrame(Y=phe2))        

    end
end

























# TODO: Simulate genetic correlations using combinations of same parameters







# OK, lets try
male_groups = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"1"),:]);
male_groups[:, "Age"] = rownumber.(eachrow(male_groups));
male_groups[:, "FE"] = -2.2025*male_groups[:, "Age"] .+ -0.0485;

female_groups = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"2"),:]);
female_groups[:, "Age"] = rownumber.(eachrow(female_groups));
female_groups[:, "FE"] = -2.2025*female_groups[:, "Age"];
new_fe = vcat(male_groups, female_groups);
calc_pe = calculate_pre(new_fe[!, "FE"]);

prev_df_test


r_fixed = rand(size(fixed_groups_u, 1))
fixed_groups_u[!, "FE"] = rand(Uniform(-4, 2), size(fixed_groups_u, 1));
fixed_groups_w_fe = leftjoin(fixed_groups, fixed_groups_u, on=:Strata);
pre=calculate_pre(fixed_groups_w_fe[!, "FE"])





# Choose prevalences
prevs_to_sim = DataFrame(
    Prevalence = [0.001, 0.01, 0.1, 0.2], 
    LowerU = [-1000, -100, -10, -8],
    UpperU = [1, 1, 1, 2]);
h_to_sim = [0.05, 0.1, 0.2, 0.4, 0.5, 0.8];


# OK, lets simulate from three prevalences
info_list=[]	
for i in 1:100
    fixed_groups_u[!, "FE"] = rand(Normal(-2, 1), size(fixed_groups_u, 1));
    fixed_groups_w_fe = leftjoin(fixed_groups, fixed_groups_u, on=:Strata);
    push!(info_list,calculate_pre(fixed_groups_w_fe[!, "FE"]))
end
mean(info_list)





info_list=[]	
for i in 1:100
    fixed_groups_u[!, "FE"] = rand(Uniform(-3, -2), size(fixed_groups_u, 1));
    fixed_groups_w_fe = leftjoin(fixed_groups, fixed_groups_u, on=:Strata);
    push!(info_list,calculate_pre(fixed_groups_w_fe[!, "FE"]))
    phe=make_pheno(var=[he],fixed=fixed_groups_w_fe[!, "FE"],chol_list=[A_chol])
end
mean(info_list)

sum(cdf.(Normal(),fixed_groups_w_fe[!, "FE"]))

fixed_fixed=sample([-3,-1.73],size(A,1));
CSV.write(out_folder*"fixed.csv",DataFrame(f=fixed_fixed))

info_list=[]	
for (i,he) in enumerate([i*5/100 for i=0:20])
    #fixed,pre=realistic_fixed(pre=true)
    fixed=fixed_fixed;

    pre=calculate_pre(fixed);
    phe=make_pheno(var=[he],fixed=fixed,chol_list=[A_chol])
    CSV.write(out_folder*string(i)*".csv",DataFrame(y=phe))

    push!(info_list,[i,he,pre])

end		
info=TMA(info_list)
info=DataFrame(index=info[1],A=info[2],pre=info[3])
CSV.write(out_folder*"info.csv",info)




A_chol = MatrixMarket.mmread(in_folder*"A_chol.mtx");
A = MatrixMarket.mmread(in_folder*"A.mtx");
Sib_chol = MatrixMarket.mmread(in_folder*"Sibling_chol.mtx");
Sib = MatrixMarket.mmread(in_folder*"Sibling.mtx");
Spouse_chol = MatrixMarket.mmread(in_folder*"Spouse_chol.mtx");
Spouse = MatrixMarket.mmread(in_folder*"Spouse.mtx");
Fam1_chol = MatrixMarket.mmread(in_folder*"Fam1_chol.mtx");
Fam1 = MatrixMarket.mmread(in_folder*"Fam1.mtx");
Fam2_chol = MatrixMarket.mmread(in_folder*"Fam2_chol.mtx");
Fam2 = MatrixMarket.mmread(in_folder*"Fam2.mtx");


# Fixed effect groups
fixed_groups = CSV.read(joinpath(in_folder, "static_covars.tsv"), DataFrame; delim='\t');
fixed_groups_u = unique(fixed_groups, "Strata");


prev_df_test = DataFrame(age=[], sex=[], prev=[]);
for i in 1:1000
    #age_hyper = -1
    #sex_hyper = -0.1
    age_hyper = rand(Uniform(-3, 0), 1)[1];
    sex_hyper = rand(Uniform(-0.1, 0.1), 1)[1];
    
    male_groups = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"1"),:], "Strata");
    male_groups[:, "Age"] = rownumber.(eachrow(male_groups));
    male_groups[:, "FE"] = age_hyper*male_groups[:, "Age"] .+ sex_hyper;

    female_groups = sort!(fixed_groups_u[startswith.(fixed_groups_u[!, :Strata],"2"),:], "Strata");
    female_groups[:, "Age"] = rownumber.(eachrow(female_groups));
    female_groups[:, "FE"] = age_hyper*female_groups[:, "Age"];

    new_fe = vcat(male_groups, female_groups);
    fixed_groups_w_fe = leftjoin(fixed_groups, new_fe[:, ["Strata", "FE"]], on=:Strata);
    calc_pe = calculate_pre(fixed_groups_w_fe[!, "FE"]);
    push!(prev_df_test, [age_hyper, sex_hyper, calc_pe]);
end
sort!(prev_df_test, "sex")
sort!(prev_df_test, "prev")

# OK, lets select parameters for some specific prevalences
prev_df_test[(prev_df_test.prev .> 0.001) .& (prev_df_test.prev .< 0.0011), :]
prev_df_test[(prev_df_test.prev .> 0.01) .& (prev_df_test.prev .< 0.012), :]
prev_df_test[(prev_df_test.prev .> 0.1) .& (prev_df_test.prev .< 0.12), :]
prev_df_test[(prev_df_test.prev .> 0.2) .& (prev_df_test.prev .< 0.21), :]