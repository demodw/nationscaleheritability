# nationscaleheritability
Calculation of heritability on nation-scale  data


# How to run ...

1. Get a pedigree and some strata for each individual

2. Calculate the A matrix for the pedigree: `make_A.R`

3. Simulate a phenotype: `julia HeritabilitySimulate.jl`

4. Run the estimation `julia HeritabilityEstimate.jl`

5. Done


