# nationscaleheritability
Calculation of heritability in nation-scale  data

# TODO

* Add anonymised pedigree (Joleen)
* Write documentation for how to estimate heritability / genetic correlations (Joleen / David)


# How to run ...

1. Get a pedigree and some strata for each individual

2. Calculate the A matrix for the pedigree: `make_A.R`

3. Simulate a phenotype: `julia HeritabilitySimulate.jl`

4. Run the estimation `julia HeritabilityEstimate.jl`

5. Done


