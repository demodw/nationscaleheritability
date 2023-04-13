__precompile__()
module HeritabilityModule
	using DataFrames
	using Roots
	using DataStructures
	using SparseArrays
	using LinearAlgebra
	using CSV
	using Distributions
	using MatrixMarket
	using Cubature
	using Optim

	export calculate_cov, preprocessing

	# Functions for simulating data
	"""
	Simulate phenotype for univariate analysis of heritability

	var: list of variances to simulate

	chol_list: list of covariances matrixes in their Cholesky form

	fixed: list of fixed effects for eacn individual
	
	output: vector of simulated phenotype
	"""
	function make_pheno(;var,chol_list,fixed)
		if !(length(var) == length(chol_list))
			println("error")
			return
			end
		d=Normal()
		
		ep = zeros(length(fixed))
		l = zeros(length(fixed))
	
		for i in 1:length(chol_list)
			ep. = rand(d,length(fixed))
			l.+ = chol_list[i]*ep*sqrt(var[i])
		end
		
		ep. = rand(d,length(fixed))
		l.+ = sqrt(1-sum(var))*sparse(I,length(fixed),length(fixed))*ep
		
		l = l.+fixed
		phe = 1*(l.>0)
		return phe
	end

	function make_gc_phe(;var1,var2,chol,cor_list,fixed1,fixed2)
		#all lists contain residual variance. So, the last enry in chol_list should be Id
		N = length(fixed1)*2
		cov = [x*sqrt(y)*sqrt(z) for (x,y,z) in zip(cor_list,var1,var2)]
		l = zeros(N)
		ep = zeros(N)
	
		for i in 1:length(chol)
			covar_mat = [[var1[i],cov[i]] [cov[i],var2[i]]]
			covar_chol = cholesky(covar_mat).L
			l.+ = kron(sparse(covar_chol),chol[i])*rand(Normal(),N)
		end
	
		l.+ = vcat(fixed1,fixed2)
		phe = 1*(l.>0)
		phe1 = phe[1:Int64(N/2)]
		phe2 = phe[Int64(N/2)+1:N]
		return (phe1, phe2)
	end
	

	
	"""
	Format data to long format
	
	mat_list: A list of the variance matrices (without the identity matrix)

	output: preprocessed long format data
	"""
    function preprocessing(;mat_list)
		#create a sparse matrix and add identity matrix?
		non_zeros = sparse(I, size(mat_list[1],1), size(mat_list[1],1)) 
		for mat in mat_list
			non_zeros += mat
		end	
		
		mat_list = append!(mat_list,[sparse(I,size(mat_list[1],1),size(mat_list[1],1))])
		
		col = col_values(non_zeros)
		long = DataFrame(row=non_zeros.rowval,col=col)
		long = long[long.row.>=long.col,:]
		for (index,mat) in enumerate(mat_list)
			col = col_values(mat)
			colname = string(index)
			tmp = DataFrame(row=mat.rowval,col=col,value=round.(mat.nzval,digits=5))
			rename!(tmp,:value=>colname)
			long = leftjoin(long,tmp,on=[:row,:col])
		end
		
		#remember to set missing =0
		long = coalesce.(long,0) 
		long = hcat(long[:,[:row,:col]],Tuple.(eachrow(select(long,Not([:row,:col])))))
		rename!(long,:x1=>:mat_val)

		return long
	end


    """
	Helper function for preprocessing
                
    mat: a variance matrix from mat_list
    
    output: ??

	"""
    function col_values(mat)
		col = zeros(length(mat.nzval))
		column = 1
		colptr = mat.colptr
		for j in 1:length(col)
			if colptr[column+1] > j
				col[j] = column
			else
				col[j] = column+1
				column += 1
			end
		end
		return col
	end


    """
	compute the fixed effects for each subgroup J using equation (4) from the supplement
                
    phe: the phenotype as a vector of 0 and 1
    fixed: the fixed effects groups as a vector of strings. length(fixed) equals length(phe).

    output: estimated fixed effects

	"""
	function calculate_fixed(phe,fixed)
		float_fixed = zeros(Float64, length(fixed))
		for f in unique(fixed)
			pre = sum(phe[fixed .==f])/length(fixed[fixed .==f])
			float_fixed[fixed .==f] .= quantile.(Normal(),pre)	
		end

		return float_fixed
	end


	# The functions below are all helper functions for calculating the loss (calculate integrals, expected values, densities)
	"""
	gives E(Y1Y2) = expected value of correlation between two phenotypes
	
	c: as defined defined in the appendix to the supplement
	fixed: f1 and f2, the fixed effects, as defined in the appendix to the supplement

	output: E(Y1,Y2)
	"""
	function expected_value(;c,fixed)
		if fixed[2] == -Inf
			return 0
		end

		#find maximum of density using multidimensional h-adaptive integration
		#is E(Y1Y2) from appendix
		return hcubature(y->density(y,fix=fixed,c=c),[0,0],[1,1],abstol=1e-5)[1] 
	end


    """
	helper function for the gradient (= partial derivative with respect to p_i of E(Y1Y2)), as defined in the appendix to the supplement
        is this the differentiated normal pdf?

    x: value at which gradient is computed    
	"""
    function phi_diff(x)
		d = -x*exp(-x^2/2)/sqrt(2*pi)
		return d
	end


    """
	helper function to function calculate_gradient_integral_pair
        gives the gradient at a specific x (= partial derivative with respect to p_i (= prevalence of phenotype in group i??) of E(Y1Y2)), see appendix to the supplement
	x: value at which gradient is computed
    f: f1 and f2, the fixed effects, as defined in the appendix to the supplement
    c: as defined defined in the appendix to the supplement

    output: gradient
	"""
    function density_grad(x; f=fixed,c=c) 
		r = (f[1]*c*(x[2]-1)-f[2]*(x[2]-1)-x[2])/((1-c^2)^(3/2)*(x[2]-1)) *
		pdf(Normal(),-f[2]+x[2]/(1-x[2])) *
		phi_diff(-(f[1]+c*(-f[2]+x[2]/(1-x[2])))/sqrt(1-c^2)+x[1]/(1-x[1])) * #is this the normal pdf differentiated wrt. something?
		1/(1-x[1])^2/(1-x[2])^2 
	end


    """
	helper function for calculating E(Y1, Y2)
        gives the density of E(Y1,Y2)

	y: x1 and x2 in in E(Y1,Y2) as defined in the appendix to the supplement
    fix: f1 and f2, the fixed effects, as defined in the appendix to the supplement
    c: as defined defined in the appendix to the supplement

    output: density value
	"""
	function density(y;fix=fixed,c=c)
	        r = pdf(Normal(),-fix[2]+y[2]/(1-y[2])) *
	         pdf(Normal(),-(fix[1]+c*(-fix[2]+y[2]/(1-y[2])))/(1-c^2)^0.5+y[1]/(1-y[1])) *
	             1/(1-y[1])^2/(1-y[2])^2 
	        return r
	end


    """
	helper function to calculate loss
        
    mv: ??
    p: ??
    i: ??
    gc: true if genetic covariance is calculated

    output: ??
	"""
	function calculate_factor(mv,p,i,gc)
		if gc == true
			return -mv[i]
		else
			factor = -sum(p*mv[i].-p.*mv)/sum(p)^2
			return factor
		end
	end


    """
	helper function to calculate loss:
        computes the gradient (= partial derivative with respect to p_i of E(Y1Y2)) using multidimensional h-adaptive integration
        
    c: as defined defined in the appendix to the supplement
    fixed: f1 and f2, the fixed effects, as defined in the appendix to the supplement

    output: ??
	"""
    function calculate_gradient_integral_pair(;fixed,c)
		if fixed[2] == -Inf
			return 0
		end
		return hcubature(y->density_grad(y,f=fixed,c=c),[0,0],[1,1],abstol=1e-5)[1] 
	end


    """
	calculate loss as defined in equation (5) from the supplement:
                
    F: phenotype??
    G: fixed effects??
    parameters:
    to_cal: 
    gc: true if genetic covariance is calculated

    output: the loss as defined in equation (5) from the supplement

	"""
	function calculate_loss(F,G,parameters,to_cal;gc=false)
		var_par = parameters	
		if gc == false
			# reparametrization of variance components to theta
			var_par = [x/sum(parameters) for x in parameters] 
		else
			var_par=[var_par...,0]
		end
		
		to_cal = transform(to_cal,[:mat_val]=> ByRow(x-> dot(x,var_par))=> :c)
		to_cal = transform(to_cal,[:fixed_low,:fixed_high,:c]=>ByRow((fl,fh,c)->expected_value(c=c,fixed=(fh,fl)))=>:E)
		to_cal = transform(to_cal, [:prod,:gsize,:E]=>ByRow((p,n,e)->p*(1-e)^2+(n-p)*e^2)=>:loss)

		if G != nothing		
			for i=1:length(parameters)
				to_cal = transform(to_cal, [:mat_val]=>ByRow((mv)->calculate_factor(mv,parameters,i,gc))=>:factor)
				to_cal = transform(to_cal, [:fixed_low,:fixed_high,:c,:factor]=>ByRow((fl,fh,c,fac)->fac*calculate_gradient_integral_pair(fixed=(fh,fl),c=c))=>:gradient_pair)
				to_cal = transform(to_cal, [:gradient_pair,:E,:prod,:gsize]=>ByRow((gp,E,p,n)->p*(gp*2*(E-1))+(n-p)*(gp*2*(E-0)))=>:gradient_group)
				G[i] = sum(to_cal.gradient_group)
			end
		end

		if F != nothing
			loss = sum(to_cal.loss)
			return loss
		end
	end


	"""
	phe: the phenotype as a vector of 0 and 1

	fixed: the fixed effects groups as a vector of strings. length(fixed) equals length(phe)

	mat_list: A list of the variance matrices (without the identity matrix). Does not have to be provided if long is provided.

	long: A DataFrame that contains information about the variance matrices. Can be calculated via the function preprocessing(mat_list). If calculate_cov  is to be called on several phenotypes with the same mat_list, it is recommended to calculate long beforehand and pass this instead of mat_list.

	phe2: A second phenoptype if correlations are to be calculated.

	var1: If correlations are to be calculated, we need the variance component of each phenotype. If these have already been estimated they can be passed as two  lists. var1 takes the variance components for phe1. length(var1)=length(mat_list)+1. var1 is the output of calculate_cov(phe=phe1, ... ).

	var2: Variance components of phe2.

	output: If phe2 is not passed, the output is a list of variance components corresponding to the matrices in mat_list. length(calculate_cov(...phe2=nothing...))=length(mat_list)+1. If phe2 is passed, the output is a list of correlations. The correlations in the residual matrix are not currently calculated as these seem to be highly imprecise. So length(calculate_cov(...phe2=phe2...))=length(mat_list).
	"""
	function calculate_cov(;phe,fixed,mat_list=nothing,long=nothing,phe2=nothing, var1=nothing,var2=nothing)
		gc = false #what does gc indicate??
		N = length(fixed)

		# create long form
		if isnothing(long)
			long = preprocessing(mat_list=copy(mat_list))
		end

		# remove rows if they are the same person
		long = long[long.row.!=long.col,:]

		if !isnothing(phe2) # if we want to compute covariance between two phenotypes 
			gc = true
			
			if isnothing(var1) 
				var1 = calculate_cov(phe=phe,fixed=fixed,mat_list=mat_list,long=long) #calculate variance components phenotype 1 
			end
		
			if isnothing(var2) 
				var2 = calculate_cov(phe=phe2,fixed=fixed,mat_list=mat_list,long=long) #calculate variance components phenotype 2 
			end
			
			long2 = transform(long,[:row,:col]=>ByRow((x,y)-> (y,x))=>[:row,:col])
			long = vcat(long,long2)
			long = transform(long,[:row]=> ByRow(x->x+N)=>[:row])

			fixed1 = fixed.*"1"
			fixed2 = fixed.*"2"
			fixed = vcat(fixed1,fixed2)

			phe = vcat(phe,phe2)
		end

		fixed = calculate_fixed(phe,fixed) # estimate fixed effects
		fixed = DataFrame(fixed=fixed,rowid=1:length(fixed))

		long = leftjoin(long,fixed,on=:row => :rowid)
		rename!(long,:fixed=>:fixed_row)
		long = leftjoin(long,fixed,on=:col => :rowid)
		rename!(long,:fixed=>:fixed_col)
		long = transform(long,[:fixed_row,:fixed_col]=> ByRow((x,y)->if x<y; x; else; y ; end)=> :fixed_low)	
		long = transform(long,[:fixed_row,:fixed_col]=> ByRow((x,y)->if x>y; x; else; y ; end)=> :fixed_high)	
	 	long = select(long,Not([:fixed_row,:fixed_col]))

		
		##find unique entries to be calculated
		phe = DataFrame(phe=phe,rowid=1:length(phe))
		long = leftjoin(long,phe,on=:row => :rowid)
		rename!(long,:phe=>:phe_row)	
		long = leftjoin(long,phe,on=:col => :rowid)
		rename!(long,:phe=>:phe_col)
		
		to_cal = combine(groupby(long,[:fixed_low,:fixed_high,:mat_val]),[:phe_row,:phe_col]=>dot=>:prod,[:phe_row]=>length=>:gsize)

		#number of matrices
		nm = length(first(long.mat_val))

		if gc == false
			lower = [1.0 for x in 1:nm]
			upper = [100.0 for x in 1:nm]
			initial = [50.0 for x in 1:nm]
			res = optimize(Optim.only_fg!((F,G,x)->calculate_loss(F,G,x,to_cal,gc=false)),lower,upper,initial,Fminbox(LBFGS()), Optim.Options(x_tol=1e-12, f_tol=1e-12, g_tol=1e-12))
		else
			lower = [-sqrt(var1[i])*sqrt(var2[i]) for i in 1:(nm-1)] #lower bound for covariance 
			upper = [sqrt(var1[i])*sqrt(var2[i]) for i in 1:(nm-1)] #upper bound for covariance
			initial = [rand(Uniform(lower[i], upper[i]), 1)[1] for i in 1:(nm-1)] #initialize optimizer with random uniform starting values between upper and lower bound
			res = optimize(Optim.only_fg!((F,G,x)->calculate_loss(F,G,x,to_cal,gc=true)),lower,upper,initial,Fminbox(LBFGS()))
		end

		if gc == false
			esti = Optim.minimizer(res)./(sum(Optim.minimizer(res)))
			return esti
		else
			cov=Optim.minimizer(res) # minimize loss function to compute the genetic covariances
			cor=[cov[i]/sqrt(var1[i]*var2[i]) for i in 1:length(cov)] #compute genetic correlations using equation (6) from the supplement
			return cor, cov, var1, var2
		end
	end
end