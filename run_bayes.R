library('tidyverse')
library('rstanarm')
library('betareg')
library('data.table')
library('rstan')
library('bayestestR')
library('optparse')
library('foreach')


parser <- OptionParser()
parser <- add_option(parser, 
	c('--indir'), 
	type="character",
	default="/data/workdata/222190/heritability/data/article_data/extended_all/simulated_phenotypes_A_estimates/")
parser <- add_option(parser, 
	c('--outfile'), 
	type="character",
	default="/data/workdata/222190/heritability/data/article_data/extended_all/sim_h2.tsv")
parser <- add_option(parser, 
	c('--ncpu'), 
	type="integer",
	default=39)
parser <- add_option(parser, 
	c('--rg'), 
	type="logical",
	action='store_true',
	default=FALSE)
arguments <- parse_args(parser)


reg <- function(data_vector, sm, rg=FALSE) {
	if (isTRUE(rg)) {
		data_vector[data_vector >= 0.999] <- 0.99
		data_vector[data_vector <= -0.999] <- -0.99
		data_vector <- (data_vector + 1) / 2
		# This is a prior focused on zero
		#prior <- rstanarm::student_t(7, 0, 1)
		# This is a prior focused on one, for bidirectional associations
		prior <- rstanarm::student_t(7, 2, 1)
	} else {
		prior <- rstanarm::student_t(7, -2, 1)
	}

	fit_dwe <- rstanarm::stan_betareg(y ~ 1, link='logit', link.phi='log', cores=1, refresh=0, 
		chains=4, iter=10000,
		data=data.frame(y=data_vector), prior_intercept=prior)
	foo2 <- sum(as.data.frame(rstan::get_sampler_params(fit_dwe$stanfit, inc_warmup=FALSE))$divergent__)
	X <- dplyr::tibble(data.frame(plogis(summary(fit_dwe, probs=c(0.025, 0.5, 0.975), pars=c('(Intercept)'), digits=3))))
	if (isTRUE(rg)) {
		X <- X %>%
			mutate(X2.5. = (X2.5.*2)-1, X50.=(X50.*2)-1, X97.5.=(X97.5.*2)-1)
	}
	
	X <- X[, c('X2.5.', 'X50.', 'X97.5.', 'n_eff', 'Rhat')] #dplyr::`%>%` select(-c(mean, mcse, sd))
	X$divergent <- foo2
	X$model <- 'dwe'
	X$Rhat <- qlogis(X$Rhat)

	X$rope1 <- bayestestR::rope(fit_dwe, range=c(-Inf, qlogis(0.01)), ci=0.95, ci_method='ETI')$ROPE_Percentage
	X$rope2 <- bayestestR::rope(fit_dwe, range=c(-Inf, qlogis(0.02)), ci=0.95, ci_method='HDI')$ROPE_Percentage
	X$mass1 <- sum(as.data.frame(fit_dwe)$`(Intercept)` <= qlogis(0.01)) / 20000
	X$mass2 <- sum(as.data.frame(fit_dwe)$`(Intercept)` <= qlogis(0.02)) / 20000

	return(X)
}

all.o <- list.files(path=arguments$indir, pattern='^estimate_*')
all_e <- NULL
for (i in all.o) {
	i_model <- stringr::str_split(tools::file_path_sans_ext(i), pattern='_')[[1]][3]
	f <- fread(file.path(arguments$indir, i))
	f$model <- i_model
	all_e <- bind_rows(all_e, f)
}
# Save it
fwrite(all_e, file.path(arguments$indir, 'all_estimates.tsv'), sep='\t')



cl <- parallel::makeCluster(arguments$ncpu)
doParallel::registerDoParallel(cl)

if (isTRUE(arguments$rg)) {
	cols <- c('A')
} else {
	cols <- c('A', 'D', 'Sib', 'Sp', 'F1', 'F2')
}

r <- NULL
r <- foreach (i=1:length(all.o), .combine=rbind, .packages=c('dplyr')) %dopar% {
	print(i)
	f <- all.o[i]
	df <- data.table::fread(file.path(arguments$indir, f))

	# Model
	model <- stringr::str_split(tools::file_path_sans_ext(f), pattern='_')[[1]][3]

	# Loop through phenotypes
	rtmp <- dplyr::tibble()
	for (p in unique(df$Phenotype)) {
		df_p <- df[df$Phenotype == p,]
		if (any(df_p == -1)) {
			next
		}
		# Regress
		for (c in cols) {
			if (c %in% colnames(df_p)) {
				tmp_vec <- round(df_p[[c]], 5)
				if (length(unique(tmp_vec)) == 1) {
					next
				}
				df_A <- try(reg(tmp_vec, m1, arguments$rg), silent=TRUE)
				if (class(df_A) == 'try-error') {
					next
				}
				df_A$term <- c
				df_A$phenotype <- p
				df_A$components <- model
				rtmp <- dplyr::bind_rows(rtmp, df_A)
			}
		}
	}
	rtmp
}

# Save it
r %>% write_tsv(arguments$outfile)

parallel::stopCluster(cl)