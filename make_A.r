library('tidyverse')
library('nadiv')
library('optparse')

option_list <- list(
  make_option(c('--ped'), type='character', help='Path to pedigree file',
    default='K:/c2_exchange/heritability_2023-04-20/pedigree.tsv'),
  make_option(c('--out'), type='character', help='Path to output matrix A',
  default='K:/c2_exchange/heritability_2023-04-20/')
  )
opt = parse_args(OptionParser(option_list=option_list))

fam = read_tsv(opt$ped) %>%
  mutate(
    missing_parent = ifelse(MOTHER_PERSON_ID == 0 | FATHER_PERSON_ID == 0, TRUE, FALSE),
    FATHER_PERSON_ID = ifelse(missing_parent == TRUE, 0, FATHER_PERSON_ID),
    MOTHER_PERSON_ID = ifelse(missing_parent == TRUE, 0, MOTHER_PERSON_ID),
    )

# Make sure pedigree is complete
fam2 <- prepPed(data.frame(fam %>% select(PERSON_ID, MOTHER_PERSON_ID, FATHER_PERSON_ID)))

# Make the A matrix. Can take quite some time
A <- nadiv::makeA(fam2)

# Keep only individuals from input
A_idx <- fam %>% select(PERSON_ID) %>% mutate(Y=1) %>%
  left_join(fam2 %>%
    mutate(
      Aidx = 1:n()
    ),
    by=c('PERSON_ID'='PERSON_ID')
  ) %>%
  filter(!is.na(Y)) %>%
  select(PERSON_ID, Aidx)
A <- A[A_idx$Aidx, A_idx$Aidx]
writeMM(A, file.path(opt$out, 'A.mtx'))

chol = Cholesky(A)
L = as(chol,"Matrix")
P = as(chol,"pMatrix")
c = t(P)%*%L
writeMM(c, file.path(opt$out, 'A_chol.mtx'))

