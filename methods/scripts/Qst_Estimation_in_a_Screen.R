# Qst Estimation code for quick-pasting into a screen.

# Load tidyverse
library(tidyverse)
library(ggridges)
library(brms)
library(AGHmatrix)

# multi-core loo to speed it up. May not be good for multiple runs
options(mc.cores = parallel::detectCores())

# Read in the data
dat0 <- read.table("~/Dropbox/oil_envr/data/trait/oil_bioclim_ancestry_N153.txt",T,'\t')

# Read in the blocking effects
block <- read.table("~/Dropbox/oil_envr/data/block_effects/blocks.txt",T,'\t')

# Add the blocking effects
dat <- plyr::join(x=dat0, y=block, by = "ind_code", type = "left")

# Read in the snps as diploid, biallelic data
snps <- read.table("~/Dropbox/oil_envr/data/snps/snps.custom.genpop.txt",T,'\t')

# Make the overlap object
l <- list(snps$ind_code, dat$ind_code)
names(l) <- c("genotypes", "traits")
overlap <- gplots::venn(l)

# List of inds to keep
keep <- attr(overlap,"intersections")$`genotypes:traits`

# Filter for snps
snps1 <- snps %>% 
  filter(ind_code %in% keep) %>%
  arrange(ind_code)

# filter for phenotypes
dat1 <- dat %>% 
  filter(ind_code %in% keep) %>%
  arrange(ind_code)

# Build a relatedness matrix

# Convert the genotypes into 012 format
snps_agh <- snps1 %>%
  select(-ind_code) %>%
  mutate_all(
    funs(case_when(
      . == 11 ~ 0,
      . == 12 ~ 1,
      . == 21 ~ 1,
      . == 22 ~ 2
    ))
  )

# Replace missing sites with -9
snps_agh[is.na(snps_agh)] <- -9

# Computing the additive relationship matrix based on VanRaden 2008
G <- Gmatrix(SNPmatrix = as.matrix(snps_agh), missingValue=-9, 
                  maf=0.05, method="VanRaden")

# Set the row names 
row.names(G) <- snps1$ind_code

# Fit the model you want