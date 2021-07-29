# Qst-Fst scratch

Here is where I'll work out the details of the Qst-Fst analyses I'm going to conduct. This script will be dynamic and change until I'm content with the results.

Do we observe evidence for spatially explicit selection on fatty acid methyl ester (FAME) profiles from wild sunflower oil? To answer this question we'll use the Qst-Fst analysis to determine if there is evidence for non-neutral genotypic variation among populations. Qst will be estimated using a mixed effects model with the form: 

Y ~ 1 + deme + (1|gr(Ind, A = G))

Qst = Va_deme / Va_deme + Va_Ind

The additive genetic variance for these data can be estimated by providing a variance-covariance matrix defining relatedness between individuals. The relatedness matrix is derived from either a pedigree of the breeding history, or, if a pedigree is not available, from a set of genetic loci. In this case we have 247 genotypes from sequenced with the Golden Gate platform, and can use those loci to create a genomic relatedness matrix (G). These sequences were originally genotyped by McAssey et al. (2016) and may contain a few loci with allele frequencies suggestive of selection. I went to remove these loci, but I had a hard time connecting the locus names from the drayad folder to the eight loci under selection mentioned in the text. I'll use all the sequences I have for now.

I have trait data for only 140 individuals and sequence data for 286. It's important to remember the values in each row of the trait matrix are mean values of families. The seeds were created from random mating of fathers within a population to a known mother, then matured, collected, pooled, and crushed to extract oil. The mean family oil value is reported in the trait data. The G matrix allows us to get the 1/4 of the additive genotypic variance from each mother. 

Here is an outline of the steps in the Qst estimation protocol.

Step 1. Subset genotypes
Step 2. Build a relatedness matrix from the SNPs
Step 3. Build the animal model & record goodness of fit tests
Step 4. Estimate Qst
Step 5. Estimate Fst on a per locus basis
Step 6. Simulate Qst values from each locus
Step 7. Compare simulated Qst to observed value.





(Fst will be estimated on a per locus basis, and the Qst-Fst comparison will use the method suggested by Whitlock & Guillame (2008). Brielfy, that method of perfroms the Qst-Fst comparison by generating an empirical distribution of Qst values estimated from the genetic loci, and then compares the simulated distribution to the observed Qst estimate. )

#### Read in the data

```{r, eval = TRUE}
# Read in the data
dat0 <- read.table("~/Dropbox/oil_envr/data/trait/oil_bioclim_ancestry_N153.txt",T,'\t')

# Read in the blocking effects
block <- read.table("~/Dropbox/oil_envr/data/block_effects/blocks.txt",T,'\t')

# Add the blocking effects
dat <- plyr::join(x=dat0, y=block, by = "ind_code", type = "left")

# Read in the snps as diploid, biallelic data
snps <- read.table("~/Dropbox/oil_envr/data/snps/snps.custom.genpop.txt",T,'\t')

```

View the overlap of genotypes with trait data.

```{R}
l <- list(snps$ind_code, dat$ind_code)
names(l) <- c("genotypes", "traits")
overlap <- gplots::venn(l)
```

There are 140 individuals with trait and sequence data, 146 with only sequence data, and 13 with only trait data.

### Subset individuals

We have genotype data for 286 individuals and trait data for 153, and 140 inds have both snps and trait data. Subset genotypes and phenotypes for the common set.

```{R, eval=T, message=F}
# Load tidyverse
library(tidyverse)
library(ggridges)

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

```

#### Build a relatedness matrix

Build a kinship matrix using the `Gmatrix` function from AGHmatrix. For details see `https://pubmed.ncbi.nlm.nih.gov/27902800/`. The first issue is filter the SNPs for the sites under selection according to McAseey. The Fst statistic is ideally built from neutral sites. Since, Ed already conducted an analyses with the SNPs to identify sites under selection, we should remove the 8 sites that show evidence of selection. See McAssey et al. 2016 for more details on the selection scan.

The second issue is I need to convert the data from it's current format to 0, 1, 2, which is what Gmatrix expects. This can be done with tidyverse.

```{r, eval = T}
# SNPs under selection to remove.
## It's not clear how to connect the results from McAssey to 
## the locus names as I have them.
## I attempted to find the loci with the HA genome browswer to no avail.

## For now, just use all the snps. 
## If you don't get significance from the Qst-Fst differences, you can
## talk to ed and see if he can help you figure out which 8 loci to 
## remove. Otherwise, continue

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
```

I'll compute the genomic relatedness matrix based on VanRaden (2009) (which seems to be more popular - not sure why yet). Figure out how this works. 

```{r, eval = T,  message=F}
library(AGHmatrix)

# Computing the additive relationship matrix based on VanRaden 2008
G <- Gmatrix(SNPmatrix = as.matrix(snps_agh), missingValue=-9, 
                  maf=0.05, method="VanRaden")

# Set the row names 
row.names(G) <- snps1$ind_code

# View the kinship matrix & save it
# The names need fixing up, (or turn them off & use population labels.)
heatmap(G)

# View a histogram of the self-variance components (the diagonal)
h0 <- ggplot(data.frame(V1=diag(G)), aes(x=V1))  + 
  geom_histogram(bins = 25, color="black", fill="grey", lwd=0.3) +
  theme_minimal() + 
  labs(x="Self-variances of A")

# View a histogram of the variance estimates.
h1 <- ggplot(data.frame(V1=G[upper.tri(G)]), aes(x=V1))  + 
  geom_histogram(bins = 25, color="black", fill="grey", lwd=0.3) +
  theme_minimal() + 
  labs(x="Covariances of A")

# Absolute value of covariance estimates
h2 <- ggplot(data.frame(V1=abs(G[upper.tri(G)])), aes(x=V1))  + 
  geom_histogram(bins = 25, color="black", fill="grey", lwd=0.3) +
  theme_minimal() + 
  labs(x="Absolute value of covariances of A")

# plot them together & save
pout <- cowplot::plot_grid(h0, h1, h2, ncol = 1, nrow = 3)
pout

#setwd("/Users/kf/Dropbox/oil_envr/results/kinship/figs")
#cowplot::save_plot(
#  filename = "Variance-Covariance_histograms_N140_v1.pdf",
#  plot = pout,
#  ncol = 1,
#  nrow = 1,
#  base_height = 6,
#  base_asp = 0.55,
#  base_width = NULL
#  )

# Combine the kinship matrix and the histograms into a single plot.
# You'll need to remake the kinship matrix with a ggplot2 function to combine
p0 <- cowplot::as_grob(gplots::heatmap.2(G, dendrogram = 'none', 
                        Rowv=TRUE, Colv=TRUE, trace='none',
                        key = FALSE))


# I think I need to make it in ggplot2 in order for it to work...
# This is relatively unimportant, work on it later

```

#### Fit the animal model

Now we can fit the animal model. Here are descriptions of the various model components:

- unsat: % unsaturated fatty acids methyl esters (FAMEs) measured from pressed seed oil. 
- block: blocking effects of greenhouse experiment. Levels = 1:4
- deme: genetic populations (i.e. 'demes') estimated from SNP data with structure. Levels = north, south
- ind_code: individual ID of mother plants
  - G: genomic relatedness matrix estimated from SNP loci.

I have a few options for fitting the parameters of the model as population or group effects. Clearly the ind_code will be fit as a grouping effect in order to estimate the additive variance. I have the choice of also fitting the block and deme effects as group or population effects. The rationale behind fitting block and deme as group effects is that the data naturally fit into groups: the data were blocked as they were collected and the data are drawn from two distinct genetic demes. However, I'm not certain it's necessary to do so. I believe deme has a fixed mean, and should be treated as a fixed effect, but the intercept of the northern and southern demes can also be estimated if fit as a grouping effects. Block makes sense to fit as a grouping effect, but I'm not particularly interested in the intercept of block, just that it's effects are accounted for in the model. 

In light of this uncertainty, I'll fit three or four iterations of the basic model and compare their goodness-of-fit with a few approaches.  

##### Fit 1: ind_code as grouping effect

```{r, eval=FALSE}
library(brms)
qst_mod0 <- brm(unsat ~ block + deme + (1|gr(ind_code, cov = G)),
                data = dat1,
                data2 = list(G = G),
                family = "gaussian",
                control = list(max_treedepth = 15, adapt_delta = 0.99),
                seed = 2398, chains = 4L, iter = 6000) %>%
            add_criterion("loo", reloo = TRUE)
                
save(qst_mod0, file="~/Dropbox/oil_envr/results/brms/models/qst/qst_mod0.brms")
```

##### Fit 2: block + ind_code as grouping effects

```{r, eval=FALSE}
qst_mod1 <- brm(unsat ~ (1|block) + deme + (1|gr(ind_code, cov = G)),
                data = dat1,
                data2 = list(G = G),
                family = "gaussian",
                control = list(max_treedepth = 15, adapt_delta = 0.99),
                seed = 2398, chains = 4L, iter = 6000) %>%
            add_criterion("loo", reloo = TRUE)
                
save(qst_mod1, file="~/Dropbox/oil_envr/results/brms/models/qst/qst_mod1.brms")
```

##### Fit 3: deme + block + ind_code as grouping effects

```{r, eval=FALSE}
qst_mod2 <- brm(unsat ~ (1|block) + (1|deme) + (1|gr(ind_code, cov = G)),
                data = dat1,
                data2 = list(G = G),
                family = "gaussian",
                control = list(max_treedepth = 15, adapt_delta = 0.99),
                seed = 2398, chains = 4L, iter = 6000) %>%
            add_criterion("loo", reloo = TRUE)

save(qst_mod2, file="~/Dropbox/oil_envr/results/brms/models/qst/qst_mod2.brms")
```

##### Fit 4: deme correlated to ind_code

Here's an alternative way to specify the model. This model accounts for the correlation between the additive genetic effects within individuals (i.e. their allele frequencies), and the fact that each deme is an evolutionary breeding population. I quite like this model. `ind_code` is treated as a grouping effect, while `deme` is a population effect. I should be able to get the conditional effects of deme from this model.

```{r, eval=F, message=F}
qst_mod3 <- brm(unsat ~ block + (1 + deme|gr(ind_code, cov = G)),
                data = dat1,
                data2 = list(G = G),
                family = "gaussian",
                control = list(max_treedepth = 15, adapt_delta = 0.99),
                seed = 2398, chains = 4L, iter = 6000) %>%
            add_criterion("loo", reloo = TRUE)
                
save(qst_mod3, file="~/Dropbox/oil_envr/results/brms/models/qst/qst_mod3.brms")
```

##### Fit 5. block as group and deme correlated to ind_code

```{r, eval=F, message=F}
qst_mod4 <- brm(unsat ~ (1|block) + (1 + deme|gr(ind_code, cov = G)),
                data = dat1,
                data2 = list(G = G),
                family = "gaussian",
                control = list(max_treedepth = 15, adapt_delta = 0.99),
                seed = 2398, chains = 4L, iter = 6000) %>%
            add_criterion("loo", reloo = TRUE)
                
save(qst_mod4, file="~/Dropbox/oil_envr/results/brms/models/qst/qst_mod4.brms")


```

##### Compare model fits

Read in the models if they aren't in memory.

```{R, eval = T, echo=F}
load("/media/two_tb_A/Dropbox/oil_envr/results/brms/models/qst/qst_mod0.brms")
load("/media/two_tb_A/Dropbox/oil_envr/results/brms/models/qst/qst_mod1.brms")
load("/media/two_tb_A/Dropbox/oil_envr/results/brms/models/qst/qst_mod2.brms")
load("/media/two_tb_A/Dropbox/oil_envr/results/brms/models/qst/qst_mod3.brms")
load("/media/two_tb_A/Dropbox/oil_envr/results/brms/models/qst/qst_mod4.brms")

# Compare best fitt with loo
loo_compare(qst_mod0, qst_mod1, qst_mod2, qst_mod3, qst_mod4)
```

##### Compare the Q_{ST} estimates from each model

```{r}
# Estimate qst
out <- pop_mod4 %>% 
  brms::posterior_samples(.vars=c(2:4,6), .funs=funs(.^2)) %>%
  select("b_Intercept","b_block", "sd_ind_code__Intercept", "sd_ind_code__demesouth", 
  "cor_ind_code__Intercept__demesouth", "sigma")
  
qst= out$sd_ind_code__demesouth / (out$sd_ind_code__demesouth + out$b_block + 2 * out$sd_ind_code__Intercept)
summary(qst)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.004464 0.409064 0.449374 0.434583 0.478720 0.641416

# hmmm, very different distribtuion, but the median is about the same.
```



```{r}
#Estimate qst with this model
out <- qst_mod1 %>% brms::posterior_samples(.vars=2:5, .funs=funs(.^2))

qst = out$sd_deme__Intercept / (out$sd_deme__Intercept + out$sd_block__Intercept + 2 * out$sd_ind_code__Intercept)

# Median is 0.47. This is lower by 0.1 than my first model, but it still is pretty high,
# higher than the Fst estimate of 0.21.
summary(qst)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.05048 0.36540 0.47206 0.47548 0.58497 0.88076 
```

#### Estimate Qst

Qst is the proportion of additive genetic variance among populations to the total phenotypic variation. Since we fit a bayesian model, we have posterior samples to work with and can estimate the median and confidence intervals of the median. Pull out the posterior sample and square them. These are the variance components 

```{r, eval=T}
out <- qst_mod %>%
  brms::posterior_samples(.vars=2:4, .funs=funs(.^2))
```

Calculate Qst from the variance component for the between population effects (`b_demesouth`) and the sum of the individual variance components from the model. We are mulitplying the additive genetic variance component by 2 since we can only ... (see iseki to answer this).

```{r, eval=T}
qst = out$b_demesouth / 
  (out$b_demesouth + out$b_block + 2 * out$sd_ind_code__Intercept)

qst = out$b_demesouth / 
  (out$b_demesouth + 2 * out$sd_ind_code__Intercept)
```
Print the median & interval estimates

```{r, eval=T}
tidybayes::median_qi(qst)
```

Plot the distribution

```{r, eval=T}
p1 <- ggplot(data.frame(Qst=qst), aes(x=Qst))  + 
  geom_histogram(bins = 25, color="black", fill="grey", lwd=0.3) +
  theme_minimal() + 
  labs(x=expression(paste(italic(Q)[ST])))
p1

# This isn't working.
# Keep it moving and come back to this if you want it.
#ggplot(data.frame(Qst=qst), aes(x=Qst, fill = factor(stat(quantile)))) +
#  stat_density_ridges(
#    geom = "density_ridges_gradient",
#    calc_ecdf = TRUE,
#    quantiles = c(0.025, 0.975)
#  ) +
#  scale_fill_manual(
#    name = "Probability", values = c("#FF0000A0", "#A0A0A0A0", "#0000FFA0"),
#    labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
#  )

#setwd("~/Dropbox/oil_envr/results/qst/figs")
#cowplot::save_plot(
#  filename = "qst_values.qst_mod0.brms_6x3.pdf",
#  plot = p1,
#  ncol = 1,
#  nrow = 1,
#  base_height = 3,
#  base_asp = 2,
#  base_width = NULL
#  )

```

### Qst-Fst test

Calculate the Qst-Fst statistic for the posterior samples. Plot them, if 95% of samples don't overlap zero, you can conclude confidently that the trait is under spatial selection.

```{r, eval = T}
# Calculate the statistic
qst_fst <- data.frame(qst = qst) %>%
  mutate(stat = qst - 0.21)

# View the Qst-Fst estimates
p2 <- ggplot(qst_fst, aes(x=stat))  + 
  geom_histogram(bins = 25, color="black", fill="grey", lwd=0.3) +
  theme_minimal() + 
  labs(x=expression(paste(italic(Q)[ST] - italic(F)[ST])))
p2

```

How can I put a probability number on the question, "What's the probability that the true value of Qst-Fst is less than 0.1?"

### Estimate h2

Narrow-sense heritabiity estimate.

```{r, eval=T}
Va = out$sd_ind_code__Intercept
Vp = out$sd_ind_code__Intercept + out$sigma

tidybayes::median_qi(Va/Vp)
```


## Fst estimate

Ed previously estimated Fst to be 0.21 using all markers. He also performed an analysis where he identified 8 SNPs that show evidence of selection. I will estimate Fst with all SNPs and with the 8 non-neutral sites removed. I need to find a bayesian method for estimating Weir & Cockerham's Fst.

```{r, eval=FALSE}
# Load the SNP data

# Estimate Fst as a two population model with multiple families.
# In this model deme = population and family = breeding population




```


## Compare Qst to Fst

How to do this? I think it's simply a matter of comparing the two values and if they overlap. You can take the different and the 95% CI's of the difference to, if the CI's of the differences don't overlap, you can take that as evidence of support.