#########################################################################################################
#########################################################################################################
## Compare my mediational gcomp code to CMAverse given multiple dependent mediators 
#########################################################################################################
#########################################################################################################
rm(list=ls())
library(data.table)
library(survey)
library(mice)
library(boot)
library(MASS)
library(snow)
library(miceadds)
library(mitools)
library(ggplot2)
library(parallel)
in_dir <- 'C:/Users/ncgra/Dropbox/Penn/repos/home_value'
repo <- 'C:/Users/ncgra/Dropbox/Penn/repos/multmed_gcomp'
sim_dir <- 'C:/Users/ncgra/Documents/cmaverse_sims/' ## Directory for saving lots of tiny intermediate results (these are batched out automatically, so will not take up a lot of memory).
source(paste0(repo,"/functions.R"))
library(CMAverse)

############################################################################################################################################################################################
## Compare my code vs. CMAverse in the context of two ORDERED, DEPENDENT mediators.
## See if my separate effects via each mediator add up to CMAverse joint mediated effect through the set of mediators. 
## CMAverse example: https://bs1125.github.io/CMAverse/articles/post_exposure_confounding.html
############################################################################################################################################################################################
expit <- function(x) exp(x)/(1+exp(x))
n <- 10000
C1 <- rnorm(n, mean = 1, sd = 0.1)
C2 <- rbinom(n, 1, 0.6)
A <- rbinom(n, 1, expit(0.2 + 0.5*C1 + 0.1*C2))
L <- rnorm(n, mean = 1 + A - C1 - 0.5*C2, sd = 0.5)
M1 <- rbinom(n, 1, expit(1 + 2*A - L + 1.5*C1 + 0.8*C2))
M2 <- rbinom(n, 1, expit(1 + 20*A - L + 1.1*C1 + 2.3*C2))
Y <- rnorm(n, -3 - 0.4*A - 1.2*M1 + 0.5*A*M1 - 0.9*M2 + 1.2*A*M2 - 0.5*L + 0.3*C1 - 0.6*C2)
data <- data.frame(A, M1, M2, Y, C1, C2, L)
data <- as.data.table(data)

## CMAverse package
res_gformula <- cmest(data = data, model = "gformula", outcome = "Y", exposure = "A",
                      mediator = c("M1",'M2'), basec = c("C1", "C2"), postc = "L", EMint = TRUE,
                      mreg = list("logistic","logistic"), yreg = "linear", postcreg = list("linear"),
                      astar = 0, a = 1, mval = list(0,0), 
                      estimation = "imputation", inference = "bootstrap", nboot = 100)
summary(res_gformula)

## My gformula package
data <- as.data.table(data)
path_cb <- fread(paste0(repo, '/cmaverse_dag2.csv'))
## Fit all models.
data[, pweight_1 := 1]
data[, id := 1:.N]
inputs <- multmed_glms(path_cb, df=data)
for(obj in names(inputs)) assign(obj, inputs[[obj]])
## G-formula setup parameters
outcome_vars <- c('Y')
interaction_vars <- 'A'
treatment_var <- 'A'
duration_vars=NULL
ever_vars=NULL
cumcount_vars=NULL
dummy_vars <- NULL
ordinal_vars <- NULL
ordinal_levels <- NULL
ordinal_refs <- NULL
## Parameters for simulations.
total_sim_count <- 1 ## Number of gformula simulations (via parametric bootstraps on fitted GLMS).
mc_reps <- 50 ## Number of data replicates (to remove MC error in individual response prediction within simulations; only if you have binomial/multinomial mediators/outcomes).
sim_dir <- 'C:/Users/ncgra/Documents/cmaverse_sims/' ## Directory for saving lots of tiny intermediate results (these are batched out automatically, so will not take up a lot of memory).
use_mean_betas <- T ## Set to TRUE if you only want quick point estimates of all effects (with total_sim_count==1) or set to FALSE for lots of simulations with parametric bootstraps to get standard errors on effects.
use_np_boot <- F ## Set to TRUE if you want non-parametric bootstrap. This should be combined with use_mean_betas <- T. A parametric bootstrap would be use_mean_betas <- F & use_np_boot <- F.
do_decomp <- TRUE ## Set to TRUE to calculate isolated expected values necessary for effect decomposition via multmed_decomp() below, or just natural/intervention courses for total effects.
decomp_type <- '4way' ## Effect decomposition type ('2way' for NDE, NDE or '4way' for CDE, PAIs, PIEs).
parallelize <- F ## Whether to parallelize gformula simulations or run serially.
windows <- TRUE ## OS if parallelizing (TRUE=Windows or FALSE=Linux).
processors <- 4 ## How many parallel threads if parallelizing?
## Specify natural courses. Default would just be no changes to the data, but you could subset or something here if you wanted.
natural_courses <- list(
  function(d) {
    return(d) 
  }
)
names(natural_courses) <- 'natural_course'
## Specify intervention rules and natural rules as a list of functions to apply to the data. 
intervention_rules <- list(
  function(d,...) {
    d[, A := 0]
    return(d)
  },
  function(d,...) {
    d[, A := 1]
    return(d)
  }
)
names(intervention_rules) <- c('control_course','treatment_course')
## Specify effect pathways for decomposition (mutually exclusive and exhaustive, should include all variables in "update_vars" in path_cb).
decomp_paths <- list(
  list(paths=list('L'='L',
                  'M1'='M1',
                  'M2'='M2'),
       outcomes='Y')
)
names(decomp_paths) <- 'treatment_course'
control_course <- 'control_course'
## Estimate g-formulas
multmed_gformula(data=data,
                 tc_vars=tc_vars,
                 models=models,
                 sim_dir=sim_dir,
                 natural_courses=natural_courses,
                 intervention_rules=intervention_rules,
                 path_cb=path_cb,
                 treatment_var=treatment_var,
                 control_course=control_course,
                 total_sim_count = total_sim_count,
                 mc_replicates = mc_reps,
                 decomp_paths = decomp_paths,
                 parallelize=parallelize,
                 processors=processors,
                 windows=windows,
                 mean_betas = use_mean_betas,
                 decomp = do_decomp,
                 decomp_type = decomp_type,
                 dummy_vars=dummy_vars,
                 ordinal_levels=ordinal_levels,
                 ordinal_vars=ordinal_vars,
                 ordinal_refs=ordinal_refs,
                 duration_vars=duration_vars,
                 ever_vars=ever_vars,
                 cumcount_vars=cumcount_vars)
## Post-process all g-formula results
lapply(names(intervention_rules), post_process_course, total_sims=total_sim_count, sim_dir=sim_dir, decomp_paths=decomp_paths)
## Plot results
o <- outcome_vars
effect_table <- multmed_decomp(intervention_course = names(decomp_paths)[1],
                               compare_course = control_course,
                               decomp_paths = decomp_paths,
                               total_sims = total_sim_count,
                               sim_dir = sim_dir,
                               decomp_type = decomp_type,
                               outcome_var = o)
draw_table <- effect_table[[2]]
effect_table <- effect_table[[1]]
effect_table
## Compare fits to CMAverse for 2 M and 1 L
effect_table[, type := 'medgformula']
total_pai <- effect_table[effect %in% c('PAI_M1','PAI_M2'), list(mean=sum(mean)), by=c('type')] ## Sum my separate PAI via M1/M2 to see if it adds up to CMAverse joint PAI via M1/M2. 
total_pai[, effect := 'PAI_M']
cma_effects <- data.table(effect=names(res_gformula$effect.pe),
                          mean=res_gformula$effect.pe,
                          se=res_gformula$effect.se,
                          upper=res_gformula$effect.ci.high,
                          lower=res_gformula$effect.ci.low)
cma_intref <- rnorm(10000, cma_effects[effect=='rintref',mean], cma_effects[effect=='rintref',se]) ## Just for quick comparison of summed effects
cma_intmed <- rnorm(10000, cma_effects[effect=='rintmed',mean], cma_effects[effect=='rintmed',se]) ## Just for quick comparison of summed effects
cma_pai <- cma_intref + cma_intmed
cma_pai <- data.table(effect='PAI_M',
                      mean=mean(cma_pai),
                      upper=quantile(cma_pai,0.975),
                      lower=quantile(cma_pai,0.025))
cma_effects <- rbind(cma_effects, cma_pai, fill=T)
cma_effects[effect=='rpnie', effect := 'PIE_M']
cma_effects[effect=='te', effect := 'ATE']
cma_effects[effect=='cde', effect := 'CDE']
cma_effects <- cma_effects[effect %in% c('ATE','CDE','PAI_M','PIE_M'),]
cma_effects[, type := 'CMAverse']
all_effects <- rbindlist(list(effect_table,cma_effects,total_pai),fill=T)

## Plot my effects and CMAverse effects together; sum of separate PAI via M1/M2 adds up to CMAverse joint mediated PAI via M1/M2. 
## CDE will be different because I also treat all "L" as mediators rather than post-treatment confounders, whereas CMAverse treats those pathways as part of CDE. 
ggplot(data=all_effects) + 
  geom_hline(yintercept=0,color='black') +
  geom_linerange(aes(x=effect,
                     ymin=lower,
                     ymax=upper,
                     group=type),
                 size=2, color='black', position = position_dodge(width=0.4)) +
  geom_point(aes(x=effect,
                 y=mean,
                 fill=type),
             shape=21,
             size=5, stroke=1.5, position = position_dodge(width=0.4)) +
  labs(y=paste0('Decomposed effect estimate'),x='',fill='Method') +
  guides(color=F) + 
  coord_flip() +
  scale_fill_viridis_d() + 
  theme_minimal() +
  theme(strip.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, margin = margin(r=10)),
        axis.title.x = element_text(size = 12, margin = margin(t=10)),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))
