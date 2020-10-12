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
in_dir <- 'C:/Users/ncgra/Dropbox/Penn/repos/multmed_gcomp/'
repo <- 'C:/Users/ncgra/Dropbox/Penn/repos/multmed_gcomp/'
source(paste0(repo,"/functions.R"))

#########################################################################################################
#########################################################################################################
## (1) Load data/DAG and fit all GLMs 
#########################################################################################################
#########################################################################################################

## Load data and path_cb (DAG).
all_gens <- readRDS(paste0(in_dir,'/data.RData'))
path_cb <- fread(paste0(in_dir, '/dag.csv'))

## Fit all models.
data <- copy(all_gens) ## data must contain variables 'id' and 'pweight' (weights can just all be 1 if you want un-weighted estimates)
inputs <- multmed_glms(path_cb, df_wide=data)
for(obj in names(inputs)) assign(obj, inputs[[obj]])

#########################################################################################################
#########################################################################################################
## (2) Calculate g-formulas for mediation analyses
#########################################################################################################
#########################################################################################################

## G-formula setup parameters
outcome_vars <- c('g3_home_value_percentile')
interaction_vars <- 'g1_white'
treatment_var <- 'g1_white'
duration_vars=NULL
ever_vars=NULL
cumcount_vars=NULL
dummy_vars <- NULL
ordinal_vars <- c('g2_edu_cat','g3_edu_cat')
ordinal_levels <- list(c('less_hs','hs','some_college','college','higher'),
                       c('less_hs','hs','some_college','college','higher'))
ordinal_refs <- c('less_hs','less_hs')
names(ordinal_levels) <- ordinal_vars
names(ordinal_refs) <- ordinal_vars

## Parameters for simulations.
total_sim_count <- 1 ## Number of gformula simulations (via parametric bootstraps on fitted GLMS).
mc_reps <- 50 ## Number of data replicates (to remove MC error in individual response prediction within simulations; only if you have binomial/multinomial mediators/outcomes).
sim_dir <- 'C:/Users/ncgra/Dropbox/Penn/papers/wealth_black_white/sims/' ## Directory for saving lots of tiny intermediate results (these are batched out automatically, so will not take up a lot of memory).
use_mean_betas <- T ## Set to TRUE if you only want quick point estimates of all effects (with total_sim_count==1) or set to FALSE for lots of simulations with parametric bootstraps to get standard errors on effects.
do_decomp <- TRUE ## Set to TRUE to calculate isolated expected values necessary for effect decomposition via multmed_decomp() below, or just natural/intervention courses for total effects.
decomp_type <- '4way' ## Effect decomposition type ('2way' for NDE, NDE or '4way' for CDE, PAIs, PIEs).
parallelize <- F ## Whether to parallelize gformula simulations or run serially.
windows <- TRUE ## OS if parallelizing (TRUE=Windows or FALSE=Linux).
processors <- 8 ## How many parallel threads if parallelizing?

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
    ## Counterfactual: treated as Black
    d[, g1_white := 0]
    return(d)
  },
  function(d,...) {
    ## Counterfactual: treated as white
    d[, g1_white := 1]
    return(d)
  }
)
names(intervention_rules) <- c('control_course','treatment_course')
## Specify effect pathways for decomposition (mutually exclusive and exhaustive, should include all variables in "update_vars" in path_cb).
decomp_paths <- list(
  list(paths=c('g1_house_value_percentile','g2_edu_cat','g2_house_value_percentile','g3_edu_cat'),
       outcomes=c('g3_house_value_percentile'))
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

#########################################################################################################
#########################################################################################################
## (3) Summarize results 
#########################################################################################################
#########################################################################################################

## Decomposed effect table
effect_table <- multmed_decomp(intervention_course = names(decomp_paths)[1],
                               compare_course = control_course,
                               decomp_paths = decomp_paths,
                               total_sims = total_sim_count,
                               sim_dir = sim_dir,
                               decomp_type = decomp_type)
draw_table <- effect_table[[2]]
effect_table <- effect_table[[1]]
effect_table
## Custom clean names for plotting effects
clean_names <- data.table(effect=c("CDE",
                                   "PAI_g1_house_value_percentile","PAI_g2_edu_cat","PAI_g2_house_value_percentile","PAI_g3_edu_cat",
                                   "PIE_g1_house_value_percentile","PIE_g2_edu_cat","PIE_g2_house_value_percentile","PIE_g3_edu_cat"),
                          effect_clean=c("CDE",
                                         "G1 House value percentile (PAI)","G2 Edu attainment (PAI)","G2 House value percentile (PAI)","G3 Edu attainment (PAI)",
                                         "G1 House value percentile (PIE)","G2 Edu attainment (PIE)","G2 House value percentile (PIE)","G3 Edu attainment (PIE)"))
## Effect plot
png(paste0(in_dir,'/multmed_figure.png'),height=11,width=11, units='in',res=600)
multmed_effect_plot(effect_table=effect_table,
                    clean_names=clean_names)
dev.off()
