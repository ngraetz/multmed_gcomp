multmed_glms <- function(path_cb, df, ordinal_vars = NULL) {
  df_wide <- as.data.table(df)
  max_step <- max(as.numeric(gsub('age_','',names(df_wide)[grepl('age_',names(df_wide))])))
  ## Check if imputed list. If not, make list for everything else below.
  if(!('list' %in% class(df_wide))) df_wide <- list(list(df_wide))
  ## Make formulas and set up model families.
  make_formula <- function(v) {
    f <- paste0(v,' ~ ',paste0(path_cb[update_vars==v, gsub(',','+',tc_vars)],'+',
                               path_cb[update_vars==v, gsub(',','+',tv_vars)]))
    as.formula(f)
  }
  formulas <- lapply(unique(path_cb[, update_vars]), make_formula)
  names(formulas) <- unique(path_cb[, update_vars])
  families <- lapply(path_cb[, family], get)
  names(families) <- path_cb[, family]
  ## Fit full set of GLMs.
  all_models <- list()
  for(i in 1:length(formulas)) {
    f <- names(formulas)[i]
    f_cb <- path_cb[i,]
    message(paste0('Fitting model for ',f))
    ## IF AGE-SPECIFIC MODEL
    if(!grepl('_',f_cb[, age])) {
      ## Non-censored sample for this step (but weights are longitudinal, correcting for censoring)
      grab_period <- function(d,t) {
        t_d <- d[!is.na(get(paste0('pweight_',t))),]
        return(t_d)
      }
      t <- f_cb[, age]
      t_DF_imp_list <- df_wide
      t_DF_imp_list[[1]] <- lapply(t_DF_imp_list[[1]], grab_period, t)
      ## Survey design for this step
      ## NEED TO FIX TO ALLOW MULTIPLE IMPUTATIONS 
      ff.design <- svydesign(id=~id, weights=as.formula(paste0('~pweight_',t)), data=t_DF_imp_list[[1]][[1]])
      ## Models for this step
      t_models <- gfit.init.survey.imputed(formulas[i], families[i], survey_design=ff.design, imputed=F)
      all_models <- c(all_models,t_models)
    }
    ## IF AGE-POOLED MODEL
    if(grepl('_',f_cb[, age])) {
      ## Repeat the wide dataset over age
      rep_data <- function(d,ages) {
        tv_vars_reshape <- attr(terms(as.formula(paste0('~',f_cb[,tv_vars]))), 'term.labels')
        tv_vars_reshape <- unique(unlist(strsplit(tv_vars_reshape,':')))
        tv_vars_reshape_lag <- tv_vars_reshape[grepl('l[.]',tv_vars_reshape)]
        tv_vars_reshape <- gsub('l[.]','',tv_vars_reshape_lag)
        d <- d[order(id,step)]
        d[, (paste0('l.',tv_vars_reshape)) := shift(.SD), .SDcols=tv_vars_reshape, by='id']
        for(v in names(d)[!(names(d) %in% ordinal_vars)]) d[, (v) := as.numeric(get(v))]
        return(d)
      }
      make_df_long <- function(d) {
        tc_vars <- unlist(strsplit(f_cb[,tc_vars],','))
        tc_vars <- unique(unlist(strsplit(tc_vars,'[*]')))
        df_long <- melt(df, id.vars = c('id',tc_vars))
        for(s in max_step:1) {
          df_long[grepl(paste0('_',s),variable), step := s]
          df_long[grepl(paste0('_',s),variable), variable := gsub(paste0('_',s),'',variable)]
        }
        df_long <- dcast(df_long, as.formula(paste0(paste(c('id','step',tc_vars),collapse='+'),'~ variable')), value.var = 'value')
        return(df_long)
      }
      df_long <- suppressWarnings(lapply(df_wide[[1]], make_df_long))
      df_long <- list(df_long)
      t_DF_imp_list <- df_long
      t_DF_imp_list[[1]] <- lapply(t_DF_imp_list[[1]], rep_data, ages)
      ## Fit pooled models
      #test <- imputationList(t_DF_imp_list[[1]])
      ## NEED TO FIX TO ALLOW MULTIPLE IMPUTATIONS 
      ff.design <- svydesign(id=~id, weights=~pweight, data=t_DF_imp_list[[1]][[1]])
      models <- gfit.init.survey.imputed(formulas[i], families[i], survey_design=ff.design, imputed=F)
      all_models <- c(all_models,models)
    }
  }
  for(i in 1:length(families)) all_models[[i]]$family$family <- names(families[i])
  names(all_models) <- unlist(lapply(formulas, function(x) all.vars(x)[1]))
  ## Return a pooled survey design object 
  ff.design.pooled <- svydesign(id=~id, weights=as.formula(paste0('~pweight_',t)), data=df)
  ## Return variable lists
  tc_vars <- unique(unlist(lapply(path_cb[, tc_vars], function(x) unlist(strsplit(x,',')))))
  tc_vars <- unique(unlist(strsplit(tc_vars,'[*]')))
  tv_vars <- unique(unlist(lapply(path_cb[, tv_vars], function(x) unlist(strsplit(x,',')))))
  tv_vars <- unique(unlist(strsplit(tv_vars,'[*]')))
  update_vars <- unique(unlist(lapply(path_cb[, update_vars], function(x) unlist(strsplit(x,',')))))
  update_vars <- unique(unlist(strsplit(update_vars,'[*]')))
  ## Final output list
  return(list(models=all_models, ff.design.pooled=ff.design.pooled, tc_vars=tc_vars, tv_vars=tv_vars, update_vars=update_vars))
}

gfit.init.survey.imputed <- function(formulas, families, survey_design, data=NULL, kwargs=NULL, sub=NULL, imputed=T){
  if(class(formulas) != "list"){
    stop("formulas argument must be a list of formulas.")
  }
  if(class(families) != "list"){
    stop("families argument must be a list of family objects.")
  }
  parLengths <- sapply(list(formulas, families), length)
  if(length(unique(parLengths)) != 1){
    stop("families, formulas, and functions arguments must be same length.")
  }
  lapply(1:length(formulas), function(i){
    ## Fit model function call (svyglm) over imputation list
    if(names(families[i]) %in% c('gaussian','quasibinomial','poisson')) {
      if(imputed) {
        model_list <- with(survey_design, svyglm(formula=formulas[[i]],
                                                 family=families[[i]],
                                                 design=survey_design))
      }
      if(!imputed) {
        m <- svyglm(formula=formulas[[i]],
                    family=families[[i]],
                    design=survey_design)
        model_list = list(m)
      }
    }
    if(names(families[i]) %in% 'svyolr') {
      if(imputed) {
        model_list <- with(survey_design, svyolr(formula=formulas[[i]],
                                                 design=survey_design))
      }
      if(!imputed) {
        m <- svyolr(formula=formulas[[i]], design=survey_design)
        model_list = list(m)
      }
    }
    ## Extract fitted betas/vcovs and pool
    beta_list <- lapply(model_list, FUN = function(x){coef(x)})
    vcov_list <- lapply(model_list, FUN = function(x){vcov(x, complete=T)})
    ## vcov() is more complicated for ordinal models, see MASS:::vcov.polr()
    if(names(families[i]) %in% 'svyolr') {
      get_ordinal_vcov <- function(x){
        class(x) <- 'polr'
        return(vcov(x))
      }
      vcov_list <- lapply(model_list, get_ordinal_vcov)
    }
    pooled_model <- miceadds::pool_mi(qhat = beta_list, u = vcov_list)
    return(pooled_model)
  })
}

multmed_gformula <- function(repo,
                             data,
                             tc_vars,
                             models,
                             sim_dir,
                             natural_courses,
                             intervention_rules,
                             path_cb,
                             control_course,
                             ## Options
                             parallelize=F,
                             processors=1,
                             windows=T,
                             total_sim_count = 1,
                             mc_replicates = 50,
                             decomp_paths = NULL,        
                             steps=NULL,
                             treatment_var = NULL,
                             mean_betas = TRUE,
                             decomp = TRUE,
                             decomp_type = '2way',
                             dummy_vars=NULL,
                             ordinal_levels=NULL,
                             ordinal_vars=NULL,
                             ordinal_refs=NULL,
                             duration_vars=NULL,
                             ever_vars=NULL,
                             cumcount_vars=NULL) {
  ptm <- proc.time()
  if(!parallelize) {
    ## Calculate all g-formulas serially on local machine
    message(paste0('Running ',total_sim_count,' draws locally'))
    message(paste0('\nInterventions selected: ', names(intervention_rules)))
    message(paste0('Decomp type: ', decomp_type))
    g <- lapply(1:total_sim_count, run_gformula_bw_cluster,
                data=data,
                tc_vars=tc_vars,
                models=models,
                sim_dir=sim_dir,
                natural_courses=natural_courses,
                intervention_rules=intervention_rules,
                path_cb=path_cb,
                mc_replicates = mc_reps,
                decomp = decomp,
                decomp_compare_df = control_course, ## Must be name in intervention_rules or natural_courses
                decomp_paths = decomp_paths,                          
                treatment_var = treatment_var,
                mean_betas = use_mean_betas,
                decomp_type = decomp_type,
                dummy_vars=dummy_vars, ordinal_levels=ordinal_levels, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs,
                duration_vars=duration_vars,ever_vars=ever_vars, cumcount_vars=cumcount_vars)
  }
  if(parallelize) {
    message(paste0('Running ',total_sim_count,' draws in parallel using ', processors, ' processors'))
    message(paste0('\nInterventions selected: ', names(intervention_rules)))
    if(windows) {
      ## Initialize parallel computing cluster
      cluster <- makeCluster(processors)
      clusterExport(cluster, list('repo','path_cb','models','tc_vars'))
      clusterEvalQ(cluster, library(data.table))
      clusterEvalQ(cluster, library(pillar))
      clusterEvalQ(cluster, library(survey))
      clusterEvalQ(cluster, library(mice))
      clusterEvalQ(cluster, library(boot))
      clusterEvalQ(cluster, library(MASS))
      clusterEvalQ(cluster, library(miceadds))
      clusterEvalQ(cluster, library(mitools))
      clusterEvalQ(cluster, source(paste0(repo,"/functions_Feb2020.R")))
      ## Distribute g-formula jobs
      clusterApplyLB(cluster, 1:total_sim_count, run_gformula_bw_cluster,
                     data=data,
                     tc_vars=tc_vars,
                     models=models,
                     sim_dir=sim_dir,
                     natural_courses=natural_courses,
                     intervention_rules=intervention_rules,
                     path_cb=path_cb,
                     mc_replicates = mc_reps,
                     decomp = decomp,
                     decomp_compare_df = control_course, ## Must be name in intervention_rules or natural_courses
                     decomp_paths = decomp_paths,                          
                     treatment_var = treatment_var,
                     mean_betas = use_mean_betas,
                     decomp_type = decomp_type,
                     dummy_vars=dummy_vars, ordinal_levels=ordinal_levels, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs,
                     duration_vars=duration_vars,ever_vars=ever_vars, cumcount_vars=cumcount_vars)
      ## Stop cluster
      stopCluster(cluster)
    }
    if(!windows) {
      mclapply(1:total_sim_count, run_gformula_bw_cluster,
               data=data,
               tc_vars=tc_vars,
               models=models,
               sim_dir=sim_dir,
               natural_courses=natural_courses,
               intervention_rules=intervention_rules,
               path_cb=path_cb,
               mc_replicates = mc_reps,
               decomp = decomp,
               decomp_compare_df = control_course, ## Must be name in intervention_rules or natural_courses
               decomp_paths = decomp_paths,                          
               treatment_var = treatment_var,
               mean_betas = use_mean_betas,
               decomp_type = decomp_type,
               dummy_vars=dummy_vars, ordinal_levels=ordinal_levels, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs,
               duration_vars=duration_vars,ever_vars=ever_vars, cumcount_vars=cumcount_vars,
               mc.cores = processors)
    }
  }
  return(proc.time() - ptm)
}

run_gformula_bw_cluster <- function(b, data, outcome_vars=NULL, tc_vars, models, binary_cols=NULL, sim_dir, path_cb, mc_replicates,
                                    decomp, decomp_compare_df, decomp_paths, treatment_var, natural_courses, intervention_rules, intervention_subset = NULL,
                                    age_interaction_vars = NULL, duration_vars = NULL, cumcount_vars = NULL,
                                    ever_vars = NULL, dummy_vars = NULL, ordinal_vars = NULL, ordinal_levels = NULL,ordinal_refs=NULL,
                                    mean_betas = FALSE, special_rules = NULL, special_vars = NULL, extra_save_vars = NULL,
                                    course_draw_micro = NULL, interactive_decomp = NULL, decomp_type='4way') {
  ## Update output dir so that sims are split into groups of 250 for saving in any particular folder.
  ## Having thousands of tiny files in a single directory makes things uber slow.
  message(paste0('\nSample ', b))
  batch <- round(b/250)
  dir.create(sim_dir)
  ## Do not double-count special rules variables in the "effect_paths", even though we use this to decompose total effects.
  update_vars <- path_cb[, update_vars]
  binary_cols <- paste0(binary_cols,'1')
  ## If using non-parametric bootstrap, do that now and refit models.
  # if(use_np_boot) {
  #   
  # }
  # Multiply dataset X times to reduce MC error from individual-level prediction.
  # Update "id" variable to be "id+mc" to avoid lagging incorrectly within "id" in functions below.
  # Order by id, age.
  mc_data <- function(mc, d) {
    d <- copy(data)
    ## Add MC id
    d[, id := paste0(id,'_',mc)]
    d[, mc := mc]
    return(d)
  }
  data <- as.data.table(data)
  if(sum(c('id','pweight_1') %in% names(data))!=2) stop('Data must contain variable id, pweight_1')
  data <- rbindlist(lapply(1:mc_replicates, mc_data, d=data))
  ## Make batches in case we need to parallelize really long datasets
  mc_list <- lapply(seq(0,2000,50), function(i) c((i+1):(i+50)))
  mc_list <- mc_list[unlist(lapply(mc_list,function(i) min(i) <= mc_replicates))]
  # Create an unweighted, nationally representative black and white dataset using the weighted data.
  # We already accounted for the weighting in the models, so we now have nationally representative models anyway.
  # This makes sampling in interventions (e.g. drawing from the white natural course) more straightforward.
  data <- data[sample(.N,replace=T,prob=pweight_1)]
  ## Need to update ids to stay unique after sampling with replacement
  data[, new_id := 1:.N, by=c('id')]
  data[, id := paste0(id,'_',new_id)]
  data[, new_id := NULL]
  # Need to stochasitcally draw from predicted probabilities of binary variables pooled across imputed datasets.
  # draw_vars <- binary_cols
  # for(v in ordinal_vars) draw_vars <- draw_vars[!(grepl(v,draw_vars))]
  # data[, (draw_vars) := lapply(.SD, function(x) {rbinom(nrow(data), 1, x)}), .SDcols = draw_vars]
  # Sample betas
  if(mean_betas) betas_draw <- lapply(models, function(m) coef(m)) ## For testing points estimates only
  if(!mean_betas) betas_draw <- lapply(models, function(m) mvrnorm(1, mu = coef(m), Sigma = vcov(m)))
  names(betas_draw) <- names(models)
  # Estimate natural course(s).
  for(n in names(natural_courses)) {
    message(paste0('   Simulating ', n, ' course...'))
    natural_function <- natural_courses[[n]]
    int_data <- copy(data)
    natural_course <- rbindlist(lapply(mc_list, gformula_sim, data=int_data, tc_vars=tc_vars, models=models, 
                                       betas_draw=betas_draw, effect_cb=path_cb,
                                       dummy_vars=dummy_vars, ordinal_levels=ordinal_levels, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs,
                                       duration_vars=duration_vars,
                                       ever_vars=ever_vars, cumcount_vars=cumcount_vars))
    natural_course[, sim := b]
    factor_vars <- names(natural_course)[grep('factor|character', sapply(natural_course, class))]
    collapse_vars <- names(natural_course)[!(names(natural_course) %in% c('id','age','actual_age','pweight',factor_vars))]
    natural_agg <- natural_course[, lapply(.SD, mean), .SDcols=collapse_vars, by='sim']
    natural_agg[, name := n]
    dir.create(paste0(sim_dir, '/',n,'/'), showWarnings=F)
    dir.create(paste0(sim_dir, '/',n,'/', batch),showWarnings=F)
    saveRDS(natural_agg, paste0(sim_dir,'/',n,'/',batch,'/sim',b,'.RDS'))
    assign(paste0(n,'_micro'), natural_course)
    rm(natural_course)
  }
  # Estimate intervention courses
  for(n in names(intervention_rules)) {
    message(paste0('   Simulating ', n, ' course...'))
    intervention_function <- intervention_rules[[n]]
    int_data <- copy(data)
    micro <- NULL
    if(!is.null(course_draw_micro)) micro <- get(course_draw_micro)
    treatment_course_micro <- rbindlist(lapply(mc_list, gformula_sim, data=int_data, tc_vars=tc_vars, models=models, 
                                               betas_draw=betas_draw, effect_cb=path_cb,
                                               dummy_vars=dummy_vars, ordinal_levels=ordinal_levels, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs,
                                               duration_vars=duration_vars,
                                               ever_vars=ever_vars, cumcount_vars=cumcount_vars,
                                               intervention_rules=intervention_function, micro=micro))
    treatment_course_micro[, sim := b]
    factor_vars <- names(treatment_course_micro)[grep('factor|character', sapply(treatment_course_micro, class))]
    collapse_vars <- names(treatment_course_micro)[!(names(treatment_course_micro) %in% c('id','age','actual_age','pweight',factor_vars))]
    collapse_vars <- collapse_vars[!grepl(':',collapse_vars)]
    intervention_agg <- treatment_course_micro[, lapply(.SD, mean), .SDcols=collapse_vars, by='sim']
    intervention_agg[, name := n]
    dir.create(paste0(sim_dir, '/',n,'/'), showWarnings=F)
    dir.create(paste0(sim_dir, '/',n,'/', batch),showWarnings=F)
    saveRDS(intervention_agg, paste0(sim_dir,'/',n,'/',batch,'/sim',b,'.RDS'))
    if(n==decomp_compare_df) assign('control_course_micro', treatment_course_micro)
    ## Decomposition of RATE for this intervention into CDE and one set of [INT_REF,INT_MED,PIE] for each mediator.
    if(decomp & n!=decomp_compare_df) {
      ## Read microdata simulants for the two courses we are drawing from to decompose effects.
      o <- decomp_paths[[n]][['outcomes']]
      decomp_cb <- path_cb[update_vars %in% o,]
      #message(paste0('  ATE: ',round(mean(treatment_course_micro[,get(o)])-mean(control_course_micro[,get(o)]),2)))
      message(decomp_type)
      if(decomp_type=='4way') {
        ## Using predicted values from the two treatment regimes, re-estimate E[Y] under variants of the g-formula necessary
        ## for calculating the CDE and [INT_REF,INT_MED,PIE] for each mediator.
        ## This requires 4 + 2 * (# Mediators) simulations. 
        ## 1) Calculate all four expected values of the outcome that do not rely on a specific mediator.
        ##  Y1 and Y2 are the treated/control courses where we fix all mediators at their reference values (0 for now, but could be provided).
        for(i in 1:2) {
          if(i==1) {
            med_eff_vars <- treatment_var
            med_other_vars <- decomp_paths[[n]][['paths']]
            ref_vars <- decomp_paths[[n]][['paths']]
          }
          if(i==2) {
            med_eff_vars <- NULL
            med_other_vars <- c(treatment_var,decomp_paths[[n]][['paths']])   
            ref_vars <- decomp_paths[[n]][['paths']]
          }
          message(paste0('   Effect decomposition: Y',i,'...'))
          course <- rbindlist(lapply(mc_list, gformula_sim, data=int_data, tc_vars=tc_vars, models=models, 
                                     dummy_vars=dummy_vars, ordinal_levels=ordinal_levels, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs,
                                     duration_vars=duration_vars,
                                     ever_vars=ever_vars, cumcount_vars=cumcount_vars,
                                     betas_draw=betas_draw, effect_cb=decomp_cb,
                                     natural_DF=control_course_micro, intervention_DF=treatment_course_micro,
                                     natural_vars=med_other_vars, intervention_var=med_eff_vars,
                                     reference_vars=ref_vars))
          course[, sim := b]
          factor_vars <- names(course)[grep('factor|character', sapply(course, class))]
          collapse_vars <- names(course)[!(names(course) %in% c('id','age','actual_age','pweight',factor_vars))]
          course_agg <- course[, lapply(.SD, mean), .SDcols=collapse_vars, by='sim']
          course_agg[, name := paste0('y',i)]
          dir.create(paste0(sim_dir, '/',n,'/effects'),showWarnings=F)
          dir.create(paste0(sim_dir, '/',n,'/effects/y',i),showWarnings=F)
          dir.create(paste0(sim_dir, '/',n,'/effects/y',i, '/', batch),showWarnings=F)
          saveRDS(course_agg, paste0(sim_dir, '/',n,'/effects/y',i, '/', batch,'/sim',b,'.RDS'))
        }
        ## 2) Calculate the two expected values (Y5, Y6) that rely on each specific mediator.
        ##  Y3 is treated, target mediator at control value, all other mediators at reference values (0).
        ##  Y4 is control, target mediator at control value, all other mediators at reference values (0).
        ##  Y5 is treated, target mediator at treated value, all other mediators at reference values (0).
        ##  Y6 is control, target mediator at treated value, all other mediators at reference values (0).
        for(path in names(decomp_paths[[n]][['paths']])) {
          all_paths <- unlist(decomp_paths[[n]][['paths']])
          these_intervention_vars <- unlist(decomp_paths[[n]][['paths']][path])
          #these_natural_vars <- all_paths[!grepl(path,all_paths)]
          these_natural_vars <- all_paths
          for(v in these_intervention_vars) these_natural_vars <- these_natural_vars[!grepl(v,these_natural_vars)]
          for(i in 3:6) {
            if(i==3) {
              med_eff_vars <- treatment_var
              med_other_vars <- c(these_intervention_vars,these_natural_vars)
              ref_vars <- these_natural_vars
            }
            if(i==4) {
              med_eff_vars <- NULL
              med_other_vars <- c(these_intervention_vars,these_natural_vars,treatment_var)  
              ref_vars <- these_natural_vars
            }
            if(i==5) {
              med_eff_vars <- c(these_intervention_vars,treatment_var)
              med_other_vars <- these_natural_vars
              ref_vars <- these_natural_vars
            }
            if(i==6) {
              med_eff_vars <- these_intervention_vars
              med_other_vars <- c(these_natural_vars,treatment_var)  
              ref_vars <- these_natural_vars
            }
            message(paste0('   Effect decomposition: Y',i,' via ', path,'...'))
            course <- rbindlist(lapply(mc_list, gformula_sim, data=int_data, tc_vars=tc_vars, models=models, 
                                       betas_draw=betas_draw, effect_cb=decomp_cb,
                                       dummy_vars=dummy_vars, ordinal_levels=ordinal_levels, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs,
                                       duration_vars=duration_vars,
                                       ever_vars=ever_vars, cumcount_vars=cumcount_vars,
                                       natural_DF=control_course_micro, intervention_DF=treatment_course_micro,
                                       natural_vars=med_other_vars, intervention_var=med_eff_vars,
                                       reference_vars=ref_vars))
            course[, sim := b]
            factor_vars <- names(course)[grep('factor|character', sapply(course, class))]
            collapse_vars <- names(course)[!(names(course) %in% c('id','age','actual_age','pweight',factor_vars))]
            course_agg <- course[, lapply(.SD, mean), .SDcols=collapse_vars, by='sim']
            course_agg[, name := path]
            dir.create(paste0(sim_dir, '/',n,'/effects'),showWarnings=F)
            dir.create(paste0(sim_dir, '/',n,'/effects/y',i,'_', path),showWarnings=F)
            dir.create(paste0(sim_dir, '/',n,'/effects/y',i,'_', path, '/', batch),showWarnings=F)
            saveRDS(course_agg, paste0(sim_dir, '/',n,'/effects/y',i,'_', path, '/', batch,'/sim',b,'.RDS'))
          }
        }
      }
      if(decomp_type=='2way') {
        for(path in c(treatment_var,names(decomp_paths[[n]][['paths']]))) {
          all_paths <- c(treatment_var,unlist(decomp_paths[[n]][['paths']]))
          these_intervention_vars <- all_paths[grepl(path,all_paths)]
          these_natural_vars <- all_paths[!grepl(path,all_paths)]
          message(paste0('   Effect decomposition (2-way): ',path,'...'))
          course <- rbindlist(lapply(mc_list, gformula_sim, data=int_data, tc_vars=tc_vars, models=models, 
                                     betas_draw=betas_draw, effect_cb=decomp_cb,
                                     dummy_vars=dummy_vars, ordinal_levels=ordinal_levels, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs,
                                     duration_vars=duration_vars,
                                     ever_vars=ever_vars, cumcount_vars=cumcount_vars,
                                     natural_DF=control_course_micro, intervention_DF=treatment_course_micro,
                                     natural_vars=these_natural_vars, intervention_var=these_intervention_vars,
                                     reference_vars=NULL))
          course[, sim := b]
          factor_vars <- names(course)[grep('factor|character', sapply(course, class))]
          collapse_vars <- names(course)[!(names(course) %in% c('id','age','actual_age','pweight',factor_vars))]
          course_agg <- course[, lapply(.SD, mean), .SDcols=collapse_vars, by='sim']
          course_agg[, name := path]
          dir.create(paste0(sim_dir, '/',n,'/effects'),showWarnings=F)
          dir.create(paste0(sim_dir, '/',n,'/effects/', path),showWarnings=F)
          dir.create(paste0(sim_dir, '/',n,'/effects/', path, '/', batch),showWarnings=F)
          saveRDS(course_agg, paste0(sim_dir, '/',n,'/effects/',path, '/', batch,'/sim',b,'.RDS'))
          message(paste0('      Effect path estimate: ',round(course_agg[,get(o)]-mean(control_course_micro[,get(o)]),2)))
        }
      }
    }
  }
  message('Sample finished!')
  return(NULL)
}

gformula_sim <- function(this_mc,
                         data,
                         tc_vars,
                         models,
                         betas_draw,
                         effect_cb,
                         ## Arguments for effect decomposition
                         natural_vars=NULL, 
                         intervention_var=NULL,
                         natural_DF=NULL,
                         intervention_DF=NULL,
                         reference_vars=NULL,
                         ## Extra arguments for special cases/exceptions.
                         duration_vars=NULL,
                         cumcount_vars=NULL,
                         ever_vars=NULL,
                         dummy_vars=NULL,
                         ordinal_vars=NULL,
                         ordinal_levels=NULL,
                         ordinal_refs=NULL,
                         interaction_vars=NULL,
                         intervention_rules=NULL,
                         special_rules=NULL, 
                         special_vars=NULL,
                         special_var_names=NULL,
                         exp_vars=NULL,
                         paths=NULL,
                         micro=NULL,
                         decomp_main=NULL) {
  ## Setup
  ptm <- proc.time()
  simDF <- copy(data[mc %in% this_mc,])
  ## We have to assert the intervention rules at t=0 as well as in updating below.
  if(!is.null(intervention_rules)) simDF <- intervention_rules(simDF,micro=micro)
  ## We have to assert decomp rules for effect calculations at t=0 as well as in updating below.
  if(!is.null(natural_DF)) simDF <- effect_rules_fast(simDF, natural_vars, intervention_var, natural_DF, intervention_DF, reference_vars, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs, ordinal_levels=ordinal_levels)
  ## Reset duration-weighted variables to their t=1 values just in case the intervention affects these values (all because of pooled imputation).
  simDF <- update_dummy_interaction(simDF, dummy_vars, interaction_vars, effect_cb, ordinal_levels, tc_vars)
  ## Add intercept variable for manual predicting with draws from our models.
  simDF[, ('(Intercept)') := 1]
  ## Simulate forward (wide) predicting all post-treatment variables (mediators or treatment-induced mediator-outcome confounders).
  message('      Beginning gformula simulation')
  ## Update one variable at a time in the sequential order they are listed in the effects codebook.
  all_ages <- as.numeric(unlist(effect_cb[,tstrsplit(age,'_')]))
  all_ages <- all_ages[!is.na(all_ages)]
  steps <- min(all_ages):max(all_ages)
  for(s in steps) {
    message(paste0('         ',s))
    this_step_vars <- c()
    for(v in effect_cb[, update_vars]) {
      var_s <- effect_cb[update_vars==v, age]
      if(grepl('_',var_s)) var_s <- unlist(strsplit(var_s,'_'))[1]:unlist(strsplit(var_s,'_'))[2]
      if(s %in% var_s) this_step_vars <- c(this_step_vars,v)
    }
    for(v in this_step_vars) {
      #message(paste0('         ',v))
      age_cb <- copy(effect_cb[update_vars==v,])
      y <- age_cb[, age]
      pool <- 0
      upDF <- copy(simDF)
      all_new_vars <- c()
      all_grab_vars <- c()
      ## If pooled model (one model fit across repeated observations), grab lagged variables.
      if(grepl('_',age_cb[update_vars==v, age])) {
        pool <- 1
        #lag_vars <- unique(unlist(lapply(age_cb[update_vars==v, tv_vars], function(x) unlist(strsplit(x,',')))))
        #grab_vars_current <- lag_vars[!grepl('l[.]',lag_vars)]
        #lag_vars <- lag_vars[grepl('l[.]',lag_vars)]
        lag_vars <- attr(terms(as.formula(paste0('~',age_cb[update_vars==v,tv_vars]))), 'term.labels')
        lag_vars <- unique(unlist(strsplit(lag_vars,':')))
        grab_vars_current <- lag_vars[!grepl('l[.]',lag_vars)]
        grab_vars_current <- grab_vars_current[!(grab_vars_current %in% interaction_vars)]
        lag_vars <- lag_vars[grepl('l[.]',lag_vars)]
        for(l in lag_vars) {
          clean_name <- gsub('as[.]factor[(]','',l)
          clean_name <- gsub('[)]','',clean_name)
          clean_name <- gsub('l[.]','',clean_name)
          grab_vars <- names(upDF)[grepl(paste0(clean_name,'_',s-1),names(upDF))]
          grab_vars <- grab_vars[!grepl('c[.]|d[.]|e[.]',grab_vars)]
          for(g in grab_vars) {
            g_new <- gsub(paste0('_',s-1),'',g)
            if(grepl('factor',l)) g_new <- gsub('as.factor[(]','as.factor(l.',g_new)
            if(!grepl('factor',l)) g_new <- paste0('l.',g_new)
            upDF[, (g_new) := get(g)]
            all_new_vars <- c(all_new_vars,g_new)
            all_grab_vars <- c(all_grab_vars,g)
          }
        }
        for(g in grab_vars_current) {
          current <- paste0(g,'_',s)
          if(grepl('factor',current)) upDF[, (gsub(paste0('_',s),'',current)) := get(current)]
          if(!grepl('factor',current)) upDF[, (g) := get(current)]
          all_new_vars <- c(all_new_vars,g)
          all_grab_vars <- c(all_grab_vars,g)
        }
        #upDF <- update_dummy_interaction(upDF, 1, dummy_vars, interaction_vars, effect_cb, ordinal_levels, tc_vars)
        upDF <- update_dummy_interaction(upDF, dummy_vars, interaction_vars, effect_cb, ordinal_levels, tc_vars)
      }
      ## Subset upDF to ONLY necessary inputs for this update variable from simDF (according to variable names from parameter sample).
      input_vars <- c('id',interaction_vars,names(betas_draw[[v]])[!grepl('[|]',names(betas_draw[[v]]))])
      upDF <- upDF[, unique(input_vars), with=F]
      ## 1) Predict new variable from models. If pooled model, rename prediction to "_[age]" at the end.
      upDF <- natural_rules_fast(upDF, s, v, models, betas_draw)
      for(test_v in input_vars) if(dim(upDF[is.na(get(test_v)),])[1]>0) stop(paste0('Missing values in ',test_v))
      if(pool==1) {
        if(!(paste0(v,'_',s) %in% ordinal_vars)) setnames(upDF, v, paste0(v,'_',s))
        v <- paste0(v,'_',s)
        for(delete in all_new_vars) suppressWarnings(upDF[, (delete) := NULL])
      }
      ## 2) Replace predictions with either values for effect decomp or intervention rules (for time-varying treatment regimes), if necessary.
      if(!is.null(natural_DF)) {
        if(v %in% natural_vars) upDF <- effect_rules_fast(upDF, step_natural_vars=v, step_intervention_var=NULL, natural_DF, intervention_DF, reference_vars, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs, ordinal_levels=ordinal_levels) 
        if(v %in% intervention_var) upDF <- effect_rules_fast(upDF, step_natural_vars=NULL, step_intervention_var=v, natural_DF, intervention_DF, reference_vars, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs, ordinal_levels=ordinal_levels) 
      }
      if(!is.null(intervention_rules)) upDF <- intervention_rules(upDF,micro)
      ## 3) Update manual interaction terms for future predictions.
      #if(v %in% names(upDF)) upDF <- update_dummy_interaction(upDF, s, ifelse(v %in% dummy_vars,v,NA), interaction_vars, effect_cb, ordinal_levels, tc_vars)
      upDF <- update_dummy_interaction(upDF, dummy_vars, interaction_vars, effect_cb, ordinal_levels, tc_vars)
      ## 4) Update for next prediction/draw steps if it is an ordinal variable.
      # if(sum(grepl(v,ordinal_vars))!=0) {
      #   for(val in unique(upDF[, get(v)])) {
      #     upDF[, (paste0(v,val)) := NULL]
      #     upDF[, (paste0(v,val)) := ifelse(get(v)==val,1,0)]
      #   }
      # }
      ## 5) Tack on only new updated variable to the full wide simDF.
      new_vars <- names(upDF)[grepl(gsub(paste0('_',s),'',v),names(upDF))]
      for(x in new_vars) if(x %in% names(simDF)) simDF <- simDF[, (x) := NULL]
      simDF <- cbind(simDF, upDF[, new_vars, with=F])
      ## 6) Add updated duration-weighted variable to simDF if necessary
      if(gsub(paste0('_',s),'',v) %in% c(duration_vars,ever_vars,cumcount_vars)) {
        simDF <- update_duration_ever(simDF, s, ordinal_vars, gsub(paste0('_',s),'',v), ever_vars, cumcount_vars)
        if(!is.null(natural_DF)) {
          for(time in c('d.','c.','e.')) {
            if(paste0(time,v) %in% natural_vars) simDF <- effect_rules_fast(simDF, step_natural_vars=paste0(time,v), step_intervention_var=NULL, natural_DF, intervention_DF, reference_vars, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs, ordinal_levels=ordinal_levels)
            if(paste0(time,v) %in% intervention_var) simDF <- effect_rules_fast(simDF, step_natural_vars=NULL, step_intervention_var=paste0(time,v), natural_DF, intervention_DF, reference_vars, ordinal_vars=ordinal_vars, ordinal_refs=ordinal_refs, ordinal_levels=ordinal_levels)
          }
        }
      }
      for(test_v in new_vars) if(dim(simDF[is.na(get(test_v)),])[1]>0) stop('Missing values in new prediction; data should be complete at t=1.')
      if(!is.null(intervention_rules)) simDF <- intervention_rules(simDF,micro)
    }
  }
  print(proc.time() - ptm)
  return(simDF)
}

effect_rules_fast <- function(d, step_natural_vars, step_intervention_var, natural_DF, intervention_DF, reference_vars=NULL, ordinal_vars=NULL, ordinal_refs=NULL, ordinal_levels=NULL) {
  #message(paste(ordinal_vars,collapse = ' '))
  ## Replace all variable values with those from the natural course. 
  if(!is.null(step_natural_vars)) {
    for(v in step_natural_vars) {
      ## If normal/binomial.
      if(!(v %in% ordinal_vars)) {
        d[, (v) := NULL]
        d <- cbind(d, natural_DF[, (v), with=F])
      }
      ## If ordinal
      if(v %in% ordinal_vars) {
        d[, (paste0(v,ordinal_levels[[v]])) := NULL]
        d <- cbind(d, natural_DF[, (paste0(v,ordinal_levels[[v]])), with=F])
      }
    }
  }  
  ## Replace effect path variable with values from the intervention course.
  if(!is.null(step_intervention_var)) {
    for(v in step_intervention_var) {
      ## If normal/binomial.
      if(!(v %in% ordinal_vars)) {
        d[, (v) := NULL]
        d <- cbind(d, intervention_DF[, (v), with=F])
      }
      ## If ordinal
      if(v %in% ordinal_vars) {
        d[, (paste0(v,ordinal_levels[[v]])) := NULL]
        d <- cbind(d, intervention_DF[, (paste0(v,ordinal_levels[[v]])), with=F])
      }
    }
  }  
  ## Replace any reference variables with 0s (though could be any fixed value at which the CDE is evaluated).
  for(v in reference_vars) {
    ## If normal/binomial.
    if(!(v %in% ordinal_vars)) {
      d[, (v) := 0]
    }
    ## If ordinal.
    if(v %in% ordinal_vars) {
      ## Grab reference for this variable, create dummy=1, and set all other dummies=0.
      ordinal_ref <- ordinal_refs[[v]]
      for(level in ordinal_levels[[v]][ordinal_levels[[v]]!=ordinal_ref]) d[, (paste0(v,level)) := 0]
      d[, (paste0(v,ordinal_ref)) := 1]
    }
  }
  return(d)
}

natural_rules_fast <- function(d, s, update_vars, models, betas) {
  for(v in update_vars) d <- simPredict_fast(s, v, d, models, betas)
  return(d)
}

simPredict_fast <- function(s, v, DF, models, betas){
  model_index <- match(v, names(models))
  model_ <- models[[model_index]]
  betas_ <- betas[[model_index]]
  if(model_$family$family %in% c('binomial','quasibinomial')) {
    # newDF <- as.data.table(DF[age==max(age), ])
    newDF <- copy(DF)
    beta_names <- names(betas_)
    newDF <- newDF[, beta_names, with=F]
    setcolorder(newDF, beta_names)
    predicted_probs <- inv.logit(as.numeric(betas_ %*% t(as.matrix(newDF))))
    # sim <- as.integer(rbinom(nrow(newDF), 1, predicted_probs))
    # DF[, (v) := predicted_probs]
    DF[, (v) := as.integer(rbinom(nrow(newDF), 1, predicted_probs))]
  }
  if(model_$family$family == 'poisson') {
    newDF <- copy(DF)
    beta_names <- names(betas_)
    newDF <- newDF[, beta_names, with=F]
    setcolorder(newDF, beta_names)
    sim <- exp(as.numeric(betas_ %*% t(as.matrix(newDF))))
    DF[, (v) := sim]
  }
  if(model_$family$family == 'gaussian') {
    # newDF <- as.data.table(DF[age==max(age), ])
    newDF <- copy(DF)
    beta_names <- names(betas_)
    newDF <- newDF[, beta_names, with=F]
    setcolorder(newDF, beta_names)
    sim <- as.numeric(betas_ %*% t(as.matrix(newDF)))
    ## MANUALLY ADD LINEAR PROBABILITY MODEL EXCEPTIONS HERE (CODEBOOK LATER) - "sim" is the linear model predicted probabilities, so need to draw 1/0s.
    if(v %in% c('college_3','in_school_2','employed_2','env')) {
      sim[sim>1] <- 1
      sim[sim<0] <- 0
      sim <- as.integer(rbinom(nrow(newDF), 1, sim))
    }
    DF[, (v) := sim]
  }
  if(model_$family$family == 'svyolr') {
    ints <- betas_[grep('[|]', names(betas_))]
    beta_vals <- betas_[!(names(betas_) %in% names(ints))]
    beta_names <- names(beta_vals)
    newDF <- copy(DF)
    newDF <- newDF[, beta_names, with=F]
    setcolorder(newDF, beta_names)
    linear_part <- data.table(linear=as.numeric(beta_vals %*% t(as.matrix(newDF))))
    all_cat_names <- unlist(lapply(strsplit(names(ints),'[|]'), function(x) x[1]))
    linear_part[, (all_cat_names[1]) := (1 / (1 + exp(-(ints[1]-linear))))]
    for(int in 2:length(ints)) {
      cat <- all_cat_names[int]
      linear_part[, (cat) := (1 / (1 + exp(-(ints[int]-linear)))) - (1 / (1 + exp(-(ints[int-1]-linear))))]
    }
    last_cat <- strsplit(names(ints)[length(ints)],'[|]')[[1]][2]
    linear_part[, (last_cat) := 1 - (1 / (1 + exp(-(ints[length(ints)]-linear))))]
    linear_part[, linear := NULL]
    # sim <- as.character(rMultinom(as.matrix(linear_part),1))
    for(val in names(linear_part)) { 
      # DF[, (paste0('as.factor(',v,'_',unique(DF[,age]),')',val)) := linear_part[, get(val)]]
      # DF[, (paste0(v,'_',s,val)) := linear_part[, get(val)]]
      DF[, (paste0(v,val)) := linear_part[, get(val)]]
    }
  }
  return(DF)
}

rMultinom <- function(probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  ran <- matrix(lev[1], ncol=m, nrow=n)
  z <- apply(probs, 1, sum)
  if(any(abs(z-1) > 1e-05)) stop('do not sum to 1')
  U <- apply(probs, 1, cumsum)
  for(i in 1:m) {
    un <- rep(runif(n), rep(k,n))
    ran[, i] <- lev[1 + apply(un > U, 2, sum)]
  }
  ran
}

update_dummy_interaction <- function(d, dummy_vars, interaction_vars, path_cb, ordinal_levels, tc_vars) {
  # for(v in dummy_vars) {
  #  for(val in ordinal_levels[[v]]) d[, (paste0(v,val)) := ifelse(get(v)==val, 1, 0)]
  # }
  ## Infer all necessary interactions based on model coefficient names.
  get_coef_names <- function(m) {
    c <- names(coef(m))
    c <- c[grepl('[:]',c)]
    return(c)
  }
  all_ints <- unique(unlist(lapply(models,get_coef_names)))
  ## Update all interactions relevant for this time step.
  for(i in all_ints) {
    i2 <- unlist(strsplit(i,'[:]'))
    if(i2[1] %in% names(d) & i2[2] %in% names(d)) d[, (i) := get(i2[1]) * get(i2[2])]
  }
  return(d)
}

update_duration_ever <- function(d, a, ordinal_vars, duration_vars, ever_vars, cumcount_vars) {
  for(v in duration_vars) {
    if(paste0(v,'_',a) %in% ordinal_vars) {
      for(this_age in 1:a) {
        for(w in 1:5) d[, (paste0(v,'_num_',w,'_',this_age)) := get(paste0(v,'_',this_age,v,w)) * w]
      }
      d[, (paste0('d.',v,'_',a)) := rowSums(.SD, na.rm=T), .SDcols = grep(paste0(v,'_num_'),names(d))]
      d[, (paste0('d.',v,'_',a)) := get(paste0('d.',v,'_',a)) / a]
      for(extra in names(d)[grep(paste0(v,'_num_'),names(d))]) d[, (extra) := NULL]
    }
    if(!(paste0(v,'_',a) %in% ordinal_vars)) {
      previous_vars <- names(d)[grepl(paste(paste0(v,'_',1:a),collapse='|'),names(d))]
      previous_vars <- previous_vars[!grepl('d[.]',previous_vars)]
      previous_vars <- previous_vars[!(grepl('as[.]factor',previous_vars))]
      d[, (paste0('d.',v,'_',a)) := rowSums(.SD, na.rm=T), .SDcols = previous_vars]
      d[, (paste0('d.',v,'_',a)) := get(paste0('d.',v,'_',a)) / a]
    }
  }
  for(v in cumcount_vars) {
    previous_vars <- names(d)[grepl(paste(paste0(v,'_',1:a),collapse='|'),names(d))]
    previous_vars <- previous_vars[!grepl('c[.]',previous_vars)]
    previous_vars <- previous_vars[!(grepl('as[.]factor',previous_vars))]
    d[, (paste0('c.',v,'_',a)) := rowSums(.SD, na.rm=T), .SDcols = previous_vars]
  }
  for(v in ever_vars) {
    previous_vars <- names(d)[grepl(paste(paste0(v,'_',1:a),collapse='|'),names(d))]
    previous_vars <- previous_vars[!grepl('e[.]',previous_vars)]
    d[, (paste0('e.',v,'_',a)) := rowSums(.SD, na.rm=T), .SDcols = previous_vars]
    d[, (paste0('e.',v,'_',a)) := ifelse(get(paste0('e.',v,'_',a))==0,0,1)]
  }
  return(d)
}

post_process_course <- function(course, total_sims, sim_dir, decomp_paths, decomp_type='4way', treatment_var=NULL) {
  message(paste0('Processing ', course, ' course across ', total_sims, ' simulations...'))
  read_course_sim <- function(x, course) {
    batch <- round(x/250)
    readRDS(paste0(sim_dir,'/',course,'/',batch,'/sim',x,'.RDS'))
  }
  all_course_sims <- rbindlist(lapply(1:total_sims,read_course_sim,course))
  saveRDS(all_course_sims, paste0(sim_dir,'/', course, '_', total_sims, '.RDS'))
  ## Also process any direct/indirect courses for this ATE (just check any nested folders).
  if(dir.exists(paste0(sim_dir,'/',course,'/effects/'))) {
    if(decomp_type=='4way') {
      paths <- c('y1','y2',
                 paste0('y3_',names(decomp_paths[[course]]$paths)),
                 paste0('y4_',names(decomp_paths[[course]]$paths)),
                 paste0('y5_',names(decomp_paths[[course]]$paths)),
                 paste0('y6_',names(decomp_paths[[course]]$paths)))
      for(path in paths) {
        read_path_sim <- function(x, course, path) {
          batch <- round(x/250)
          readRDS(paste0(sim_dir,'/',course,'/effects/',path,'/',batch,'/sim',x,'.RDS'))
        }
        all_path_sims <- rbindlist(lapply(1:total_sims,read_path_sim,course,path))
        saveRDS(all_path_sims, paste0(sim_dir,'/', course, '_', path, '_', total_sims, '.RDS'))
      }
    }
    if(decomp_type=='2way') {
      for(path in c(treatment_var,names(decomp_paths[[course]][['paths']]))) {
        read_path_sim <- function(x, course, path) {
          batch <- round(x/250)
          readRDS(paste0(sim_dir,'/',course,'/effects/',path,'/',batch,'/sim',x,'.RDS'))
        }
        all_path_sims <- rbindlist(lapply(1:total_sims,read_path_sim,course,path))
        saveRDS(all_path_sims, paste0(sim_dir,'/', course, '_', path, '_', total_sims, '.RDS'))
      }
    }
  }
}

multmed_decomp <- function(intervention_course, compare_course, total_sims, sim_dir, decomp_paths,
                           decomp_type='4way', outcome_var=NULL, treatment_var=NULL, scale='add') {
  if(is.null(outcome_var)) outcome_var <- decomp_paths[[intervention_course]]$outcomes
  paths <- names(decomp_paths$treatment_course$paths)
  compare_course <- readRDS(paste0(sim_dir, '/', compare_course, '_', total_sims, '.RDS'))
  setnames(compare_course, outcome_var, 'compare_outcome')
  compare_course <- compare_course[, c('sim','compare_outcome')]
  ## Merge on intervention course and calculate ATE
  int_course <- readRDS(paste0(sim_dir, '/', intervention_course, '_', total_sims, '.RDS'))
  setnames(int_course, outcome_var, 'int_outcome')
  int_course <- int_course[, c('sim','int_outcome')]
  effect_table <- merge(compare_course, int_course, by=c('sim'))
  ##  Y1 and Y2 are the treated/control courses where we fix all mediators at their reference values (0 for now, but could be provided).
  ##  Y3 is treated, target mediator at control value, all other mediators at reference values (0).
  ##  Y4 is control, target mediator at control value, all other mediators at reference values (0).
  ##  Y5 is treated, target mediator at treated value, all other mediators at reference values (0).
  ##  Y6 is control, target mediator at treated value, all other mediators at reference values (0).
  if(scale=='add') {
    ## ATE p-values
    effect_table[, ATE := int_outcome - compare_outcome]
    if(effect_table[, mean(ATE)]>=0) effect_table[, ATE_pval := ifelse(ATE<0, 1, 0)]
    if(effect_table[, mean(ATE)]<0) effect_table[, ATE_pval := ifelse(ATE>0, 1, 0)]
    ## Merge on all expected values of Y to calculate four-way decomposition of ATE.
    if(decomp_type=='4way') {
      for(path in c('y1','y2',paste0('y3_',paths),paste0('y4_',paths),paste0('y5_',paths),paste0('y6_',paths))) {
        path_course <- readRDS(paste0(sim_dir,'/', intervention_course, '_', path, '_', total_sims, '.RDS'))
        setnames(path_course, outcome_var, path)
        path_course <- path_course[, c('sim',path), with=F]
        effect_table <- merge(effect_table, path_course, by=c('sim'))
      }
      ## CDE
      effect_table[, CDE := y1-y2]
      if(effect_table[, mean(CDE)]>=0) effect_table[, CDE_pval := ifelse(CDE<0, 1, 0)]
      if(effect_table[, mean(CDE)]<0) effect_table[, CDE_pval := ifelse(CDE>0, 1, 0)]
      ## PAI and PIE for each mediator
      # for(path in paths) {
      #   effect_table[, (paste0('PAI_',path)) := get(paste0('y5_',path)) - get(paste0('y6_',path)) - y1 + y2]
      #   effect_table[, (paste0('PIE_',path)) := get(paste0('y6_',path)) - get(paste0('y4_',path))]
      # }
      # draw_table <- copy(effect_table)
      ## Make p-values for each decomposition effect.
      # for(path in paths) {
      #   ## PAI
      #   if(effect_table[, mean(get(paste0('PAI_',path)))]>=0) effect_table[, (paste0('PAI_',path,'_pval')) := ifelse(get(paste0('PAI_',path))<0, 1, 0)]
      #   if(effect_table[, mean(get(paste0('PAI_',path)))]<0) effect_table[, (paste0('PAI_',path,'_pval')) := ifelse(get(paste0('PAI_',path))>0, 1, 0)]
      #   ## PIE
      #   if(effect_table[, mean(get(paste0('PIE_',path)))]>=0) effect_table[, (paste0('PIE_',path,'_pval')) := ifelse(get(paste0('PIE_',path))<0, 1, 0)]
      #   if(effect_table[, mean(get(paste0('PIE_',path)))]<0) effect_table[, (paste0('PIE_',path,'_pval')) := ifelse(get(paste0('PIE_',path))>0, 1, 0)]
      # }
      ## For each mediator: PNDE (Y3-Y4), TNIE (Y5-Y3), PNIE (Y6-Y4)
      for(path in paths) {
        effect_table[, (paste0('PNDE_',path)) := get(paste0('y3_',path)) - get(paste0('y4_',path))]
        effect_table[, (paste0('TNIE_',path)) := get(paste0('y5_',path)) - get(paste0('y3_',path))]
        effect_table[, (paste0('PNIE_',path)) := get(paste0('y6_',path)) - get(paste0('y4_',path))]
        effect_table[, (paste0('INTref_',path)) := get(paste0('PNDE_',path)) - CDE]
        effect_table[, (paste0('INTmed_',path)) := get(paste0('TNIE_',path)) - get(paste0('PNIE_',path))]
        effect_table[, (paste0('PAI_',path)) := get(paste0('INTref_',path)) + get(paste0('INTmed_',path))]
        ## Make p-values
        for(v in c('PNDE_','TNIE_','PNIE_','INTref_','INTmed_','PAI_')) {
          if(effect_table[, mean(get(paste0(v,path)))]>=0) effect_table[, (paste0(v,path,'_pval')) := ifelse(get(paste0(v,path))<0, 1, 0)]
          if(effect_table[, mean(get(paste0(v,path)))]<0) effect_table[, (paste0(v,path,'_pval')) := ifelse(get(paste0(v,path))>0, 1, 0)]
        }
      }
      draw_table <- copy(effect_table)
      ## Summarize point estimates (mean of effect estimates), 95% (quantile of effect estimates), and p-value (mean) by collapsing over rows (sims/bootstraps)
      means <- effect_table[, lapply(.SD,mean), .SDcols=c('ATE','CDE',
                                                          paste0('PNDE_',paths),
                                                          paste0('TNIE_',paths),
                                                          paste0('PNIE_',paths),
                                                          paste0('INTref_',paths),
                                                          paste0('INTmed_',paths),
                                                          paste0('PAI_',paths))]
      means <- suppressWarnings(melt(means,value.name='mean',variable.name='effect'))
      lowers <- effect_table[, lapply(.SD,quantile,probs=0.025), .SDcols=c('ATE','CDE',
                                                                           paste0('PNDE_',paths),
                                                                           paste0('TNIE_',paths),
                                                                           paste0('PNIE_',paths),
                                                                           paste0('INTref_',paths),
                                                                           paste0('INTmed_',paths),
                                                                           paste0('PAI_',paths))]
      lowers <- suppressWarnings(melt(lowers,value.name='lower',variable.name='effect'))
      uppers <- effect_table[, lapply(.SD,quantile,probs=0.975), .SDcols=c('ATE','CDE',
                                                                           paste0('PNDE_',paths),
                                                                           paste0('TNIE_',paths),
                                                                           paste0('PNIE_',paths),
                                                                           paste0('INTref_',paths),
                                                                           paste0('INTmed_',paths),
                                                                           paste0('PAI_',paths))]
      uppers <- suppressWarnings(melt(uppers,value.name='upper',variable.name='effect'))
      pvalues <- effect_table[, lapply(.SD,mean), .SDcols=c('ATE_pval','CDE_pval',
                                                            paste0('PNDE_',paths,'_pval'),
                                                            paste0('TNIE_',paths,'_pval'),
                                                            paste0('PNIE_',paths,'_pval'),
                                                            paste0('INTref_',paths,'_pval'),
                                                            paste0('INTmed_',paths,'_pval'),
                                                            paste0('PAI_',paths,'_pval'))]
      pvalue <- suppressWarnings(melt(pvalues,value.name='pvalue',variable.name='effect'))
      pvalue[, effect := gsub('_pval','',effect)]
      full_table <- Reduce(merge,list(means,lowers,uppers,pvalue))
    }
    if(decomp_type=='2way') {
      paths <- c(treatment_var,names(decomp_paths[[intervention_course]]$paths))
      for(path in paths) {
        path_course <- readRDS(paste0(sim_dir,'/', intervention_course, '_', path, '_', total_sims, '.RDS'))
        path_course[, (path) := NULL]
        setnames(path_course, outcome_var, path)
        path_course <- path_course[, c('sim',path), with=F]
        effect_table <- merge(effect_table, path_course, by=c('sim'))
        effect_table[, (paste0('NIE_',path)) := get(path) - compare_outcome]
        if(effect_table[, mean(get(paste0('NIE_',path)))]>=0) effect_table[, (paste0('NIE_',path,'_pval')) := ifelse(get(paste0('NIE_',path))<0, 1, 0)]
        if(effect_table[, mean(get(paste0('NIE_',path)))]<0) effect_table[, (paste0('NIE_',path,'_pval')) := ifelse(get(paste0('NIE_',path))>0, 1, 0)]
      }
      draw_table <- copy(effect_table)
      means <- effect_table[, lapply(.SD,mean), .SDcols=c('ATE',paste0('NIE_',paths))]
      means <- suppressWarnings(melt(means,value.name='mean',variable.name='effect'))
      lowers <- effect_table[, lapply(.SD,quantile,probs=0.025), .SDcols=c('ATE',paste0('NIE_',paths))]
      lowers <- suppressWarnings(melt(lowers,value.name='lower',variable.name='effect'))
      uppers <- effect_table[, lapply(.SD,quantile,probs=0.975), .SDcols=c('ATE',paste0('NIE_',paths))]
      uppers <- suppressWarnings(melt(uppers,value.name='upper',variable.name='effect'))
      pvalues <- effect_table[, lapply(.SD,mean), .SDcols=c('ATE_pval',paste0('NIE_',paths,'_pval'))]
      pvalue <- suppressWarnings(melt(pvalues,value.name='pvalue',variable.name='effect'))
      pvalue[, effect := gsub('_pval','',effect)]
      ses <- effect_table[, lapply(.SD,sd), .SDcols=c('ATE',paste0('NIE_',paths))]
      ses <- suppressWarnings(melt(ses,value.name='se',variable.name='effect'))
      full_table <- Reduce(merge,list(means,lowers,uppers,pvalue,ses))
      ATE <- full_table[effect=='ATE', mean]
      full_table[, percent_ATE := mean / ATE]
      full_table[grepl(treatment_var,effect), effect := gsub('NIE','NDE',effect)]
    }
    ## Make sure effects from decomposition add up to ATE
    message(paste0('ATE          = ',sum(full_table[effect=='ATE',mean])))
    message(paste0('CDE+PAI+PNIE = ',sum(full_table[grepl('CDE|PAI_|PNIE_',effect),mean])))
    out <- list(full_table,draw_table)
  }
  if(scale=='mult') {
    ## ATE p-values
    effect_table[, ATE := int_outcome / compare_outcome]
    if(effect_table[, mean(ATE)]>=0) effect_table[, ATE_pval := ifelse(ATE<0, 1, 0)]
    if(effect_table[, mean(ATE)]<0) effect_table[, ATE_pval := ifelse(ATE>0, 1, 0)]
    ## Merge on all expected values of Y to calculate four-way decomposition of ATE.
    if(decomp_type=='4way') {
      for(path in c('y1','y2',paste0('y3_',paths),paste0('y4_',paths),paste0('y5_',paths),paste0('y6_',paths))) {
        path_course <- readRDS(paste0(sim_dir,'/', intervention_course, '_', path, '_', total_sims, '.RDS'))
        setnames(path_course, outcome_var, path)
        path_course <- path_course[, c('sim',path), with=F]
        effect_table <- merge(effect_table, path_course, by=c('sim'))
      }
      ## CDE
      effect_table[, CDE := y1 / y2]
      if(effect_table[, mean(CDE)]>=0) effect_table[, CDE_pval := ifelse(CDE<0, 1, 0)]
      if(effect_table[, mean(CDE)]<0) effect_table[, CDE_pval := ifelse(CDE>0, 1, 0)]
      effect_table[, erCDE := (y1 - y2) / compare_outcome]
      if(effect_table[, mean(erCDE)]>=0) effect_table[, erCDE_pval := ifelse(CDE<0, 1, 0)]
      if(effect_table[, mean(erCDE)]<0) effect_table[, erCDE_pval := ifelse(CDE>0, 1, 0)]
      ## PAI and PIE for each mediator
      # for(path in paths) {
      #   effect_table[, (paste0('PAI_',path)) := get(paste0('y5_',path)) - get(paste0('y6_',path)) - y1 + y2]
      #   effect_table[, (paste0('PIE_',path)) := get(paste0('y6_',path)) - get(paste0('y4_',path))]
      # }
      # draw_table <- copy(effect_table)
      ## Make p-values for each decomposition effect.
      # for(path in paths) {
      #   ## PAI
      #   if(effect_table[, mean(get(paste0('PAI_',path)))]>=0) effect_table[, (paste0('PAI_',path,'_pval')) := ifelse(get(paste0('PAI_',path))<0, 1, 0)]
      #   if(effect_table[, mean(get(paste0('PAI_',path)))]<0) effect_table[, (paste0('PAI_',path,'_pval')) := ifelse(get(paste0('PAI_',path))>0, 1, 0)]
      #   ## PIE
      #   if(effect_table[, mean(get(paste0('PIE_',path)))]>=0) effect_table[, (paste0('PIE_',path,'_pval')) := ifelse(get(paste0('PIE_',path))<0, 1, 0)]
      #   if(effect_table[, mean(get(paste0('PIE_',path)))]<0) effect_table[, (paste0('PIE_',path,'_pval')) := ifelse(get(paste0('PIE_',path))>0, 1, 0)]
      # }
      ## For each mediator: PNDE (Y3-Y4), TNIE (Y5-Y3), PNIE (Y6-Y4)
      for(path in paths) {
        effect_table[, (paste0('PNDE_',path)) := get(paste0('y3_',path)) / get(paste0('y4_',path))]
        effect_table[, (paste0('TNIE_',path)) := get(paste0('y5_',path)) / get(paste0('y3_',path))]
        effect_table[, (paste0('PNIE_',path)) := get(paste0('y6_',path)) / get(paste0('y4_',path))]
        effect_table[, (paste0('INTref_',path)) := get(paste0('PNDE_',path)) - 1 - erCDE]
        effect_table[, (paste0('INTmed_',path)) := get(paste0('TNIE_',path)) * get(paste0('PNDE_',path)) - get(paste0('PNDE_',path)) - get(paste0('PNIE_',path)) + 1]
        effect_table[, (paste0('PAI_',path)) := get(paste0('INTref_',path)) + get(paste0('INTmed_',path))]
        ## Make p-values
        for(v in c('PNDE_','TNIE_','PNIE_','INTref_','INTmed_','PAI_')) {
          if(effect_table[, mean(get(paste0(v,path)))]>=0) effect_table[, (paste0(v,path,'_pval')) := ifelse(get(paste0(v,path))<0, 1, 0)]
          if(effect_table[, mean(get(paste0(v,path)))]<0) effect_table[, (paste0(v,path,'_pval')) := ifelse(get(paste0(v,path))>0, 1, 0)]
        }
      }
      draw_table <- copy(effect_table)
      ## Summarize point estimates (mean of effect estimates), 95% (quantile of effect estimates), and p-value (mean) by collapsing over rows (sims/bootstraps)
      means <- effect_table[, lapply(.SD,mean), .SDcols=c('ATE','CDE',
                                                          paste0('PNDE_',paths),
                                                          paste0('TNIE_',paths),
                                                          paste0('PNIE_',paths),
                                                          paste0('INTref_',paths),
                                                          paste0('INTmed_',paths),
                                                          paste0('PAI_',paths))]
      means <- suppressWarnings(melt(means,value.name='mean',variable.name='effect'))
      lowers <- effect_table[, lapply(.SD,quantile,probs=0.025), .SDcols=c('ATE','CDE',
                                                                           paste0('PNDE_',paths),
                                                                           paste0('TNIE_',paths),
                                                                           paste0('PNIE_',paths),
                                                                           paste0('INTref_',paths),
                                                                           paste0('INTmed_',paths),
                                                                           paste0('PAI_',paths))]
      lowers <- suppressWarnings(melt(lowers,value.name='lower',variable.name='effect'))
      uppers <- effect_table[, lapply(.SD,quantile,probs=0.975), .SDcols=c('ATE','CDE',
                                                                           paste0('PNDE_',paths),
                                                                           paste0('TNIE_',paths),
                                                                           paste0('PNIE_',paths),
                                                                           paste0('INTref_',paths),
                                                                           paste0('INTmed_',paths),
                                                                           paste0('PAI_',paths))]
      uppers <- suppressWarnings(melt(uppers,value.name='upper',variable.name='effect'))
      pvalues <- effect_table[, lapply(.SD,mean), .SDcols=c('ATE_pval','CDE_pval',
                                                            paste0('PNDE_',paths,'_pval'),
                                                            paste0('TNIE_',paths,'_pval'),
                                                            paste0('PNIE_',paths,'_pval'),
                                                            paste0('INTref_',paths,'_pval'),
                                                            paste0('INTmed_',paths,'_pval'),
                                                            paste0('PAI_',paths,'_pval'))]
      pvalue <- suppressWarnings(melt(pvalues,value.name='pvalue',variable.name='effect'))
      pvalue[, effect := gsub('_pval','',effect)]
      full_table <- Reduce(merge,list(means,lowers,uppers,pvalue))
    }
    if(decomp_type=='2way') {
      paths <- c(treatment_var,names(decomp_paths[[intervention_course]]$paths))
      for(path in paths) {
        path_course <- readRDS(paste0(sim_dir,'/', intervention_course, '_', path, '_', total_sims, '.RDS'))
        path_course[, (path) := NULL]
        setnames(path_course, outcome_var, path)
        path_course <- path_course[, c('sim',path), with=F]
        effect_table <- merge(effect_table, path_course, by=c('sim'))
        effect_table[, (paste0('NIE_',path)) := get(path) - compare_outcome]
        if(effect_table[, mean(get(paste0('NIE_',path)))]>=0) effect_table[, (paste0('NIE_',path,'_pval')) := ifelse(get(paste0('NIE_',path))<0, 1, 0)]
        if(effect_table[, mean(get(paste0('NIE_',path)))]<0) effect_table[, (paste0('NIE_',path,'_pval')) := ifelse(get(paste0('NIE_',path))>0, 1, 0)]
      }
      draw_table <- copy(effect_table)
      means <- effect_table[, lapply(.SD,mean), .SDcols=c('ATE',paste0('NIE_',paths))]
      means <- suppressWarnings(melt(means,value.name='mean',variable.name='effect'))
      lowers <- effect_table[, lapply(.SD,quantile,probs=0.025), .SDcols=c('ATE',paste0('NIE_',paths))]
      lowers <- suppressWarnings(melt(lowers,value.name='lower',variable.name='effect'))
      uppers <- effect_table[, lapply(.SD,quantile,probs=0.975), .SDcols=c('ATE',paste0('NIE_',paths))]
      uppers <- suppressWarnings(melt(uppers,value.name='upper',variable.name='effect'))
      pvalues <- effect_table[, lapply(.SD,mean), .SDcols=c('ATE_pval',paste0('NIE_',paths,'_pval'))]
      pvalue <- suppressWarnings(melt(pvalues,value.name='pvalue',variable.name='effect'))
      pvalue[, effect := gsub('_pval','',effect)]
      ses <- effect_table[, lapply(.SD,sd), .SDcols=c('ATE',paste0('NIE_',paths))]
      ses <- suppressWarnings(melt(ses,value.name='se',variable.name='effect'))
      full_table <- Reduce(merge,list(means,lowers,uppers,pvalue,ses))
      ATE <- full_table[effect=='ATE', mean]
      full_table[, percent_ATE := mean / ATE]
      full_table[grepl(treatment_var,effect), effect := gsub('NIE','NDE',effect)]
    }
    ## Make sure effects from decomposition add up to ATE
    message(paste0('ATE          = ',sum(full_table[effect=='ATE',mean])))
    message(paste0('CDE+PAI+PNIE = ',sum(full_table[grepl('CDE|PAI_|PNIE_',effect),mean])))
    out <- list(full_table,draw_table)
  }
  return(out)
}

multmed_effect_plot <- function(effect_table, clean_names) {
  for(e in c('ATE','CDE','PAI','PNIE')) effect_table[grepl(e,effect), effect_type := e]
  effect_table[, effect_type := factor(effect_type, levels=c('CDE','PNIE','PAI'))]
  effect_table <- merge(effect_table, clean_names, by='effect', all.x=T)
  effect_table[, effect_clean := factor(effect_clean, levels=clean_names[,unique(rev(effect_clean))])]
  ATE <- effect_table[effect=='ATE',mean]
  gg <- ggplot(data=effect_table[effect_type!='ATE' & !is.na(effect_clean),]) + 
    geom_hline(yintercept=0,color='black') +
    geom_segment(aes(x=effect_clean,
                     xend=effect_clean,
                     y=lower,
                     yend=upper,
                      color=effect_type),
                  size=2, color='black') +
    geom_point(aes(x=effect_clean,
                   y=mean,
                   fill=effect_type),
               shape=21,
               size=8, stroke=1.5) +
    labs(y=paste0('Decomposed effect estimate'),x='',fill='Effect\ntype') +
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
  message(paste0('% Explained by mediators: ', 100-round(effect_table[effect=='CDE',mean/ATE*100])))
  return(gg)
}
