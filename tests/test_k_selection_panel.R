# ============================================================
# Script: test_k_selection_panel.R
# Purpose: Verify fixed panel coverage and model-specific adoption criteria
# Author: Andrew Walther
# Created: 2026-10-16
# Dependencies: tests/test_helpers.R; panel summary functions
# Run: Rscript -e 'source("tests/test_helpers.R"); source("tests/test_k_selection_panel.R"); report_results("panel")'
# ============================================================
if(!exists('run_test'))source('tests/test_helpers.R')
# Read definitions without running the experiment or its CLI.
for(e in parse('results/benchmark_sim/run_k_selection_panel.R'))
  if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.call(e[[3]])&&
     identical(e[[3]][[1]],as.name('function')))eval(e,.GlobalEnv)
run_test('KPanel-T1: every model/rank/noise has all prescribed seeds and arms', {
  g<-k_panel_grid()
  assert_equal(nrow(g),800L)
  assert_equal(sum(g$model=='single'),320L)
  assert_equal(sum(g$model=='multimodal'),480L)
  assert_equal(length(unique(g$seed)),20L)
  assert_equal(length(unique(g$K)),4L)
  assert_true(setequal(unique(g$arm[g$model=='single']),c('reference','protected')))
  assert_true(setequal(unique(g$arm[g$model=='multimodal']),c('reference','selected','protected')))
  assert_true(all(vapply(split(g,list(g$model,g$arm,g$K,g$noise),drop=TRUE),
    function(x)length(unique(x$seed))==20L,logical(1))))
})
run_test('KPanel-T2: use each model rule; do not mistake unavailable MM ELBO for failure', {
  g<-expand.grid(seed=1:20,model=c('single','multimodal'),arm=c('reference','protected'),
    stringsAsFactors=FALSE)
  g$K<-7L;g$noise<-.25;g$error<-NA_character_;g$recovered<-TRUE
  g$false_active<-.2;g$K_final<-3L;g$converged<-TRUE
  g$elbo_converged<-ifelse(g$model=='single',TRUE,NA);g$outputs_converged<-TRUE
  g$c_validation<-.7
  s<-k_panel_summary(g);p<-s[s$arm=='protected',]
  assert_true(all(p$recovery_pass & p$false_positive_pass & p$convergence_pass & p$rank_pass))
  bad<-g;bad$converged[bad$model=='single' & bad$arm=='protected' & bad$seed==1]<-FALSE
  s<-k_panel_summary(bad)
  assert_false(s$convergence_pass[s$model=='single' & s$arm=='protected'])
  bad<-g;bad$recovered[bad$arm=='protected' & bad$seed<=5]<-FALSE
  assert_true(all(!k_panel_summary(bad)$recovery_pass[k_panel_summary(bad)$arm=='protected']))
  bad<-g;bad$false_active[bad$arm=='protected']<-.8
  assert_true(all(!k_panel_summary(bad)$false_positive_pass[k_panel_summary(bad)$arm=='protected']))
  bad<-g;bad$K_final[bad$arm=='protected']<-4L
  assert_true(all(!k_panel_summary(bad)$rank_pass[k_panel_summary(bad)$arm=='protected']))
  bad<-g[-1,];s<-k_panel_summary(bad)
  assert_false(s$complete[s$model==g$model[1] & s$arm==g$arm[1]])
})
run_test('KPanel-T3: entry point writes the fixed design without changing inherited bindings', {
  dest<-tempfile('panel_entry_')
  # Zero workers stops before any fit. The design must already have been saved,
  # and an inherited locked graphics::grid must not intercept configuration.
  err<-tryCatch(run_k_selection_panel(0L,output_dir=dest),error=function(e)e)
  assert_true(inherits(err,'error'))
  assert_true(file.exists(file.path(dest,'design.csv')),
    msg='entry reached worker setup rather than failing on a locked binding')
  design<-read.csv(file.path(dest,'design.csv'))
  assert_equal(nrow(design),800L)
  assert_false(grepl('locked binding',conditionMessage(err),fixed=TRUE))
  unlink(dest,recursive=TRUE)
})
run_test('KPanel-T4: a real unpruned fit retains metrics when pruning history is absent', {
  load_panel_code()
  env<-new.env(parent=.GlobalEnv)
  env$grid<-data.frame(seed=42L,noise=.25,K=3L,model='single',arm='reference')
  env$cfg<-list(cavi=list(max_iter=2L,tol=1e-5))
  env$out_dir<-tempfile('panel_real_cell_');dir.create(file.path(env$out_dir,'cells'),recursive=TRUE)
  cell<-panel_cell;environment(cell)<-env
  result<-cell(1L)
  assert_true(is.na(result$error),msg=paste('real fit reporting error:',result$error))
  assert_equal(result$pruning_caps,0L)
  assert_true(is.finite(result$K_final) && result$K_final==3L)
  assert_true(file.exists(file.path(env$out_dir,'cells','cell_0001.csv')))
  unlink(env$out_dir,recursive=TRUE)
})
run_test('KPanel-T5: posterior-state refitting is forwarded only when requested', {
  load_panel_code()
  env<-new.env(parent=.GlobalEnv)
  env$grid<-data.frame(seed=42L,noise=.25,K=3L,model='single',arm='protected')
  env$cfg<-list(cavi=list(max_iter=2L,tol=1e-5))
  env$out_dir<-tempfile('panel_posterior_');dir.create(file.path(env$out_dir,'cells'),recursive=TRUE)
  original<-fit_cox_on_yf;seen<-logical()
  # Test the adapter separately from posterior-state fitting, covered by its
  # numerical tests. Removing the observed flag keeps this a tiny old-API fit.
  env$fit_cox_on_yf<-function(...) {
    a<-list(...);seen<<-c(seen,isTRUE(a$prune_preserve_posterior))
    a$prune_preserve_posterior<-NULL;do.call(original,a)
  }
  cell<-panel_cell;environment(cell)<-env
  off<-cell(1L);on<-cell(1L,posterior_start=TRUE)
  assert_true(is.na(off$error) && is.na(on$error),msg='both adapter paths complete')
  assert_true(any(seen) && !seen[1],msg='posterior state is opt-in at the fitter boundary')
  unlink(env$out_dir,recursive=TRUE)
})
run_test('KPanel-T6: posterior panel refuses an existing attempt before any fit', {
  source('results/benchmark_sim/run_posterior_k_panel.R')
  p<-tempfile('posterior_panel_exists_');dir.create(p)
  writeLines('preserve',file.path(p,'checkpoint'))
  e<-tryCatch(run_posterior_k_panel(out_dir=p),error=function(e)e)
  assert_true(inherits(e,'error'),msg='existing output refused')
  assert_true(identical(readLines(file.path(p,'checkpoint')),'preserve'),msg='old checkpoint untouched')
  unlink(p,recursive=TRUE)
})
