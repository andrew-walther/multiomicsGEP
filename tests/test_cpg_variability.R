# ============================================================
# Script: test_cpg_variability.R
# Purpose: Verify observed-value CpG moments, descriptive filters and plot counts.
# Author: Codex (reviewed by Andrew Walther)
# Created: 2026-10-15
# Dependencies: base R, existing test helpers, ggplot2
# ============================================================
if (!exists("run_test")) source("tests/test_helpers.R")
if (file.exists("results/multimodal_real/summarize_cpg_variability.R"))
  source("results/multimodal_real/summarize_cpg_variability.R")

run_test("CpG-T1: means and sample variances match independent row calculations", {
  x <- rbind(a=c(0,.2,.5,1), b=c(.1,NA,.4,.8))
  z <- cpg_variability_table(x)
  assert_near(z$mean_beta, apply(x,1,mean,na.rm=TRUE), tol=1e-14)
  assert_near(z$var_beta, apply(x,1,var,na.rm=TRUE), tol=1e-14)
  assert_near(z$sd_beta, apply(x,1,sd,na.rm=TRUE), tol=1e-14)
  assert_equal(z$probe, rownames(x))
  assert_true(identical(z$n_observed,c(4L,3L)))
})
run_test("CpG-T2: constant and missing probes retain explicit undefined SDs", {
  x <- rbind(zero=rep(0,4),constant=rep(.1,4),one=rep(1,4),
    empty=rep(NA_real_,4),singleton=c(.5,NA,NA,NA))
  z <- cpg_variability_table(x)
  assert_true(all(z$sd_beta[1:3]==0) && all(z$constant[1:3]))
  assert_true(all(is.na(z$sd_beta[4:5])) && is.na(z$mean_beta[4]))
  assert_true(is.na(z$keep_mean[4]) && is.na(z$keep_proportion[4]))
  assert_true(z$n_observed[5]==1L && z$mean_beta[5]==.5 && z$n_missing[5]==3L)
})
run_test("CpG-T3: mean filter includes .1 and .9, extreme proportions use strict boundaries", {
  x <- rbind(low=rep(.1,20),high=rep(.9,20),below=rep(.099,20),above=rep(.901,20))
  z <- cpg_variability_table(x)
  assert_equal(z$keep_mean,c(TRUE,TRUE,FALSE,FALSE))
  assert_true(all(z$prop_below_01[1:2]==0) && all(z$prop_above_09[1:2]==0))
  assert_equal(z$keep_proportion,c(TRUE,TRUE,FALSE,FALSE))
})
run_test("CpG-T4: mean and 95-percent proportion filters are not conflated", {
  x <- rbind(a=c(rep(.099,19),1),b=c(rep(.89,10),rep(.93,10)))
  z <- cpg_variability_table(x)
  assert_equal(z$keep_mean,c(TRUE,FALSE))
  assert_equal(z$keep_proportion,c(FALSE,TRUE))
  assert_true(z$prop_below_01[1]==.95 && z$prop_outside[2]==.5)
})
run_test("CpG-T5: arcsine SD is calculated on transformed observations", {
  x <- rbind(a=c(0,.25,.75,1),b=c(.2,NA,.6,.9))
  z <- cpg_variability_table(x)
  assert_near(z$sd_asin,apply(asin(2*x-1),1,sd,na.rm=TRUE),tol=1e-14)
  assert_true(any(abs(z$sd_asin-z$sd_beta)>.1))
})
run_test("CpG-T6: block size changes memory use, not values or probe order", {
  x <- matrix(seq(0,1,length.out=35),7,5,dimnames=list(paste0("cg",1:7),NULL))
  x[2,3] <- NA
  assert_equal(cpg_variability_table(x,1L),cpg_variability_table(x,4L))
})
run_test("CpG-T7: invalid beta values and ambiguous probe identities fail", {
  x <- rbind(a=c(.2,.4),b=c(.5,.8))
  cpg_variability_table(x)
  for (bad in list(x+1,replace(x,1,Inf),unname(x),x[c(1,1),],matrix(numeric(),0,2))) {
    e <- tryCatch(cpg_variability_table(bad),error=identity)
    assert_true(inherits(e,"error"))
  }
  assert_true(inherits(tryCatch(cpg_variability_table(x,0),error=identity),"error"))
})
run_test("CpG-T8: summary counts include missing SD and screened-set overlap", {
  x <- rbind(a=rep(.2,4),b=c(0,0,0,1),c=rep(NA_real_,4),d=rep(.95,4))
  z <- cpg_variability_table(x)
  s <- summarize_cpg_filters(z,c("a","b"))
  all <- s[s$subset=="All available",]
  assert_true(all$n_probes==4L && all$n_sd_beta==3L && all$n_sd_missing==1L &&
    all$n_constant==2L && all$n_screened==2L)
  both <- s[s$subset=="Both filters",]
  assert_true(both$n_probes==2L && both$n_screened==2L)
  assert_true(inherits(tryCatch(summarize_cpg_filters(z,"absent"),error=identity),"error"))
})
run_test("CpG-T9: plot bin counts account for every eligible finite probe", {
  x <- rbind(a=c(.1,.2,.3,.4),b=c(.2,.3,.5,.7),c=rep(NA_real_,4))
  z <- cpg_variability_table(x);s <- summarize_cpg_filters(z,"a")
  g <- cpg_variability_plots(z,"a",4L)
  h <- ggplot2::ggplot_build(g$sd_histogram)$data[[1]]
  wanted <- sum(s$n_sd_beta[s$subset %in%
    c("All available","Mean filter","Proportion filter","Current screened set")])*2L
  assert_true(sum(h$count)==wanted,msg="all finite beta and arcsine observations counted")
  m <- ggplot2::ggplot_build(g$mean_sd)$data[[1]]
  assert_true(sum(m$count)==2L,msg="mean/SD figure excludes only explicitly undefined moments")
})
run_test("CpG-T10: patient alignment uses names and preserves every raw probe", {
  raw <- matrix(seq(.1,.9,length.out=12),3,4,
    dimnames=list(c("a","b","c"),c("p1","p2","p3","p4")))
  cached <- t(raw[c("a","c"),c("p4","p2"),drop=FALSE])
  got <- align_cpg_training(raw,cached)
  assert_equal(got,raw[,c("p4","p2"),drop=FALSE])
  missing <- cached;rownames(missing)[1]<-"missing"
  assert_true(inherits(tryCatch(align_cpg_training(raw,missing),error=identity),"error"))
  duplicate <- cached;rownames(duplicate)<-c("p2","p2")
  assert_true(inherits(tryCatch(align_cpg_training(raw,duplicate),error=identity),"error"))
})
run_test("CpG-T11: descriptive driver preserves inputs and refuses overwrite", {
  root <- tempfile("cpg_fixture_");dir.create(root)
  raw <- matrix(seq(.1,.9,length.out=12),3,4,
    dimnames=list(c("a","b","c"),c("p1","p2","p3","p4")))
  cached <- t(raw[c("a","c"),c("p4","p2"),drop=FALSE])
  rp<-file.path(root,"raw.rds");tp<-file.path(root,"training.rds")
  bp<-file.path(root,"benchmark.rds");out<-file.path(root,"output")
  saveRDS(list(meth=raw),rp);saveRDS(list(Y=list(methylation=cached)),tp)
  saveRDS(list(joint=list(EF=list(methylation=matrix(1,2,1)),
    training_spec=list(feature_names=list(methylation=c("a","c")),
      subject_ids=rownames(cached)))),bp)
  before <- tools::md5sum(c(rp,tp,bp))
  run_cpg_variability(rp,tp,bp,out)
  summary <- read.csv(file.path(out,"filter_summary.csv"))
  assert_true(summary$n_probes[summary$subset=="All available"]==3L)
  assert_true(file.exists(file.path(out,"cpg_sd_histogram.png")) &&
    nrow(read.csv(file.path(out,"warnings.csv")))==0L)
  assert_equal(tools::md5sum(c(rp,tp,bp)),before)
  assert_true(inherits(tryCatch(run_cpg_variability(rp,tp,bp,out),error=identity),"error"))
  unlink(root,recursive=TRUE)
})
if (sys.nframe()==0L) {
  report_results("CpG variability")
  if (.test_env$failed>0L) quit(status=1L)
}
