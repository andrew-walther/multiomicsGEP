# ============================================================
# Script: run_fixed_q_loading_case.R
# Purpose: Check one fixed-q prior change against the saved profiled Single cell021.
# Created: 2026-10-16
# Dependencies: existing panel modules and base R numerical integration
# ============================================================
# Task2-9: commit the prospective design/code before running once, serially.
# One simulated fit only. Original seed20261026/noise.25,n400/120,p50,K7,
# profileTRUE,w_S1,Normal beta,no supervised start/pruning,300 cap/ELBO rule.
# No patient/external data, alternate cells, automatic retry or adoption panel.
source("results/benchmark_sim/run_round2_fixed_prior.R")

#' Add the sole analytical change to the saved profiled controls
#' @return Existing controls plus loading_prior_update="fixed_q_em".
#' @examples
#' fixed_q_loading_case_controls()
fixed_q_loading_case_controls <- function() {
  a <- fixed_prior_controls("single","profiled")
  a$loading_prior_update <- "fixed_q_em"
  a
}

#' Integrate a point-Laplace posterior without inverse-Mills approximations
#'
#' Each half-line uses u=|F|sqrt(A), kernel exp(h*u-u²/2). Recenter at
#' positive h; rescale negative h to avoid cancellation. Gaussian integration
#' omits only tails below machine-underflow scale. Tolerance1e-11 controls
#' quadrature only, never fitting or scientific adoption.
#' @param A,B Finite scalar precision A>0 and linear signal.
#' @param prior Incoming point-Laplace prior with pi>0 and finite positive rate.
#' @return Integrated mean, second, mean_abs, prob and log_ml.
#' @examples
#' fixed_q_integrated_posterior(4,2,list(pi=.4,rate=1.3))
fixed_q_integrated_posterior <- function(A,B,prior) {
  stopifnot(length(A)==1L,is.finite(A),A>0,length(B)==1L,is.finite(B),
    prior$pi>0,prior$pi<=1,is.finite(prior$rate),prior$rate>0)
  s <- 1/sqrt(A)
  pieces <- lapply(c(B-prior$rate,-B-prior$rate)*s,function(h) {
    scale <- if(h<0)max(1,-h) else 1
    lower <- if(h<0)0 else -min(h,sqrt(-2*log(.Machine$double.xmin)))
    v <- vapply(0:2,function(r) {
      fn <- if(h<0)function(t)t^r*exp(h*t/scale-t^2/(2*scale^2)) else
        function(t)t^r*exp(-t^2/2)
      integrate(fn,lower,0,rel.tol=1e-11)$value+integrate(fn,0,Inf,rel.tol=1e-11)$value
    },numeric(1))
    if(h<0)list(logJ=log(v[1]/scale),mean=s*v[2]/v[1]/scale,
      second=s^2*v[3]/v[1]/scale^2) else
      list(logJ=h^2/2+log(v[1]),mean=s*(h+v[2]/v[1]),
        second=s^2*(h^2+2*h*v[2]/v[1]+v[3]/v[1]))
  })
  ll <- c(if(prior$pi<1)log1p(-prior$pi)+dnorm(B/A,0,s,log=TRUE) else -Inf,
    log(prior$pi)+log(prior$rate/2)-log(2*pi)/2-(B*s)^2/2+
      vapply(pieces,`[[`,numeric(1),"logJ"))
  log_ml <- max(ll)+log(sum(exp(ll-max(ll))));w <- exp(ll-log_ml)
  list(mean=w[2]*pieces[[1]]$mean-w[3]*pieces[[2]]$mean,
    second=w[2]*pieces[[1]]$second+w[3]*pieces[[2]]$second,
    mean_abs=w[2]*pieces[[1]]$mean+w[3]*pieces[[2]]$mean,prob=w[2]+w[3],log_ml=log_ml)
}

#' Report one fit under the original planted-program and prediction rules
#' @param fit Single fit with posterior moments and native history.
#' @param tb Original simulation test bed, including truth and validation data.
#' @return One row retaining convergence, C, recovery/sign, false programs and norms.
#' @examples
#' # fixed_q_case_metrics(fit,tb)
#' @seealso [fixed_prior_cell()]
fixed_q_case_metrics <- function(fit,tb) {
  risk <- predict_cox_on_yf(tb$Yva,fit$EF,fit$EBeta,fit$EF_norms)$risk_scores
  rr <- suppressWarnings(cor(tb$F_true[,3],fit$EF));rr[!is.finite(rr)] <- 0
  k <- which.max(abs(rr));rho <- rr[k]
  sd_beta <- sqrt(pmax(fit$EBeta2-fit$EBeta^2,0))
  active <- sd_beta>0 & abs(fit$EBeta)>=1.96*sd_beta
  data.frame(iterations=fit$history$n_iter,converged=fit$history$converged,
    validation_C=frozen_reverse_cindex(risk,tb$tva,tb$sva),risk_sd=stats::sd(risk),
    cor_planted=abs(rho),recovered=abs(rho)>=.7 && sign(rho*fit$EBeta[k])==1,
    n_active=sum(active),false_active=sum(active[setdiff(seq_along(active),if(abs(rho)>=.7)k else integer())]),
    max_loading_norm=max(sqrt(colSums(fit$EF^2))),max_activity_norm=max(sqrt(colSums(fit$EL^2))))
}

if(sys.nframe()==0L) {
  if(length(commandArgs(TRUE)))stop("This registered case accepts no alternate cells or controls.")
  out <- "results/benchmark_sim/outputs/round2/fixed_q_seed26"
  detailed <- "data/round2_fixed_q_seed26.rds"
  if(dir.exists(out)||file.exists(detailed))stop("Fresh output required; no retry or overwrite.")
  reference_path <- "results/benchmark_sim/outputs/round2/profiled/cells/cell_021.rds"
  sha <- system2("shasum",c("-a","256",reference_path),stdout=TRUE)
  if(length(sha)!=1L || !startsWith(sha,"5b83bf02f24603e9be1ebc337b7cb3f9beba9bf53d7c2a0458885d429132b820 "))
    stop("Registered saved reference hash mismatch.")
  load_panel_code();source("code/concordance_ci.R")
  paths <- c(list.files("code",pattern="\\.(R|cpp)$",full.names=TRUE),
    "config/globals.yml",reference_path,
    list.files("results/benchmark_sim",pattern="\\.R$",full.names=TRUE),
    "results/multimodal_real/run_round2_collapse_diagnosis.R")
  hashes <- tools::md5sum(paths);stopifnot(!anyNA(hashes))
  dir.create(out,recursive=TRUE)
  write.csv(data.frame(path=paths,md5=unname(hashes)),file.path(out,"input_hashes.csv"),row.names=FALSE)
  tb <- make_testbed(20261026L,.25);reference <- readRDS(reference_path)
  original <- multimodal_yfb_update_F_mk;rows <- list();selected <- list()
  #' Observe prior gains and four predetermined scalar quadrature cases
  #' @param ... Original loading-update arguments.
  #' @return Original update, unchanged.
  #' @examples
  #' # observed(...)
  observed <- function(...) {
    frame <- parent.frame();args <- list(...);g <- args[[14]]
    ans <- original(...);d <- ans$details;f <- loading_prior_fixed_q(d,g)
    stopifnot(identical(f$prior,ans$prior),identical(f$kl,ans$kl))
    rows[[length(rows)+1L]] <<- data.frame(iteration=get("iter",frame),factor=get("k",frame),
      prior_gain=f$prior_gain,kl=ans$kl,old_pi=g$pi,new_pi=ans$prior$pi,
      old_rate=g$rate,new_rate=ans$prior$rate)
    if(g$pi>0 && all(d$A>0)) {
      hp <- (d$B-g$rate)/sqrt(d$A);hm <- (-d$B-g$rate)/sqrt(d$A)
      score <- pmin(abs(hp+8),abs(hm+8))
      candidates <- list(first=1L,last=length(d$A),near_tail=which.min(score),
        largest_signal=which.max(abs(d$B)/sqrt(d$A)))
      for(nm in names(candidates)) {
        j <- candidates[[nm]];value <- if(nm=="near_tail")score[j] else abs(d$B[j])/sqrt(d$A[j])
        keep <- is.null(selected[[nm]]) || nm=="last" ||
          (nm=="near_tail" && value<selected[[nm]]$score) ||
          (nm=="largest_signal" && value>selected[[nm]]$score)
        if(keep)selected[[nm]] <<- list(A=d$A[j],B=d$B[j],prior=g,mean=d$mean[j],
          second=d$second[j],log_ml=d$log_ml[j],score=value,iteration=get("iter",frame),
          factor=get("k",frame),coordinate=j)
      }
    }
    ans
  }
  multimodal_yfb_update_F_mk <- observed
  warnings <- character();started <- proc.time()[["elapsed"]];set.seed(20261026L)
  fit <- tryCatch(withCallingHandlers(do.call(fit_cox_on_yf,c(list(tb$Ytr,tb$ttr,tb$str),
    fixed_q_loading_case_controls())),warning=function(w)warnings<<-c(warnings,conditionMessage(w))),
    error=function(e){writeLines(conditionMessage(e),file.path(out,"error.txt"));stop(e)},
    finally={multimodal_yfb_update_F_mk<-original
      write.csv(data.frame(message=warnings),file.path(out,"warnings.csv"),row.names=FALSE)
      saveRDS(list(prior_updates=rows,density_cases=selected),file.path(out,"telemetry.rds"))})
  saveRDS(fit,detailed)
  if(!identical(hashes,tools::md5sum(paths)))stop("Inputs changed during execution.")
  d <- do.call(rbind,rows)
  write.csv(data.frame(updates=nrow(d),negative_gains=sum(d$prior_gain < -1e-8),
    nonfinite=sum(!is.finite(d$prior_gain)|!is.finite(d$kl)),min_gain=min(d$prior_gain),
    min_kl=min(d$kl),warnings=length(warnings),runtime_seconds=proc.time()[["elapsed"]]-started),
    file.path(out,"prior_summary.csv"),row.names=FALSE)
  metrics <- rbind(reference=fixed_q_case_metrics(reference,tb),candidate=fixed_q_case_metrics(fit,tb))
  metrics$arm <- rownames(metrics);rownames(metrics)<-NULL
  write.csv(metrics,file.path(out,"comparison.csv"),row.names=FALSE)
  checks <- do.call(rbind,lapply(names(selected),function(nm) {
    z <- selected[[nm]];q <- fixed_q_integrated_posterior(z$A,z$B,z$prior)
    data.frame(case=nm,iteration=z$iteration,factor=z$factor,coordinate=z$coordinate,A=z$A,B=z$B,
      mean_error=z$mean-q$mean,second_error=z$second-q$second,normalizer_error=z$log_ml-q$log_ml)
  }))
  write.csv(checks,file.path(out,"density_checks.csv"),row.names=FALSE)
  invalid <- any(d$prior_gain < -1e-8 | !is.finite(d$prior_gain) | !is.finite(d$kl))
  writeLines(c("One simulation is not adoption. Keep caps/density limitations; stop without retry.",
    "Four selected scalar checks do not certify every density or a common implemented ELBO.",
    if(invalid)"INVALID local prior-step check; interpretation blocked." else
      "Finite local prior gains checked; scientific and full-objective gates remain unpassed."),
    file.path(out,"interpretation.txt"))
  if(invalid)stop("Local prior-step check failed; no retry.")
}
