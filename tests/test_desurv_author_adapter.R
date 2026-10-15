# ============================================================
# Script: test_desurv_author_adapter.R
# Purpose: Verify the pinned DeSurv pipeline adapter against author APIs.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: DeSurv in the isolated Round2 library, existing test helpers
# ============================================================
if(!exists("run_test"))source("tests/test_helpers.R")
.libPaths(c(normalizePath("data/round2_R_library"),.libPaths()))
if(file.exists("results/benchmark_sim/desurv_author_adapter.R"))
  source("results/benchmark_sim/desurv_author_adapter.R")
run_test("DeSurvAuthor-T1: training union and filtering match the authors' exact preprocessing", {
  mk <- function(ex,ds)list(ex=ex,featInfo=rownames(ex),
    sampInfo=data.frame(ID=colnames(ex),time=seq_len(ncol(ex))+1,
      event=rep(c(1,0),length.out=ncol(ex)),dataset=ds,keep=1))
  a <- matrix(c(1,7,2,3,9,4,5,11,6,7,13,8),3,4,
    dimnames=list(c("g1","g1","g2"),paste0("a",1:4)))
  b <- matrix(seq_len(12),3,4,dimnames=list(c("g1","g2","g3"),paste0("b",1:4)))
  x <- list(A=mk(a,"A"),B=mk(b,"B"));x$B$sampInfo$keep[4] <- 0
  got <- desurv_author_prepare(x)
  merged <- matrix(0,3,7,dimnames=list(c("g1","g2","g3"),c(colnames(a),colnames(b)[1:3])))
  merged[c("g1","g1","g2"),1:4] <- a
  merged[c("g1","g2","g3"),5:7] <- b[,1:3]
  expected <- DeSurv::preprocess_data(merged,c(2:5,2:4),c(1,0,1,0,1,0,1),
    c(rep("A",4),rep("B",3)),ngene=3000L,method_trans_train="rank",verbose=FALSE)
  assert_true(identical(got$ex,expected$ex),msg="same union, duplicate assignment, sample filter and within-sample ranks")
  assert_true(all(got$ex>=0) && identical(colnames(got$ex),colnames(merged)),
    msg="nonnegative method input preserves the exact patients and order")
})
run_test("DeSurvAuthor-T2: validation zero-filling precedes the documented rank transformation", {
  ex <- matrix(c(4,2,8,6),2,2,dimnames=list(c("g1","g2"),c("s1","s2")))
  raw <- list(ex=ex,featInfo=rownames(ex),sampInfo=data.frame(ID=colnames(ex),
    time=c(1,2),event=c(1,0),dataset="test",keep=1))
  got <- desurv_author_prepare(list(test=raw),genes=c("g1","g2","missing"))
  expected <- DeSurv::preprocess_data(rbind(ex,missing=c(0,0)),c(1,2),c(1,0),
    c("test","test"),genes=c("g1","g2","missing"),method_trans_train="rank",verbose=FALSE)
  assert_true(identical(got$ex,expected$ex),msg="missing genes are zero-filled before the authors' per-sample ranks")
})
run_test("DeSurvAuthor-T3: the paper's top-gene union projection uses fixed Cox orientation", {
  W <- matrix(c(3,1,2,0,1,4,0,2),4,2,dimnames=list(paste0("g",1:4),c("P1","P2")))
  X <- matrix(seq_len(12),4,3,dimnames=list(rownames(W),paste0("s",1:3)))
  fit <- list(W=W,beta=c(-2,1))
  genes <- unique(unlist(DeSurv::desurv_get_top_genes(W,1L)$top_genes))
  expected <- as.vector(t(X[genes,,drop=FALSE])%*%W[genes,,drop=FALSE]%*%fit$beta)
  assert_near(desurv_author_risk(fit,X,1L),expected,tol=0,msg="exact union-of-top-genes X' W beta, including negative coefficients")
})
