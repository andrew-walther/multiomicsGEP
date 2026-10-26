# Small simulated data sets shared by the tests
sim_dd <- function(n = 120, p = c(expression = 50, methylation = 50), K = 3, beta = NULL, seed = 1, ...) {
  sim <- ssbmf_simulate(n = n, p = p, K = K, beta = beta, seed = seed, ...)
  list(sim = sim, dd = ssbmf_data(sim$X, sim$time, sim$status))
}

# A state after the setup pass and a few sweeps, for testing single updates
state_after <- function(dd, K, sweeps = 3, control = ssbmf_control()) {
  d <- ssbmf:::.prep_data(dd, control)
  st <- ssbmf:::.init_state(d, K)
  st <- ssbmf:::.update_mu(st, d); st <- ssbmf:::.update_tau(st, d)
  for (k in seq_len(K)) st <- ssbmf:::.update_L(st, d, k, safeguard = FALSE)
  for (m in seq_len(d$M)) for (k in seq_len(K)) st <- ssbmf:::.update_F(st, d, m, k, safeguard = FALSE)
  st <- ssbmf:::.objective(st, d)
  for (i in seq_len(sweeps)) {
    st <- ssbmf:::.update_mu(st, d); st <- ssbmf:::.update_tau(st, d)
    for (k in seq_len(K)) st <- ssbmf:::.update_L(st, d, k)
    for (m in seq_len(d$M)) for (k in seq_len(K)) st <- ssbmf:::.update_F(st, d, m, k)
    st <- ssbmf:::.update_beta(st, d)
  }
  list(st = st, d = d)
}
