############################################################
# Simulation study:
# Proposed category-preserving Bernoulli swap
# vs Data Swapping, PRAM, and Randomized Response
############################################################

library(dplyr)
library(tidyr)
library(purrr)

set.seed(2026)

############################################################
# 1. Generate categorical microdata
############################################################

generate_data <- function(N = 300, g = 3, K = 4,
                          scenario = "balanced") {
  
  if (scenario == "balanced") {
    p_A <- rep(1 / g, g)
    p_Y_given_A <- matrix(1 / K, nrow = g, ncol = K)
  }
  
  if (scenario == "moderate_imbalance") {
    p_A <- c(0.50, 0.30, 0.20)[1:g]
    p_A <- p_A / sum(p_A)
    
    p_Y_given_A <- matrix(1 / K, nrow = g, ncol = K)
    p_Y_given_A[1, ] <- c(0.40, 0.30, 0.20, 0.10)[1:K]
    p_Y_given_A[2, ] <- c(0.25, 0.35, 0.25, 0.15)[1:K]
    p_Y_given_A[3, ] <- c(0.15, 0.25, 0.30, 0.30)[1:K]
    p_Y_given_A <- p_Y_given_A / rowSums(p_Y_given_A)
  }
  
  if (scenario == "severe_imbalance") {
    p_A <- c(0.70, 0.20, 0.10)[1:g]
    p_A <- p_A / sum(p_A)
    
    p_Y_given_A <- matrix(1 / K, nrow = g, ncol = K)
    p_Y_given_A[1, ] <- c(0.65, 0.20, 0.10, 0.05)[1:K]
    p_Y_given_A[2, ] <- c(0.20, 0.50, 0.20, 0.10)[1:K]
    p_Y_given_A[3, ] <- c(0.05, 0.15, 0.30, 0.50)[1:K]
    p_Y_given_A <- p_Y_given_A / rowSums(p_Y_given_A)
  }
  
  if (scenario == "sparse") {
    p_A <- rep(1 / g, g)
    
    p_Y_given_A <- matrix(1 / K, nrow = g, ncol = K)
    p_Y_given_A[1, ] <- c(0.80, 0.10, 0.08, 0.02)[1:K]
    p_Y_given_A[2, ] <- c(0.10, 0.80, 0.08, 0.02)[1:K]
    p_Y_given_A[3, ] <- c(0.10, 0.08, 0.80, 0.02)[1:K]
    p_Y_given_A <- p_Y_given_A / rowSums(p_Y_given_A)
  }
  
  if (scenario == "non_swappable") {
    p_A <- rep(1 / g, g)
    
    p_Y_given_A <- matrix(1 / K, nrow = g, ncol = K)
    p_Y_given_A[1, ] <- c(0.85, 0.10, 0.05, 0.00)[1:K]
    p_Y_given_A[2, ] <- c(0.00, 0.80, 0.15, 0.05)[1:K]
    p_Y_given_A[3, ] <- c(0.00, 0.05, 0.15, 0.80)[1:K]
    p_Y_given_A <- p_Y_given_A / rowSums(p_Y_given_A)
  }
  
  A <- sample(1:g, size = N, replace = TRUE, prob = p_A)
  
  Y <- numeric(N)
  for (u in 1:N) {
    Y[u] <- sample(1:K, size = 1, prob = p_Y_given_A[A[u], ])
  }
  
  data.frame(
    id = 1:N,
    A = A,
    Y = Y
  )
}

############################################################
# 2. Contingency table
############################################################

make_table <- function(dat, g, K) {
  tab <- table(
    factor(dat$A, levels = 1:g),
    factor(dat$Y, levels = 1:K)
  )
  as.matrix(tab)
}

############################################################
# 3. Proposed category-preserving Bernoulli swap
############################################################

proposed_swap <- function(dat, pi = 0.7, T = 100) {
  
  dat_star <- dat
  
  for (t in 1:T) {
    
    k <- sample(unique(dat_star$Y), size = 1)
    idx_k <- which(dat_star$Y == k)
    
    if (length(idx_k) < 2) next
    
    pair <- sample(idx_k, size = 2, replace = FALSE)
    u <- pair[1]
    v <- pair[2]
    
    if (dat_star$A[u] != dat_star$A[v]) {
      if (runif(1) < pi) {
        temp <- dat_star$A[u]
        dat_star$A[u] <- dat_star$A[v]
        dat_star$A[v] <- temp
      }
    }
  }
  
  dat_star
}

############################################################
# 4. Classical data swapping
############################################################

data_swap <- function(dat, swap_rate = 0.20) {
  
  dat_star <- dat
  N <- nrow(dat_star)
  n_swap <- floor(N * swap_rate / 2)
  
  for (s in 1:n_swap) {
    pair <- sample(1:N, size = 2, replace = FALSE)
    u <- pair[1]
    v <- pair[2]
    
    temp <- dat_star$A[u]
    dat_star$A[u] <- dat_star$A[v]
    dat_star$A[v] <- temp
  }
  
  dat_star
}

############################################################
# 5. PRAM
############################################################

pram <- function(dat, g, p_stay = 0.80) {
  
  dat_star <- dat
  
  P <- matrix((1 - p_stay) / (g - 1), nrow = g, ncol = g)
  diag(P) <- p_stay
  
  dat_star$A <- sapply(dat$A, function(a) {
    sample(1:g, size = 1, prob = P[a, ])
  })
  
  dat_star
}

############################################################
# 6. Randomized response label reassignment
############################################################

randomized_response <- function(dat, g, rr_prob = 0.20) {
  
  dat_star <- dat
  
  dat_star$A <- sapply(dat$A, function(a) {
    if (runif(1) < rr_prob) {
      sample(setdiff(1:g, a), size = 1)
    } else {
      a
    }
  })
  
  dat_star
}

############################################################
# 7. Utility metrics
############################################################

utility_metrics <- function(C, C_star) {
  
  suppressWarnings({
    chisq_orig <- chisq.test(C, correct = FALSE)
    chisq_star <- chisq.test(C_star, correct = FALSE)
  })
  
  data.frame(
    delta_chisq = abs(as.numeric(chisq_star$statistic) -
                        as.numeric(chisq_orig$statistic)),
    delta_p = abs(chisq_star$p.value - chisq_orig$p.value)
  )
}

############################################################
# 8. Non-swappable category proportion
############################################################

non_swappable_rate <- function(C) {
  
  n_k <- colSums(C)
  
  non_swappable <- apply(C, 2, function(x) {
    sum(x > 0) <= 1
  })
  
  sum(n_k[non_swappable]) / sum(n_k)
}

############################################################
# 9. Disclosure-risk diagnostics
############################################################

disclosure_metrics <- function(dat, release_fun, B = 200, g, K, ...) {
  
  N <- nrow(dat)
  A_orig <- dat$A
  Y_orig <- dat$Y
  
  A_star_mat <- matrix(NA, nrow = N, ncol = B)
  
  for (b in 1:B) {
    released <- release_fun(dat, ...)
    A_star_mat[, b] <- released$A
  }
  
  H_u <- numeric(N)
  S_u <- numeric(N)
  
  for (u in 1:N) {
    
    p_hat <- table(factor(A_star_mat[u, ], levels = 1:g)) / B
    
    H_u[u] <- -sum(ifelse(p_hat > 0, p_hat * log(p_hat), 0))
    S_u[u] <- mean(A_star_mat[u, ] == A_orig[u])
  }
  
  C <- make_table(dat, g, K)
  swappable_category <- apply(C, 2, function(x) sum(x > 0) >= 2)
  is_swappable_subject <- swappable_category[Y_orig]
  
  data.frame(
    H_bar = mean(H_u),
    S_bar = mean(S_u),
    H_bar_swappable = mean(H_u[is_swappable_subject]),
    S_bar_swappable = mean(S_u[is_swappable_subject]),
    H_bar_nonswappable = mean(H_u[!is_swappable_subject]),
    S_bar_nonswappable = mean(S_u[!is_swappable_subject])
  )
}

############################################################
# 10. Full-mixing benchmark for proposed method
############################################################

full_mixing_benchmark <- function(dat, g, K) {
  
  C <- make_table(dat, g, K)
  N <- nrow(dat)
  
  H_ref_u <- numeric(N)
  S_ref_u <- numeric(N)
  
  for (u in 1:N) {
    
    k <- dat$Y[u]
    a <- dat$A[u]
    
    probs <- C[, k] / sum(C[, k])
    
    H_ref_u[u] <- -sum(ifelse(probs > 0, probs * log(probs), 0))
    S_ref_u[u] <- probs[a]
  }
  
  data.frame(
    H_ref = mean(H_ref_u),
    S_ref = mean(S_ref_u)
  )
}

############################################################
# 11. Run one simulation replicate
############################################################

run_one <- function(scenario,
                    N = 300,
                    g = 3,
                    K = 4,
                    pi = 0.7,
                    T = 100,
                    B = 200) {
  
  dat <- generate_data(N = N, g = g, K = K, scenario = scenario)
  C <- make_table(dat, g, K)
  r_NS <- non_swappable_rate(C)
  ref <- full_mixing_benchmark(dat, g, K)
  
  methods <- list(
    Proposed = function(d) proposed_swap(d, pi = pi, T = T),
    DataSwap = function(d) data_swap(d, swap_rate = 0.20),
    PRAM = function(d) pram(d, g = g, p_stay = 0.80),
    RandomizedResponse = function(d) randomized_response(d, g = g, rr_prob = 0.20)
  )
  
  out <- lapply(names(methods), function(m) {
    
    release_fun <- methods[[m]]
    
    released_once <- release_fun(dat)
    C_star <- make_table(released_once, g, K)
    
    util <- utility_metrics(C, C_star)
    
    disc <- disclosure_metrics(
      dat = dat,
      release_fun = release_fun,
      B = B,
      g = g,
      K = K
    )
    
    data.frame(
      scenario = scenario,
      method = m,
      pi = ifelse(m == "Proposed", pi, NA),
      T = ifelse(m == "Proposed", T, NA),
      r_NS = r_NS,
      util,
      disc,
      ref
    )
  })
  
  bind_rows(out)
}

############################################################
# 12. Run full simulation study
############################################################

scenarios <- c(
  "balanced",
  "moderate_imbalance",
  "severe_imbalance",
  "sparse",
  "non_swappable"
)

pi_values <- c(0.5, 0.7, 0.9)
T_values <- c(10, 50, 100, 500, 1000)

n_rep <- 100
B <- 200

results <- list()
counter <- 1

for (sc in scenarios) {
  for (pi in pi_values) {
    for (T in T_values) {
      for (rep in 1:n_rep) {
        
        results[[counter]] <- run_one(
          scenario = sc,
          N = 300,
          g = 3,
          K = 4,
          pi = pi,
          T = T,
          B = B
        )
        
        counter <- counter + 1
      }
    }
  }
}

sim_results <- bind_rows(results)

############################################################
# 13. Summarize results
############################################################

summary_results <- sim_results %>%
  group_by(scenario, method, pi, T) %>%
  summarise(
    mean_delta_chisq = mean(delta_chisq, na.rm = TRUE),
    mean_delta_p = mean(delta_p, na.rm = TRUE),
    mean_H = mean(H_bar, na.rm = TRUE),
    mean_S = mean(S_bar, na.rm = TRUE),
    mean_H_swappable = mean(H_bar_swappable, na.rm = TRUE),
    mean_S_swappable = mean(S_bar_swappable, na.rm = TRUE),
    mean_r_NS = mean(r_NS, na.rm = TRUE),
    mean_H_ref = mean(H_ref, na.rm = TRUE),
    mean_S_ref = mean(S_ref, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_results)

############################################################
# 14. Select proposal budget T*
############################################################

T_selection <- summary_results %>%
  filter(method == "Proposed") %>%
  mutate(
    entropy_ratio = mean_H / mean_H_ref,
    retention_gap = abs(mean_S - mean_S_ref)
  ) %>%
  group_by(scenario, pi) %>%
  summarise(
    T_entropy_95 = min(T[entropy_ratio >= 0.95], na.rm = TRUE),
    T_retention_005 = min(T[retention_gap < 0.05], na.rm = TRUE),
    .groups = "drop"
  )

print(T_selection)

############################################################
# 15. Optional: save results
############################################################

write.csv(
  sim_results,
  "simulation_raw_results.csv",
  row.names = FALSE
)

write.csv(
  summary_results,
  "simulation_summary_results.csv",
  row.names = FALSE
)

write.csv(
  T_selection,
  "simulation_T_selection.csv",
  row.names = FALSE
)
