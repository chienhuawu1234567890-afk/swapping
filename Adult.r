############################################################
## Real Data Illustration: Adult Census Income Dataset
## Category-Preserving Label Swaps for Categorical Outcomes
##
## A = sensitive categorical label, e.g., occupation
## Y = binary income outcome, >50K vs <=50K
############################################################

rm(list = ls())
set.seed(12345)

############################################################
## 0. Packages
############################################################

pkgs <- c("dplyr", "readr", "stringr")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = TRUE)
}

############################################################
## 1. Download and load Adult data
############################################################

adult_url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data"

col_names <- c(
  "age",
  "workclass",
  "fnlwgt",
  "education",
  "education_num",
  "marital_status",
  "occupation",
  "relationship",
  "race",
  "sex",
  "capital_gain",
  "capital_loss",
  "hours_per_week",
  "native_country",
  "income"
)

adult_raw <- read_csv(
  adult_url,
  col_names = col_names,
  na = c("?", " ?"),
  trim_ws = TRUE,
  show_col_types = FALSE
)

adult_raw <- adult_raw %>%
  na.omit()

############################################################
## 2. Construct analysis dataset
############################################################
## Sensitive label A: occupation
## Outcome Y: income category

analysis_df <- adult_raw %>%
  mutate(
    id = row_number(),
    A = factor(occupation),
    Y = factor(income, levels = c("<=50K", ">50K"))
  ) %>%
  select(id, A, Y, everything())

stopifnot(is.data.frame(analysis_df))
stopifnot("A" %in% names(analysis_df))
stopifnot("Y" %in% names(analysis_df))

cat("\nSample size:\n")
print(nrow(analysis_df))

cat("\nOriginal occupation by income table:\n")
print(table(analysis_df$A, analysis_df$Y))

############################################################
## 3. Safe Fisher test
############################################################

safe_fisher <- function(tab, B = 10000) {
  fisher.test(tab, simulate.p.value = TRUE, B = B)
}

############################################################
## 4. Original inference
############################################################

C_original <- table(analysis_df$A, analysis_df$Y)

chi_original <- chisq.test(C_original)
fisher_original <- safe_fisher(C_original, B = 10000)

cat("\nOriginal Pearson chi-square test:\n")
print(chi_original)

cat("\nOriginal simulated Fisher test:\n")
print(fisher_original)

############################################################
## 5. Proposed category-preserving Bernoulli swap
############################################################

category_preserving_swap <- function(data, A_col, Y_col,
                                     pi_swap = 0.8,
                                     T_swap = 1000) {
  
  data_star <- data
  A_star <- as.character(data_star[[A_col]])
  Y <- data_star[[Y_col]]
  
  categories <- unique(Y)
  
  for (t in seq_len(T_swap)) {
    
    k <- sample(categories, size = 1)
    idx_k <- which(Y == k)
    
    if (length(idx_k) < 2) next
    
    pair <- sample(idx_k, size = 2, replace = FALSE)
    u <- pair[1]
    v <- pair[2]
    
    if (A_star[u] != A_star[v] && runif(1) < pi_swap) {
      temp <- A_star[u]
      A_star[u] <- A_star[v]
      A_star[v] <- temp
    }
  }
  
  data_star$A_star <- factor(A_star, levels = levels(data[[A_col]]))
  return(data_star)
}

############################################################
## 6. Existing comparison methods
############################################################

no_swap <- function(data, A_col) {
  
  data_star <- data
  
  data_star$A_star <- factor(
    data_star[[A_col]],
    levels = levels(data[[A_col]])
  )
  
  return(data_star)
}

random_within_category_swap <- function(data, A_col, Y_col,
                                        T_swap = 1000) {
  
  data_star <- data
  A_star <- as.character(data_star[[A_col]])
  Y <- data_star[[Y_col]]
  
  categories <- unique(Y)
  
  for (t in seq_len(T_swap)) {
    
    k <- sample(categories, size = 1)
    idx_k <- which(Y == k)
    
    if (length(idx_k) < 2) next
    
    pair <- sample(idx_k, size = 2, replace = FALSE)
    
    temp <- A_star[pair[1]]
    A_star[pair[1]] <- A_star[pair[2]]
    A_star[pair[2]] <- temp
  }
  
  data_star$A_star <- factor(A_star, levels = levels(data[[A_col]]))
  return(data_star)
}

global_random_swap <- function(data, A_col,
                               pi_swap = 0.8,
                               T_swap = 1000) {
  
  data_star <- data
  A_star <- as.character(data_star[[A_col]])
  N <- nrow(data_star)
  
  for (t in seq_len(T_swap)) {
    
    pair <- sample(seq_len(N), size = 2, replace = FALSE)
    
    if (A_star[pair[1]] != A_star[pair[2]] && runif(1) < pi_swap) {
      temp <- A_star[pair[1]]
      A_star[pair[1]] <- A_star[pair[2]]
      A_star[pair[2]] <- temp
    }
  }
  
  data_star$A_star <- factor(A_star, levels = levels(data[[A_col]]))
  return(data_star)
}

pram_randomized_response <- function(data, A_col,
                                     p_stay = 0.7) {
  
  data_star <- data
  groups <- levels(data[[A_col]])
  A_old <- as.character(data[[A_col]])
  A_star <- A_old
  
  for (i in seq_along(A_old)) {
    
    if (runif(1) > p_stay) {
      other_groups <- setdiff(groups, A_old[i])
      A_star[i] <- sample(other_groups, size = 1)
    }
  }
  
  data_star$A_star <- factor(A_star, levels = groups)
  return(data_star)
}

full_permutation <- function(data, A_col) {
  
  data_star <- data
  
  A_star <- sample(
    as.character(data[[A_col]]),
    replace = FALSE
  )
  
  data_star$A_star <- factor(A_star, levels = levels(data[[A_col]]))
  return(data_star)
}

############################################################
## 7. Utility evaluation
############################################################

evaluate_utility <- function(C_original, C_release, method_name) {
  
  chi_original <- chisq.test(C_original)
  chi_release  <- chisq.test(C_release)
  
  fisher_original <- safe_fisher(C_original, B = 10000)
  fisher_release  <- safe_fisher(C_release, B = 10000)
  
  data.frame(
    Method = method_name,
    Chi_Original = unname(chi_original$statistic),
    Chi_Released = unname(chi_release$statistic),
    Chi_Absolute_Change =
      abs(unname(chi_release$statistic) -
            unname(chi_original$statistic)),
    Chi_P_Original = chi_original$p.value,
    Chi_P_Released = chi_release$p.value,
    Chi_P_Absolute_Change =
      abs(chi_release$p.value - chi_original$p.value),
    Fisher_P_Original = fisher_original$p.value,
    Fisher_P_Released = fisher_release$p.value,
    Fisher_P_Absolute_Change =
      abs(fisher_release$p.value - fisher_original$p.value),
    Table_Preserved =
      identical(as.matrix(C_original), as.matrix(C_release))
  )
}

############################################################
## 8. Disclosure diagnostics
############################################################

compute_disclosure_diagnostics <- function(data, release_function,
                                           method_name,
                                           B = 200, ...) {
  
  stopifnot(is.data.frame(data))
  stopifnot("A" %in% names(data))
  
  A_original <- data$A
  groups <- levels(A_original)
  N <- nrow(data)
  g <- length(groups)
  
  release_mat <- matrix(NA, nrow = N, ncol = B)
  
  for (b in seq_len(B)) {
    
    data_b <- release_function(data = data, ...)
    
    if (!"A_star" %in% names(data_b)) {
      stop("Release function must return a column named A_star.")
    }
    
    release_mat[, b] <- as.character(data_b$A_star)
  }
  
  entropy_u <- numeric(N)
  retention_u <- numeric(N)
  reachable_u <- numeric(N)
  
  for (u in seq_len(N)) {
    
    probs <- table(
      factor(release_mat[u, ], levels = groups)
    ) / B
    
    probs <- as.numeric(probs)
    
    entropy_u[u] <-
      -sum(ifelse(probs > 0, probs * log(probs), 0))
    
    retention_u[u] <-
      mean(release_mat[u, ] == as.character(A_original[u]))
    
    reachable_u[u] <- sum(probs > 0)
  }
  
  data.frame(
    Method = method_name,
    B = B,
    Average_Entropy = mean(entropy_u),
    Normalized_Entropy = mean(entropy_u) / log(g),
    Average_Retention_Probability = mean(retention_u),
    Label_Change_Rate = 1 - mean(retention_u),
    Identification_Probability = mean(retention_u),
    Reachable_Diversity = mean(reachable_u)
  )
}

############################################################
## 9. Run all methods
############################################################

pi_swap <- 0.8
T_swap <- 1000
B_release <- 200

method_list <- list(
  
  Original_No_Swap = list(
    fun = no_swap,
    args = list(A_col = "A")
  ),
  
  Proposed_Category_Preserving_Swap = list(
    fun = category_preserving_swap,
    args = list(
      A_col = "A",
      Y_col = "Y",
      pi_swap = pi_swap,
      T_swap = T_swap
    )
  ),
  
  Random_Within_Outcome_Swap = list(
    fun = random_within_category_swap,
    args = list(
      A_col = "A",
      Y_col = "Y",
      T_swap = T_swap
    )
  ),
  
  Global_Random_Swap = list(
    fun = global_random_swap,
    args = list(
      A_col = "A",
      pi_swap = pi_swap,
      T_swap = T_swap
    )
  ),
  
  PRAM_Randomized_Response = list(
    fun = pram_randomized_response,
    args = list(
      A_col = "A",
      p_stay = 0.7
    )
  ),
  
  Full_Permutation = list(
    fun = full_permutation,
    args = list(
      A_col = "A"
    )
  )
)

utility_results <- data.frame()
disclosure_results <- data.frame()
released_tables <- list()
released_data <- list()

for (method_name in names(method_list)) {
  
  cat("\nRunning method:", method_name, "\n")
  
  method <- method_list[[method_name]]
  
  data_release <- do.call(
    method$fun,
    c(list(data = analysis_df), method$args)
  )
  
  C_release <- table(data_release$A_star, data_release$Y)
  
  released_tables[[method_name]] <- C_release
  released_data[[method_name]] <- data_release
  
  utility_now <- evaluate_utility(
    C_original = C_original,
    C_release = C_release,
    method_name = method_name
  )
  
  disclosure_now <- do.call(
    compute_disclosure_diagnostics,
    c(
      list(
        data = analysis_df,
        release_function = method$fun,
        method_name = method_name,
        B = B_release
      ),
      method$args
    )
  )
  
  utility_results <- rbind(utility_results, utility_now)
  disclosure_results <- rbind(disclosure_results, disclosure_now)
}

cat("\nUtility results:\n")
print(utility_results)

cat("\nDisclosure results:\n")
print(disclosure_results)

############################################################
## 10. Proposal-budget comparison
############################################################

T_grid <- c(10, 50, 100, 500, 1000, 2000)
T_results <- data.frame()

for (T_now in T_grid) {
  
  cat("\nRunning proposal budget T =", T_now, "\n")
  
  diag_now <- compute_disclosure_diagnostics(
    data = analysis_df,
    release_function = category_preserving_swap,
    method_name = paste0("Proposed_T_", T_now),
    B = B_release,
    A_col = "A",
    Y_col = "Y",
    pi_swap = pi_swap,
    T_swap = T_now
  )
  
  T_results <- rbind(T_results, diag_now)
}

cat("\nProposal-budget comparison:\n")
print(T_results)

############################################################
## 11. Save outputs
############################################################

write.csv(
  analysis_df,
  "Adult_analysis_data.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame.matrix(C_original),
  "Adult_original_table.csv"
)

write.csv(
  utility_results,
  "Adult_utility_method_comparison.csv",
  row.names = FALSE
)

write.csv(
  disclosure_results,
  "Adult_disclosure_method_comparison.csv",
  row.names = FALSE
)

write.csv(
  T_results,
  "Adult_T_budget_results.csv",
  row.names = FALSE
)

for (method_name in names(released_tables)) {
  
  write.csv(
    as.data.frame.matrix(released_tables[[method_name]]),
    paste0("Adult_table_", method_name, ".csv")
  )
  
  write.csv(
    released_data[[method_name]],
    paste0("Adult_data_", method_name, ".csv"),
    row.names = FALSE
  )
}

cat("\nAll Adult dataset results saved successfully.\n")
