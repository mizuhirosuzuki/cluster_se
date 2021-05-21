packages <- c(
  "tidyverse",
  "estimatr",
  "lfe",
  "fastmatch",
  "magrittr"
)

pacman::p_load(packages, character.only = TRUE)

# ---------------------------------------- #

sim_func <- function(
  cluster_sample,
  cluster_fe,
  hetero_te,
  cluster_assignment
) {
  
  # Number of people in each cluster (population, not sample)
  NC <- 10000
  # Number of simulations
  num_sim <- 1000
  # Sample fraction (completely random, not cluster sampling)
  s_frac <- 0.01
  # Sample cluster
  s_cluster <- 100
  
  # Number of clusters in the population
  if (cluster_sample == TRUE) {
    MC <- 10000
  } else {
    MC <- 100
  }
  
  # Treatment assignment probability (can vary across clusters)
  mu <- 1/2 # mean assignment
  if (cluster_assignment == TRUE) {
    max_sigma <- 1/4 # variability of assignment probability across clusters
  } else {
    max_sigma <- 0 # variability of assignment probability across clusters
  }
  
  df <- tibble(
    prob = rep(
      runif(n = MC, min = mu - 2 * max_sigma, max = mu + 2 * max_sigma), 
      each = NC
      ),
    # Cluster index
    C = rep(seq(MC), each = NC),
    # random treatment 
    W = rbinom(MC * NC, 1, prob),
    # error term
    epsilon = rnorm(NC * MC),
    # Cluster fixed effects
    Cfe = (cluster_fe == TRUE) * rep(rnorm(n = MC), each = NC),
    # treatment effect (can be cluster-specific)
    tauC = (hetero_te == TRUE) * ((C <= MC / 2) * 2 - 1),
    # outcome
    Y = tauC * W + Cfe + epsilon
    ) %>% 
    select(Y, C, W)
  
  set.seed(123)
  pval_robust <- rep(0, num_sim)
  pval_cluster <- rep(0, num_sim)
  pval_robust_fe <- rep(0, num_sim)
  pval_cluster_fe <- rep(0, num_sim)
  for (i in seq(num_sim)) {
    
    sample_cluster <- sample(seq(MC), s_cluster, replace = FALSE)
    df_sample <- df[df$C %in% sample_cluster,] %>% 
      slice_sample(prop = s_frac)
    
    res <- felm(
      Y ~ -1 + W | 0 | 0 | C, 
      data = df_sample
      )
    pval_cluster[i] <- res$cpval["W"]
    
    res <- felm(
      Y ~ -1 + W | 0 | 0 | 0, 
      data = df_sample
      )
    pval_robust[i] <- res$rpval["W"]
    
    res <- felm(
      Y ~ -1 + W | C | 0 | C, 
      data = df_sample
      )
    pval_cluster_fe[i] <- res$cpval["W"]
    
    res <- felm(
      Y ~ -1 + W | C | 0 | 0, 
      data = df_sample
      )
    pval_robust_fe[i] <- res$rpval["W"]
    
  }
  
  return(
    list(
      mean(pval_cluster < 0.05),
      mean(pval_robust < 0.05),
      mean(pval_cluster_fe < 0.05),
      mean(pval_robust_fe < 0.05)
    )
  )

}

# Clusters not sampled, no cluster FE, not heterogeneous TE, assignment no within-cluster correlation

res_1 <- sim_func(
  cluster_sample = FALSE,
  cluster_fe = FALSE,
  hetero_te = FALSE,
  cluster_assignment = FALSE
  )

print(res_1[[1]])
print(res_1[[2]])
print(res_1[[3]])
print(res_1[[4]])

# Clusters not sampled, cluster FE, not heterogeneous TE, assignment no within-cluster correlation

res_2 <- sim_func(
  cluster_sample = FALSE,
  cluster_fe = TRUE,
  hetero_te = FALSE,
  cluster_assignment = FALSE
  )

print(res_2[[1]])
print(res_2[[2]])
print(res_2[[3]])
print(res_2[[4]])

# Clusters not sampled, no cluster FE, heterogeneous TE, assignment no within-cluster correlation

res_3 <- sim_func(
  cluster_sample = FALSE,
  cluster_fe = FALSE,
  hetero_te = TRUE,
  cluster_assignment = FALSE
  )

print(res_3[[1]])
print(res_3[[2]])
print(res_3[[3]])
print(res_3[[4]])

# Clusters not sampled, no cluster FE, no heterogeneous TE, assignment within-cluster correlation

res_4 <- sim_func(
  cluster_sample = FALSE,
  cluster_fe = FALSE,
  hetero_te = FALSE,
  cluster_assignment = TRUE
  )

print(res_4[[1]])
print(res_4[[2]])
print(res_4[[3]])
print(res_4[[4]])

# Clusters not sampled, cluster FE, heterogeneous TE, assignment no within-cluster correlation

res_5 <- sim_func(
  cluster_sample = FALSE,
  cluster_fe = TRUE,
  hetero_te = TRUE,
  cluster_assignment = FALSE
  )

print(res_5[[1]])
print(res_5[[2]])
print(res_5[[3]])
print(res_5[[4]])

# Clusters not sampled, cluster FE, not heterogeneous TE, assignment within-cluster correlation

res_6 <- sim_func(
  cluster_sample = FALSE,
  cluster_fe = TRUE,
  hetero_te = FALSE,
  cluster_assignment = TRUE
  )

print(res_6[[1]])
print(res_6[[2]])
print(res_6[[3]])
print(res_6[[4]])

# Clusters not sampled, no cluster FE, heterogeneous TE, assignment within-cluster correlation

res_7 <- sim_func(
  cluster_sample = FALSE,
  cluster_fe = FALSE,
  hetero_te = TRUE,
  cluster_assignment = TRUE
  )

print(res_7[[1]])
print(res_7[[2]])
print(res_7[[3]])
print(res_7[[4]])

# Clusters not sampled, cluster FE, heterogeneous TE, assignment within-cluster correlation

res_8 <- sim_func(
  cluster_sample = FALSE,
  cluster_fe = TRUE,
  hetero_te = TRUE,
  cluster_assignment = TRUE
  )

print(res_8[[1]])
print(res_8[[2]])
print(res_8[[3]])
print(res_8[[4]])

# Clusters sampled, no cluster FE, no heterogeneous TE, assignment no within-cluster correlation

res_9 <- sim_func(
  cluster_sample = TRUE,
  cluster_fe = FALSE,
  hetero_te = FALSE,
  cluster_assignment = FALSE
  )

print(res_9[[1]])
print(res_9[[2]])
print(res_9[[3]])
print(res_9[[4]])

# Clusters sampled, cluster FE, no heterogeneous TE, assignment no within-cluster correlation

res_10 <- sim_func(
  cluster_sample = TRUE,
  cluster_fe = TRUE,
  hetero_te = FALSE,
  cluster_assignment = FALSE
  )

print(res_10[[1]])
print(res_10[[2]])
print(res_10[[3]])
print(res_10[[4]])

# Clusters sampled, no cluster FE, heterogeneous TE, assignment no within-cluster correlation

pmt <- Sys.time()
res_11 <- sim_func(
  cluster_sample = TRUE,
  cluster_fe = FALSE,
  hetero_te = TRUE,
  cluster_assignment = FALSE
  )

print(res_11[[1]])
print(res_11[[2]])
print(res_11[[3]])
print(res_11[[4]])

Sys.time() - pmt

# Clusters sampled, no cluster FE, no heterogeneous TE, assignment within-cluster correlation

res_12 <- sim_func(
  cluster_sample = TRUE,
  cluster_fe = FALSE,
  hetero_te = FALSE,
  cluster_assignment = TRUE
  )

print(res_12[[1]])
print(res_12[[2]])
print(res_12[[3]])
print(res_12[[4]])

# Clusters sampled, cluster FE, heterogeneous TE, assignment no within-cluster correlation

res_13 <- sim_func(
  cluster_sample = TRUE,
  cluster_fe = TRUE,
  hetero_te = TRUE,
  cluster_assignment = FALSE
  )

print(res_13[[1]])
print(res_13[[2]])
print(res_13[[3]])
print(res_13[[4]])

# Clusters sampled, cluster FE, no heterogeneous TE, assignment within-cluster correlation

res_14 <- sim_func(
  cluster_sample = TRUE,
  cluster_fe = TRUE,
  hetero_te = FALSE,
  cluster_assignment = TRUE
  )

print(res_14[[1]])
print(res_14[[2]])
print(res_14[[3]])
print(res_14[[4]])

# Clusters sampled, no cluster FE, heterogeneous TE, assignment within-cluster correlation

res_15 <- sim_func(
  cluster_sample = TRUE,
  cluster_fe = FALSE,
  hetero_te = TRUE,
  cluster_assignment = TRUE
  )

print(res_15[[1]])
print(res_15[[2]])
print(res_15[[3]])
print(res_15[[4]])

# Clusters sampled, cluster FE, heterogeneous TE, assignment within-cluster correlation

res_16 <- sim_func(
  cluster_sample = TRUE,
  cluster_fe = TRUE,
  hetero_te = TRUE,
  cluster_assignment = TRUE
  )

print(res_16[[1]])
print(res_16[[2]])
print(res_16[[3]])
print(res_16[[4]])

out_mat <- rbind(
  unlist(res_1),
  unlist(res_2),
  unlist(res_3),
  unlist(res_4),
  unlist(res_5),
  unlist(res_6),
  unlist(res_7),
  unlist(res_8),
  unlist(res_9),
  unlist(res_10),
  unlist(res_11),
  unlist(res_12),
  unlist(res_13),
  unlist(res_14),
  unlist(res_15),
  unlist(res_16)
)

colnames(out_mat) <- c(
  "Cluster-robust SE w/o FE",
  "Hetero-robust SE w/o FE",
  "Cluster-robust SE w/ FE",
  "Hetero-robust SE w/ FE"
)

out_mat <- as_tibble(out_mat) %>% 
  cbind(DGP = c(
    "Clusters not sampled, no cluster FE, not heterogeneous TE, assignment no within-cluster correlation",
    "Clusters not sampled, cluster FE, not heterogeneous TE, assignment no within-cluster correlation",
    "Clusters not sampled, no cluster FE, heterogeneous TE, assignment no within-cluster correlation",
    "Clusters not sampled, no cluster FE, no heterogeneous TE, assignment within-cluster correlation",
    "Clusters not sampled, cluster FE, heterogeneous TE, assignment no within-cluster correlation",
    "Clusters not sampled, cluster FE, not heterogeneous TE, assignment within-cluster correlation",
    "Clusters not sampled, no cluster FE, heterogeneous TE, assignment within-cluster correlation",
    "Clusters not sampled, cluster FE, heterogeneous TE, assignment within-cluster correlation",
    "Clusters sampled, no cluster FE, no heterogeneous TE, assignment no within-cluster correlation",
    "Clusters sampled, cluster FE, no heterogeneous TE, assignment no within-cluster correlation",
    "Clusters sampled, no cluster FE, heterogeneous TE, assignment no within-cluster correlation",
    "Clusters sampled, no cluster FE, no heterogeneous TE, assignment within-cluster correlation",
    "Clusters sampled, cluster FE, heterogeneous TE, assignment no within-cluster correlation",
    "Clusters sampled, cluster FE, no heterogeneous TE, assignment within-cluster correlation",
    "Clusters sampled, no cluster FE, heterogeneous TE, assignment within-cluster correlation",
    "Clusters sampled, cluster FE, heterogeneous TE, assignment within-cluster correlation"
    )) %>% 
  relocate(DGP)

write_csv(out_mat, "cluster_se_shiny/out_mat.csv")
