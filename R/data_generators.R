



simulate_clustered_data <- function(treatment_arm_cluster_size = 15, treatment_arm_clusters = 50,
                                    control_arm_cluster_size = 5, control_arm_clusters = 50,
                                    txt_eff = .5, ranef_sd = 1, covar1_coef = 1, seed = NA,
                                    informative_cluster_size = F){
  if(!is.na(seed)){ set.seed(seed) }
  total_units_txt <- treatment_arm_cluster_size*treatment_arm_clusters
  total_units_ctr <- control_arm_cluster_size*control_arm_clusters
  total_units <- total_units_txt + total_units_ctr
  unit_id <- 1:(total_units)
  cluster_id_txt <- rep(1:treatment_arm_clusters, each = treatment_arm_cluster_size)
  cluster_id_ctr <- rep((max(cluster_id_txt)+1):(max(cluster_id_txt) + control_arm_clusters), each = control_arm_cluster_size)
  cluster_id <- c(cluster_id_txt, cluster_id_ctr)
  total_clusters <- max(cluster_id)
  cluster_lengths <- c(rep(treatment_arm_cluster_size, times = treatment_arm_clusters),
                       rep(control_arm_cluster_size, times = control_arm_clusters))
  cluster_cov <- sd(cluster_lengths) / mean(cluster_lengths)
  N_j <- rep(cluster_lengths, times = cluster_lengths)
  
  A <- c(rep(1, times = total_units_txt), rep(0, times = total_units_ctr))
  
  if(ranef_sd == 0){
    UE <- 0
  } else {
    UEs <- exp(rnorm(n = total_clusters, mean = 0, sd = ranef_sd))
    UE <- rep(UEs, times = cluster_lengths)
  }
  UY <- runif(n = total_units, min = 0, max = 2)
  covar0 <- rnorm(n = total_units, mean = UE, sd = .5)
  covar1 <- rbinom(n = total_units, size = 1, prob = .2)
  covar2 <- rnorm(n = total_units, mean = 0, sd = .5)
  covar3 <- rnorm(n = total_units, mean = 0, sd = .5)
  
  Y0 <- covar1_coef*covar1 + covar0^2 + covar0*covar1 + UY + UE
  Y1 <- covar1_coef*covar1 + covar0^2 + covar0*covar1 + UY + UE + txt_eff
  if(informative_cluster_size){
    Y1 <- Y1 + txt_eff*N_j/total_units
  }
  Y <- ifelse(A == 1, Y1, Y0)
  dat <- cbind.data.frame(unit_id, cluster_id, A,
                          UY, covar0, covar1, covar2, covar3, UE,
                          txt_eff, Y0, Y1, Y,
                          N_j = N_j) %>% dplyr::mutate(id = cluster_id, alpha = 1, U = 1)
  cluster_means <- dat %>% dplyr::group_by(cluster_id, A) %>% dplyr::summarise(cluster_mean0 = mean(Y0),
                                                                               cluster_mean1 = mean(Y1),
                                                                               cluster_meanY = mean(Y))
  ind_weight <- total_clusters / total_units
  ind_weighted_cluster_means <- dat %>% dplyr::group_by(cluster_id, A) %>%
    dplyr::summarise(ind_weighted_cluster_mean0 = sum(Y0)*ind_weight,
                     ind_weighted_cluster_mean1 = sum(Y1)*ind_weight,
                     ind_weighted_cluster_meanY = sum(Y)*ind_weight)
  cluster_mean_diff <- mean(cluster_means$cluster_mean1) - mean(cluster_means$cluster_mean0)
  cluster_mean_irr <- mean(cluster_means$cluster_mean1) / mean(cluster_means$cluster_mean0)
  ind_weighted_cluster_mean_diff <- mean(ind_weighted_cluster_means$ind_weighted_cluster_mean1) - mean(ind_weighted_cluster_means$ind_weighted_cluster_mean0)
  ind_weighted_cluster_mean_irr <- mean(ind_weighted_cluster_means$ind_weighted_cluster_mean1) / mean(ind_weighted_cluster_means$ind_weighted_cluster_mean0)
  
  return(cbind.data.frame(dat,
                          covar1_coef = covar1_coef,
                          true_rd_input = txt_eff,
                          true_covar_coef_input = covar1_coef,
                          ranef_sd_input = ranef_sd,
                          cluster_cov = cluster_cov,
                          sample_mean0 = mean(Y0),
                          sample_mean1 = mean(Y1),
                          unit_mean_diff = mean(Y1) - mean(Y0),
                          unit_irr = mean(Y1) / mean(Y0),
                          cluster_mean_diff = cluster_mean_diff,
                          cluster_mean_irr = cluster_mean_irr,
                          ind_weighted_cluster_mean_diff = ind_weighted_cluster_mean_diff,
                          ind_weighted_cluster_mean_irr = ind_weighted_cluster_mean_irr
  ))
}
