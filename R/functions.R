#-------------------------------------------------------------------------------
#                                 Functions
#-------------------------------------------------------------------------------

Read_Gen_Data <- function(fileName, dimVals, groupName1, groupName2) {
  
  # This function reads data and prepares the data list for stan.
  # The stimulus dimension column must be named "x", responses as "y", group as
  # "group", subject ID as "subj"
  
  # Note that the "x" column must match the "dimVals" argument in this function
  # in terms of the range, and the position of the CS+ should be 0.
  
  # Changing the dimVals parameter may also require changing the parameters of 
  # the prior distributions in the model string.
  
  data <- read.csv(fileName, header = TRUE)
  
  # groupNames <- unique(data$group)
  groupNames <- c(groupName1, groupName2)
  
  out <- vector("list", 2)
  
  for (i in 1:length(groupNames)) {
    subset_data <- data %>% 
      dplyr::filter(group == groupNames[i]) %>%
      arrange(subj, x)
    out[[i]] <- list(subj = subset_data[["subj"]],
                     responses = matrix(subset_data[["y"]], ncol = length(dimVals), 
                                        byrow = TRUE),
                     nSubj = n_distinct(subset_data[["subj"]]),
                     nStim = length(dimVals),
                     xs = dimVals)
  }
  
  return(list(out, groupNames))
}
#_______________________________________________________________________________

Run_Model <- function(dataList, modelName, modelFile, params) {
  
  stanfit <- stan(file = modelFile,
                  data = dataList, 
                  pars = params,
                  iter = n_iter, 
                  warmup = n_burnin, 
                  thin = n_thin, 
                  chains = n_chains, 
                  init = "random",
                  algorithm = "NUTS",
                  cores = parallel::detectCores())
  
  diag <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  samples <- rstan::extract(stanfit)
  
  summary <- rstan::summary(stanfit, probs = c(0.025, 0.50, 0.975))$summary
  write.csv(summary, file = paste0(file_name_root, modelName, "-summary.csv"), row.names = TRUE)
  
  # calculate waic
  loglik <- loo::extract_log_lik(stanfit)
  waic <- loo::waic(loglik)
  
  # output
  out <- list(stanfit, diag, samples, summary, waic)
  names(out) <- c("stanfit", "diag", "samples", "summary", "waic")
  return(out)
}
#_______________________________________________________________________________

Posterior_Preds <- function(samples, responses, modelName, nSubj, subjList, 
                            nStim, summary, nRow = nRow, figMult) {
  
  # This function plots the posterior predictives overlayed on the empirical 
  # gradients facetting by subject
  
  post_preds <- matrix(NA, nrow = nSubj*nStim*n_samp, ncol = 4)
  post_preds <- as.data.frame(post_preds)
  colnames(post_preds) <- c("subj", "dim", "samp", "pred")
  post_preds$subj <- rep(1:nSubj, each = n_samp*nStim)
  post_preds$dim <- rep(rep(1:nStim, each = n_samp), times = nSubj)
  post_preds$samp <- rep(1:n_samp, times = nSubj*nStim)
  for (subj in 1:nSubj) {
    for (dim in 1:nStim) {
      temp <- sample(samples$predR[, subj, dim], size = n_samp) 
      for (samp in 1:n_samp) {
        start <- (subj-1) * nStim * n_samp+(dim-1) * n_samp + 1
        end <- (subj-1) * nStim * n_samp + (dim-1) * n_samp + n_samp
        post_preds$pred[start:end] <- temp
      }
    }
  }
  
  # add responses
  post_preds$response <- NA
  post_preds$response <- rep(responses, each = n_samp)
  
  # add mean posterior estimates
  label <- rep(NA, nSubj)
  for (i in 1:nSubj) {
    label[i] <- paste0("M: ", round(summary$mean[i], 2),
                       " SD +: ", round(summary$mean[nSubj + i], 2),
                       " SD -: ", round(summary$mean[nSubj*2 + i]),
                       " H: ", round(summary$mean[nSubj*3 + i]))
  }
  post_preds$label <- rep(label, each = n_samp * nStim)
  
  # plot empirical gradients with posterior samples overlayed
  grad_layers <- list(
    geom_line(stat = "identity", size = 1.25),
    labs(title = "", x = "dimension", y = "responding"),
    scale_colour_manual(values = fig_cols),
    geom_vline(xintercept = 6, linetype = "dotted", colour = "black"), 
    scale_x_continuous(limits = c(1, nStim), breaks = c(1, 6, nStim),
                       labels = c(min(dim_vals), 0, max(dim_vals))),
    scale_y_continuous(limits = c(0, 120), breaks = c(0, 50, 100)),
    theme_classic(),
    theme(panel.background = element_rect(colour = "black", size = 0.5, 
                                          linetype = "solid", fill = NA),
          axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
          axis.line.y = element_line(colour="black", size=0.5, linetype="solid"),
          legend.position="none")
  )
  
  fig <- ggplot(post_preds, aes(y = response, x = dim)) + 
    grad_layers +
    geom_text(data = post_preds, mapping = aes(label = label, x = 6, y = 115),  
              size = 2, colour = "black")
  ggsave(file = paste0(file_name_root, modelName, "-gradients", graph_file_type), 
         plot = fig + facet_wrap(~ subj, nrow = nRow),
         width = gg_width*figMult, height = gg_height*figMult, 
         # width = 25, height = 20,
         units = "cm", dpi = dpi)
  
  fig <- ggplot(post_preds, aes(y = response, x = dim)) + 
    grad_layers +
    geom_point(aes(y = pred, x = dim), colour = scat_col, size = scat_size, 
               shape = scat_shape) +
    geom_text(data = post_preds, mapping = aes(label = label, x = 6, y = 115),  
              size = 2, colour = "black")
  ggsave(file = paste0(file_name_root, modelName, "-postpreds", graph_file_type), 
         plot = fig + facet_wrap(~ subj, nrow = nRow),
         width = gg_width*figMult, height = gg_height*figMult, 
         # width = 25, height = 20,
         units = "cm", dpi = dpi)
}
#_______________________________________________________________________________

Plot_Densities <- function(samples1, samples2, groupName1, groupName2, graphName,
                           paramNames, dimVals) {
  
  # This function plots the densities of each augmented Gaussian parameter for 2 groups
  
  estimates <- as.data.frame(cbind(c(as.vector(samples1$M_group), as.vector(samples2$M_group)),
                                   c(as.vector(samples1$SDMinus_group), as.vector(samples2$SDMinus_group)),
                                   c(as.vector(samples1$SDPlus_group), as.vector(samples2$SDPlus_group)),
                                   c(as.vector(samples1$height_group), as.vector(samples2$height_group))))
  colnames(estimates) <- paramNames
  estimates$group <- c(rep(groupName1, length(as.vector(samples1$M_group))),
                       rep(groupName2, length(as.vector(samples2$M_group))))
  estimates$group <- as.factor(estimates$group)
  estimates$group <- fct_relevel(estimates$group, groupName1, groupName2)
  
  density_layers <- list(geom_density(alpha = .25),
                         scale_fill_manual(values = density_cols),
                         theme_classic())
  
  M_fig <- ggplot(estimates, aes(M, fill = group)) + 
    density_layers +
    guides(fill = FALSE) +
    scale_x_continuous(limits = c(min(dimVals), max(dimVals)), breaks = dimVals) +
    # geom_vline(xintercept = mean(samples1$M_group), linetype = "solid", colour = density_cols[1],
    #            size = 1) +
    # geom_vline(xintercept = mean(samples2$M_group), linetype = "solid", colour = density_cols[2],
    #            size = 1) +
    ggtitle("a) Mean")
  
  SDMinus_fig <- ggplot(estimates, aes(SDMinus, fill = group)) + 
    density_layers +
    guides(fill = FALSE) +
    scale_x_continuous(limits = c(0, 1)) +
    # geom_vline(xintercept = mean(samples1$SDMinus_group), linetype = "solid", colour = density_cols[1],
    #            size = 1) +
    # geom_vline(xintercept = mean(samples2$SDMinus_group), linetype = "solid", colour = density_cols[2],
    #            size = 1) +
    ggtitle("b) Width -")
  
  SDPlus_fig <- ggplot(estimates, aes(SDPlus, fill = group)) + 
    density_layers +
    guides(fill = FALSE) +
    scale_x_continuous(limits = c(0, 1)) +
    # geom_vline(xintercept = mean(samples1$SDPlus_group), linetype = "solid", colour = density_cols[1],
    #            size = 1) +
    # geom_vline(xintercept = mean(samples2$SDPlus_group), linetype = "solid", colour = density_cols[2],
    #            size = 1) +
    ggtitle("c) Width +")
  
  height_fig <- ggplot(estimates, aes(height, fill = group)) + 
    density_layers +
    scale_x_continuous(limits = c(40, 100), breaks = seq(40, 100, 10)) +
    # geom_vline(xintercept = mean(samples1$height_group), linetype = "solid", colour = density_cols[1],
    #            size = 1) +
    # geom_vline(xintercept = mean(samples2$height_group), linetype = "solid", colour = density_cols[2],
    #            size = 1) +
    ggtitle("d) Height") +
    theme(legend.position="none")
  
  fig_panel <- grid.arrange(M_fig + theme(legend.title = element_blank()), 
                            SDMinus_fig + theme(legend.title = element_blank()), 
                            SDPlus_fig + theme(legend.title = element_blank()), 
                            height_fig + theme(legend.title = element_blank()), nrow = 1)
  ggsave(paste0(file_name_root, graphName, "density", graph_file_type), fig_panel, 
         "jpeg", height = gg_height*.8, width = gg_width*3, units = "cm", dpi = dpi)
  
}

#_______________________________________________________________________________

Run_Analysis <- function(fileName, dimVals, nRow = c(6,6), figMult, graphName,
                         paramNames, groupName1, groupName2, modelFile, params) {
  
  # Function that reads data, runs the analysis, and saves the output
  # Note that groupName1 and groupName2 must match those in the data files
  
  # 1. read data
  out <- Read_Gen_Data(fileName, dimVals, groupName1, groupName2)
  data_list_1 <- out[[1]][[1]]
  data_list_2 <- out[[1]][[2]]
  
  # 2. fit models for each group
  mcmc_out_1 <- Run_Model(data_list_1, modelName = groupName1, modelFile, params = params)
  samples_1 <- mcmc_out_1[["samples"]]
  mcmc_out_2 <- Run_Model(data_list_2, modelName = groupName2, modelFile, params = params)
  samples_2 <- mcmc_out_2[["samples"]]
  
  # 3. plot posterior predictives
  Posterior_Preds(samples = samples_1, responses = as.vector(t(data_list_1$responses)), 
                  modelName = groupName1, nSubj = data_list_1$nSubj, 
                  nStim = data_list_1$nStim, summary = as.data.frame(mcmc_out_1$summary), 
                  nRow = nRow[1], figMult = figMult)
  Posterior_Preds(samples = samples_2, responses = as.vector(t(data_list_2$responses)), 
                  modelName = groupName2, nSubj = data_list_2$nSubj, 
                  nStim = data_list_2$nStim, summary = as.data.frame(mcmc_out_2$summary), 
                  nRow = nRow[2], figMult = figMult)
  
  # 4. compare posteriors between groups
  Plot_Densities(samples_1, samples_2, groupName1 = groupName1, groupName2 = groupName2, 
                 graphName = graphName, paramNames = paramNames, dimVals = dimVals)
  
  # 5. calculate HDIs for posterior estimates
  hdi_m1 <- bayestestR::hdi(as.vector(samples_1$M_group), ci = hdi_limit) 
  hdi_m2 <- bayestestR::hdi(as.vector(samples_2$M_group), ci = hdi_limit) 
  hdi_sdplus1 <- bayestestR::hdi(as.vector(samples_1$SDPlus_group), ci = hdi_limit) 
  hdi_sdplus2 <- bayestestR::hdi(as.vector(samples_2$SDPlus_group), ci = hdi_limit)
  hdi_sdminus1 <- bayestestR::hdi(as.vector(samples_1$SDMinus_group), ci = hdi_limit) 
  hdi_sdminus2 <- bayestestR::hdi(as.vector(samples_2$SDMinus_group), ci = hdi_limit)
  hdi_h1 <- bayestestR::hdi(as.vector(samples_1$height_group), ci = hdi_limit) 
  hdi_h2 <- bayestestR::hdi(as.vector(samples_2$height_group), ci = hdi_limit)
  temp <- as.data.frame(t(c(hdi_m1$CI_low, hdi_m1$CI_high, hdi_m2$CI_low, hdi_m2$CI_high, 
                            hdi_sdplus1$CI_low, hdi_sdplus1$CI_high, hdi_sdplus2$CI_low, hdi_sdplus2$CI_high,
                            hdi_sdminus1$CI_low, hdi_sdminus1$CI_high, hdi_sdminus2$CI_low, hdi_sdminus2$CI_high,
                            hdi_h1$CI_low, hdi_h1$CI_high, hdi_h2$CI_low, hdi_h2$CI_high)))
  colnames(temp) <- paste0(rep(c("m", "sdplus", "sdminus", "h"), each = 4),
                           rep(c("1", "1", "2", "2"), times = 4), rep(c("_low", "_high"), times = 8))
  write_csv(temp, paste0(file_name_root, "HDIs.csv"))
  
  # 6. calculate HDIs and ROPE for group difference scores
  # M
  diff_m <- as.vector(samples_1$M_group) - c(as.vector(samples_2$M_group))
  hdi_m <- bayestestR::hdi(diff_m, ci = hdi_limit) 
  rope_m <- bayestestR::rope(diff_m) 
  # W+
  diff_wplus <- as.vector(samples_1$SDPlus_group) - c(as.vector(samples_2$SDPlus_group))
  hdi_wplus <- bayestestR::hdi(diff_wplus, ci = hdi_limit) 
  rope_wplus <- bayestestR::rope(diff_wplus) 
  # W-
  diff_wminus <- as.vector(samples_1$SDMinus_group) - c(as.vector(samples_2$SDMinus_group))
  hdi_wminus <- bayestestR::hdi(diff_wminus, ci = hdi_limit) 
  rope_wminus <- bayestestR::rope(diff_wminus)
  # H
  diff_h <- as.vector(samples_1$height_group) - c(as.vector(samples_2$height_group))
  hdi_h <- bayestestR::hdi(diff_h, ci = hdi_limit) 
  rope_h <- bayestestR::rope(diff_h)
  
  temp <- as.data.frame(t(c(hdi_m$CI_low, hdi_m$CI_high, 
                            hdi_wplus$CI_low, hdi_wplus$CI_high, 
                            hdi_wminus$CI_low, hdi_wminus$CI_high,
                            hdi_h$CI_low, hdi_h$CI_high)))
  colnames(temp) <- paste0(rep(c("m", "sdplus", "sdminus", "h"), each = 2),
                           rep(c("_low", "_high"), times = 4))
  write_csv(temp, paste0(file_name_root, "group_diff_HDIs.csv"))
  
  out <- list(data_list_1, data_list_2, mcmc_out_1, mcmc_out_2, 
              samples_1, samples_2, diff_m, hdi_m, rope_m, diff_wplus, 
              hdi_wplus, rope_wplus, diff_wminus, hdi_wminus, rope_wminus, 
              diff_h, hdi_h, rope_h)
  names(out) <- c("data_list_1", "data_list_2", "mcmc_out_1", "mcmc_out_2", 
                  "samples_1", "samples_2", "diff_m", "hdi_m", "rope_m", 
                  "diff_wplus", "hdi_wplus", "rope_wplus", "diff_wminus", 
                  "hdi_wminus", "rope_wminus", "diff_h", "hdi_h", "rope_h")
  return(out)
}
#_______________________________________________________________________________