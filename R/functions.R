#-------------------------------------------------------------------------------
#                                 Functions
#-------------------------------------------------------------------------------

Simulate_Data <- function(params, nSubj, dimVals = dim_vals, fName = "toy_data") {
  # This function simulates augmented Gaussian gradients for two groups of N=nSubj
  # Note
  # - params is a list of length 2 (number of groups), and each item consists of 
  # a list: M, W1, W2, H, and noise are vectors of length nSubj
  
  toyData <- list()
  
  for (g in 1:2) {
    
    simData <- matrix(nrow = nSubj, ncol = length(dimVals))
    M <- params[[g]]$M
    WM <- params[[g]]$WM
    WP <- params[[g]]$WP
    H <- params[[g]]$H
    noise <- params[[g]]$noise
    
    for (i in 1:nSubj) {
      gausLeft <- H[i] * exp(1)^-(((dimVals-M[i])^2) / (2 * WM[i]^2)) + 
        rnorm(length(dimVals), 0, noise)
      gausRight <- H[i] * exp(1)^-(((dimVals-M[i])^2) / (2 * WP[i]^2)) + 
        rnorm(length(dimVals), 0, noise)
      mIdx <- which(round(dimVals,1) == round(M[i], 1))
      if (mIdx == 1) {
        simData[i,] <- c(gausLeft[1], gausRight[(mIdx+1):length(dimVals)])
      } else {
        simData[i,] <- c(gausLeft[1:(mIdx-1)], gausRight[mIdx:length(dimVals)])
      }
    }
    subj <- rep(((g-1)*nSubj + 1):(nSubj*g), each = length(dimVals))
    x <- rep(dimVals, times = nSubj)
    y <- as.vector(t(simData))
    y[y > 100] <- 100
    y[y < 0] <- 0
    group <- rep(paste0("group", g), nSubj)
    toyData[[g]] <- data.frame(cbind(subj, group, x, y))
  }
  
  out <- rbind(toyData[[1]], toyData[[2]])
  out <- out %>%
    mutate(x = as.numeric(as.character(x)),
           y = as.numeric(as.character(y)),
           subj = as.numeric(as.character(subj)))
  
  layers <- list(
    geom_line(alpha = .5, colour = "black"),
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey"),
    scale_x_continuous(limits = c(min(dimVals), max(dimVals)),
                       breaks = c(min(dimVals), 0, max(dimVals)),
                       labels = c("", "CS+", "")),
    scale_y_continuous(limits = c(0, 105), breaks = c(0, 50, 100)),
    theme_classic(),
    facet_wrap(~ subj, nrow = round(sqrt(nSubj)))
  )
  
  fig1 <- ggplot(filter(out, group == "group1"), aes(x = x, y = y)) +
    layers + ggtitle("Group 1")
  fig2 <- ggplot(filter(out, group == "group2"), aes(x = x, y = y)) +
    layers + ggtitle("Group 2")
  ggsave(file = paste0("output/", fName, graph_file_type), 
         plot = fig1 + fig2, width = gg_width*3, height = gg_height*1.5, 
         units = "cm", dpi = dpi)
  write_csv(out, paste0("data/", fName ,".csv"))
  return(out)
}

#_______________________________________________________________________________
Read_Gen_Data <- function(fileName, dimVals, groupName1, groupName2) {
  
  # This function reads long data, prepares the data list for stan, and plots 
  # the gradients
  # Note
  # - Responses should be labelled as "y", subject ID as "subj", group names 
  # labelled as "group" (and should match groupName1 and groupName2)
  # - The "x" column does not have to match the dimVals argument but the dimVals 
  # argument should match the order of the "x" column, and the position of the 
  # CS+ should be 0 (see demo1.csv)
  
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
  
  layers <- list(
    geom_line(size = 1.5),
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey"),
    scale_x_continuous(limits = c(min(dimVals), max(dimVals)),
                       breaks = dimVals,
                       labels = c(min(dimVals), "", "", "", "", "CS+", 
                                  "", "", "", "", max(dimVals))),
    scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 10),
                       labels = c(0, "", "", "", "", 50, "", "", "", "", 100)),
    scale_colour_manual(values = density_cols),
    labs(x = "stimulus", y = "mean response"),
    theme_classic()
    # theme(legend.position = c(0.1, 0.8))
  )
  
  gradients <- data %>%
    mutate(group = fct_relevel(group, levels = c(groupName1, groupName2))) %>%
    group_by(group, x) %>%
    summarise(y = mean(y)) %>%
    mutate(x = dimVals)
  
  fig <- ggplot(gradients, aes(x = x, y = y, group = group, colour = group)) +
    layers
  
  return(list(out, groupNames, fig))
}
#_______________________________________________________________________________

Run_Model <- function(dataList, modelName, modelFile, params) {
  
  # This function runs the model in stan
  
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
  
  # save diagnostics
  diag <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  
  # save samples
  samples <- rstan::extract(stanfit)
  
  # save summary file
  summary <- rstan::summary(stanfit, probs = c(0.025, 0.50, 0.975))$summary
  write.csv(summary, file = paste0(file_name_root, modelName, "-summary.csv"), row.names = TRUE)
  
  # calculate waic
  loglik <- extract_log_lik(stanfit)
  waic <- waic(loglik)
  loo <- loo(loglik)
  
  # output
  out <- list(stanfit, diag, samples, summary, waic, loo)
  names(out) <- c("stanfit", "diag", "samples", "summary", "waic", "loo")
  return(out)
}
#_______________________________________________________________________________

Posterior_Preds <- function(samples, responses, modelName, nSubj, subjList, 
                            nStim, summary, nRow = nRow, figMult, labels) {
  
  # This function plots the posterior predictives for each subject overlayed on 
  # the empirical gradients
  # - Optional labels: display the mean of the posterior for each parameter
  
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
                       " SD -: ", round(summary$mean[nSubj*2 + i], 2),
                       " H: ", round(summary$mean[nSubj*3 + i]))
  }
  post_preds$label <- rep(label, each = n_samp * nStim)
  
  # figure layers
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
  
  # plot gradients alone
  if (labels == TRUE) {
    fig <- ggplot(post_preds, aes(y = response, x = dim)) + 
      grad_layers +
      geom_text(data = post_preds, mapping = aes(label = label, x = 6, y = 115),  
                     size = 2.5, colour = "black")
  } else if (labels == FALSE) {
    fig <- ggplot(post_preds, aes(y = response, x = dim)) + 
      grad_layers
  }
  ggsave(file = paste0(file_name_root, modelName, "-gradients", graph_file_type), 
         plot = fig + facet_wrap(~ subj, nrow = nRow),
         width = gg_width*figMult, height = gg_height*figMult, 
         units = "cm", dpi = dpi)
  
  # plot gradients with posterior predictives
  if (labels == TRUE) {
    fig <- ggplot(post_preds, aes(y = response, x = dim)) + 
      grad_layers +
      geom_point(aes(y = pred, x = dim), colour = scat_col, size = scat_size, 
                 shape = scat_shape) +
      geom_text(data = post_preds, mapping = aes(label = label, x = 6, y = 115),  
                size = 2.5, colour = "black")
  } else if (labels == FALSE) {
    fig <- ggplot(post_preds, aes(y = response, x = dim)) + 
      grad_layers +
      geom_point(aes(y = pred, x = dim), colour = scat_col, size = scat_size, 
                 shape = scat_shape)
  }
  ggsave(file = paste0(file_name_root, modelName, "-postpreds", graph_file_type), 
         plot = fig + facet_wrap(~ subj, nrow = nRow),
         width = gg_width*figMult, height = gg_height*figMult, 
         units = "cm", dpi = dpi)
}
#_______________________________________________________________________________

Plot_Densities <- function(samples1, samples2, groupName1, groupName2, graphName,
                           paramNames, dimVals) {
  
  # This function produces a multi-panelled figure with:
  # a) The mean empirical gradients for each group
  # b) Posterior distribution for the Mean parameter for each group
  # c) Posterior distribution for the Height parameter for each group
  # d) Posterior distribution for the Width- parameter for each group
  # e) Posterior distribution for the Width+ parameter for each group
  
  estimates <- as.data.frame(
    cbind(
      c(as.vector(samples1$M_group), as.vector(samples2$M_group)),
      c(as.vector(samples1$SDPlus_group), as.vector(samples2$SDPlus_group)),
      c(as.vector(samples1$SDMinus_group), as.vector(samples2$SDMinus_group)),
      c(as.vector(samples1$height_group), as.vector(samples2$height_group))))
  colnames(estimates) <- paramNames
  estimates$group <- c(rep(groupName1, length(as.vector(samples1$M_group))),
                       rep(groupName2, length(as.vector(samples2$M_group))))
  estimates$group <- as.factor(estimates$group)
  estimates$group <- fct_relevel(estimates$group, groupName1, groupName2)
  
  density_layers <- list(geom_density(alpha = .25),
                         scale_fill_manual(values = density_cols),
                         theme_classic())
  
  M_fig <- ggplot(estimates, aes(M_group, fill = group)) + 
    density_layers +
    scale_x_continuous(limits = c(min(dimVals), max(dimVals)), breaks = dimVals, 
                       labels = c(min(dimVals), "", "", "", "", 0, "", "", "", "", max(dimVals))) +
    theme(axis.title.x = element_blank()) +
    ggtitle("b) Mean")
  
  height_fig <- ggplot(estimates, aes(height_group, fill = group)) + 
    density_layers +
    scale_x_continuous(limits = c(40, 100), breaks = seq(40, 100, 10)) +
    theme(axis.title.x = element_blank()) +
    ggtitle("c) Height")
  
  SDMinus_fig <- ggplot(estimates, aes(SDMinus_group, fill = group)) + 
    density_layers +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    theme(axis.title.x = element_blank()) +
    ggtitle("d) Width -")
  
  SDPlus_fig <- ggplot(estimates, aes(SDPlus_group, fill = group)) + 
    density_layers +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    theme(axis.title.x = element_blank()) +
    ggtitle("e) Width +")
  
  fig_panel <- M_fig + height_fig + SDMinus_fig + SDPlus_fig +
    plot_layout(nrow = 2, byrow = TRUE)
}

#_______________________________________________________________________________
Fit_Aug_Gaussian <- function(fileName, modelFile, groupName1, groupName2, graphName, 
                             dimVals, params, groupParams, 
                             ropeLow, ropeHigh, ropeLowDiffs, ropeHighDiffs,
                             hdiLim, nRow = c(5,5), figMult = 2, labels = TRUE) {
  
  # Master function that reads the data file, runs the analysis, and generates summary output
  
  # Note
  # - groupName1 and groupName2 must match those in the data files
  # - all ROPE parameters must be of length(groupParams) and specified in the 
  # same order as groupParams
  # - ropeLow and ropeHigh refer to custom ROPE limits for each parameter
  # - ropeLowDiffs and ropeHighDiffs refer to custom ROPE limits for the group
  # difference in each parameter
  # - also calculates standard ROPE limits (+- 0.1xSD) for each group difference
  
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
                  nRow = nRow[1], figMult = figMult, labels = labels)
  Posterior_Preds(samples = samples_2, responses = as.vector(t(data_list_2$responses)), 
                  modelName = groupName2, nSubj = data_list_2$nSubj, 
                  nStim = data_list_2$nStim, summary = as.data.frame(mcmc_out_2$summary), 
                  nRow = nRow[2], figMult = figMult, labels = labels)
  
  # 4. plot gradients + posteriors
  density_fig <- Plot_Densities(samples_1, samples_2, groupName1 = groupName1, 
                                groupName2 = groupName2, graphName = graphName, 
                                paramNames = groupParams, dimVals = dimVals)
  
  fig_panel <- out[[3]] + ggtitle("a) Generalisation gradients ") +
    density_fig +
    plot_layout(widths = c(1, 1.75), guides = "collect")
  ggsave(paste0(file_name_root, graphName, "density", graph_file_type), fig_panel, 
         "jpeg", height = gg_height, width = gg_width*2.5, units = "cm", dpi = dpi)
  
  # 5. HDIS and posterior summary stats for the group parameters 
  temp <- data.frame(group = rep(c(groupName1, groupName2), each = length(groupParams)),
                     param = rep(groupParams, times = 2),
                     hdi_lim = hdiLim,
                     hdi_low = rep(NA, length(groupParams)*2),
                     hdi_high = rep(NA, length(groupParams)*2),
                     p_dir = rep(NA, length(groupParams)*2),
                     rope_low = rep(ropeLow, 2),
                     rope_high = rep(ropeHigh, 2),
                     prop_rope = rep(NA, length(groupParams)*2))
  for (i in 1:length(groupParams)) {
    temp$hdi_low[i] <- hdi(as.vector(samples_1[[groupParams[i]]]), ci = hdiLim)$CI_low
    temp$hdi_low[length(groupParams) + i] <- hdi(as.vector(samples_2[[groupParams[i]]]), ci = hdiLim)$CI_low
    temp$hdi_high[i] <- hdi(as.vector(samples_1[[groupParams[i]]]), ci = hdiLim)$CI_high
    temp$hdi_high[length(groupParams) + i] <- hdi(as.vector(samples_2[[groupParams[i]]]), ci = hdiLim)$CI_high
    temp$p_dir[i] <- p_direction(as.vector(samples_1[[groupParams[i]]]), ci = hdiLim)
    temp$p_dir[length(groupParams) + i] <- p_direction(as.vector(samples_2[[groupParams[i]]]), ci = hdiLim)
    # custom ROPE
    temp$prop_rope[i] <- rope(as.vector(samples_1[[groupParams[i]]]), ci = hdiLim,
                              range = c(ropeLow[i], ropeHigh[i]))$ROPE_Percentage
    temp$prop_rope[length(groupParams) + i] <- rope(as.vector(samples_2[[groupParams[i]]]), ci = hdiLim,
                                                    range = c(ropeLow[i], ropeHigh[i]))$ROPE_Percentage
  }
  write_csv(format(temp, scientific = FALSE), paste0(file_name_root, "hdis.csv"))
  
  # 6. HDIS and posterior summary stats for the group differences 
  diff_m <- as.vector(samples_1$M_group) - c(as.vector(samples_2$M_group))
  diff_wplus <- as.vector(samples_1$SDPlus_group) - c(as.vector(samples_2$SDPlus_group))
  diff_wminus <- as.vector(samples_1$SDMinus_group) - c(as.vector(samples_2$SDMinus_group))
  diff_h <- as.vector(samples_1$height_group) - c(as.vector(samples_2$height_group))
  samples_diffs <- as.data.frame(cbind(diff_m, diff_wplus, diff_wminus, diff_h))
  colnames(samples_diffs) <- groupParams
  temp <- data.frame(param = groupParams,
                     hdi_lim = hdiLim,
                     hdi_low = rep(NA, length(groupParams)),
                     hdi_high = rep(NA, length(groupParams)),
                     p_dir = rep(NA, length(groupParams)),
                     rope_low = ropeLowDiffs, # custom ROPE limits
                     rope_high = ropeHighDiffs, 
                     prop_rope = rep(NA, length(groupParams)),
                     rope_low_stand = rep(NA, length(groupParams)), # standardised ROPE limits
                     rope_high_stand = rep(NA, length(groupParams)),
                     prop_rope_stand = rep(NA, length(groupParams)))
  for (i in 1:length(groupParams)) {
    temp$hdi_low[i] <- hdi(as.vector(samples_diffs[[groupParams[i]]]), ci = hdiLim)$CI_low
    temp$hdi_high[i] <- hdi(as.vector(samples_diffs[[groupParams[i]]]), ci = hdiLim)$CI_high
    temp$p_dir[i] <- p_direction(as.vector(samples_diffs[[groupParams[i]]]), ci = hdiLim)
    # custom ROPE
    temp$prop_rope[i] <- rope(as.vector(samples_diffs[[groupParams[i]]]), ci = hdiLim,
                                    range = c(temp$rope_low[i], temp$rope_high[i]))$ROPE_Percentage
    # standardised ROPE
    temp$rope_low_stand[i] <- -0.1 * sd(samples_diffs[[groupParams[i]]])
    temp$rope_high_stand[i] <- +0.1 * sd(samples_diffs[[groupParams[i]]])
    temp$prop_rope_stand[i] <- rope(as.vector(samples_diffs[[groupParams[i]]]), ci = hdiLim,
                              range = c(temp$rope_low_stand[i], temp$rope_high_stand[i]))$ROPE_Percentage
  }
  write_csv(format(temp, scientific = FALSE), paste0(file_name_root, "group_diff_hdis.csv"))
  
  # 7. predictive accuracy measures
  loo_1 <- mcmc_out_1[["loo"]]
  loo_2 <- mcmc_out_2[["loo"]]
  waic_1 <- mcmc_out_1[["waic"]]
  waic_2 <- mcmc_out_2[["waic"]]
  temp <- data.frame(group = c(groupName1, groupName2),
                     elpd_waic = c(waic_1$estimates[1,1], waic_2$estimates[1,1]),
                     p_waic = c(waic_1$estimates[2,1], waic_2$estimates[2,1]),
                     waic = c(waic_1$estimates[3,1], waic_2$estimates[3,1]))
  write_csv(temp, paste0(file_name_root, "waics.csv"))
  
  # return output list
  out <- list(data_list_1, data_list_2, mcmc_out_1, mcmc_out_2, waic_1, waic_2,
              loo_1, loo_2)
  names(out) <- c("data_list_1", "data_list_2", "mcmc_out_1", "mcmc_out_2", 
                  "waic_1", "waic_2", "loo_1", "loo_2")
  return(out)
}

#_______________________________________________________________________________
Plot_Param_Recovery <- function(params, samples1, samples2, nSubj, fName) {
  
  # This function plots simulated parameters (params) against recovered 
  # parameters (samples1 and samples2)
  # - params is a list of length 2 (number of groups), and each item consists of 
  # a list: M, W1, W2, H, and noise are vectors of length nSubj
  
  H1 <- WM1 <- WP1 <- M1 <- H2 <- WM2 <- WP2 <- M2 <- rep(NA, nSubj)
  
  for (i in 1:nSubj) {
    M1[i] <- mean(samples1$M[,i])
    H1[i] <- mean(samples1$height[,i])
    WM1[i] <- mean(samples1$SDMinus[,i])
    WP1[i] <- mean(samples1$SDPlus[,i])
    M2[i] <- mean(samples2$M[,i])
    H2[i] <- mean(samples2$height[,i])
    WM2[i] <- mean(samples2$SDMinus[,i])
    WP2[i] <- mean(samples2$SDPlus[,i])
  }
  
  temp1 <- data.frame(param = rep(c("M", "H", "WM", "WP"), each = nSubj),
                      simulated = c(params[[1]]$M, params[[1]]$H,
                                    params[[1]]$WM, params[[1]]$WP),
                      recovered = c(M1, H1, WM1, WP1))
  temp2 <- data.frame(param = rep(c("M", "H", "WM", "WP"), each = nSubj),
                      simulated = c(params[[2]]$M, params[[2]]$H,
                                    params[[2]]$WM, params[[2]]$WP),
                      recovered = c(M2, H2, WM2, WP2))
  
  layers <- list(
    geom_point(),
    geom_abline(),
    theme_classic(),
    facet_wrap(~param, scales = "free")
  )
  
  fig1 <- ggplot(temp1, aes(x = simulated, y = recovered)) + layers +
    ggtitle("Group 1")
  fig2 <- ggplot(temp2, aes(x = simulated, y = recovered)) + layers +
    ggtitle("Group 2")
  
  fig_panel <- fig1 + fig2
  ggsave(file = paste0("output/",fName, "-param_recovery", graph_file_type), 
         plot = fig_panel, width = gg_width*2.5, height = gg_height*1.5, 
         units = "cm", dpi = dpi)
  write_csv(rbind(temp1, temp2), paste0("output/", fName, "-recovered_params.csv"))
}
#_______________________________________________________________________________
Compare_Gaussians <- function(mod1_loo, mod2_loo, mod1_waic, mod2_waic) {
  
  # This function compares standard vs. augmented Gaussian models using waic 
  # and loo
  
  aug_v_norm <- list(
    loo_compare(mod1_loo, mod2_loo),
    loo_compare(mod1_waic, mod2_waic)
  )
  names(aug_v_norm) <- c("loo_comp", "waic_comp")
  return(aug_v_norm)
}
#_______________________________________________________________________________