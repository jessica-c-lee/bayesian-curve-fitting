#-------------------------------------------------------------------------------
#                                 Functions
#-------------------------------------------------------------------------------

Read_Gen_Data <- function(fileName, dimVals, groupName1, groupName2) {
  
  # This function reads data and prepares the data list to be inputted to stan.
  # The stimulus dimension column must be named "x", responses as "y", group as
  # "group", subject ID as "subj"
  
  # Note that the "x" column must match the "dimVals" argument in this function. 
  # For example, if you have 11 stimuli (S1-S11) and the CS+ is the middle 
  # stimulus (S6), then the position of the CS+ should be 0 in the xs argument, 
  # and range between -.5 to +.5.
  
  # Note that the range of the xs argument can be changed if needed, but then 
  # the parameters of the prior distributions must also be changed accordingly
  # in stan.
  
  data <- read.csv(fileName, header = TRUE)
  
  # groupNames <- unique(data$group)
  groupNames <- c(groupName1, groupName2)
  
  out <- vector("list", 2)
  
  for (i in 1:length(groupNames)) {
    subset_data <- data %>% 
      filter(group == groupNames[i]) %>%
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

Run_Aug_Gaussian_Mod <- function(dataList, modelName, modelFile) {
  
  stanfit <- stan(file = modelFile,
                  data = dataList, 
                  pars = c("M", "SDPlus", "SDMinus", "height", "noise", "predR", "log_lik"),
                  iter = n_iter, 
                  warmup = n_burnin, 
                  thin = n_thin, 
                  chains = n_chains, 
                  init = "random",
                  algorithm = "NUTS",
                  cores = 1)
  
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
  
  # This function samples from the posterior for each individual and plots the  
  # samples overlayed on the empirical gradients
  
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
  # return(fig)
}

#_______________________________________________________________________________

Compare_Groups <- function(samples1, samples2, groupName1, groupName2, graphName,
                           paramNames, dimVals) {
  
  # This function plots the densities of each Gaussian parameter for 2 groups
  
  estimates <- as.data.frame(cbind(c(as.vector(samples1$M), as.vector(samples2$M)),
                                   c(as.vector(samples1$SDMinus), as.vector(samples2$SDMinus)),
                                   c(as.vector(samples1$SDPlus), as.vector(samples2$SDPlus)),
                                   c(as.vector(samples1$height), as.vector(samples2$height))))
  colnames(estimates) <- paramNames
  estimates$group <- c(rep(groupName1, length(as.vector(samples1$M))),
                       rep(groupName2, length(as.vector(samples2$M))))
  estimates$group <- as.factor(estimates$group)
  estimates$group <- fct_relevel(estimates$group, groupName1, groupName2)
  
  density_layers <- list(geom_density(alpha = .25),
                         scale_fill_manual(values = density_cols),
                         theme_classic())
  
  M_fig <- ggplot(estimates, aes(M, fill = group)) + 
    density_layers +
    guides(fill = FALSE) +
    scale_x_continuous(limits = c(min(dimVals), max(dimVals)), breaks = dimVals) +
    geom_vline(xintercept = mean(samples1$M), linetype = "solid", colour = density_cols[1],
               size = 1.5) +
    geom_vline(xintercept = mean(samples2$M), linetype = "solid", colour = density_cols[2],
               size = 1.5) +
    ggtitle("a) Mean")
  
  SDMinus_fig <- ggplot(estimates, aes(SDMinus, fill = group)) + 
    density_layers +
    guides(fill = FALSE) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_vline(xintercept = mean(samples1$SDMinus), linetype = "solid", colour = density_cols[1],
               size = 1.5) +
    geom_vline(xintercept = mean(samples2$SDMinus), linetype = "solid", colour = density_cols[2],
               size = 1.5) +
    ggtitle("b) Width -")
  
  SDPlus_fig <- ggplot(estimates, aes(SDPlus, fill = group)) + 
    density_layers +
    guides(fill = FALSE) +
    scale_x_continuous(limits = c(0, 1)) +
    geom_vline(xintercept = mean(samples1$SDPlus), linetype = "solid", colour = density_cols[1],
               size = 2) +
    geom_vline(xintercept = mean(samples2$SDPlus), linetype = "solid", colour = density_cols[2],
               size = 2) +
    ggtitle("c) Width +")
  
  height_fig <- ggplot(estimates, aes(height, fill = group)) + 
    density_layers +
    scale_x_continuous(limits = c(40, 100), breaks = seq(40, 100, 10)) +
    geom_vline(xintercept = mean(samples1$height), linetype = "solid", colour = density_cols[1],
               size = 1.5) +
    geom_vline(xintercept = mean(samples2$height), linetype = "solid", colour = density_cols[2],
               size = 1.5) +
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
                         paramNames, groupName1, groupName2, modelFile) {
  
  # Handler that calls functions to run analysis and save output
  # Note that groupName1 and groupName2 must match those in the data files
  
  # 1. read data
  out <- Read_Gen_Data(fileName, dimVals, groupName1, groupName2)
  data_list_1 <- out[[1]][[1]]
  data_list_2 <- out[[1]][[2]]
  # groupNames <- out[2][[1]]
  
  # 2. fit models for each group
  mcmc_out_1 <- Run_Aug_Gaussian_Mod(data_list_1, modelName = groupName1, modelFile)
  samples_1 <- mcmc_out_1[["samples"]]
  mcmc_out_2 <- Run_Aug_Gaussian_Mod(data_list_2, modelName = groupName2, modelFile)
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
  Compare_Groups(samples_1, samples_2, groupName1 = groupName1, 
                 groupName2 = groupName2, graphName = graphName,  
                 paramNames = paramNames, dimVals = dimVals)
  
  # 5. calculate HDI and ROPE for difference scores
  # M
  diff_m <- as.vector(samples_1$M) - c(as.vector(samples_2$M))
  hdi_m <- bayestestR::hdi(diff_m, ci = hdi_limit) 
  rope_m <- bayestestR::rope(diff_m) 
  # W+
  diff_wplus <- as.vector(samples_1$SDPlus) - c(as.vector(samples_2$SDPlus))
  hdi_wplus <- bayestestR::hdi(diff_wplus, ci = hdi_limit) 
  rope_wplus <- bayestestR::rope(diff_wplus) 
  # W-
  diff_wminus <- as.vector(samples_1$SDMinus) - c(as.vector(samples_2$SDMinus))
  hdi_wminus <- bayestestR::hdi(diff_wminus, ci = hdi_limit) 
  rope_wminus <- bayestestR::rope(diff_wminus)
  # H
  diff_h <- as.vector(samples_1$height) - c(as.vector(samples_2$height))
  hdi_h <- bayestestR::hdi(diff_h, ci = hdi_limit) 
  rope_h <- bayestestR::rope(diff_h)
  
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