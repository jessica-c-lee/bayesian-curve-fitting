#-------------------------------------------------------------------------------
#                                 Functions
#-------------------------------------------------------------------------------

Simulate_Data <- function(M, SD, simHeight, noise, nSubj, dimVals) {
  # This function simulates Gaussian gradients
  # M, SD, simHeight are all vectors of length nSubj
  # dimVals must range between -.5 and +.5
  
  simData <- matrix(nrow = nSubj, ncol = length(dimVals))
  for (i in 1:nSubj) {
    simData[i,] <- simHeight[i] * exp(1)^-(((dimVals-M[i])^2) / (2 * SD[i]^2)) + rnorm(length(dimVals), 0, noise)
  }
  subj <- rep(1:nSubj, each = length(dimVals))
  x <- rep(dimVals, times = nSubj)
  y <- as.vector(t(simData))
  y[y > 100] <- 100
  y[y < 0] <- 0
  toyData <- data.frame(cbind(subj, x, y))
  fig <- ggplot(toyData, aes(x = x, y = y)) +
    geom_line(alpha = .5) + 
    geom_vline(xintercept = 0, linetype = "dotted", colour = "black") +
    scale_x_continuous(limits = c(min(dimVals), max(dimVals)), 
                       breaks = c(min(dimVals), 0, max(dimVals))) +
    scale_y_continuous(limits = c(0, 120), breaks = c(0, 50, 100)) +
    theme_classic() +
    facet_wrap(~ subj, nrow = n_row) 
  out <- list(toyData, fig)
  names(out) <- c("data", "fig")
  return(out)
}

#_______________________________________________________________________________

Read_Demo <- function(fileName, dimVals) {
  
  # This function reads data in the format specified in the demo_data.csv file and 
  # prepares the data list to be inputted to stan. 
  
  # Note that the "x" column must match the "dimVals" argument in this function. 
  # For example, if you have 11 stimuli (S1-S11) and the CS+ is the middle 
  # stimulus (S6), then the position of the CS+ should be 0 in the xs argument, 
  # and range between -.5 to +.5.
  # If you have a 2-choice task and have 5 test stimuli (S0-S4) then dimVals 
  # should start from 0 and range between 0 and +.5 (e.g., 0, +.125, +.25, +.375, +.5). 
  # Note that the range of the xs argument can be changed if needed, but then 
  # the parameters of the prior distributions must also be changed accordingly
  # in stan.
  
  data <- read.csv(fileName, header = TRUE)
  data <- arrange(data, subj, x)
  subj <- data[["subj"]]
  x <- data[["x"]]
  y <- data[["y"]]
  nSubj <- length(unique(subj))
  nStim <- length(unique(x))
  responses <- matrix(data[["y"]], ncol = nStim, byrow = TRUE)
  data_list <- list(subj = subj,
                    responses = responses,
                    nSubj = nSubj,
                    nStim = nStim,
                    xs = dimVals)
  # fig <- data %>% group_by(x) %>%
  #   summarise(mean = mean(y)) %>%
  #   ggplot(aes(y = mean, x = x)) + 
  #   geom_line() + theme_classic() 
  # ggsave(filename = paste0(file_name_root, "meangradient", graph_file_type), fig,
  #        height = gg_height, width = gg_width, units = "cm")
  return(list(data_list, nSubj, nStim, y))
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

Read_Subgroups <- function(fileName, dimVals) {
  
  # This function reads data in the format specified in the NSW02-Data.csv file and 
  # prepares the data list to be inputted to stan. 
  
  # Note that the "x" column must match the "dimVals" argument in this function. 
  # For example, if you have 11 stimuli (S1-S11) and the CS+ is the middle 
  # stimulus (S6), then the position of the CS+ should be 0 in the xs argument, 
  # and range between -.5 to +.5.
  
  # Note that the range of the xs argument can be changed if needed, but then 
  # the parameters of the prior distributions must also be changed accordingly
  # in stan.
  
  data <- read.csv(fileName, header = TRUE)
  
  # similarity
  data_sim <- data %>% 
    filter(group == "diff") %>%
    filter(rule == "Similarity") %>%
    arrange(subj, dimVal)
  subj_sim <- data_sim[["subj"]]
  x_sim <- data_sim[["dimVal"]]
  y_sim <- data_sim[["response"]]
  nSubj_sim <- length(unique(subj_sim))
  nStim_sim <- length(unique(x_sim))
  responses_sim <- matrix(data_sim[["response"]], ncol = nStim_sim, byrow = TRUE)
  data_list_sim <- list(subj = subj_sim,
                        responses = responses_sim,
                        nSubj = nSubj_sim,
                        nStim = nStim_sim,
                        xs = dimVals)
  
  # linear
  data_lin <- data %>% 
    filter(group == "diff") %>%
    filter(rule == "Linear") %>%
    arrange(subj, dimVal)
  subj_lin <- data_lin[["subj"]]
  x_lin <- data_lin[["dimVal"]]
  y_lin <- data_lin[["response"]]
  nSubj_lin <- length(unique(subj_lin))
  nStim_lin <- length(unique(x_lin))
  responses_lin <- matrix(data_lin[["response"]], ncol = nStim_lin, byrow = TRUE)
  data_list_lin <- list(subj = subj_lin,
                        responses = responses_lin,
                        nSubj = nSubj_lin,
                        nStim = nStim_lin,
                        xs = dimVals)
  
  fig <- data %>% 
    filter(group == "diff") %>%
    filter(rule == "Similarity" | rule == "Linear") %>%
    group_by(dimVal, rule) %>%
    summarise(mean = mean(response)) %>%
    ggplot(aes(y = mean, x = dimVal, group = rule, colour = rule)) + 
    geom_line() + theme_classic() 
  ggsave(filename = paste0(file_name_root, "subgroups", graph_file_type), fig,
         height = gg_height, width = gg_width*1.25, units = "cm")
  return(list(data_list_sim, nSubj_sim, nStim_sim, y_sim, 
              data_list_lin, nSubj_lin, nStim_lin, y_lin))
}
#_______________________________________________________________________________

Compare_Groups <- function(samples1, samples2, groupName1, groupName2, graphName,
                           paramNames, dimVals) {
  
  # This function plots the densities of each Gaussian parameter for each of
  # 2 groups
  
  # estimates <- data.frame(group = factor(), param = factor(), sample = numeric())
  # temp <- list()
  # for (i in 1:length(paramNames)) {
  #   param <- paramNames[i]
  #   temp1 <- as.vector(samples1[[param]])
  #   temp2 <- as.vector(samples2[[param]])
  #   temp[i] <- list(list(samples = c(temp1, temp2),
  #                        group = c(rep(groupName1, length(as.vector(temp1))),
  #                                  rep(groupName2, length(as.vector(temp2)))),
  #                        param=param))
  # temp[i] <- list(cbind(c(temp1, temp2),
  #                  c(rep(groupName1, length(as.vector(temp1))),
  #                    rep(groupName2, length(as.vector(temp2)))),
  #                  param))
  # }
  # str(temp)
  # estimates <- map(temp, bind_rows)
  
  # for (i in 1:length(paramNames)) {
  #   estimates$group <- c(temp[i][[1]]$group)
  #   estimates$param <- c(temp[i][[1]]$param)
  #   estimates$sample <- c(temp[i][[1]]$samples)
  #   aa <- cbind(c(temp[i][[1]]$group),
  #               c(temp[i][[1]]$param),
  #               c(temp[i][[1]]$samples))
  # }
  # colnames(estimates) <- paramNames
  # estimates$group <- c(rep(groupName1, length(as.vector(samples1$M))), 
  #                      rep(groupName2, length(as.vector(samples2$M))))
  # estimates$group <- fct_relevel(estimates$group, groupName1, groupName2)
  
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


# Compare_Groups <- function(samples1, samples2, groupName1, groupName2, graphName,
#                            paramNames = c("M", "SD", "height")) {
#   
#   # This function plots the densities of each Gaussian parameter for each of
#   # 2 groups
#   
#   estimates <- as.data.frame(cbind(c(as.vector(samples1$M), as.vector(samples2$M)),
#                                    c(as.vector(samples1$SD), as.vector(samples2$SD)),
#                                    c(as.vector(samples1$height), as.vector(samples2$height))))
#   colnames(estimates) <- paramNames
#   estimates$group <- c(rep(groupName1, length(as.vector(samples1$M))), 
#                        rep(groupName2, length(as.vector(samples2$M))))
#   estimates$group <- fct_relevel(estimates$group, groupName1, groupName2)
#   
#   density_layers <- list(geom_density(alpha = .25),
#                          scale_fill_manual(values = density_cols),
#                          theme_classic())
#   
#   M_fig <- ggplot(estimates, aes(M, fill = group)) + 
#     density_layers +
#     guides(fill = FALSE) +
#     geom_vline(xintercept = mean(samples1$M), linetype = "solid", colour = density_cols[1],
#                size = 2) +
#     geom_vline(xintercept = mean(samples2$M), linetype = "solid", colour = density_cols[2],
#                size = 2) +
#     ggtitle("a) Mean")
#   
#   SD_fig <- ggplot(estimates, aes(SD, fill = group)) + 
#     density_layers +
#     guides(fill = FALSE) +
#     geom_vline(xintercept = mean(samples1$SD), linetype = "solid", colour = density_cols[1],
#                size = 2) +
#     geom_vline(xintercept = mean(samples2$SD), linetype = "solid", colour = density_cols[2],
#                size = 2) +
#     ggtitle("b) Width")
#   
#   height_fig <- ggplot(estimates, aes(height, fill = group)) + 
#     density_layers +
#     geom_vline(xintercept = mean(samples1$height), linetype = "solid", colour = density_cols[1],
#                size = 2) +
#     geom_vline(xintercept = mean(samples2$height), linetype = "solid", colour = density_cols[2],
#                size = 2) +
#     ggtitle("c) Height") +
#     theme(legend.position="none")
#   
#   fig_panel <- grid.arrange(M_fig + theme(legend.title = element_blank()), 
#                             SD_fig + theme(legend.title = element_blank()), 
#                             height_fig + theme(legend.title = element_blank()), nrow = 1)
#   ggsave(paste0(file_name_root, graphName, "density", graph_file_type), fig_panel, 
#          "jpeg", height = gg_height*.8, width = gg_width*2, units = "cm", dpi = dpi)
#   
# }

#_______________________________________________________________________________

Read_Groups <- function(fileName, dimVals, groupName1, groupName2, graphName) {
  
  # This function reads data in the format specified in the NSW09-Data.csv file and 
  # prepares the data list to be inputted to stan. 
  
  # Requires input for group names as strings
  
  # Note that the "x" column must match the "dimVals" argument in this function. 
  # For example, if you have 11 stimuli (S1-S11) and the CS+ is the middle 
  # stimulus (S6), then the position of the CS+ should be 0 in the xs argument, 
  # and range between -.5 to +.5.
  
  # Note that the range of the xs argument can be changed if needed, but then 
  # the parameters of the prior distributions must also be changed accordingly
  # in stan.
  
  data <- read.csv(fileName, header = TRUE)
  
  # group 1 (single cue)
  data_1 <- data %>% 
    filter(group == groupName1) %>%
    arrange(subj, dimVal)
  subj_1 <- data_1[["subj"]]
  x_1 <- data_1[["dimVal"]]
  y_1 <- data_1[["response"]]
  nSubj_1 <- length(unique(subj_1))
  nStim_1 <- length(unique(x_1))
  responses_1 <- matrix(data_1[["response"]], ncol = nStim_1, byrow = TRUE)
  data_list_1 <- list(subj = subj_1,
                      responses = responses_1,
                      nSubj = nSubj_1,
                      nStim = nStim_1,
                      xs = dimVals)
  
  # group 2 (distant neg)
  data_2<- data %>% 
    filter(group == groupName2) %>%
    arrange(subj, dimVal)
  subj_2 <- data_2[["subj"]]
  x_2 <- data_2[["dimVal"]]
  y_2 <- data_2[["response"]]
  nSubj_2 <- length(unique(subj_2))
  nStim_2 <- length(unique(x_2))
  responses_2 <- matrix(data_2[["response"]], ncol = nStim_2, byrow = TRUE)
  data_list_2 <- list(subj = subj_2,
                      responses = responses_2,
                      nSubj = nSubj_2,
                      nStim = nStim_2,
                      xs = dimVals)
  
  fig <- data %>% 
    filter(group == groupName1 | group == groupName2) %>%
    group_by(dimVal, group) %>%
    summarise(mean = mean(response)) %>%
    ggplot(aes(y = mean, x = dimVal, group = group, colour = group)) + 
    geom_line() + theme_classic() 
  ggsave(filename = paste0(file_name_root, graphName, "gradients", graph_file_type), 
         fig, height = gg_height, width = gg_width*1.25, units = "cm")
  return(list(data_list_1, nSubj_1, nStim_1, y_1, 
              data_list_2, nSubj_2, nStim_2, y_2))
}

#_______________________________________________________________________________

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
    # subj_1 <- subset_data[["subj"]]
    # x_1 <- subset_data[["x"]]
    # y_1 <- subset_data[["y"]]
    # nSubj_1 <- length(unique(subj_1))
    # nStim_1 <- length(dimVals)
    # responses_1 <- matrix(data_1[["response"]], ncol = length(dimVals), byrow = TRUE)
    out[i] <- list(list(subj = subset_data[["subj"]],
                        responses = matrix(subset_data[["y"]], ncol = length(dimVals), byrow = TRUE),
                        nSubj = length(unique(subset_data[["subj"]])),
                        nStim = length(dimVals),
                        xs = dimVals))
  }
  
  # # group 1: single
  # data_1 <- data %>% 
  #   filter(group == "single") %>%
  #   arrange(subj, dimVal)
  # subj_1 <- data_1[["subj"]]
  # x_1 <- data_1[["dimVal"]]
  # y_1 <- data_1[["response"]]
  # nSubj_1 <- length(unique(subj_1))
  # nStim_1 <- length(dimVals)
  # responses_1 <- matrix(data_1[["response"]], ncol = length(dimVals), byrow = TRUE)
  # data_list_1 <- list(subj = subj_1,
  #                       responses = responses_1,
  #                       nSubj = nSubj_1,
  #                       nStim = nStim_1,
  #                       xs = dimVals)
  # 
  # # group 2: differential
  # data_2 <- data %>% 
  #   filter(group == "diff") %>%
  #   arrange(subj, dimVal)
  # subj_2 <- data_2[["subj"]]
  # x_2 <- data_2[["dimVal"]]
  # y_2 <- data_2[["response"]]
  # nSubj_2 <- length(unique(subj_2))
  # nStim_2 <- length(dimVals)
  # responses_2 <- matrix(data_2[["response"]], ncol = length(dimVals), byrow = TRUE)
  # data_list_2 <- list(subj = subj_2,
  #                     responses = responses_2,
  #                     nSubj = nSubj_2,
  #                     nStim = nStim_2,
  #                     xs = dimVals)
  # 
  # return(list(data_list_1, nSubj_1, nStim_1, y_1, 
  #             data_list_2, nSubj_2, nStim_2, y_2))
  
  return(list(out, groupNames))
}

#_______________________________________________________________________________

Run_Analysis <- function(fileName, dimVals, nRow = c(6,6), figMult, graphName,
                         paramNames, groupName1, groupName2) {
  
  # handler that calls functions to run analysis and save output
  # note that groupName1 and groupName2 must match those in the data files
  
  # 1. read data
  out <- Read_Gen_Data(fileName, dimVals, groupName1, groupName2)
  data_list_1 <- out[1][[1]][1][[1]]
  data_list_2 <- out[1][[1]][2][[1]]
  # groupNames <- out[2][[1]]
  
  # 2. fit models for each group
  mcmc_out_1 <- Run_GaussianA_Mod(data_list_1, modelName = groupName1)
  samples_1 <- mcmc_out_1[["samples"]]
  mcmc_out_2 <- Run_GaussianA_Mod(data_list_2, modelName = groupName2)
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
}
#_______________________________________________________________________________
