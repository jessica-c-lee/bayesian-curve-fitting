#-------------------------------------------------------------------------------
#                                BUILD MODEL STRING
#-------------------------------------------------------------------------------

#_______________________________________________________________________________

Build_Model <- function(funcs, data, params, tParams, likelihood, prior, quants) {
  # This function inserts strings in each block in the stan skeleton model
  
  modString <- paste0(read_lines("./models/skeleton.stan"), collapse = "\n")
  modString <- str_replace_all(modString, "INS_FUNCS", funcs)
  modString <- str_replace_all(modString, "INS_DATA", data)
  modString <- str_replace_all(modString, "INS_PARAMS", params)
  modString <- str_replace_all(modString, "INS_T_PARAMS", tParams)
  modString <- str_replace_all(modString, "INS_LIKELIHOOD", likelihood)
  modString <- str_replace_all(modString, "INS_PRIOR", prior)
  modString <- str_replace_all(modString, "INS_QUANTS", quants)
}
#_______________________________________________________________________________
# write augmented gaussian model
write(Build_Model(funcs = mod_funcs, 
                  data = paste0(mod_data, "\n", mod_data_3ex), 
                  params = mod_params_nonhier, 
                  tParams = paste0(mod_tparams_nonhier, "\n", subj_loop_open, "\n", stim_loop_open), 
                  exfuncs = paste0(mod_exfuncs_nonhier, "\n", mod_expred_av3_nonhier, "\n", loop_close, "\n", loop_close), 
                  likelihood = paste0(subj_loop_open, "\n", stim_loop_open, "\n", mod_like_singav_nonhier, "\n", loop_close), 
                  prior = paste0(mod_priors_unif_nonhier, "\n", loop_close), 
                  quants = mod_genq_av3_nonhier), 
      file = paste0("models/", "av3-nh-u", ".stan"))
#_______________________________________________________________________________