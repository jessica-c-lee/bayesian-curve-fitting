data {
  INS_DATA
}

parameters {
   INS_PARAMS
}

transformed parameters {
   INS_T_PARAMS
}

model {
   INS_LIKELIHOOD
   INS_PRIOR
}

generated quantities {
   INS_QUANTS
}
