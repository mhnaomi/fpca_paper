my_packages <- c("MASS", "nlme", "lme4", "lmerTest")
suppressMessages(lapply(my_packages, library, character.only = TRUE))

b_sigma = 1
b_tau_j = c(0, 1)

nsim <- 1000

# Define function for calculating covariance matrix
calcSigma <- function(sigma,X1,X2) {
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      Sigma[i,j] <- exp(-sigma*abs(X1[i]-X2[j]))
    }
  }
  return(Sigma)
}


for(j in 1:2){
  
  b_tau = b_tau_j[j]
  
  F_test <- rep(0, nsim)
  F_crit <- rep(0, nsim)
  t_test <- rep(0, nsim)
  t_p <- rep(0, nsim)
  
for (i in 1:nsim){
  print(i)
  T <- seq(0, 1, by = 0.01)
  G <- 2
  L <- 500
  J <- 5
  V <- 10
  
  G_ind <- rep(1:G, each = L*J*V)
  L_ind <- rep(1:(G*L), each = J*V)
  J_ind <- rep(1:(G*L*J), each = V)
  V_ind <- seq(1, G*L*J*V)
  
  # Generate overall mean function
  mu <- 2*sin(2*pi*T)
  
  # Generate fixed-effect of group
  tau <- matrix(0, ncol = G, nrow = length(T))
  for (g in 1:G){
    tau[, g] <- (b_tau*g) * exp(-(T - 0.4))
  }
  
  # Calculate the covariance matrix
  sigma <- calcSigma(b_sigma, T, T)
  
  # Generate matrices for random effects
  K_U <- 0.8*sigma
  U_values <- matrix(mvrnorm(G*L, rep(0, length(T)), K_U), ncol = G*L, byrow = FALSE)
  
  K_W <- 0.5*sigma
  W_values = matrix(mvrnorm(G*L*J, rep(0, length(T)), K_W), 101, ncol = G*L*J, byrow = FALSE)
  
  K_e <- 0.2*sigma
  e_values <- matrix(mvrnorm(G*L*J*V, rep(0, length(T)), K_e), ncol = G*L*J*V)

  sim_data <- replicate(G*L*J*V, mu) +
    tau[, rep(1:ncol(tau), each = L*J*V)] +
    U_values[, rep(1:ncol(U_values), each = J*V)] +
    W_values[, rep(1:ncol(W_values), each = V)] +
    e_values
  sim_data <- t(sim_data)
  
  ## test1 - hierarchical HT
  y_bar <- apply(sim_data, 2, mean)
  y_g_bar <- apply(sim_data, 2, function(x) {by(x, G_ind, mean)})
  y_gl_bar <- apply(sim_data, 2, function(x) {by(x, L_ind, mean)})
  
  est <- 0
  for (gg in 1:G) {
    est <- est +
      (t(y_gl_bar[((gg-1) * L + 1):(gg * L), ])- replicate(L, y_g_bar[gg, ])) %*%
      (y_gl_bar[((gg-1) * L + 1):(gg * L), ] - t(replicate(L, y_g_bar[gg, ])))
  }
  eigen_vals <- eigen(est/(2*(L-1)))$values
  ratio <- sum(eigen_vals)^2/sum(eigen_vals^2)
  
  SS_group <- sum(apply(y_g_bar, 1, function(x) return((x - y_bar)^2*L*J*V)))
  SS_lesion <- sum((y_gl_bar - y_g_bar[rep(1:G, each = L), ])^2*J*V)
  F_test[i] <- (SS_group/(G-1))/(SS_lesion/(G*(L-1)))
  F_crit[i] <- qf(.95, df1 = ratio, df2 = ratio*G*(L-1))
  
  ## test2 - Mixed effect model
  sim_pca <- prcomp(sim_data, scale. = F)
  pca_data <- data.frame(pc1 = sim_pca$x[,1], group = G_ind, lesion = L_ind, layer = J_ind)
  fit <- summary(lmer(pc1~ factor(group) + (1|lesion) + (1|layer)  , data = pca_data))
  t_test[i] <- fit$coefficients[2,4]
  t_p[i] <- (fit$coefficients[2,5] < 0.05)
}

re <- list(btau = b_tau, 
           bsigma = b_sigma, 
           ftest = F_test,
           fcrit = F_crit,
           ttest = t_test,
           tp = t_p)

filename =paste0("function1_btau", b_tau, "_bsigma", b_sigma, ".rds")
save(re, file = filename)
}
