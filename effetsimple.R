install.packages("C:\\Users\\think\\Desktop\\DOBAD_1.0.6.tar.gz", repos = NULL, type = "source")
library(DOBAD)

# parameters
Lamb <- .3; mu <- .4; nu <- 1

# similuate birth depth processes
N <- 30; t <- 10
X <- list() # initialisation
X_t0 = sample(6:15, N, replace = TRUE) # start points
for (i in 1:N){
  X <- c(X, birth.death.simulant(t=t, X0=X_t0[i], lambda=Lamb, mu=mu, nu=nu)) # simulation
}

# likelihood for init parameter
ll_init <- numeric(N)
for (i in 1:N){
  latent_variables <- sufficient_stats(X[[i]])
  U_k <- latent_variables[[1]]; D_k <- latent_variables[[2]]; T_k <- latent_variables[[3]]
  len <- length(U_k_temp)
  k_vector <- seq_along(U_k)
  ll_init[i] <- sum(U_k * log(k_vector*0.3089843 + 0.1006474)) + sum(D_k * log(0.1030675)) - sum(T_k * (k_vector*(0.3089843+0.1030675) + 0.1006474))
}
mean(ll_init)



# k_bar for full process
maximum <- 0
minimum <- 100
for (i in 1:N){
  if (max(X[[i]]@states)>maximum){maximum <- max(X[[i]]@states)}
  if (min(X[[i]]@states)<minimum){minimum <- min(X[[i]]@states)}
}

for (i in 1:N){
  plot(X[[i]]@times, X[[i]]@states, type = "l", xlab = "Time", ylab = "States", main = "BDMC States over Time")
  legend("topleft", legend = c(round(Lambs[i],digits = 3), round(mus[i], digits = 3), round(nus[i], digits = 3)), col = c("blue", "red", "green"), lwd = 2)
}
# get discret observations : 4 observation per process
discret_times <- c(0, .21,.62,.73, 1.44, 1.95, 3.56, 4.17)
X_discret <- lapply(X, function(x) getPartialData(discret_times, x))

sufficient_stats <- function(bd_process){
  # access the states and corresponding times of the bd_process
  states <- bd_process@states
  times <- bd_process@times
  
  # get the maximum state value
  max_k <- max(states)
  
  # initialize empty vectors to store U_k, D_k, and T_k for each state k
  U_k <- numeric(max_k)
  D_k <- numeric(max_k)
  T_k <- numeric(max_k)
  
  # iterate through each state
  for (k in unique(states)) {
    # find indices where the state equals k
    indices <- which(states == k)
    
    # calculate the number of up steps (births) and down steps (deaths)
    jumps <- c(diff(states), 0)   # treat the "out of index" issue
    U_k[k] <- sum(jumps[indices] > 0)
    D_k[k] <- sum(jumps[indices] < 0)
    
    # calculate the total time spent in state k
    times_spent <- c(diff(times),0)
    T_k[k] <- sum(times_spent[indices])
  }
  
  return (list(U_k, D_k, T_k))
}


# Psi
Psi <- function(lambda, mu, nu, k_max) {
  k <- 1:k_max
  f <- numeric(k_max * 3)
  f[1:(k_max * 3)] <- c(log(k * lambda + nu), rep(log(mu),k_max) , -k * (lambda + mu) + nu)
  return(f)
}

# critère
cri <- function(lamb, mu, eta, k_bar){
  (eta/k_bar + lamb)/mu
}


# Stochastic Proximal Gradient Algorithms
SGD <- function(X_discret, learning_rate = 0.01, epochs = 100, batch_size = 1, N_latent_sim = 10, k_max = 50) {
  n <- length(X_discret)
  
  # initialize parameters
  Lamb <- runif(1); mu <- runif(1); nu = runif(1)
  cat("initialisation : ", Lamb,mu,nu)
  
  # Iterate over epochs
  for (epoch in 1:epochs) {
    
    # Shuffle the data
    indices <- sample(n)
    X_discret <- X_discret[indices]
    
    # Iterate over batches
    for (i in seq(1, N, by = batch_size)) {
      # Extract batch
      X_batch_discret <- X_discret[i:min(i + batch_size - 1, n)]
      # initialize latent variables matrix for the batch
      batch_S_Z <- matrix(0, nrow = 3, ncol = batch_size)
      ll <- numeric(N)
      for (batch_ind in 1:batch_size){
        # conditional simulation
        X_batch_cond <- sim.condBD(N=N_latent_sim, bd.PO=X_batch_discret[[batch_ind]], L=Lamb, m=mu, nu=nu) # latent variable x
        
        # compute U_k, D_k, T_k by the algorithm in "Stochastic Proximal Gradient Algorithms for Penalized Mixed Models"
        U_k_mat <- matrix(0, nrow = k_max, ncol = N_latent_sim)
        D_k_mat <- matrix(0, nrow = k_max, ncol = N_latent_sim)
        T_k_mat <- matrix(0, nrow = k_max, ncol = N_latent_sim)
        pi_z <- numeric(N_latent_sim)
        
        for (j in 1:N_latent_sim) {
          # get sufficient statistics
          latent_variables <- sufficient_stats(X_batch_cond[[j]])
          U_k_temp <- latent_variables[[1]]; D_k_temp <- latent_variables[[2]]; T_k_temp <- latent_variables[[3]]
          len <- length(U_k_temp)
          
          # (non normalized) measure of latent variable x
          psi <- Psi(Lamb, mu, nu, len)
          pi_z[j] <- c(U_k_temp, D_k_temp, T_k_temp) %*% psi
          
          U_k_mat[1:len,j] <- U_k_temp
          D_k_mat[1:len,j] <- D_k_temp
          T_k_mat[1:len,j] <- T_k_temp
        }
      
        # integral of S(x) wrt the measure pi_theta(x)
        pi_z <- exp(pi_z) / sum(exp(pi_z))
        U_k <- U_k_mat %*% pi_z
        D_k <- D_k_mat %*% pi_z
        T_k <- T_k_mat %*% pi_z
      
        # compute gradient (by jacobien matrix of psi(theta))
        k_vector <- seq_along(U_k)
        k_max <- max(k_vector); k_min <- min(k_vector)
        batch_S_Z[1,batch_ind]  <- sum(U_k * (k_vector / (k_vector * Lamb + nu)) - k_vector * T_k)
        batch_S_Z[2,batch_ind] <- sum(D_k / mu - k_vector * T_k)
        batch_S_Z[3,batch_ind] <- sum(U_k / (k_vector * Lamb + nu) - T_k)
        ll[batch_ind] <- sum(U_k * log(k_vector*Lamb + nu)) + sum(D_k * log(mu)) - sum(T_k * (k_vector*(Lamb+mu) + nu))
      }
      vraisamblence <- mean(ll)
      G <- rowMeans(batch_S_Z)
      # update parameters
      Lamb <- Lamb + learning_rate * G[1]
      mu <- mu + learning_rate * G[2]
      nu <- nu + learning_rate * G[3]
    }
    c <- cri(Lamb, mu, nu, (k_max+k_min)/2)
    # Print for monitoring
    cat("Epoch:", epoch, " batch:", i," Parameters:", Lamb, mu, nu, " vraisamblance", vraisamblence," critere:",c, "\n")
  }
  
  return(list(Lamb, mu, nu))
}


parameters <- SGD(X_discret, learning_rate = 0.01, epochs = 15, batch_size = 30, N_latent_sim = 1)

