install.packages("DOBAD_1.0.6.tar.gz", repos = NULL, type = "source")
install.packages("truncnorm")
install.packages('TukeyRegion')
library(DOBAD)
library(truncnorm)
library(TukeyRegion)
# parameters
Lamb <- .4; mu <- .6; nu <- .6

# similuate birth depth processes
N <- 10; t <- 10
X <- list() # initialisation
X_t0 = sample(6:15, N, replace = TRUE) # start points

Lambs <- rtruncnorm(N, a = 0, b = Lamb*2, mean = Lamb, sd = 0.02)
mus <- rtruncnorm(N, a = 0, b = mu*2, mean = mu, sd = 0.02)
nus <- rtruncnorm(N, a = 0, b = nu*2, mean = nu, sd = 0.02)



for (i in 1:N){
  X <- c(X, birth.death.simulant(t=t, X0=X_t0[i], lambda=Lambs[i], mu=mus[i], nu=nus[i])) # simulation
}


# visualize a generated process
for (i in 1:N){
  plot(X[[i]]@times, X[[i]]@states, type = "l", xlab = "Time", ylab = "States", main = "BDMC States over Time")
  legend("topleft", legend = c(round(Lambs[i],digits = 3), round(mus[i], digits = 3), round(nus[i], digits = 3)), col = c("blue", "red", "green"), lwd = 2)
}
# get discret observations : 9 observation per process
discret_times <- c(0, 1, 2, 2.5, 3.4, 4.3, 5.2, 6.1, 7.2, 8, 9)
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
  
  return (c(U_k, D_k, T_k))
}

# S
S <- function(lambda_i, mu_i, nu_i, latent_variables, len) {
  S_0 <- c(lambda_i, mu_i, nu_i)
  k <- 1:len
  S_k <-  latent_variables * c(log(k * lambda_i + nu_i), log(k * mu_i) , -k * (lambda_i + mu_i) + nu_i)
  return(c(S_0, S_k))
}

# Psi
Psi <- function(lambda, mu, nu, len) {
  return(c(lambda, mu, nu, rep(1, len*3)))
}

# k_bar
maximum <- 0
minimum <- 100
for (i in 1:N){
  if (max(X_discret[[i]]@states)>maximum){maximum <- max(X_discret[[i]]@states)}
  if (min(X_discret[[i]]@states)<minimum){minimum <- min(X_discret[[i]]@states)}
}

# Stochastic Proximal Gradient Algorithms
SGD <- function(X_discret, learning_rate = 0.01, epochs = 100, batch_size = 1, N_latent_params = 10, N_latent_x = 2) {
  n <- length(X_discret)
  
  # initialize parameters
  Lamb <- 0.3; mu <- 0.3; nu = 0.3
  cat("Initialization parameters:", Lamb, mu, nu, "\n")
  
  # Iterate over epochs
  for (epoch in 1:epochs) {
    
    # Shuffle the data
    indices <- sample(n)
    X_discret <- X_discret[indices]
    
    # Iterate over batches
    for (i in seq(1, N, by = batch_size)) {
      # Extract batch
      X_batch_discret <- X_discret[i:min(i + batch_size - 1, n)]
      
      # generate individual parameters
      Lamb_i <- rtruncnorm(N_latent_params, a = 0, b = 1, mean = Lamb, sd = 0.2)
      mu_i <- rtruncnorm(N_latent_params, a = 0, b = 1, mean = mu, sd = 0.2)
      nu_i <- rtruncnorm(N_latent_params, a = 0, b = 1, mean = nu, sd = 0.2)
      cat("Lamb_i : ",Lamb_i, '\n')
      cat("mu_i : ",mu_i, '\n')
      cat("nu_i : ",nu_i, '\n')
      # initialize latent variables matrix for the batch
      batch_S_Z <- matrix(0, nrow = 3, ncol = batch_size)
      
      for (batch_ind in 1:batch_size){
        # initialize latent variables matrix
        S_z_mat <- matrix(0, nrow = 3, ncol = N_latent_params * N_latent_x)
        
        # initialize vector of measures of latent variables
        pi_z <- numeric(N_latent_params * N_latent_x)
        for (j1 in seq(1, N_latent_params)) {
          # conditional simulation
          X_batch_cond <- sim.condBD(N=N_latent_x, bd.PO=X_batch_discret[[batch_ind]], L=Lamb_i[j1], m=mu_i[j1], nu=nu_i[j1]) # latent variable x
          
          for (j2 in 1:N_latent_x) {
            # get sufficient statistics
            latent_variables <- sufficient_stats(X_batch_cond[[j2]])
            len <- as.integer(length(latent_variables)/3)
            
            S_z <- S(Lamb_i[j1], mu_i[j1], nu_i[j1], latent_variables, len)
            psi <- Psi(Lamb, mu, nu, len)
            S_z_mat[1:3, (j1-1)*N_latent_x + j2] <- S_z[1:3]
            
            # (non normalized) measure of latent variable z
            pi_z[(j1-1)*N_latent_x + j2] <- S_z %*% psi
  
          }
        }
        
        # integral of S(z) wrt the measure pi_theta(z)
        #pi_z <- exp(pi_z) / sum(exp(pi_z))
        #S_Z <- S_z_mat %*% pi_z
        
        S_Z <- S_z_mat[,which.max(pi_z)]
        cat("S_Z : ",S_Z[1:3], "\n")
        batch_S_Z[1:3,batch_ind] <- S_Z[1:3]
      }
      
      batch_mean_S <- rowMeans(batch_S_Z)
      cat("batch_mean : ",batch_mean_S, "\n")
      
      # compute gradient (by jacobien matrix of psi(theta) and phi(theta))
      
      gradient_Lamb <- batch_mean_S[1]-Lamb
      gradient_mu <- batch_mean_S[2]-mu
      gradient_nu <- batch_mean_S[3]-nu
      
      # update parameters
      Lamb <- Lamb + learning_rate * gradient_Lamb
      mu <- mu + learning_rate * gradient_mu
      nu <- nu + learning_rate * gradient_nu
      
      # Print the loss for monitoring
      cat("Epoch:", epoch, " Individual",i, " Parameters:", Lamb, mu, nu, "\n")
    }
    
    # Print the loss for monitoring
    cat("Epoch:", epoch, " Parameters:", Lamb, mu, nu, "\n")
  }
  
  return(list(Lamb, mu, nu))
}

parameters <- SGD(X_discret, learning_rate = 0.1, epochs = 3, batch_size = 10, N_latent_params = 15, N_latent_x = 1)


