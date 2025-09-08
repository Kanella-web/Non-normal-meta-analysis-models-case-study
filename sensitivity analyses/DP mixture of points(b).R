pre_term_data_58 = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\pre_term_data58.csv")

### DP mixture of points(b) ###
if (requireNamespace("neojags", quietly = TRUE)){
  neojags::load.neojagsmodule()
} 
#> module neojags loaded
if (requireNamespace("neojags", quietly = TRUE)){
  library(rjags)
} 

## impliment Fernandez and Steel skew normal distribution for the base distribution 
library(R2jags)
set.seed(358)

model <- "
model {
  for (i in 1:ns) { 
      m[i] ~ dnorm(0,.0001)
    
    # Mean outcome for group 1
      y1[i] ~ dnorm(mu1[i], prec1[i])
      mu1[i] <- m[i]
      
      # Mean outcome for group 2
      y2[i] ~ dnorm(mu2[i], prec2[i])
      mu2[i] <- m[i] + delta12[i]*pooled_sd[i]
      
    
      delta12[i] <- theta[Z[i]]
    
      Z[i] ~ dcat(p[]) #Z is an integer variable
      
      # Precision parameters for group 1 and 2 based on observed standard deviations
      prec1[i] <- 1 / sd1[i]*sd1[i]
      prec2[i] <- 1 / sd2[i]*sd2[i]

      # Calculate pooled standard deviation
      pooled_sd[i] <- sqrt(((n1[i] - 1)*pow(sd1[i], 2) + (n2[i] - 1)*pow(sd2[i], 2)) / (n1[i] + n2[i] - 2))
  }

  # Constructive DP prior
  # stick-breaking prior
  p[1] <- r[1]
  for (j in 2:(N-1)) {
    p[j] <- r[j] * (1 - r[j-1]) * p[j-1] / r[j-1]
  }
  for (k in 1:(N-1)) {
    r[k] ~ dbeta(1, alpha) T(0, 0.99)
  }
  # assumption to ensure sum p[] is 1 Ishwaran truncation
  ps <- sum(p[1:(N-1)])
  for (k in N:N) {
    p[k] <- 1 - ps
  }

  # Base distribution with normal prior
  for (k in 1:N) {
    
    theta[k] ~ dfskew.t(xi, omega1, nu, shape)
  }

   # Priors
   
      xi ~ dnorm(0, .0001)
      omega1 <- 1 / omega_sqr
      omega_sqr <- omega * omega
      omega ~ dunif(0,10)
      shape ~ dnorm(0.1,5)
      nu ~ dexp(0.10) T(2.5, )
      ## shape ~ dgamma(0.5, 0.318) a prior they propose but it is similae to unif(0,5) we applied

  # concentration parameter prior
  alpha ~ dunif(0.3, 5)
  
  # Random effects distribution mean
  for (i in 1:N) {
    meancl[i] <- p[i] * theta[i]
  }
  poptrue <- sum(meancl[])  # E[X]

  # Random effects distribution variance
  for (i in 1:N) {
    mom2[i] <- p[i] * theta[i] * theta[i]  # E[X2]
  }
  mom2.true <- sum(mom2[])
  var.true <- mom2.true - (poptrue * poptrue)  # E[X2] - E[X]^2

  # Summary statistics
  for (i in 1:ns) {
    for (j in 1:N) {
      SC[i, j] <- equals(j, Z[i])
    }
  }

  # Total clusters K
  for (j in 1:N) {
    cl[j] <- step(sum(SC[, j]) - 1)
  }
  K <- sum(cl[])

  for (i in 1:ns) {
    rho[i] <- (exp(delta12[i]) / (1 + exp(delta12[i])))
  }

  for (i in 1:ns) {
    for (j in 1:ns) {
      equalsmatrix[i, j] <- equals(rho[i], rho[j])
    }
    equalsres[i] <- sum(equalsmatrix[i, ])
  }
  
  pred ~ dcat(p[1:N]) ##First randomly assigning the new study to one of the mixture components, according to the estimated weights p.
  delta_new = theta[pred] ##Then drawing the new study's effect size from the random effects distribution.
  
}
"

inits1 = list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 3143)
inits2 = list(".RNG.name" = "base::Wichmann-Hill",".RNG.seed" = 3144)

dati1  <- list(ns = nrow(pre_term_data_58),
              y1 = pre_term_data_58$mean_FT,
              y2 = pre_term_data_58$mean_EPT.VPT,
              sd1 = pre_term_data_58$sd_FT,
              sd2 = pre_term_data_58$sd_EPT.VPT,
              n1 = pre_term_data_58$n_FT,
              n2 = pre_term_data_58$n_EPT.VPT,
              N = 26)
#pi =3.141593)

model <- jags.model(textConnection(model), inits = list(inits1, inits2),data = dati1, n.chains=2)

samplesv <- coda.samples(model, variable.names =c(
  "theta",
  "xi", ## mu of the Normal base  distribution
  "omega", ## tau of the Normal base distribution
  "shape",
  "nu",
  # "base_mu",
  # "base_tau2",
  "poptrue",   ## overall mu
  "var.true", ## between-study variance
  "delta12", #### study-specific effects
  "K",       ### total number of clusters 
  "p",      ## weights of the process 
  "alpha",  ## concentration parameter
  "SC"   ,  ## probability of each cluster assignment
  "delta_new",
  "pred"
),                      
n.iter = 50000,
n.burnin = 10000)

samples_summary <- summary(samplesv)


##### to check the convergence ###########
library(coda)

# Define the parameters of interest
parameters_of_interest <- c("poptrue", "var.true", "alpha", "xi", "omega", "shape", "nu", "delta_new")

# Subset the samplesv object
subset_samples <- samplesv[, parameters_of_interest]

# Convert to a matrix
samples_matrix <- as.matrix(subset_samples)

# Set up a 4x4 plotting layout
par(mfrow = c(3, 3))  # Adjust rows and columns as needed

# Loop through each parameter and plot its trace
for (i in 1:ncol(samples_matrix)) {
  plot(samples_matrix[, i], type = "l",
       main = colnames(samples_matrix)[i],  # Use parameter name as title
       xlab = "Iteration", ylab = "Value")
}

# Reset plotting layout
par(mfrow = c(1, 1))

# Extract delta12[1:65] parameter names programmatically
##### the first 16 delta12 due to lack of space 
delta12_params <- paste0("delta12[", 1:16, "]")

# Subset samplesv for only delta12 parameters
delta12_samples <- samplesv[, delta12_params]

# Convert to a matrix
delta12_matrix <- as.matrix(delta12_samples)

# Set up the plotting area (e.g., 8x8 grid for 65 parameters)
par(mfrow = c(4, 4))
# Loop through each delta12 parameter and plot
for (i in 1:ncol(delta12_matrix)) {
  plot(delta12_matrix[, i], type = "l",
       main = colnames(delta12_matrix)[i],  # Parameter name as title
       xlab = "Iteration", ylab = "Value")
}

# Reset plotting layout
par(mfrow = c(1, 1))


# Convert the summary statistics into a data frame
DPpresults_U_U_26_sk_t_base1 <- as.data.frame(samples_summary$quantiles)
DPpresults_U_U_26_sk_t_base2 <- as.data.frame(samples_summary$statistics)
DPpresults_U_U_26_sk_t_base <- cbind.data.frame(DPpresults_U_U_26_sk_t_base2,DPpresults_U_U_26_sk_t_base1 )

DPp_results_U_U_26_sk_t_base <-round(DPpresults_U_U_26_sk_t_base, digits=2) 
View(DPp_results_U_U_26_sk_t_base)
View(DPpresults_U_U_26_sk_t_base)

# Probability that poptrue < 0 (use raw samples if available)
if ("poptrue" %in% colnames(as.matrix(samplesv))) {
  Pr_mu_DPp_U_U26_sk_t_base <- mean(as.matrix(samplesv)[, "poptrue"] < 0)
} else {
  stop("The parameter 'poptrue' is not found in the samples.")
}
Pr_mu_DPp_U_U26_sk_t_base


# Probability that poptrue < 0 (use raw samples if available)
if ("delta_new" %in% colnames(as.matrix(samplesv))) {
  Pr_delta_new_DPp_U_U26_sk_t_base <- mean(as.matrix(samplesv)[, "delta_new"] < 0)
} else {
  stop("The parameter 'delta_new' is not found in the samples.")
}
Pr_delta_new_DPp_U_U26_sk_t_base
########## CLUSTER ASSIGNEMENT ##########
DPpresults_U_U_26_sk_t_base$ind <- colnames(as.matrix(samplesv))

############# extraction of the parameters of interest ###########
f55sc_DPp_U_U26_sk_t_base <- grepl("SC", row.names(DPpresults_U_U_26_sk_t_base))
m55sc_DPp_U_U26_sk_t_base <- DPpresults_U_U_26_sk_t_base[f55sc_DPp_U_U26_sk_t_base,]
SC_norm_mix_base <- m55sc_DPp_U_U26_sk_t_base$mean

#### BASED ON THE max(SC[i,j]) ####
##### Distribution of data points to clusters ###

extra_col = rownames(m55sc_DPp_U_U26_sk_t_base)
m55sc_DPp_U_U26_sk_t_base$names = extra_col

prob_list <- data.frame(Column10 = m55sc_DPp_U_U26_sk_t_base[, 10], Column1 = m55sc_DPp_U_U26_sk_t_base[, 1])

split_dataframe <- function(df, chunk_size) {
  split(df, rep(1:ceiling(nrow(df) / chunk_size), each = chunk_size, length.out = nrow(df)))
}

split_prob_list <- split_dataframe(prob_list, nrow(pre_term_data_58))

split_prob_list

# Function to find the maximum probabilities for data points to be distributed into clusters ####
find_max_values_and_indices <- function(split_list) {
  num_rows <- nrow(split_list[[1]])
  num_sublists <- length(split_list)
  
  max_values <- numeric(num_rows)
  max_indices <- integer(num_rows)
  
  for (i in 1:num_rows) {
    values <- sapply(split_list, function(x) x[i, 2])
    max_value <- max(values, na.rm = TRUE)
    max_index <- which(values == max_value)[1]  
    max_values[i] <- max_value
    max_indices[i] <- max_index
  }
  
  return(list(max_values = max_values, max_indices = max_indices))
}

result <- find_max_values_and_indices(split_prob_list)
max_values <- result$max_values #### probabilities of cluster assignment
max_indices <- result$max_indices #### cluster assignments 

write.csv(max_values, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\max_cluster_prob[i,j]_pred.csv",row.names=FALSE )  
write.csv(max_indices, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\clusters_pred.csv",row.names=FALSE )  

### clusters with data points ###
k <- unique(max_indices)

positions_DPp_U_U26_sk_t_base <- vector("list", length(k))
names(positions_DPp_U_U26_sk_t_base) <- k

for(value in k) {
  positions_DPp_U_U26_sk_t_base[[as.character(value)]] <- which(max_indices == value)
}

cluster1 = positions_DPp_U_U26_sk_t_base[["1"]] 
cluster2 = positions_DPp_U_U26_sk_t_base[["2"]] 
cluster5 = positions_DPp_U_U26_sk_t_base[["5"]] 
cluster10 = positions_DPp_U_U26_sk_t_base[["10"]] 


##### resutls are saved here as the order of clusters presented in the main article
write.csv(cluster1, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\cluster1_pred.csv",row.names=FALSE )  
write.csv(cluster2, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\cluster2_pred.csv",row.names=FALSE )  
write.csv(cluster5, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\cluster5_pred.csv",row.names=FALSE )  
write.csv(cluster10, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\cluster10_pred.csv",row.names=FALSE )  


fdd <- grepl("delta12", row.names(DPpresults_U_U_26_sk_t_base))
mdd <- DPpresults_U_U_26_sk_t_base[fdd,]
rel_effDPp_U_U26_sk_t_base <- mdd$`50%`
LB_rel_effDPp_U_U26_sk_t_base <- mdd$`2.5%`
UB_rel_effDPp_U_U26_sk_t_base <- mdd$`97.5%`
sd_rel_effDPp_U_U26_sk_t_base <- mdd$SD

fdd1 <- grepl("theta", row.names(DPpresults_U_U_26_sk_t_base))
mdd1 <- DPpresults_U_U_26_sk_t_base[fdd1,]
theta_DPp_U_U26_sk_t_base <- mdd1$`50%`

plot(density(theta_DPp_U_U26_sk_t_base))

#### DENSITY PLOT OF THE RELATIVE EFFECTS #######
plot(density(rel_effDPp_U_U26_sk_t_base))

#### Mean of each cluster ######
cl1 = rel_effDPp_U_U26_sk_t_base[cluster1]
mean_cl1 = mean(cl1)
cl2 = rel_effDPp_U_U26_sk_t_base[cluster2]
mean_cl2 = mean(cl2)
cl5 = rel_effDPp_U_U26_sk_t_base[cluster5]
mean_cl5 = mean(cl5)

#### DENSITY PLOT OF THE RELATIVE EFECTS WITH THE CLUSTER MEANS
plot(density(rel_effDPp_U_U26_sk_t_base))

clusters_means_DPp_26U_U_sk_t_base = cbind.data.frame(mean_cl1, mean_cl2, mean_cl5, rel_effDPp_U_U26_sk_t_base[cluster10])
write.csv(clusters_means_DPp_26U_U_sk_t_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\clusters_mean.csv")  

abline(v= clusters_means_DPp_26U_U_sk_t_base)
abline(v= clusters_means_DPp_26U_U_sk_t_base,
       col = c("black", "#b71c1c","#0d47a1",  "#737373" ) ,
       lwd = 2, lty = c(1,1,1,2))

# Calculate densities for both clusters 
density1 <- density(rel_effDPp_U_U26_sk_t_base[cluster1])
density2 <- density(rel_effDPp_U_U26_sk_t_base[cluster2])
density5 <- density(rel_effDPp_U_U26_sk_t_base[cluster5])

# Plot the first density
plot(density1, col = "black", lwd = 2, main = "Density Plot of Two Clusters",
     xlab = "Value", ylab = "Density", xlim= c(-2,0), ylim = c(0, 6))

# Add the second density
lines(density2, col = "#b71c1c", lwd = 2)
lines(density5, col = "#0d47a1", lwd = 2)

# Add a legend
legend("topright", legend = c("Cluster 1", "Cluster 2", "Cluster 3"), 
       col = c("black", "#b71c1c", "#0d47a1" ), lty = 1, lwd = 2)


#### performing m-a using the metafor package 
#### library(metafor)

c = rma(yi = rel_effDPp_U_U26_sk_t_base , vi = (sd_rel_effDPp_U_U26_sk_t_base)^2, method = "REML",
        slab = first_second_authors)

forest(c, cex= 0.5)

library(metafor)
###### TO VISUALIZE EACH MAIN CLUSTER ###########
###### MAIN CLUSTER 1 ########
c1 = rma(yi = rel_effDPp_U_U26_sk_t_base[cluster1] , vi = (sd_rel_effDPp_U_U26_sk_t_base[cluster1])^2, 
         method = "REML",  slab = pre_term_data_58$Study[cluster1])
##### DENSITY PLOT #######
plot(density(rel_effDPp_U_U26_sk_t_base[cluster1]))
##### FOREST PLOT #######
forest(c1)
###### MAIN CLUSTER 2 ########
c2 = rma(yi = rel_effDPp_U_U26_sk_t_base[cluster2] , vi = (sd_rel_effDPp_U_U26_sk_t_base[cluster2])^2, 
         method = "REML",  slab = pre_term_data_58$Study[cluster2])
##### DENSITY PLOT #######
plot(density(rel_effDPp_U_U26_sk_t_base[cluster2]))
##### FOREST PLOT #######
forest(c2)

###### MAIN CLUSTER 3 ########
c5 = rma(yi = rel_effDPp_U_U26_sk_t_base[cluster5] , vi = (sd_rel_effDPp_U_U26_sk_t_base[cluster5])^2, 
         method = "REML",  slab = pre_term_data_58$Study[cluster5])
##### DENSITY PLOT #######
plot(density(rel_effDPp_U_U26_sk_t_base[cluster5]))
##### FOREST PLOT #######
forest(c5)
########### SAVE THE REST OF THE RESULTS ###########

f11_U_U26 <- grepl("poptrue", row.names(DPpresults_U_U_26_sk_t_base))
m11_U_U26 <- DPpresults_U_U_26_sk_t_base[f11_U_U26,]
mu_DPp_U_U26_sk_t_base<- m11_U_U26$`50%`
LB_mu_DPp_U_U26_sk_t_base <- m11_U_U26$`2.5%`
UB_mu_DPp_U_U26_sk_t_base <- m11_U_U26$`97.5%`
precDPp_mu_U_U26_sk_t_base <- UB_mu_DPp_U_U26_sk_t_base - LB_mu_DPp_U_U26_sk_t_base


f117_U_U26 <- grepl("delta_new", row.names(DPpresults_U_U_26_sk_t_base))
m117_U_U26 <- DPpresults_U_U_26_sk_t_base[f117_U_U26,]
delta_new_DPp_U_U26_sk_t_base<- m117_U_U26$`50%`
LB_delta_new_DPp_U_U26_sk_t_base <- m117_U_U26$`2.5%`
UB_delta_new_DPp_U_U26_sk_t_base <- m117_U_U26$`97.5%`
precDPp_delta_new_U_U26_sk_t_base <- UB_delta_new_DPp_U_U26_sk_t_base - LB_delta_new_DPp_U_U26_sk_t_base


f22_U_U26 <- grepl("var.true", row.names(DPpresults_U_U_26_sk_t_base))
m22_U_U26 <- DPpresults_U_U_26_sk_t_base[f22_U_U26,]
tau2_DPp_U_U26_sk_t_base <- m22_U_U26$`50%`
LB_tau2_DPp_U_U26_sk_t_base <- m22_U_U26$`2.5%`
UB_tau2_DPp_U_U26_sk_t_base <- m22_U_U26$`97.5%`
precDPp_tau2_U_U26_sk_t_base <- UB_tau2_DPp_U_U26_sk_t_base -  LB_tau2_DPp_U_U26_sk_t_base


f33_U_U26 <- grepl("xi", row.names(DPpresults_U_U_26_sk_t_base))
m33_U_U26 <- DPpresults_U_U_26_sk_t_base[f33_U_U26,]
xi_DPp_U_U26_sk_t_base <- m33_U_U26$`50%`
LB_xi_DPp_U_U26_sk_t_base <- m33_U_U26$`2.5%`
UB_xi_DPp_U_U26_sk_t_base <- m33_U_U26$`97.5%`
prec_xi_DPp_U_U26_sk_t_base <- UB_xi_DPp_U_U26_sk_t_base - LB_xi_DPp_U_U26_sk_t_base

f44_U_U26 <- grepl("omega", row.names(DPpresults_U_U_26_sk_t_base))
m44_U_U26 <- DPpresults_U_U_26_sk_t_base[f44_U_U26,]
omega_DPp_U_U26_sk_t_base <- m44_U_U26$`50%`
omega2_DPp_U_U26_sk_t_base <- (m44_U_U26$`50%`)^2
LB_omega2_DPp_U_U26_sk_t_base <- (m44_U_U26$`2.5%`)^2
UB_omega2_DPp_U_U26_sk_t_base <- (m44_U_U26$`97.5%`)^2
prec_omega2_DPp_U_U26_sk_t_base <- UB_omega2_DPp_U_U26_sk_t_base - LB_omega2_DPp_U_U26_sk_t_base

f55_U_U26 <- grepl("shape", row.names(DPpresults_U_U_26_sk_t_base))
m55_U_U26 <- DPpresults_U_U_26_sk_t_base[f55_U_U26,]
shape_DPp_U_U26_sk_t_base <- m55_U_U26$`50%`
LB_shape_DPp_U_U26_sk_t_base <- m55_U_U26$`2.5%`
UB_shape_DPp_U_U26_sk_t_base <- m55_U_U26$`97.5%`
prec_shape_DPp_U_U26_sk_t_base <- UB_shape_DPp_U_U26_sk_t_base - LB_shape_DPp_U_U26_sk_t_base

f551_U_U26 <- grepl("nu", row.names(DPpresults_U_U_26_sk_t_base))
m551_U_U26 <- DPpresults_U_U_26_sk_t_base[f551_U_U26,]
df_DPp_U_U26_sk_t_base <- m551_U_U26$`50%`
LB_df_DPp_U_U26_sk_t_base <- m551_U_U26$`2.5%`
UB_df_DPp_U_U26_sk_t_base <- m551_U_U26$`97.5%`
prec_df_DPp_U_U26_sk_t_base <- UB_df_DPp_U_U26_sk_t_base - LB_df_DPp_U_U26_sk_t_base

# f111_U_U26 <- grepl("base_mu", row.names(DPpresults_U_U_26_sk_t_base))
# m111_U_U26 <- DPpresults_U_U_26_sk_t_base[f111_U_U26,]
# base_mu_DPp_U_U26_sk_t_base<- m111_U_U26$`50%`
# LB_base_mu_DPp_U_U26_sk_t_base <- m111_U_U26$`2.5%`
# UB_base_mu_DPp_U_U26_sk_t_base <- m111_U_U26$`97.5%`
# precDPp_base_mu_U_U26_sk_t_base <- UB_base_mu_DPp_U_U26_sk_t_base - LB_base_mu_DPp_U_U26_sk_t_base
# 
# f222_U_U26 <- grepl("base_tau2", row.names(DPpresults_U_U_26_sk_t_base))
# m222_U_U26 <- DPpresults_U_U_26_sk_t_base[f222_U_U26,]
# base_tau2_DPp_U_U26_sk_t_base <- m222_U_U26$`50%`
# LB_base_tau2_DPp_U_U26_sk_t_base <- m222_U_U26$`2.5%`
# UB_base_tau2_DPp_U_U26_sk_t_base <- m222_U_U26$`97.5%`
# precDPp_base_tau2_U_U26_sk_t_base <- UB_base_tau2_DPp_U_U26_sk_t_base -  LB_base_tau2_DPp_U_U26_sk_t_base


f55_U_U26 <- grepl("alpha", row.names(DPpresults_U_U_26_sk_t_base))
m55_U_U26 <- DPpresults_U_U_26_sk_t_base[f55_U_U26,]
alpha_DPp_U_U26_sk_t_base <- m55_U_U26$`50%`
LB_alpha_DPp_U_U26_sk_t_base <- m55_U_U26$`2.5%`
UB_alpha_DPp_U_U26_sk_t_base <- m55_U_U26$`97.5%`
prec_alpha_DPp_U_U26_sk_t_base <- UB_alpha_DPp_U_U26_sk_t_base - LB_alpha_DPp_U_U26_sk_t_base

fcl_U_U26 <- grepl("K", row.names(DPpresults_U_U_26_sk_t_base))   
mcl_U_U26 <- DPpresults_U_U_26_sk_t_base[fcl_U_U26,]
median_K_DPp_U_U26_sk_t_base <- mcl_U_U26$`50%` 
LB_K_DPp_U_U26_sk_t_base <- mcl_U_U26$`2.5%`
UB_K_DPp_U_U26_sk_t_base <- mcl_U_U26$`97.5%`

fp_U_U26 <- grepl("p", row.names(DPpresults_U_U_26_sk_t_base))
mp_U_U26 <- DPpresults_U_U_26_sk_t_base[fp_U_U26,]
mp_U_U26 <- mp_U_U26[!grepl("pop", row.names(mp_U_U26)),]
mp_U_U26 <- mp_U_U26[!grepl("alpha", row.names(mp_U_U26)),]
pi_DPp_U_U26_sk_t_base <- mp_U_U26$mean


listaDPp_U_U26_sk_t_base <- cbind.data.frame(rel_effDPp_U_U26_sk_t_base,sd_rel_effDPp_U_U26_sk_t_base, LB_rel_effDPp_U_U26_sk_t_base,
                                             UB_rel_effDPp_U_U26_sk_t_base)


numclus_DPp_U_U26_sk_t_base <- unlist(median_K_DPp_U_U26_sk_t_base)
LB_K_DPp_U_U26_sk_t_base <- unlist(LB_K_DPp_U_U26_sk_t_base)
UB_K_DPp_U_U26_sk_t_base <- unlist(UB_K_DPp_U_U26_sk_t_base)

resDPp_U_U26_sk_t_base <- cbind.data.frame(
  #base_mu_DPp_U_U26_sk_t_base, LB_base_mu_DPp_U_U26_sk_t_base, UB_base_mu_DPp_U_U26_sk_t_base,
  #base_tau2_DPp_U_U26_sk_t_base, LB_base_tau2_DPp_U_U26_sk_t_base, UB_base_tau2_DPp_U_U26_sk_t_base,
  xi_DPp_U_U26_sk_t_base,LB_xi_DPp_U_U26_sk_t_base,UB_xi_DPp_U_U26_sk_t_base,
  omega_DPp_U_U26_sk_t_base,
  omega2_DPp_U_U26_sk_t_base,LB_omega2_DPp_U_U26_sk_t_base,UB_omega2_DPp_U_U26_sk_t_base,
  shape_DPp_U_U26_sk_t_base, LB_shape_DPp_U_U26_sk_t_base, UB_shape_DPp_U_U26_sk_t_base,
  df_DPp_U_U26_sk_t_base, LB_df_DPp_U_U26_sk_t_base, UB_df_DPp_U_U26_sk_t_base,
  mu_DPp_U_U26_sk_t_base,LB_mu_DPp_U_U26_sk_t_base, UB_mu_DPp_U_U26_sk_t_base,
  delta_new_DPp_U_U26_sk_t_base,LB_delta_new_DPp_U_U26_sk_t_base, UB_delta_new_DPp_U_U26_sk_t_base,
  tau2_DPp_U_U26_sk_t_base, LB_tau2_DPp_U_U26_sk_t_base, UB_tau2_DPp_U_U26_sk_t_base,
  #precDPp_base_mu_U_U26_sk_t_base, precDPp_base_tau2_U_U26_sk_t_base,
  prec_df_DPp_U_U26_sk_t_base,
  precDPp_mu_U_U26_sk_t_base, precDPp_delta_new_U_U26_sk_t_base,
  precDPp_tau2_U_U26_sk_t_base, prec_xi_DPp_U_U26_sk_t_base, prec_omega2_DPp_U_U26_sk_t_base,
  prec_alpha_DPp_U_U26_sk_t_base,  alpha_DPp_U_U26_sk_t_base , LB_alpha_DPp_U_U26_sk_t_base, UB_alpha_DPp_U_U26_sk_t_base,
  numclus_DPp_U_U26_sk_t_base, LB_K_DPp_U_U26_sk_t_base, UB_K_DPp_U_U26_sk_t_base,
  Pr_mu_DPp_U_U26_sk_t_base, Pr_delta_new_DPp_U_U26_sk_t_base)

extra_col = row.names(DPpresults_U_U_26_sk_t_base)
DPpresults_U_U_26_sk_t_base$extra_col = extra_col
prob_DPp_U_U26_sk_t_base = round(max_values,2)


write.csv(resDPp_U_U26_sk_t_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\res_DPp_26_U_U_sk_t_base_pred.csv",row.names=FALSE )  
write.csv(listaDPp_U_U26_sk_t_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\rel_eff_DPp_26_U_U_sk_t_base_pred.csv",row.names=FALSE )  
write.csv(prob_DPp_U_U26_sk_t_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\max_prob_cluster_DPp_26_U_U_sk_t_base_pred.csv",row.names=FALSE )  
write.csv(clusters_means_DPp_26U_U_sk_t_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\clusters_means_DPp_26_U_U_sk_t_base_pred.csv",row.names=FALSE )  
write.csv(DPpresults_U_U_26_sk_t_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\all_res_DPp_26_U_U_sk_t_base_pred.csv",row.names=FALSE )  


########## TO CREATE BEATIFULL FOREST PLOT  ###########
max_values = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\max_prob_cluster_DPp_26_U_U_sk_t_base_pred.csv")  
max_indices = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\clusters_pred.csv")  
rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\rel_eff_DPp_26_U_U_sk_t_base_pred.csv")  
es = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\res_DPp_26_U_U_sk_t_base_pred.csv")  
cl_mean = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\clusters_mean.csv")  

mean_cl1 = cl_mean$mean_cl1
mean_cl2 = cl_mean$mean_cl2
mean_cl3 = cl_mean$mean_cl5

numbers <- unlist(max_indices[[1]])

# Define a color palette
color_palette <- c("black", "#b71c1c","#0d47a1",  "#737373" ) 

# Assign colors based on the unique numbers in the vector
unique_numbers <- unique(numbers)
color_map <- setNames(color_palette[seq_along(unique_numbers)], unique_numbers)

# Create the color vector
color_vector <- color_map[as.character(numbers)]

# # View the color vector
# color_vector
# 
# 

# Unique numbers in the original vector
unique_numbers <- unique(numbers)

replacement_values <- c(16, 15, 14, 12)  

replacement_map <- setNames(replacement_values, unique_numbers)

new_numbers <- replacement_map[as.character(numbers)]

##new_numbers

######## forest plot ########
clusters = unlist(max_indices[[1]])
cl_vec <- paste0('"', clusters, '"', collapse = ", ")
cl = cat(cl_vec)
clusters = c("2", "2", "1", "1", "2", "2", "2", "1", "2", "2", "1", "2", "1", "2", 
             "2", "2", "2", "2", "1", "2", "1", "3", "2", "2", "2", "3", "2", "2", 
             "2", "2", "2", "2", "1", "3", "1", "2", "2", "3", "2", "3", "2", "2", 
             "2", "4", "2", "1", "2", "2", "1", "2", "2", "2", "3", "1", "3", "2",
             "2", "2")

pr <- unlist(max_values[[1]])
pr <- as.vector(round(pr, 2))
# 
# string_vec <- paste0('"', pr, '"', collapse = ", ")
# pr = cat(string_vec)
# pr = c("0.18", "0.18", "0.14", "0.14", "0.15", "0.17", "0.17", "0.14", "0.15", "0.14", "0.16", "0.15",
#        "0.14", "0.14", "0.14", "0.16", "0.18", "0.11", "0.18", "0.17", "0.14", "0.1", "0.08/0.08", "0.1", 
#        "0.18", "0.14", "0.15", "0.11", "0.17", "0.12", "0.18", "0.17", "0.14", "0.14", "0.1", "0.13",
#        "0.17", "0.1", "0.14", "0.17", "0.15", "0.11", "0.15", "0.11", "0.16", "0.11", "0.14", "0.18",
#        "0.16", "0.07", "0.13", "0.14", "0.18", "0.16", "0.18", "0.14", "0.11", "0.18", "0.1", "0.11",
#        "0.09", "0.17", "0.16", "0.18")

PMB = pre_term_data_58$Birth_level

matched =pre_term_data_58$Matched
matched = ifelse(matched == 1, "yes", "no")

mean_ga = pre_term_data_58$Mean..ga
mean_bw = pre_term_data_58$Mean.bw
median_by = pre_term_data_58$Median.by
quality = pre_term_data_58$Quality.
country = pre_term_data_58$Country
age = pre_term_data_58$Mean.age


add_data = cbind.data.frame(mean_ga, mean_bw ,matched, PMB,  clusters, pr) 

##Combine the effect size data with additional study info
sorted_data <- cbind.data.frame(
  rel_eff, add_data, color_vector, new_numbers, pre_term_data_58$Study
)



# Sort by PMB first, then by effect size within each color
sorted_data <- sorted_data[order(sorted_data$color_vector, sorted_data$PMB,
                                 as.numeric(trimws(sorted_data$mean_bw)),
                                 as.numeric(trimws(sorted_data$mean_ga)),
                                 #sorted_data$rel_effDPp_U_U26_normal_base,
                                 
                                 sorted_data$matched,
                                 
                                 decreasing = TRUE), ]

# Extract sorted values
sorted_rel_eff <- sorted_data$rel_effDPp_U_U26_sk_t_base
sorted_LB <- sorted_data$LB_rel_effDPp_U_U26_sk_t_base
sorted_UB <- sorted_data$UB_rel_effDPp_U_U26_sk_t_base
sorted_authors <- sorted_data$`pre_term_data_58$Study`
sorted_color <- sorted_data$color_vector
sorted_pch <- sorted_data$new_numbers
sorted_add_data <- sorted_data[, c(   "mean_ga", "mean_bw" ,"matched", "PMB",  "clusters", "pr")]

library(metafor)

# Generate the sorted forest plot
forest(x = sorted_rel_eff, 
       ci.lb = sorted_LB, 
       ci.ub = sorted_UB,
       slab = sorted_authors,
       xlim = c(-6,2), 
       psize = 2,
       cex = 0.45,
       ilab = sorted_add_data, 
       ilab.xpos =  c(-5,-4.5,-4,-3.5, 0.5, 1), 
       lwd = 1.4,
       col = sorted_color, 
       pch = sorted_pch,
       refline =  FALSE)

abline(v=c(mean_cl1, mean_cl2, mean_cl3), col = c("black", "#b71c1c","#0d47a1"), lwd=1.5)

addpoly(x= es$mu_DPp_U_U26_sk_t_base, ci.lb = es$LB_mu_DPp_U_U26_sk_t_base ,
        ci.ub = es$UB_mu_DPp_U_U26_sk_t_base , rows=-0.2)
arrows(x0 = es$LB_delta_new_DPp_U_U26_sk_t_base,
       x1 = es$UB_delta_new_DPp_U_U26_sk_t_base,
       y0 = -1.2, y1 = -1.2,
       angle = 90, code = 3, length = 0.05, lty = 1, lwd = 1.5, col =  "#D2691E")
# if you want to add the point estimate of the PI 
#points(es$delta_new_DPp_U_U26_normal_base, -1.2, pch = 15, col = "black")
abline(h=0.4, lwd=0.1, col="black", lty=1)

# Get row numbers where PMB changes
pm_levels <- sorted_data$PMB
row_positions <- which(diff(as.numeric(as.factor(pm_levels))) != 0)

# Adjust for plotting direction (forest() plots top row at high y)
n_rows <- length(pm_levels)
adjusted_rows <- n_rows - row_positions + 0.5  # +0.5 centers the line between groups


# Add horizontal lines after the forest plot
for (r in adjusted_rows) {
  abline(h = r, col = "gray60", lty = "dotted", lwd = 0.7)
}
abline(h = 7.5, col = "gray60", lty = "dotted", lwd = 0.7)




## SAVE THE FOREST PLOT IN A SPECIFIC FILE WITH 300 dpi
# Open TIFF device
tiff("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\DPMp_ST_base.tiff", 
     units = "in", width = 15, height = 9, res = 300)


# Generate the sorted forest plot
forest(x = sorted_rel_eff, 
       ci.lb = sorted_LB, 
       ci.ub = sorted_UB,
       slab = sorted_authors,
       xlim = c(-6,2), 
       psize = 2,
       cex = 0.45,
       ilab = sorted_add_data, 
       ilab.xpos =  c(-5,-4.5,-4,-3.5, 0.5, 1), 
       lwd = 1.4,
       col = sorted_color, 
       pch = sorted_pch,
       refline =  FALSE)

abline(v=c(mean_cl1, mean_cl2, mean_cl3), col = c("black", "#b71c1c","#0d47a1"), lwd=1.5)

addpoly(x= es$mu_DPp_U_U26_sk_t_base, ci.lb = es$LB_mu_DPp_U_U26_sk_t_base ,
        ci.ub = es$UB_mu_DPp_U_U26_sk_t_base , rows=-0.2)
arrows(x0 = es$LB_delta_new_DPp_U_U26_sk_t_base,
       x1 = es$UB_delta_new_DPp_U_U26_sk_t_base,
       y0 = -1.2, y1 = -1.2,
       angle = 90, code = 3, length = 0.05, lty = 1, lwd = 1.5, col =  "#D2691E")
# if you want to add the point estimate of the PI 
#points(es$delta_new_DPp_U_U26_normal_base, -1.2, pch = 15, col = "black")
abline(h=0.4, lwd=0.1, col="black", lty=1)

# Get row numbers where PMB changes
pm_levels <- sorted_data$PMB
row_positions <- which(diff(as.numeric(as.factor(pm_levels))) != 0)

# Adjust for plotting direction (forest() plots top row at high y)
n_rows <- length(pm_levels)
adjusted_rows <- n_rows - row_positions + 0.5  # +0.5 centers the line between groups


# Add horizontal lines after the forest plot
for (r in adjusted_rows) {
  abline(h = r, col = "gray60", lty = "dotted", lwd = 0.7)
}
abline(h = 7.5, col = "gray60", lty = "dotted", lwd = 0.7)
# Close the device to save the file
dev.off()







