pre_term_data_58 = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\pre_term_data58.csv")

set.seed(458)
library(R2jags)
cat("
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
  
  # Constructive DP
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

  # Baseline distribution with spike-and-mixab prior 
  ### mixab is concentrated at 1 and spike at 0
  for (k in 1:N) {
    del[k] ~ dbern(pi)  # Inclusion indicator
    theta[k] ~ dnorm((ifelse(del[k] == 1, mu11, mu22)), (ifelse(del[k] == 1, tau1p, tau2p)))  # 
  }

   # Priors
   pi ~ dbeta(1, 5)  # more probabilities of values between [0,1] closer to the component with lower variance
  
   tau1p <- 1 / tau1_sqr  
   tau1_sqr = tau1 * tau1
   tau1 ~ dunif(0,10)  
   
   #tau2p = 3/ 0.05
   
   tau2p <- 1 / tau2_sqr  
   tau2_sqr = tau2 * tau2
   tau2 ~ dunif(0,0.001)  

   mu11 ~ dnorm(0, 0.0001)   
   mu22 ~ dnorm(0, 0.0001) 
   #mu_spike = 0
   
   # mu_mixab ~ dnorm(0, t1)  # Prior for mu_mixab
  
  
  
  # concentration parameter prior
  alpha ~ dunif(0.3, 5)
  
  # t1 = 1/ t*t
  # t ~ dnorm(0, 1)T(0,)

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
  
}", file = "DPp_mix.txt")
modfile = 'DPp_mix.txt'


datiDP  <- list(ns = nrow(pre_term_data_58),
               y1 = pre_term_data_58$mean_FT,
               y2 = pre_term_data_58$mean_EPT.VPT,
               sd1 = pre_term_data_58$sd_FT,
               sd2 = pre_term_data_58$sd_EPT.VPT,
               n1 = pre_term_data_58$n_FT,
               n2 = pre_term_data_58$n_EPT.VPT,
               N = 26)

run.modelDPp26_U_U_U_mix_base = jags(
  data = datiDP ,
  inits = NULL,
  parameters.to.save = c(
    "theta",
    "poptrue",   ## overall mu
    "var.true", ## between-study variance 
    "delta12", #### study-specific effects
    "K",       ### number of clusters 
    "p",      ## weights of the process 
    "alpha",  ## concentration parameter
    "SC",    ## probability of each cluster assignment
    "mu11", ## 
    "mu22", ### 
    "tau1", ## 
    "del", #### spike and mix component assignment 
    "pi" , ### prob of mu11 and mix component assignment 
    "pred",
    "delta_new"
    #"equalsres",
    #"equalsmatrix"
  ),  
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

DPpresults_U_U_U_26_mix_base <- as.data.frame(run.modelDPp26_U_U_U_mix_base$BUGSoutput$summary) 

DPp_results_U_U_U_26_mix_base <- as.data.frame(run.modelDPp26_U_U_U_mix_base$BUGSoutput$summary) 

DPp_results_U_U_U_26_mix_base <-round(DPpresults_U_U_U_26_mix_base, digits=2) 

View(DPp_results_U_U_U_26_mix_base)

######## Pr(mu<0) #######
Pr_mu_DPp_U_U_U_26_mix_base = mean(run.modelDPp26_U_U_U_mix_base$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )
######## Pr(mu<0) #######
Pr_delta_new_DPp_U_U_U_26_mix_base = mean(run.modelDPp26_U_U_U_mix_base$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )

########## CLUSTER ASSIGNEMENT ##########
DPpresults_U_U_U_26_mix_base$ind <- row.names(DPpresults_U_U_U_26_mix_base)

############# extraction of the parameters of interest ###########
f55sc_DPp_U_U_U_26_mix_base <- grepl("SC", row.names(DPpresults_U_U_U_26_mix_base))
m55sc_DPp_U_U_U_26_mix_base <- DPpresults_U_U_U_26_mix_base[f55sc_DPp_U_U_U_26_mix_base,]
SC_norm_mix_base <- m55sc_DPp_U_U_U_26_mix_base$mean

#### BASED ON  THE max(SC[i,j]) ####
##### Distribution of data points to clusters ###

extra_col = rownames(m55sc_DPp_U_U_U_26_mix_base)
m55sc_DPp_U_U_U_26_mix_base$names = extra_col

prob_list <- data.frame(Column10 = m55sc_DPp_U_U_U_26_mix_base[, 10], Column1 = m55sc_DPp_U_U_U_26_mix_base[, 1])

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

write.csv(max_values, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\max_prob_mix.csv")  
write.csv(max_indices, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\clustering_assignment_mix.csv")  

### clusters with data points ###
k <- unique(max_indices)

positions_DPp_U_U_U_26_mix_base <- vector("list", length(k))
names(positions_DPp_U_U_U_26_mix_base) <- k

for(value in k) {
  positions_DPp_U_U_U_26_mix_base[[as.character(value)]] <- which(max_indices == value)
}

#### cluster with 22 elements ###
cluster1 = positions_DPp_U_U_U_26_mix_base[["1"]]
#### cluster with 29 elements ###
cluster4 = positions_DPp_U_U_U_26_mix_base[["4"]]

#### clusterS with 2 elementS ###
cluster3 = positions_DPp_U_U_U_26_mix_base[["3"]]
cluster5 = positions_DPp_U_U_U_26_mix_base[["5"]]
#### cluster with 1 element ###
cluster8 = positions_DPp_U_U_U_26_mix_base[["8"]]

write.csv(cluster1, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\cluster1_mix.csv",row.names=FALSE )  
write.csv(cluster4, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\cluster4_mix.csv",row.names=FALSE )  

write.csv(cluster3, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\cluster3_mix.csv",row.names=FALSE )  
write.csv(cluster5, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\cluster5_mix.csv",row.names=FALSE )  
write.csv(cluster8, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\cluster8_mix.csv",row.names=FALSE )  


#### cluster with only 1 element ###
fdd1 <- grepl("theta", row.names(DPpresults_U_U_U_26_mix_base))
mdd1 <- DPpresults_U_U_U_26_mix_base[fdd1,]
thetaDPp_U_U_U_26_mix_base <- mdd1$`50%`

plot(density(thetaDPp_U_U_U_26_mix_base))


fdd <- grepl("delta12", row.names(DPpresults_U_U_U_26_mix_base))
mdd <- DPpresults_U_U_U_26_mix_base[fdd,]
rel_effDPp_U_U_U_26_mix_base <- mdd$`50%`
LB_rel_effDPp_U_U_U_26_mix_base <- mdd$`2.5%`
UB_rel_effDPp_U_U_U_26_mix_base <- mdd$`97.5%`
sd_rel_effDPp_U_U_U_26_mix_base <- mdd$sd
Rhat_deltaDPp_U_U_U_26_mix_base <- mdd$Rhat

#### Mean of each cluster ######
cl1 = rel_effDPp_U_U_U_26_mix_base[cluster1]
mean_cl1 = mean(cl1)
cl4 = rel_effDPp_U_U_U_26_mix_base[cluster4]
mean_cl4 = mean(cl4)


cl3 = rel_effDPp_U_U_U_26_mix_base[cluster3]
mean_cl3 = mean(cl3)
cl5 = rel_effDPp_U_U_U_26_mix_base[cluster5]
mean_cl5 = mean(cl5)

clusters_means_DPp_26HNU_mix_base = c(mean_cl1, mean_cl4, mean_cl3, mean_cl5,
                                      rel_effDPp_U_U_U_26_mix_base[cluster8])


library(metafor)

###### TO VISUALIZE EACH MAIN CLUSTER ###########
###### MAIN CLUSTER 1 ########
c1 = rma(yi = rel_effDPp_U_U_U_26_mix_base[cluster1] , vi = (sd_rel_effDPp_U_U_U_26_mix_base[cluster1])^2, 
         method = "REML",  slab = pre_term_data_58$Study[cluster1])
##### DENSITY PLOT #######
plot(density(rel_effDPp_U_U_U_26_mix_base[cluster1]))
##### FOREST PLOT #######
forest(c1)
###### MAIN CLUSTER 4 ########
c4 = rma(yi = rel_effDPp_U_U_U_26_mix_base[cluster4] , vi = (sd_rel_effDPp_U_U_U_26_mix_base[cluster4])^2,
         method = "REML", slab = pre_term_data_58$Study[cluster4])
##### DENSITY PLOT #######
plot(density(rel_effDPp_U_U_U_26_mix_base[cluster4]))
##### FOREST PLOT #######
forest(c4)


###### TO VISUALISE ALL THE DENSITIES OF CLUSTERS TOGETHER ##########
# Create density objects for each cluster
density1 <- density(rel_effDPp_U_U_U_26_mix_base[cluster1])
density4 <- density(rel_effDPp_U_U_U_26_mix_base[cluster4])

# Plot the first density
plot(density1, col = "black", lwd = 2, main = "Density Plot for Clusters", 
     xlab = "Relative Efficiency", ylab = "Density", xlim = c(-2,0.2), ylim = c(0, 5))

# Add the second and third densities
lines(density4, col = "#b71c1c", lwd = 2)


# Add a legend for clarity
legend("topright", legend = c("Cluster 1", "Cluster 2"), 
       col = c( "#b71c1c", "black"), lwd = 2)


########### SAVE THE REST OF THE RESULTS ###########

f11_U_U_U_26 <- grepl("poptrue", row.names(DPpresults_U_U_U_26_mix_base))
m11_U_U_U_26 <- DPpresults_U_U_U_26_mix_base[f11_U_U_U_26,]
mu_DPp_U_U_U_26_mix_base<- m11_U_U_U_26$`50%`
LB_mu_DPp_U_U_U_26_mix_base <- m11_U_U_U_26$`2.5%`
UB_mu_DPp_U_U_U_26_mix_base <- m11_U_U_U_26$`97.5%`
Rhat_muDPp_U_U_U_26_mix_base <- m11_U_U_U_26$Rhat
precDPp_mu_U_U_U_26_mix_base <- UB_mu_DPp_U_U_U_26_mix_base - LB_mu_DPp_U_U_U_26_mix_base

f17_U_U_U_26 <- grepl("delta_new", row.names(DPpresults_U_U_U_26_mix_base))
m17_U_U_U_26 <- DPpresults_U_U_U_26_mix_base[f17_U_U_U_26,]
delta_new_DPp_U_U_U_26_mix_base<- m17_U_U_U_26$`50%`
LB_delta_new_DPp_U_U_U_26_mix_base <- m17_U_U_U_26$`2.5%`
UB_delta_new_DPp_U_U_U_26_mix_base <- m17_U_U_U_26$`97.5%`
Rhat_delta_newDPp_U_U_U_26_mix_base <- m17_U_U_U_26$Rhat
precDPp_delta_new_U_U_U_26_mix_base <- UB_delta_new_DPp_U_U_U_26_mix_base - LB_delta_new_DPp_U_U_U_26_mix_base

f22_U_U_U_26 <- grepl("var.true", row.names(DPpresults_U_U_U_26_mix_base))
m22_U_U_U_26 <- DPpresults_U_U_U_26_mix_base[f22_U_U_U_26,]
tau2_DPp_U_U_U_26_mix_base <- m22_U_U_U_26$`50%`
LB_tau2_DPp_U_U_U_26_mix_base <- m22_U_U_U_26$`2.5%`
UB_tau2_DPp_U_U_U_26_mix_base <- m22_U_U_U_26$`97.5%`
Rhat_tauDPp_U_U_U_26_mix_base <- m22_U_U_U_26$Rhat
precDPp_tau2_U_U_U_26_mix_base <- UB_tau2_DPp_U_U_U_26_mix_base -  LB_tau2_DPp_U_U_U_26_mix_base


f33_U_U_U_26 <- grepl("mu11", row.names(DPpresults_U_U_U_26_mix_base))
m33_U_U_U_26 <- DPpresults_U_U_U_26_mix_base[f33_U_U_U_26,]
mu1_mix_DPp_U_U_U_26_mix_base <- m33_U_U_U_26$`50%`
LB_mu1_mix_DPp_U_U_U_26_mix_base <- m33_U_U_U_26$`2.5%`
UB_mu1_mix_DPp_U_U_U_26_mix_base <- m33_U_U_U_26$`97.5%`
Rhat_mu1_mix_DPp_U_U_U_26_mix_base <- m33_U_U_U_26$Rhat
prec_mu1_mix_DPp_U_U_U_26_mix_base <- UB_mu1_mix_DPp_U_U_U_26_mix_base - LB_mu1_mix_DPp_U_U_U_26_mix_base

f331_U_U_U_26 <- grepl("mu22", row.names(DPpresults_U_U_U_26_mix_base))
m331_U_U_U_26 <- DPpresults_U_U_U_26_mix_base[f331_U_U_U_26,]
mu2_mix_DPp_U_U_U_26_mix_base <- m331_U_U_U_26$`50%`
LB_mu2_mix_DPp_U_U_U_26_mix_base <- m331_U_U_U_26$`2.5%`
UB_mu2_mix_DPp_U_U_U_26_mix_base <- m331_U_U_U_26$`97.5%`
Rhat_mu2_mix_DPp_U_U_U_26_mix_base <- m331_U_U_U_26$Rhat
prec_mu2_mix_DPp_U_U_U_26_mix_base <- UB_mu2_mix_DPp_U_U_U_26_mix_base - LB_mu2_mix_DPp_U_U_U_26_mix_base


f44_U_U_U_26 <- grepl("tau1", row.names(DPpresults_U_U_U_26_mix_base))
m44_U_U_U_26 <- DPpresults_U_U_U_26_mix_base[f44_U_U_U_26,]
tau1_mix_sqr_DPp_U_U_U_26_mix_base <- m44_U_U_U_26$`50%`
LB_tau1_mix_sqr_DPp_U_U_U_26_mix_base <- m44_U_U_U_26$`2.5%`
UB_tau1_mix_sqr_DPp_U_U_U_26_mix_base <- m44_U_U_U_26$`97.5%`
Rhat_tau1_mix_sqr_DPp_U_U_U_26_mix_base <- m44_U_U_U_26$Rhat
prec_tau1_mix_sqr_DPp_U_U_U_26_mix_base <- UB_tau1_mix_sqr_DPp_U_U_U_26_mix_base - LB_tau1_mix_sqr_DPp_U_U_U_26_mix_base


f55_U_U_U_26 <- grepl("alpha", row.names(DPpresults_U_U_U_26_mix_base))
m55_U_U_U_26 <- DPpresults_U_U_U_26_mix_base[f55_U_U_U_26,]
alpha_DPp_U_U_U_26_mix_base <- m55_U_U_U_26$`50%`
LB_alpha_DPp_U_U_U_26_mix_base <- m55_U_U_U_26$`2.5%`
UB_alpha_DPp_U_U_U_26_mix_base <- m55_U_U_U_26$`97.5%`
Rhat_alpha_DPp_U_U_U_26_mix_base <- m55_U_U_U_26$Rhat
prec_alpha_DPp_U_U_U_26_mix_base <- UB_alpha_DPp_U_U_U_26_mix_base - LB_alpha_DPp_U_U_U_26_mix_base

fcl_U_U_U_26 <- grepl("K", row.names(DPpresults_U_U_U_26_mix_base))   
mcl_U_U_U_26 <- DPpresults_U_U_U_26_mix_base[fcl_U_U_U_26,]
median_K_DPp_U_U_U_26_mix_base <- mcl_U_U_U_26$`50%` 
LB_K_DPp_U_U_U_26_mix_base <- mcl_U_U_U_26$`2.5%`
UB_K_DPp_U_U_U_26_mix_base <- mcl_U_U_U_26$`97.5%`

fp_U_U_U_26 <- grepl("p", row.names(DPpresults_U_U_U_26_mix_base))
mp_U_U_U_26 <- DPpresults_U_U_U_26_mix_base[fp_U_U_U_26,]
mp_U_U_U_26 <- mp_U_U_U_26[!grepl("pop", row.names(mp_U_U_U_26)),]
mp_U_U_U_26 <- mp_U_U_U_26[!grepl("alpha", row.names(mp_U_U_U_26)),]
pi_DPp_U_U_U_26_mix_base <- mp_U_U_U_26$mean


listaDPp_U_U_U_26_mix_base <- cbind.data.frame(rel_effDPp_U_U_U_26_mix_base,sd_rel_effDPp_U_U_U_26_mix_base, LB_rel_effDPp_U_U_U_26_mix_base,
                                            UB_rel_effDPp_U_U_U_26_mix_base,Rhat_deltaDPp_U_U_U_26_mix_base)


numclus_DPp_U_U_U_26_mix_base <- unlist(median_K_DPp_U_U_U_26_mix_base)
LB_K_DPp_U_U_U_26_mix_base <- unlist(LB_K_DPp_U_U_U_26_mix_base)
UB_K_DPp_U_U_U_26_mix_base <- unlist(UB_K_DPp_U_U_U_26_mix_base)

resDPp_U_U_U_26_mix_base <- cbind.data.frame(mu2_mix_DPp_U_U_U_26_mix_base,mu1_mix_DPp_U_U_U_26_mix_base,
                                          LB_mu2_mix_DPp_U_U_U_26_mix_base, LB_mu1_mix_DPp_U_U_U_26_mix_base,
                                          UB_mu2_mix_DPp_U_U_U_26_mix_base,UB_mu1_mix_DPp_U_U_U_26_mix_base,
                                          tau1_mix_sqr_DPp_U_U_U_26_mix_base,
                                          LB_tau1_mix_sqr_DPp_U_U_U_26_mix_base,LB_tau1_mix_sqr_DPp_U_U_U_26_mix_base,
                                          mu_DPp_U_U_U_26_mix_base,LB_mu_DPp_U_U_U_26_mix_base, 
                                          UB_mu_DPp_U_U_U_26_mix_base, 
                                          delta_new_DPp_U_U_U_26_mix_base,LB_delta_new_DPp_U_U_U_26_mix_base, 
                                          UB_delta_new_DPp_U_U_U_26_mix_base, 
                                          tau2_DPp_U_U_U_26_mix_base, LB_tau2_DPp_U_U_U_26_mix_base, UB_tau2_DPp_U_U_U_26_mix_base,
                                          precDPp_mu_U_U_U_26_mix_base,  precDPp_delta_new_U_U_U_26_mix_base, 
                                          precDPp_tau2_U_U_U_26_mix_base, 
                                          prec_mu2_mix_DPp_U_U_U_26_mix_base, prec_mu1_mix_DPp_U_U_U_26_mix_base, 
                                          prec_tau1_mix_sqr_DPp_U_U_U_26_mix_base,
                                          prec_alpha_DPp_U_U_U_26_mix_base,  
                                          alpha_DPp_U_U_U_26_mix_base , LB_alpha_DPp_U_U_U_26_mix_base, UB_alpha_DPp_U_U_U_26_mix_base,
                                          Rhat_muDPp_U_U_U_26_mix_base, Rhat_delta_newDPp_U_U_U_26_mix_base, 
                                          Rhat_tauDPp_U_U_U_26_mix_base, 
                                          Rhat_mu2_mix_DPp_U_U_U_26_mix_base,Rhat_mu1_mix_DPp_U_U_U_26_mix_base,
                                          Rhat_tau1_mix_sqr_DPp_U_U_U_26_mix_base, 
                                          Rhat_alpha_DPp_U_U_U_26_mix_base, 
                                          numclus_DPp_U_U_U_26_mix_base, LB_K_DPp_U_U_U_26_mix_base, UB_K_DPp_U_U_U_26_mix_base,
                                          Pr_mu_DPp_U_U_U_26_mix_base, Pr_delta_new_DPp_U_U_U_26_mix_base)

extra_col = row.names(DPpresults_U_U_U_26_mix_base)
DPpresults_U_U_U_26_mix_base$extra_col = extra_col
prob_DPp_U_U_U_26_mix_base = round(max_values,2)


write.csv(resDPp_U_U_U_26_mix_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\res_DPp_26_U_U_U_mix_base.csv",row.names=FALSE )  
write.csv(listaDPp_U_U_U_26_mix_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\rel_eff_DPp_26_U_U_U_mix_base.csv",row.names=FALSE )  
write.csv(prob_DPp_U_U_U_26_mix_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\max_prob_cluster_DPp_26_U_U_U_mix_base.csv",row.names=FALSE )  
write.csv(clusters_means_DPp_26HNU_mix_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\clusters_means_DPp_26_U_U_U_mix_base.csv",row.names=FALSE )  
write.csv(DPpresults_U_U_U_26_mix_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\all_res_DPp_26_U_U_U_mix_base.csv",row.names=FALSE )  

########## TO CREATE BEATIFULL FOREST PLOT  ###########
max_values = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\max_prob_mix.csv")  
max_indices = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\clustering_assignment_mix.csv")  
es = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\res_DPp_26_U_U_U_mix_base.csv")  
rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\rel_eff_DPp_26_U_U_U_mix_base.csv")  
mean_cl = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\clusters_means_DPp_26_U_U_U_mix_base.csv")  

mean_cl1 = mean_cl$x[1]
mean_cl2 = mean_cl$x[2]


# Original vector
numbers <- read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\clustering_assignment_mix.csv")  
numbers <- unlist(max_indices[[2]])
clusters = numbers
cl_vec <- paste0('"', clusters, '"', collapse = ", ")
cl = cat(cl_vec)
clusters  <- c("2", "2", "1", "1", "2", "2", "2", "1", "2", "2", "1", "2", "1", "2", "2",
               "2", "2", "2", "1", "2", "1", "2", "2", "2", "2", "3", "2", "2", "2", "2", 
               "2", "2", "1", "3", "1", "2", "2", "2", "2", "2", "2", "2", "2", "4", "2", 
               "1", "2", "2", "1", "2", "2", "2", "5", "1", "5", "2", "2", "2")

# Define a color palette
color_palette <- c( "black", "#b71c1c", "grey", "grey", "grey", "grey", "grey")
# Add more colors if needed

# Assign colors based on the unique numbers in the vector
unique_numbers <- unique(numbers)
color_map <- setNames(color_palette[seq_along(unique_numbers)], unique_numbers)

# Create the color vector
color_vector <- color_map[as.character(numbers)]

# View the color vector
color_vector



# Unique numbers in the original vector
unique_numbers <- unique(numbers)

# Replacement values (in order corresponding to the unique numbers)
replacement_values <- c(16, 15, 13, 14, 17, 14,17)  # Adjust the length if needed

# Create a mapping from original numbers to replacement numbers
replacement_map <- setNames(replacement_values, unique_numbers)

# Replace the original numbers with the new values
new_numbers <- replacement_map[as.character(numbers)]

# View the new vector
new_numbers

######## forest plot ########
max_values = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\max_prob_mix.csv")  

pr <- unlist(max_values[[2]])
pr <- as.vector(round(pr, 2))
# string_vec <- paste0('"', pr, '"', collapse = ", ")
# cat(string_vec)
# pr = c("0.17", "0.18", "0.19", "0.19", "0.2", "0.16", "0.18", "0.16", "0.19", "0.15", "0.18", "0.19", "0.19",
#        "0.17", "0.18", "0.16", "0.14", "0.18", "0.14/0.14", "0.18", "0.17", "0.19", "0.11", "0.14/0.13", "0.16", "0.15", 
#        "0.18", "0.18", "0.18", "0.17", "0.16", "0.18", "0.16/0.16", "0.16", "0.18", "0.13", "0.16", "0.16", "0.19", 
#        "0.17", "0.17", "0.15", "0.16", "0.14", "0.17", "0.14", "0.17", "0.17", "0.07", "0.16", "0.15", "0.19", 
#        "0.16", "0.13", "0.17", "0.18", "0.17", "0.15", "0.17", "0.13", "0.14", "0.12", "0.17", "0.17", "0.17")


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
sorted_rel_eff <- sorted_data$rel_effDPp_U_U_U_26_mix_base
sorted_LB <- sorted_data$LB_rel_effDPp_U_U_U_26_mix_base
sorted_UB <- sorted_data$UB_rel_effDPp_U_U_U_26_mix_base
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

abline(v=c(mean_cl1, mean_cl2), 
       col = c(  "black", "#b71c1c" ), lwd = 2)

addpoly(x= es$mu_DPp_U_U_U_26_mix_base, ci.lb = es$LB_mu_DPp_U_U_U_26_mix_base ,
        ci.ub = es$UB_mu_DPp_U_U_U_26_mix_base , rows=-0.2)

arrows(x0 = es$LB_delta_new_DPp_U_U_U_26_mix_base,
       x1 = es$UB_delta_new_DPp_U_U_U_26_mix_base,
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
abline(h = 53.5, col = "gray60", lty = "dotted", lwd = 0.7)

### SAVE THE FOREST PLOT IN A SPECIFIC FILE WITH 300 dpi
# Open TIFF device
tiff("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_U_26_mix_base_58\\DPMp_mix_base.tiff", 
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

abline(v=c(mean_cl1, mean_cl2), 
       col = c(  "black", "#b71c1c" ), lwd = 2)

addpoly(x= es$mu_DPp_U_U_U_26_mix_base, ci.lb = es$LB_mu_DPp_U_U_U_26_mix_base ,
        ci.ub = es$UB_mu_DPp_U_U_U_26_mix_base , rows=-0.2)

arrows(x0 = es$LB_delta_new_DPp_U_U_U_26_mix_base,
       x1 = es$UB_delta_new_DPp_U_U_U_26_mix_base,
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
abline(h = 53.5, col = "gray60", lty = "dotted", lwd = 0.7)

dev.off()





