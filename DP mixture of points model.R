
pre_term_data_58 = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\pre_term_data58.csv")

### DP mixture of points model ###
library(R2jags)
set.seed(158)
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
    
    theta[k] ~ dnorm(basemu, basetau1)  
  }

   # Priors
   
   basetau1 <- 1 / basetau_sqr  
   basetau_sqr = basetau * basetau
   basetau ~ dunif(0,10)  

   basemu ~ dnorm(0, 0.0001)   
   
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
  
}", file = "DPp_26_U_U_normal_base.txt")
modfile = 'DPp_26_U_U_normal_base.txt'

run.modelDPp26_U_U_normal_base = jags(
  dati1  <-list(ns = nrow(pre_term_data_58),
                y1 = pre_term_data_58$mean_FT,
                y2 = pre_term_data_58$mean_EPT.VPT,
                sd1 = pre_term_data_58$sd_FT,
                sd2 = pre_term_data_58$sd_EPT.VPT,
                n1 = pre_term_data_58$n_FT,
                n2 = pre_term_data_58$n_EPT.VPT,
                N = 26
  ) ,  
  inits = NULL,
  parameters.to.save = c(
    "theta",
    "basemu", ## mu of the Normal base  distribution
    "basetau", ## tau of the Normal base distribution
    "poptrue",   ## overall mu
    "var.true", ## between-study variance
    "delta12", #### study-specific effects
    "K",       ### total number of clusters 
    "p",      ## weights of the process 
    "alpha",  ## concentration parameter
    "SC",     ## probability of each cluster assignment
    "delta_new",
    "pred"
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

DPpresults_U_U_26_normal_base <- as.data.frame(run.modelDPp26_U_U_normal_base$BUGSoutput$summary) 

DPp_results_U_U_26_normal_base <- as.data.frame(run.modelDPp26_U_U_normal_base$BUGSoutput$summary) 

DPp_results_U_U_26_normal_base <-round(DPpresults_U_U_26_normal_base, digits=2) 

View(round(DPpresults_U_U_26_normal_base,2))
View(DPpresults_U_U_26_normal_base)

######## Pr(mu<0) #######
Pr_mu_DPp_U_U26_normal_base = mean(run.modelDPp26_U_U_normal_base$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )

######## Pr(mu_new<0) #######
Pr_delta_new_DPp_U_U26_normal_base = mean(run.modelDPp26_U_U_normal_base$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )

########## CLUSTER ASSIGNEMENT ##########
DPpresults_U_U_26_normal_base$ind <- row.names(DPpresults_U_U_26_normal_base)

############# extraction of the parameters of interest ###########
f55sc_DPp_U_U26_normal_base <- grepl("SC", row.names(DPpresults_U_U_26_normal_base))
m55sc_DPp_U_U26_normal_base <- DPpresults_U_U_26_normal_base[f55sc_DPp_U_U26_normal_base,]
SC_norm_mix_base <- m55sc_DPp_U_U26_normal_base$mean

#### BASED ON THE max(SC[i,j]) ####
##### Distribution of data points to clusters ###

extra_col = rownames(m55sc_DPp_U_U26_normal_base)
m55sc_DPp_U_U26_normal_base$names = extra_col

prob_list <- data.frame(Column10 = m55sc_DPp_U_U26_normal_base[, 10], Column1 = m55sc_DPp_U_U26_normal_base[, 1])

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

write.csv(max_values, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\max_prob_[i,j].csv",row.names=FALSE )
write.csv(max_indices, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\cluster_assignement.csv",row.names=FALSE )

### clusters with data points ###
k <- unique(max_indices)

positions_DPp_U_U26_normal_base <- vector("list", length(k))
names(positions_DPp_U_U26_normal_base) <- k

for(value in k) {
  positions_DPp_U_U26_normal_base[[as.character(value)]] <- which(max_indices == value)
}

#### cluster with 40 elements ###
cluster1 = positions_DPp_U_U26_normal_base[["1"]] ### cluster1 is cluster2 in the main article
#### cluster with 11 elements ###
cluster2 = positions_DPp_U_U26_normal_base[["2"]] 
## cluster with 4 elements
cluster4 = positions_DPp_U_U26_normal_base[["4"]] 

## clusters with 1 element
cluster10 = positions_DPp_U_U26_normal_base[["10"]] 
cluster5 = positions_DPp_U_U26_normal_base[["5"]] 
cluster7 = positions_DPp_U_U26_normal_base[["7"]] 

##### resutls are saved here as the order of clusters presented in the main article
write.csv(cluster1, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\cluster1.csv",row.names=FALSE )  
write.csv(cluster2, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\cluster2.csv",row.names=FALSE )  
write.csv(cluster4, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\cluster4.csv",row.names=FALSE )  

write.csv(cluster5, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\cluster5.csv",row.names=FALSE )  
write.csv(cluster7, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\cluster7.csv",row.names=FALSE )  
write.csv(cluster10, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\cluster10.csv",row.names=FALSE )  

################
fdd1 <- grepl("theta", row.names(DPpresults_U_U_26_normal_base))
mdd1 <- DPpresults_U_U_26_normal_base[fdd1,]
thetaDPp_U_U26_normal_base <- mdd1$`50%`

#### DENSITY PLOT OF THE RELATIVE EFFECTS #######
plot(density(thetaDPp_U_U26_normal_base))


fdd <- grepl("delta12", row.names(DPpresults_U_U_26_normal_base))
mdd <- DPpresults_U_U_26_normal_base[fdd,]
rel_effDPp_U_U26_normal_base <- mdd$`50%`
LB_rel_effDPp_U_U26_normal_base <- mdd$`2.5%`
UB_rel_effDPp_U_U26_normal_base <- mdd$`97.5%`
sd_rel_effDPp_U_U26_normal_base <- mdd$sd
Rhat_deltaDPp_U_U26_normal_base <- mdd$Rhat


#### DENSITY PLOT OF THE RELATIVE EFFECTS #######
plot(density(rel_effDPp_U_U26_normal_base))

#### Mean of each cluster ######
cl1 = rel_effDPp_U_U26_normal_base[cluster1]
mean_cl1 = mean(cl1)
cl2 = rel_effDPp_U_U26_normal_base[cluster2]
mean_cl2 = mean(cl2)
cl4 = rel_effDPp_U_U26_normal_base[cluster4]
mean_cl4 = mean(cl4)

##### resutls are saved here as the order of clusters presented in the main article
write.csv(mean_cl1, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\mean_cl1.csv",row.names=FALSE )  
write.csv(mean_cl2, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\mean_cl2.csv",row.names=FALSE )  
write.csv(mean_cl4, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\mean_cl4.csv",row.names=FALSE )  

#### DENSITY PLOT OF THE RELATIVE EFECTS WITH THE CLUSTER MEANS
plot(density(rel_effDPp_U_U26_normal_base))

##cluster 1 is cluster 3, black 
### cluster 2 is cluster 1, red 
### cluster 4 is cluster 2 , blue

clusters_means_DPp_26U_U_normal_base = cbind.data.frame(mean_cl1,mean_cl2,  mean_cl4,
                                                        rel_effDPp_U_U26_normal_base[cluster5],
                                                        rel_effDPp_U_U26_normal_base[cluster7],
                                                        rel_effDPp_U_U26_normal_base[cluster10])
#abline(v= clusters_means_DPp_26U_U_normal_base)
abline(v= clusters_means_DPp_26U_U_normal_base,
       col = c( "black", "#b71c1c","#0d47a1","grey","grey","grey"  ),lty =c(1,1,1,2,2,2),
       lwd = 2)
# Add a legend
## the legend will present the clusters in the order they appear
legend("topright", legend = c("mean cluster 1", "mean cluster 2","mean cluster 3", "one-study cluster"), 
       col = c("#b71c1c","#0d47a1", "black",  "#737373" ), lty =c(1,1,1,2), lwd = 2)

# Calculate densities for both clusters 
density1 <- density(rel_effDPp_U_U26_normal_base[cluster1])
density2 <- density(rel_effDPp_U_U26_normal_base[cluster2])
density4 <- density(rel_effDPp_U_U26_normal_base[cluster4])

# Plot the first density
plot(density1, col ="black", lwd = 2, main = "Density Plot of Two Clusters",
     xlab = "Value", ylab = "Density", xlim= c(-2,0), ylim = c(0,10))

# Add the second density
lines(density2, col = "#b71c1c", lwd = 2)
lines(density4, col = "#0d47a1", lwd = 2)


## the legend will present the clusters in the order they appear
legend("topright", legend = c("cluster 1", "cluster 2", "cluster 3"), 
       col = c("#b71c1c", "#0d47a1","black"  ), lty =c(1,1,1), lwd = 2)


#### performing m-a using the metafor package 
#### library(metafor)

c = rma(yi = rel_effDPp_U_U26_normal_base , vi = (sd_rel_effDPp_U_U26_normal_base)^2, method = "REML",
        slab = first_second_authors)

forest(c, cex= 0.5)

library(metafor)
###### TO VISUALIZE EACH MAIN CLUSTER ###########
###### MAIN CLUSTER 1 ########
c1 = rma(yi = rel_effDPp_U_U26_normal_base[cluster1] , vi = (sd_rel_effDPp_U_U26_normal_base[cluster1])^2, 
         method = "REML" ,slab = pre_term_data_58$Study[cluster1])
##### DENSITY PLOT #######
plot(density(rel_effDPp_U_U26_normal_base[cluster1]))
##### FOREST PLOT #######
forest(c1)
###### MAIN CLUSTER 2 ########
c2 = rma(yi = rel_effDPp_U_U26_normal_base[cluster2] , vi = (sd_rel_effDPp_U_U26_normal_base[cluster2])^2, 
         method = "REML", slab = pre_term_data_58$Study[cluster2] )
##### DENSITY PLOT #######
plot(density(rel_effDPp_U_U26_normal_base[cluster2]))
##### FOREST PLOT #######
forest(c2)

###### MAIN CLUSTER 2 ########
c4 = rma(yi = rel_effDPp_U_U26_normal_base[cluster4] , vi = (sd_rel_effDPp_U_U26_normal_base[cluster4])^2, 
         method = "REML", slab = pre_term_data_58$Study[cluster4] )
##### DENSITY PLOT #######
plot(density(rel_effDPp_U_U26_normal_base[cluster4]))
##### FOREST PLOT #######
forest(c4)


########### SAVE THE REST OF THE RESULTS ###########

f11_U_U26 <- grepl("poptrue", row.names(DPpresults_U_U_26_normal_base))
m11_U_U26 <- DPpresults_U_U_26_normal_base[f11_U_U26,]
mu_DPp_U_U26_normal_base<- m11_U_U26$`50%`
LB_mu_DPp_U_U26_normal_base <- m11_U_U26$`2.5%`
UB_mu_DPp_U_U26_normal_base <- m11_U_U26$`97.5%`
Rhat_muDPp_U_U26_normal_base <- m11_U_U26$Rhat
precDPp_mu_U_U26_normal_base <- UB_mu_DPp_U_U26_normal_base - LB_mu_DPp_U_U26_normal_base

f17_U_U26 <- grepl("delta_new", row.names(DPpresults_U_U_26_normal_base))
m17_U_U26 <- DPpresults_U_U_26_normal_base[f17_U_U26,]
delta_new_DPp_U_U26_normal_base<- m17_U_U26$`50%`
LB_delta_new_DPp_U_U26_normal_base <- m17_U_U26$`2.5%`
UB_delta_new_DPp_U_U26_normal_base <- m17_U_U26$`97.5%`
Rhat_delta_newDPp_U_U26_normal_base <- m17_U_U26$Rhat
precDPp_delta_new_U_U26_normal_base <- UB_delta_new_DPp_U_U26_normal_base - LB_delta_new_DPp_U_U26_normal_base


f22_U_U26 <- grepl("var.true", row.names(DPpresults_U_U_26_normal_base))
m22_U_U26 <- DPpresults_U_U_26_normal_base[f22_U_U26,]
tau2_DPp_U_U26_normal_base <- m22_U_U26$`50%`
LB_tau2_DPp_U_U26_normal_base <- m22_U_U26$`2.5%`
UB_tau2_DPp_U_U26_normal_base <- m22_U_U26$`97.5%`
Rhat_tauDPp_U_U26_normal_base <- m22_U_U26$Rhat
precDPp_tau2_U_U26_normal_base <- UB_tau2_DPp_U_U26_normal_base -  LB_tau2_DPp_U_U26_normal_base


f33_U_U26 <- grepl("basemu", row.names(DPpresults_U_U_26_normal_base))
m33_U_U26 <- DPpresults_U_U_26_normal_base[f33_U_U26,]
base_mu_DPp_U_U26_normal_base <- m33_U_U26$`50%`
LB_base_mu_DPp_U_U26_normal_base <- m33_U_U26$`2.5%`
UB_base_mu_DPp_U_U26_normal_base <- m33_U_U26$`97.5%`
Rhat_base_mu_DPp_U_U26_normal_base <- m33_U_U26$Rhat
prec_base_mu_DPp_U_U26_normal_base <- UB_base_mu_DPp_U_U26_normal_base - LB_base_mu_DPp_U_U26_normal_base

f44_U_U26 <- grepl("basetau", row.names(DPpresults_U_U_26_normal_base))
m44_U_U26 <- DPpresults_U_U_26_normal_base[f44_U_U26,]
base_tau_DPp_U_U26_normal_base <- m44_U_U26$`50%`
base_tau2_DPp_U_U26_normal_base <- (m44_U_U26$`50%`)^2
LB_base_tau2_DPp_U_U26_normal_base <- (m44_U_U26$`2.5%`)^2
UB_base_tau2_DPp_U_U26_normal_base <- (m44_U_U26$`97.5%`)^2
Rhat_base_tau_DPp_U_U26_normal_base <- (m44_U_U26$Rhat)^2
prec_base_tau2_DPp_U_U26_normal_base <- UB_base_tau2_DPp_U_U26_normal_base - LB_base_tau2_DPp_U_U26_normal_base

f55_U_U26 <- grepl("alpha", row.names(DPpresults_U_U_26_normal_base))
m55_U_U26 <- DPpresults_U_U_26_normal_base[f55_U_U26,]
alpha_DPp_U_U26_normal_base <- m55_U_U26$`50%`
LB_alpha_DPp_U_U26_normal_base <- m55_U_U26$`2.5%`
UB_alpha_DPp_U_U26_normal_base <- m55_U_U26$`97.5%`
Rhat_alpha_DPp_U_U26_normal_base <- m55_U_U26$Rhat
prec_alpha_DPp_U_U26_normal_base <- UB_alpha_DPp_U_U26_normal_base - LB_alpha_DPp_U_U26_normal_base

fcl_U_U26 <- grepl("K", row.names(DPpresults_U_U_26_normal_base))   
mcl_U_U26 <- DPpresults_U_U_26_normal_base[fcl_U_U26,]
median_K_DPp_U_U26_normal_base <- mcl_U_U26$`50%` 
LB_K_DPp_U_U26_normal_base <- mcl_U_U26$`2.5%`
UB_K_DPp_U_U26_normal_base <- mcl_U_U26$`97.5%`

fp_U_U26 <- grepl("p", row.names(DPpresults_U_U_26_normal_base))
mp_U_U26 <- DPpresults_U_U_26_normal_base[fp_U_U26,]
mp_U_U26 <- mp_U_U26[!grepl("pop", row.names(mp_U_U26)),]
mp_U_U26 <- mp_U_U26[!grepl("alpha", row.names(mp_U_U26)),]
pi_DPp_U_U26_normal_base <- mp_U_U26$mean


listaDPp_U_U26_normal_base <- cbind.data.frame(rel_effDPp_U_U26_normal_base,sd_rel_effDPp_U_U26_normal_base, LB_rel_effDPp_U_U26_normal_base,
                                               UB_rel_effDPp_U_U26_normal_base,Rhat_deltaDPp_U_U26_normal_base)


numclus_DPp_U_U26_normal_base <- unlist(median_K_DPp_U_U26_normal_base)
LB_K_DPp_U_U26_normal_base <- unlist(LB_K_DPp_U_U26_normal_base)
UB_K_DPp_U_U26_normal_base <- unlist(UB_K_DPp_U_U26_normal_base)

resDPp_U_U26_normal_base <- cbind.data.frame(base_mu_DPp_U_U26_normal_base,LB_base_mu_DPp_U_U26_normal_base,
                                             UB_base_mu_DPp_U_U26_normal_base,base_tau_DPp_U_U26_normal_base,base_tau2_DPp_U_U26_normal_base,
                                             LB_base_tau2_DPp_U_U26_normal_base,
                                             UB_base_tau2_DPp_U_U26_normal_base,
                                             mu_DPp_U_U26_normal_base,
                                             LB_mu_DPp_U_U26_normal_base, 
                                             UB_mu_DPp_U_U26_normal_base, tau2_DPp_U_U26_normal_base, LB_tau2_DPp_U_U26_normal_base, UB_tau2_DPp_U_U26_normal_base,
                                             delta_new_DPp_U_U26_normal_base, LB_delta_new_DPp_U_U26_normal_base, UB_delta_new_DPp_U_U26_normal_base,
                                             
                                             precDPp_mu_U_U26_normal_base,precDPp_delta_new_U_U26_normal_base,
                                             precDPp_tau2_U_U26_normal_base, prec_base_mu_DPp_U_U26_normal_base, prec_base_tau2_DPp_U_U26_normal_base,
                                             prec_alpha_DPp_U_U26_normal_base,  alpha_DPp_U_U26_normal_base , LB_alpha_DPp_U_U26_normal_base, UB_alpha_DPp_U_U26_normal_base,
                                             Rhat_muDPp_U_U26_normal_base, Rhat_delta_newDPp_U_U26_normal_base,
                                             Rhat_tauDPp_U_U26_normal_base, Rhat_base_mu_DPp_U_U26_normal_base,
                                             Rhat_base_tau_DPp_U_U26_normal_base, Rhat_alpha_DPp_U_U26_normal_base, 
                                             numclus_DPp_U_U26_normal_base, LB_K_DPp_U_U26_normal_base, UB_K_DPp_U_U26_normal_base,
                                             Pr_mu_DPp_U_U26_normal_base, Pr_delta_new_DPp_U_U26_normal_base)

extra_col = row.names(DPpresults_U_U_26_normal_base)
DPpresults_U_U_26_normal_base$extra_col = extra_col
prob_DPp_U_U26_normal_base = round(max_values,2)

write.csv(resDPp_U_U26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\res_DPp_26_U_U_normal_base.csv",row.names=FALSE )  
write.csv(listaDPp_U_U26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\rel_eff_DPp_26_U_U_normal_base.csv",row.names=FALSE )  
write.csv(prob_DPp_U_U26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\max_prob_cluster_DPp_26_U_U_normal_base.csv",row.names=FALSE )  
write.csv(clusters_means_DPp_26U_U_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\\\clusters_means_DPp_26_U_U_normal_base.csv",row.names=FALSE )  
write.csv(DPpresults_U_U_26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\all_res_DPp_26_U_U_normal_base.csv",row.names=FALSE )  


