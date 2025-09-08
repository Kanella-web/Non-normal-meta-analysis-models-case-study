pre_term_data_58 = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\pre_term_data58.csv")

### Binomial-DPp-n-(Unif/Gamma) ###
library(R2jags)
set.seed(1958)
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

  # Base distribution with normal prior
  for (k in 1:N) {
    
    theta[k] ~ dnorm(basemu, basetau1)  # normal base distribution
  }

   # Priors
   
   basetau1 <- 1 / basetau_sqr  
   basetau_sqr = basetau * basetau
   basetau ~ dunif(0,10)

   basemu ~ dnorm(0, 0.0001)  
   
  # DPP parameter prior
  alpha ~ dgamma(1,1)
  
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

  # Programming for calculating summary statistics
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
  
}", file = "DPp_n_U_G_normal_base.txt")
modfile = 'DPp_n_U_G_normal_base.txt'

run.modelDPpn_U_G_normal_base = jags(
  dati1  <-list(ns = nrow(pre_term_data_58),
                y1 = pre_term_data_58$mean_FT,
                y2 = pre_term_data_58$mean_EPT.VPT,
                sd1 = pre_term_data_58$sd_FT,
                sd2 = pre_term_data_58$sd_EPT.VPT,
                n1 = pre_term_data_58$n_FT,
                n2 = pre_term_data_58$n_EPT.VPT,
                N = nrow(pre_term_data_58)
  )  ,  
  inits = NULL,
  parameters.to.save = c(
    "basemu", ## mu of the Normal base  distribution
    "basetau", ## tau of the Normal base distribution
    "poptrue",   ## overall mu
    "var.true", ## between-study heterogeneity
    "delta12", #### study-specific effects
    "K",       ### total number of clusters 
    "p",      ## weights of the process 
    "alpha",  ## concentration parameter
    "SC" ,    ## probability of each cluster assignment
    "pred",
    "delta_new"
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

DPpresults_U_G_n_normal_base <- as.data.frame(run.modelDPpn_U_G_normal_base$BUGSoutput$summary) 

DPp_results_U_G_n_normal_base <- as.data.frame(run.modelDPpn_U_G_normal_base$BUGSoutput$summary) 

DPp_results_U_G_n_normal_base <-round(DPpresults_U_G_n_normal_base, digits=2) 

View(DPpresults_U_G_n_normal_base)

######## Pr(mu<0) #######
Pr_mu_DPp_UGn_normal_base = mean(run.modelDPpn_U_G_normal_base$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )
######## Pr(pred<0) #######
Pr_delta_new_DPp_UGn_normal_base = mean(run.modelDPpn_U_G_normal_base$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )

########## CLUSTER ASSIGNEMENT ##########
DPpresults_U_G_n_normal_base$ind <- row.names(DPpresults_U_G_n_normal_base)

############# extraction of the parameters of interest ###########
f55sc_DPp_UGn_normal_base <- grepl("SC", row.names(DPpresults_U_G_n_normal_base))
m55sc_DPp_UGn_normal_base <- DPpresults_U_G_n_normal_base[f55sc_DPp_UGn_normal_base,]
SC_norm_mix_base <- m55sc_DPp_UGn_normal_base$mean

#### BASED ON  THE max(SC[i,j]) ####
##### Distribution of data points to clusters ###

extra_col = rownames(m55sc_DPp_UGn_normal_base)
m55sc_DPp_UGn_normal_base$names = extra_col

prob_list <- data.frame(Column10 = m55sc_DPp_UGn_normal_base[, 10], Column1 = m55sc_DPp_UGn_normal_base[, 1])

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

write.csv(max_values, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\max_prob.csv")  
write.csv(max_indices, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\clustering_assignment.csv")  

### clusters with data points ###
k <- unique(max_indices)

positions_DPp_UGn_normal_base <- vector("list", length(k))
names(positions_DPp_UGn_normal_base) <- k

for(value in k) {
  positions_DPp_UGn_normal_base[[as.character(value)]] <- which(max_indices == value)
}

cluster1 = positions_DPp_UGn_normal_base[["1"]]
cluster2 = positions_DPp_UGn_normal_base[["2"]]
cluster4 = positions_DPp_UGn_normal_base[["4"]]

cluster8 = positions_DPp_UGn_normal_base[["8"]] 


write.csv(cluster1, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\cluster1.csv",row.names=FALSE )  
write.csv(cluster2, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\cluster2.csv",row.names=FALSE )  
write.csv(cluster4, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\cluster4.csv",row.names=FALSE )  

write.csv(cluster8, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\cluster8.csv",row.names=FALSE )  

fdd <- grepl("delta12", row.names(DPpresults_U_G_n_normal_base))
mdd <- DPpresults_U_G_n_normal_base[fdd,]
rel_effDPp_UGn_normal_base <- mdd$`50%`
LB_rel_effDPp_UGn_normal_base <- mdd$`2.5%`
UB_rel_effDPp_UGn_normal_base <- mdd$`97.5%`
sd_rel_effDPp_UGn_normal_base <- mdd$sd
Rhat_deltaDPp_UGn_normal_base <- mdd$Rhat


#### PLOT THE DENSITY OF THE RELATIVE EFFECTS #######
plot(density(rel_effDPp_UGn_normal_base))

#### Mean of each cluster ######
cl1 = rel_effDPp_UGn_normal_base[cluster1]
mean_cl1 = mean(cl1)
cl2 = rel_effDPp_UGn_normal_base[cluster2]
mean_cl2 = mean(cl2)
cl4 = rel_effDPp_UGn_normal_base[cluster4]
mean_cl4 = mean(cl4)


mean_cl = c(mean_cl1, mean_cl2, mean_cl4, rel_effDPp_UGn_normal_base[cluster8])
write.csv(mean_cl, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\clusters_mean.csv")  

plot(density(rel_effDPp_UGn_normal_base),lwd = 2)

abline(v= mean_cl,
       col = c("black","#0d47a1", "#b71c1c", "grey" ),
       lwd = 2)
# Add a legend
legend("topright", legend = c("mean cluster 1", "mean cluster 2", "mean cluster 3"
), 
col = c("#0d47a1", "#b71c1c", "black" ), lty = 1, lwd = 2)
# Calculate densities for both clusters
density1 <- density(rel_effDPp_UGn_normal_base[cluster1])
density2 <- density(rel_effDPp_UGn_normal_base[cluster2])
density4 <- density(rel_effDPp_UGn_normal_base[cluster4])

# Plot the first density
plot(density1, col ="black", lwd = 2, main = "Density Plot of four Clusters",
     xlab = "Value", ylab = "Density", xlim= c(-1.5,0), ylim = c(0,10) )

# Add the second density
lines(density2, col = "#b71c1c", lwd = 2)
# Add the second density
lines(density4, col =  "#0d47a1", lwd = 2)


# Add a legend
legend("topright", legend = c("Cluster 1", "Cluster 2","Cluster 3"), 
       col = c("black", "#b71c1c","#0d47a1" ), lty = 1, lwd = 2)

library(metafor)
###### TO VISUALIZE EACH MAIN CLUSTER ###########
# ###### MAIN CLUSTER 1 ########
c1 = rma(yi = rel_effDPp_UGn_normal_base[cluster1] , vi = (sd_rel_effDPp_UGn_normal_base[cluster1])^2,
         method = "REML", slab = pre_term_data_58$Study[cluster1])
##### DENSITY PLOT #######
plot(density(rel_effDPp_UGn_normal_base[cluster1]))
##### FOREST PLOT #######
forest(c1)

###### MAIN CLUSTER 2 ########
c2 = rma(yi = rel_effDPp_UGn_normal_base[cluster2] , vi = (sd_rel_effDPp_UGn_normal_base[cluster2])^2,
         method = "REML", slab = pre_term_data_58$Study[cluster2])
##### DENSITY PLOT #######
plot(density(rel_effDPp_UGn_normal_base[cluster2]))
##### FOREST PLOT #######
forest(c2)

###### MAIN CLUSTER 4 ########
c4 = rma(yi = rel_effDPp_UGn_normal_base[cluster4] , vi = (sd_rel_effDPp_UGn_normal_base[cluster4])^2, 
         method = "REML", slab = pre_term_data_58$Study[cluster4])
##### DENSITY PLOT #######
plot(density(rel_effDPp_UGn_normal_base[cluster4]))
##### FOREST PLOT #######
forest(c4)


########### SAVE THE REST OF THE RESULTS ###########

f11_UGn <- grepl("poptrue", row.names(DPpresults_U_G_n_normal_base))
m11_UGn <- DPpresults_U_G_n_normal_base[f11_UGn,]
mu_DPp_UGn_normal_base<- m11_UGn$`50%`
LB_mu_DPp_UGn_normal_base <- m11_UGn$`2.5%`
UB_mu_DPp_UGn_normal_base <- m11_UGn$`97.5%`
Rhat_muDPp_UGn_normal_base <- m11_UGn$Rhat
precDPp_mu_UGn_normal_base <- UB_mu_DPp_UGn_normal_base - LB_mu_DPp_UGn_normal_base


f17_UGn <- grepl("delta_new", row.names(DPpresults_U_G_n_normal_base))
m17_UGn <- DPpresults_U_G_n_normal_base[f17_UGn,]
delta_new_DPp_UGn_normal_base<- m17_UGn$`50%`
LB_delta_new_DPp_UGn_normal_base <- m17_UGn$`2.5%`
UB_delta_new_DPp_UGn_normal_base <- m17_UGn$`97.5%`
Rhat_delta_newDPp_UGn_normal_base <- m17_UGn$Rhat
precDPp_delta_new_UGn_normal_base <- UB_delta_new_DPp_UGn_normal_base - LB_delta_new_DPp_UGn_normal_base


f22_UGn <- grepl("var.true", row.names(DPpresults_U_G_n_normal_base))
m22_UGn <- DPpresults_U_G_n_normal_base[f22_UGn,]
tau2_DPp_UGn_normal_base <- m22_UGn$`50%`
LB_tau2_DPp_UGn_normal_base <- m22_UGn$`2.5%`
UB_tau2_DPp_UGn_normal_base <- m22_UGn$`97.5%`
Rhat_tauDPp_UGn_normal_base <- m22_UGn$Rhat
precDPp_tau2_UGn_normal_base <- UB_tau2_DPp_UGn_normal_base -  LB_tau2_DPp_UGn_normal_base


f33_UGn <- grepl("basemu", row.names(DPpresults_U_G_n_normal_base))
m33_UGn <- DPpresults_U_G_n_normal_base[f33_UGn,]
base_mu_DPp_UGn_normal_base <- m33_UGn$`50%`
LB_base_mu_DPp_UGn_normal_base <- m33_UGn$`2.5%`
UB_base_mu_DPp_UGn_normal_base <- m33_UGn$`97.5%`
Rhat_base_mu_DPp_UGn_normal_base <- m33_UGn$Rhat
prec_base_mu_DPp_UGn_normal_base <- UB_base_mu_DPp_UGn_normal_base - LB_base_mu_DPp_UGn_normal_base

f44_UGn <- grepl("basetau", row.names(DPpresults_U_G_n_normal_base))
m44_UGn <- DPpresults_U_G_n_normal_base[f44_UGn,]
base_tau_DPp_UGn_normal_base <- m44_UGn$`50%`
base_tau2_DPp_UGn_normal_base <- (m44_UGn$`50%`)^2
LB_base_tau2_DPp_UGn_normal_base <- (m44_UGn$`2.5%`)^2
UB_base_tau2_DPp_UGn_normal_base <- (m44_UGn$`97.5%`)^2
Rhat_base_tau_DPp_UGn_normal_base <- (m44_UGn$Rhat)^2
prec_base_tau2_DPp_UGn_normal_base <- UB_base_tau2_DPp_UGn_normal_base - LB_base_tau2_DPp_UGn_normal_base

f55_UGn <- grepl("alpha", row.names(DPpresults_U_G_n_normal_base))
m55_UGn <- DPpresults_U_G_n_normal_base[f55_UGn,]
alpha_DPp_UGn_normal_base <- m55_UGn$`50%`
LB_alpha_DPp_UGn_normal_base <- m55_UGn$`2.5%`
UB_alpha_DPp_UGn_normal_base <- m55_UGn$`97.5%`
Rhat_alpha_DPp_UGn_normal_base <- m55_UGn$Rhat
prec_alpha_DPp_UGn_normal_base <- UB_alpha_DPp_UGn_normal_base - LB_alpha_DPp_UGn_normal_base

fcl_UGn <- grepl("K", row.names(DPpresults_U_G_n_normal_base))   
mcl_UGn <- DPpresults_U_G_n_normal_base[fcl_UGn,]
median_K_DPp_UGn_normal_base <- mcl_UGn$`50%` 
LB_K_DPp_UGn_normal_base <- mcl_UGn$`2.5%`
UB_K_DPp_UGn_normal_base <- mcl_UGn$`97.5%`

fp_UGn <- grepl("p", row.names(DPpresults_U_G_n_normal_base))
mp_UGn <- DPpresults_U_G_n_normal_base[fp_UGn,]
mp_UGn <- mp_UGn[!grepl("pop", row.names(mp_UGn)),]
mp_UGn <- mp_UGn[!grepl("alpha", row.names(mp_UGn)),]
pi_DPp_UGn_normal_base <- mp_UGn$mean


listaDPp_UGn_normal_base <- cbind.data.frame(rel_effDPp_UGn_normal_base,sd_rel_effDPp_UGn_normal_base, LB_rel_effDPp_UGn_normal_base,
                                             UB_rel_effDPp_UGn_normal_base,Rhat_deltaDPp_UGn_normal_base)


numclus_DPp_UGn_normal_base <- unlist(median_K_DPp_UGn_normal_base)
LB_K_DPp_UGn_normal_base <- unlist(LB_K_DPp_UGn_normal_base)
UB_K_DPp_UGn_normal_base <- unlist(UB_K_DPp_UGn_normal_base)

resDPp_UGn_normal_base <- cbind.data.frame(base_mu_DPp_UGn_normal_base,LB_base_mu_DPp_UGn_normal_base,
                                           UB_base_mu_DPp_UGn_normal_base,base_tau_DPp_UGn_normal_base,base_tau2_DPp_UGn_normal_base,
                                           LB_base_tau2_DPp_UGn_normal_base,
                                           UB_base_tau2_DPp_UGn_normal_base,
                                           mu_DPp_UGn_normal_base,
                                           LB_mu_DPp_UGn_normal_base, 
                                           UB_mu_DPp_UGn_normal_base, 
                                           delta_new_DPp_UGn_normal_base,
                                           LB_delta_new_DPp_UGn_normal_base, 
                                           UB_delta_new_DPp_UGn_normal_base,
                                           tau2_DPp_UGn_normal_base, LB_tau2_DPp_UGn_normal_base, UB_tau2_DPp_UGn_normal_base,
                                           precDPp_mu_UGn_normal_base,precDPp_delta_new_UGn_normal_base,
                                           precDPp_tau2_UGn_normal_base, prec_base_mu_DPp_UGn_normal_base, prec_base_tau2_DPp_UGn_normal_base,
                                           prec_alpha_DPp_UGn_normal_base,  alpha_DPp_UGn_normal_base , LB_alpha_DPp_UGn_normal_base, UB_alpha_DPp_UGn_normal_base,
                                           Rhat_muDPp_UGn_normal_base,Rhat_delta_newDPp_UGn_normal_base,
                                           Rhat_tauDPp_UGn_normal_base, Rhat_base_mu_DPp_UGn_normal_base,
                                           Rhat_base_tau_DPp_UGn_normal_base, Rhat_alpha_DPp_UGn_normal_base, 
                                           numclus_DPp_UGn_normal_base, LB_K_DPp_UGn_normal_base, UB_K_DPp_UGn_normal_base,
                                           Pr_mu_DPp_UGn_normal_base, Pr_delta_new_DPp_UGn_normal_base)

extra_col = row.names(DPpresults_U_G_n_normal_base)
DPpresults_U_G_n_normal_base$extra_col = extra_col
prob_DPp_UGn_normal_base = round(max_values,2)


write.csv(resDPp_UGn_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\res_DPp_n_U_G_normal_base.csv",row.names=FALSE )  
write.csv(listaDPp_UGn_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\rel_eff_DPp_n_U_G_normal_base.csv",row.names=FALSE )  
write.csv(prob_DPp_UGn_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\max_prob_cluster_DPp_n_U_G_normal_base.csv",row.names=FALSE )  
write.csv(mean_cl, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\clusters_means_DPp_n_U_G_normal_base.csv",row.names=FALSE )  
write.csv(DPpresults_U_G_n_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\all_res_DPp_n_U_G_normal_base.csv",row.names=FALSE )  


########## TO CREATE BEATIFULL FOREST PLOT  ###########
max_values = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\max_prob.csv")  
max_indices = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\clustering_assignment.csv")  
es = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\res_DPp_n_U_G_normal_base.csv")  
rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\rel_eff_DPp_n_U_G_normal_base.csv")  
mean_cl = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\clusters_mean.csv")  


mean_cl1 = mean_cl$x[1]
mean_cl2 = mean_cl$x[2]
mean_cl3 = mean_cl$x[3]
mean_cl4 = mean_cl$x[4]

# Original vector
numbers <- max_indices
numbers <- unlist(max_indices$x)
clusters = numbers
cl_vec <- paste0('"', clusters, '"', collapse = ", ")
cl = cat(cl_vec)
clusters = c("1", "1", "2", "2", "1", "1", "1", "2", "1", "1", "2", "1", "2", "1", "1",
             "1", "1", "1", "2", "1", "2", "4", "1", "1", "1", "4", "1", "1", "1", "1", 
             "1", "1", "2", "4", "2", "1", "1", "4", "1", "4", "1", "4", "1", "8", "1",
             "2", "1", "1", "2", "1", "1", "1", "4", "2", "4", "1", "1", "1")
# Step 1: Convert expressions to numeric values
cluster_nums <- sapply(clusters, function(x) eval(parse(text = x)))

# Step 2: Define a named mapping
recode_map <- c("1" = 2, "2" = 1, "4" = 3, "8" = 4)

# Step 3: Apply recoding
cluster_recode <- sapply(cluster_nums, function(x) {
  # Use recoded value if it exists in the mapping, otherwise return original value
  recode_map_char <- as.character(x)
  if (recode_map_char %in% names(recode_map)) {
    recode_map[[recode_map_char]]
  } else {
    x
  }
})

# Result
cluster_recode

# Convert to a character vector (each element is a string)
cluster_recode_str_vec <- unname(as.character(cluster_recode))

# View the result
cluster_recode_str_vec

cat(cluster_recode_str_vec, sep = ", ")

# Format and print as quoted strings with commas
cat(paste0('"', cluster_recode_str_vec, '"', collapse = ", "))

clusters = c("2", "2", "1", "1", "2", "2", "2", "1", "2", "2", "1", "2", "1", "2", "2", "2", "2",
             "2", "1", "2", "1", "3", "2", "2", "2", "3", "2", "2", "2", "2", "2", "2", "1", "3", 
             "1", "2", "2", "3", "2", "3", "2", "3", "2", "4", "2", "1", "2", "2", "1", "2", "2",
             "2", "3", "1", "3", "2", "2", "2")

# Define a color palette
color_palette <- c( "black","#b71c1c","#0d47a1", "#737373"  )
# Add more colors if needed

# Assign colors based on the unique numbers in the vector
unique_numbers <- sort(unique(numbers))
color_map <- setNames(color_palette[seq_along(unique_numbers)], unique_numbers)

# Create the color vector
color_vector <- color_map[as.character(numbers)]

# View the color vector
color_vector



# Unique numbers in the original vector
unique_numbers <- unique(numbers)

# Replacement values (in order corresponding to the unique numbers)
replacement_values <- c(16, 15, 14, 12)  # Adjust the length if needed

# Create a mapping from original numbers to replacement numbers
replacement_map <- setNames(replacement_values, unique_numbers)

# Replace the original numbers with the new values
new_numbers <- replacement_map[as.character(numbers)]

# View the new vector
new_numbers

# ######## forest plot ########
pr <- unlist(max_values$x)
pr <- as.vector(round(pr, 2))

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
sorted_rel_eff <- sorted_data$rel_effDPp_UGn_normal_base
sorted_LB <- sorted_data$LB_rel_effDPp_UGn_normal_base
sorted_UB <- sorted_data$UB_rel_effDPp_UGn_normal_base
sorted_authors <- sorted_data$`pre_term_data_58$Study`
sorted_color <- sorted_data$color_vector
sorted_pch <- sorted_data$new_numbers
sorted_add_data <- sorted_data[, c(   "mean_ga", "mean_bw" ,"matched", "PMB",  "clusters", "pr")]

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

addpoly(x= es$mu_DPp_UGn_normal_base, ci.lb = es$LB_mu_DPp_UGn_normal_base ,
        ci.ub = es$UB_mu_DPp_UGn_normal_base , rows=-0.2)


arrows(x0 = es$LB_delta_new_DPp_UGn_normal_base,
       x1 = es$UB_delta_new_DPp_UGn_normal_base,
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
abline(h = 8.5, col = "gray60", lty = "dotted", lwd = 0.7)


# Add a legend
legend("topright", legend = c("Cluster 1", "Cluster 2","Cluster 3"), 
       col = c("black", "#b71c1c","#0d47a1" ), lty = 1, lwd = 2)

### SAVE THE FOREST PLOT IN A SPECIFIC FILE WITH 300 dpi
# Open TIFF device
tiff("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\DPMp_U_G_k.tiff", 
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

addpoly(x= es$mu_DPp_UGn_normal_base, ci.lb = es$LB_mu_DPp_UGn_normal_base ,
        ci.ub = es$UB_mu_DPp_UGn_normal_base , rows=-0.2)


arrows(x0 = es$LB_delta_new_DPp_UGn_normal_base,
       x1 = es$UB_delta_new_DPp_UGn_normal_base,
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
abline(h = 8.5, col = "gray60", lty = "dotted", lwd = 0.7)

dev.off()



