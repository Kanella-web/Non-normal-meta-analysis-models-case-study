pre_term_data_58 = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\pre_term_data58.csv")

### DP mixture of normal distributions model ###
library(R2jags)
set.seed(558)
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
      
     delta12[i] ~ dnorm(theta[Z[i]],tauinv[Z[i]])
    
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

  # Normal base distribution 
  
  for (k in 1:N) {
    
    theta[k] ~ dnorm(basemu, basetau1)  
    
    tau_comp[k]~dgamma(3,b)
    
    tauinv[k] = 1/ tau_comp_sqr[k]
    tau_comp_sqr[k] = tau_comp[k] * tau_comp[k] 
    
    
  }

   # Priors
   basemu~dnorm(0,0.0001)
   basetau1 <- 1/basetau_sqr
   basetau_sqr <- basetau*basetau
   basetau~dunif(0,10)
   
   b~dgamma(0.03,0.03)

   # concentration parameter prior
   alpha ~ dunif(0.3, 5)
  
  # Random effects distribution mean
  for (i in 1:N) {
    meancl[i] <- p[i] * theta[i]
  }
  overall_mean <- sum(meancl[])  # E[X]

  # Random effects distribution variance
  for (i in 1:N) {
    var1[i]<-p[i]* tau_comp_sqr[i]
    sq[i]<-theta[i]-overall_mean
    var2[i]<-p[i]*sq[i]*sq[i]   #### E[X2]
  }
  overall_variance<-sum(var1[])+sum(var2[]) ### E[X2] - E[X]2
 
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
  delta_new ~ dnorm(theta[pred],tauinv[pred]) ##Then drawing the new study's effect size from the random effects distribution.
  
}", file = "DPd_26_U_U_norm_mix_normal_base.txt")
modfile = 'DPd_26_U_U_norm_mix_normal_base.txt'

run.modelDPd26_U_U_norm_mix_normal_base = jags(
  dati1  <- list(ns = nrow(pre_term_data_58),
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
    "basemu", ## mu of the base Normal distribution
    "basetau", ## tau of the base Normal distribution
    "overall_mean",   ## overall mu
    "overall_variance", ## between-study variance
    "delta12", #### study-specific effects
    "K",       ### total number of clusters 
    "p",      ## weights of the process 
    "alpha",  ## concentration parameter
    "theta", ## mean of each cluster
    "tau_comp_sqr",  ## tau2 of each cluster
    "tau_comp",
    "Z",    ### cluster assignment 
    "SC" ,    ## the probability of each cluster assignment
    "delta_new",
    "pred"
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

DPdresults_U_U_26_norm_mix_normal_base <- as.data.frame(run.modelDPd26_U_U_norm_mix_normal_base$BUGSoutput$summary) 

DPd_results_U_U_26_norm_mix_norm_mix_norm_mix_normal_base <-round(DPdresults_U_U_26_norm_mix_normal_base, digits=2) 

View(DPd_results_U_U_26_norm_mix_norm_mix_norm_mix_normal_base)

# ##### show into how many clusters the study could be assigned in each iteration
# fdd1 <- grepl("Z", row.names(DPdresults_U_U_26_norm_mix_normal_base))
# mdd1 <- DPdresults_U_U_26_norm_mix_normal_base[fdd1,]
# clust <- mdd1$`50%`
# 
# clust1_DPd_U_U_26_norm_mix_normal_base = which(clust == 2)
# clust2_DPd_U_U_26_norm_mix_normal_base = which(clust == 4)


######## Pr(mu<0) #######
Pr_mu_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base = mean(run.modelDPd26_U_U_norm_mix_normal_base$BUGSoutput$sims.matrix[  ,"overall_mean"] < 0 )

######## Pr(mu_new<0) #######
Pr_delta_new_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base = mean(run.modelDPd26_U_U_norm_mix_normal_base$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )

# ########## CLUSTER ASSIGNEMENT ##########
# DPd_results_U_U_26_norm_mix_norm_mix_norm_mix_normal_base$ind <- row.names(DPd_results_U_U_26_norm_mix_norm_mix_norm_mix_normal_base)
# 
# ############# extraction of the parameters of interest ###########
f55sc_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base <- grepl("SC", row.names(DPd_results_U_U_26_norm_mix_norm_mix_norm_mix_normal_base))
m55sc_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base <- DPd_results_U_U_26_norm_mix_norm_mix_norm_mix_normal_base[f55sc_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base,]
SC_norm_mix_norm_mix_norm_mix_norm_mix_base <- m55sc_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base$mean
# 
# #### BASED ON  THE max(SC[i,j]) ####
# ##### Distribution of data points to clusters ###
# 
extra_col = rownames(m55sc_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base)
m55sc_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base$names = extra_col
# 
prob_list <- data.frame(Column10 = m55sc_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base[, 10], Column1 = m55sc_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base[, 1])
# 
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
max_values <- result$max_values #### the probabilities of cluster assignment
max_indices <- result$max_indices #### the cluster assignments


write.csv(max_values, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPn_U_U_26_normal_base_58\\max_prob_mix.csv",row.names=FALSE )  
write.csv(max_indices, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPn_U_U_26_normal_base_58\\clustering_assignment_mix.csv",row.names=FALSE )  

### clusters with data points ###
# Find unique indices
k <- unique(max_indices)

# Initialize the positions list
positions_DPd_U_U_26_norm_mix_normal_base <- vector("list", length(k))
names(positions_DPd_U_U_26_norm_mix_normal_base) <- k

# Fill the positions_DPd_U_U_26_norm_mix_normal_base list with indices
for(value in k) {
  positions_DPd_U_U_26_norm_mix_normal_base[[as.character(value)]] <- which(max_indices == value)
}

#### cluster with 33 elements ###
clust1_DPd_U_U_26_norm_mix_normal_base = positions_DPd_U_U_26_norm_mix_normal_base[["1"]]
clust2_DPd_U_U_26_norm_mix_normal_base = positions_DPd_U_U_26_norm_mix_normal_base[["2"]]

write.csv(clust1_DPd_U_U_26_norm_mix_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPn_U_U_26_normal_base_58\\cluster1_DPd_26_U_U_norm_mix_normal_base_pred.csv",row.names=FALSE )  
write.csv(clust2_DPd_U_U_26_norm_mix_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPn_U_U_26_normal_base_58\\cluster2_DPd_26_U_U_norm_mix_normal_base_pred.csv",row.names=FALSE )  

fdd <- grepl("delta12", row.names(DPdresults_U_U_26_norm_mix_normal_base))
mdd <- DPdresults_U_U_26_norm_mix_normal_base[fdd,]
rel_effDPd_U_U_26_norm_mix_normal_base <- mdd$`50%`
LB_rel_effDPd_U_U_26_norm_mix_normal_base <- mdd$`2.5%`
UB_rel_effDPd_U_U_26_norm_mix_normal_base <- mdd$`97.5%`
sd_rel_effDPd_U_U_26_norm_mix_normal_base <- mdd$sd
Rhat_deltaDPd_U_U_26_norm_mix_normal_base <- mdd$Rhat


cl1 = rel_effDPd_U_U_26_norm_mix_normal_base[clust1_DPd_U_U_26_norm_mix_normal_base]
mean_cl1 = mean(cl1)

cl2 = rel_effDPd_U_U_26_norm_mix_normal_base[clust2_DPd_U_U_26_norm_mix_normal_base]
mean_cl2 = mean(cl2)

mean_cl = cbind.data.frame(mean_cl1, mean_cl2)
write.csv(mean_cl, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPn_U_U_26_normal_base_58\\clusters_mean.csv",row.names=FALSE )  

#### PLOT THE DENSITY OF THE RELATIVE EFFECTS #######
plot(density(rel_effDPd_U_U_26_norm_mix_normal_base))

#### add all the cluster means
fth <- grepl("theta", row.names(DPdresults_U_U_26_norm_mix_normal_base))
mth <- DPdresults_U_U_26_norm_mix_normal_base[fth,]
clust_mean_DPd_U_U_26_norm_mix_normal_base <- mth$`50%`
abline(v= clust_mean_DPd_U_U_26_norm_mix_normal_base)


########### SAVE THE REST OF THE RESULTS ###########

f11_U_U_26 <- grepl("overall_mean", row.names(DPdresults_U_U_26_norm_mix_normal_base))
m11_U_U_26 <- DPdresults_U_U_26_norm_mix_normal_base[f11_U_U_26,]
mu_DPd_U_U_26_norm_mix_normal_base<- m11_U_U_26$`50%`
LB_mu_DPd_U_U_26_norm_mix_normal_base <- m11_U_U_26$`2.5%`
UB_mu_DPd_U_U_26_norm_mix_normal_base <- m11_U_U_26$`97.5%`
Rhat_muDPd_U_U_26_norm_mix_normal_base <- m11_U_U_26$Rhat
precDPd_mu_U_U_26_norm_mix_normal_base <- UB_mu_DPd_U_U_26_norm_mix_normal_base - LB_mu_DPd_U_U_26_norm_mix_normal_base


f17_U_U_26 <- grepl("delta_new", row.names(DPdresults_U_U_26_norm_mix_normal_base))
m17_U_U_26 <- DPdresults_U_U_26_norm_mix_normal_base[f17_U_U_26,]
delta_new_DPd_U_U_26_norm_mix_normal_base<- m17_U_U_26$`50%`
LB_delta_new_DPd_U_U_26_norm_mix_normal_base <- m17_U_U_26$`2.5%`
UB_delta_new_DPd_U_U_26_norm_mix_normal_base <- m17_U_U_26$`97.5%`
Rhat_delta_newDPd_U_U_26_norm_mix_normal_base <- m17_U_U_26$Rhat
precDPd_delta_new_U_U_26_norm_mix_normal_base <- UB_delta_new_DPd_U_U_26_norm_mix_normal_base - LB_delta_new_DPd_U_U_26_norm_mix_normal_base


f22_U_U_26 <- grepl("overall_variance", row.names(DPdresults_U_U_26_norm_mix_normal_base))
m22_U_U_26 <- DPdresults_U_U_26_norm_mix_normal_base[f22_U_U_26,]
tau2_DPd_U_U_26_norm_mix_normal_base <- m22_U_U_26$`50%`
LB_tau2_DPd_U_U_26_norm_mix_normal_base <- m22_U_U_26$`2.5%`
UB_tau2_DPd_U_U_26_norm_mix_normal_base <- m22_U_U_26$`97.5%`
Rhat_tauDPd_U_U_26_norm_mix_normal_base <- m22_U_U_26$Rhat
precDPd_tau2_U_U_26_norm_mix_normal_base <- UB_tau2_DPd_U_U_26_norm_mix_normal_base -  LB_tau2_DPd_U_U_26_norm_mix_normal_base


f33_U_U_26 <- grepl("basemu", row.names(DPdresults_U_U_26_norm_mix_normal_base))
m33_U_U_26 <- DPdresults_U_U_26_norm_mix_normal_base[f33_U_U_26,]
base_mu_DPd_U_U_26_norm_mix_normal_base <- m33_U_U_26$`50%`
LB_base_mu_DPd_U_U_26_norm_mix_normal_base <- m33_U_U_26$`2.5%`
UB_base_mu_DPd_U_U_26_norm_mix_normal_base <- m33_U_U_26$`97.5%`
Rhat_base_mu_DPd_U_U_26_norm_mix_normal_base <- m33_U_U_26$Rhat
prec_base_mu_DPd_U_U_26_norm_mix_normal_base <- UB_base_mu_DPd_U_U_26_norm_mix_normal_base - LB_base_mu_DPd_U_U_26_norm_mix_normal_base

f44_U_U_26 <- grepl("basetau", row.names(DPdresults_U_U_26_norm_mix_normal_base))
m44_U_U_26 <- DPdresults_U_U_26_norm_mix_normal_base[f44_U_U_26,]
base_tau_DPd_U_U_26_norm_mix_normal_base <- m44_U_U_26$`50%`
base_tau2_DPd_U_U_26_norm_mix_normal_base <- (m44_U_U_26$`50%`)^2
LB_base_tau2_DPd_U_U_26_norm_mix_normal_base <- (m44_U_U_26$`2.5%`)^2
UB_base_tau2_DPd_U_U_26_norm_mix_normal_base <- (m44_U_U_26$`97.5%`)^2
Rhat_base_tau_DPd_U_U_26_norm_mix_normal_base <- (m44_U_U_26$Rhat)^2
prec_base_tau2_DPd_U_U_26_norm_mix_normal_base <- UB_base_tau2_DPd_U_U_26_norm_mix_normal_base - LB_base_tau2_DPd_U_U_26_norm_mix_normal_base

f55_U_U_26 <- grepl("alpha", row.names(DPdresults_U_U_26_norm_mix_normal_base))
m55_U_U_26 <- DPdresults_U_U_26_norm_mix_normal_base[f55_U_U_26,]
alpha_DPd_U_U_26_norm_mix_normal_base <- m55_U_U_26$`50%`
LB_alpha_DPd_U_U_26_norm_mix_normal_base <- m55_U_U_26$`2.5%`
UB_alpha_DPd_U_U_26_norm_mix_normal_base <- m55_U_U_26$`97.5%`
Rhat_alpha_DPd_U_U_26_norm_mix_normal_base <- m55_U_U_26$Rhat
prec_alpha_DPd_U_U_26_norm_mix_normal_base <- UB_alpha_DPd_U_U_26_norm_mix_normal_base - LB_alpha_DPd_U_U_26_norm_mix_normal_base

fcl_U_U_26 <- grepl("K", row.names(DPdresults_U_U_26_norm_mix_normal_base))   
mcl_U_U_26 <- DPdresults_U_U_26_norm_mix_normal_base[fcl_U_U_26,]
median_K_DPd_U_U_26_norm_mix_normal_base <- mcl_U_U_26$`50%` 
LB_K_DPd_U_U_26_norm_mix_normal_base <- mcl_U_U_26$`2.5%`
UB_K_DPd_U_U_26_norm_mix_normal_base <- mcl_U_U_26$`97.5%`

fp_U_U_26 <- grepl("p", row.names(DPdresults_U_U_26_norm_mix_normal_base))
mp_U_U_26 <- DPdresults_U_U_26_norm_mix_normal_base[fp_U_U_26,]
mp_U_U_26 <- mp_U_U_26[!grepl("pop", row.names(mp_U_U_26)),]
mp_U_U_26 <- mp_U_U_26[!grepl("alpha", row.names(mp_U_U_26)),]
pi_DPd_U_U_26_norm_mix_normal_base <- mp_U_U_26$mean


listaDPd_U_U_26_norm_mix_normal_base <- cbind.data.frame(rel_effDPd_U_U_26_norm_mix_normal_base,sd_rel_effDPd_U_U_26_norm_mix_normal_base, LB_rel_effDPd_U_U_26_norm_mix_normal_base,
                                                        UB_rel_effDPd_U_U_26_norm_mix_normal_base,Rhat_deltaDPd_U_U_26_norm_mix_normal_base)


numclus_DPd_U_U_26_norm_mix_normal_base <- unlist(median_K_DPd_U_U_26_norm_mix_normal_base)
LB_K_DPd_U_U_26_norm_mix_normal_base <- unlist(LB_K_DPd_U_U_26_norm_mix_normal_base)
UB_K_DPd_U_U_26_norm_mix_normal_base <- unlist(UB_K_DPd_U_U_26_norm_mix_normal_base)

resDPd_U_U_26_norm_mix_normal_base <- cbind.data.frame(base_mu_DPd_U_U_26_norm_mix_normal_base,LB_base_mu_DPd_U_U_26_norm_mix_normal_base,
                                                      UB_base_mu_DPd_U_U_26_norm_mix_normal_base,base_tau_DPd_U_U_26_norm_mix_normal_base,base_tau2_DPd_U_U_26_norm_mix_normal_base,
                                                      LB_base_tau2_DPd_U_U_26_norm_mix_normal_base,
                                                      UB_base_tau2_DPd_U_U_26_norm_mix_normal_base,
                                                      mu_DPd_U_U_26_norm_mix_normal_base,
                                                      LB_mu_DPd_U_U_26_norm_mix_normal_base, 
                                                      UB_mu_DPd_U_U_26_norm_mix_normal_base,
                                                      delta_new_DPd_U_U_26_norm_mix_normal_base,LB_delta_new_DPd_U_U_26_norm_mix_normal_base,UB_delta_new_DPd_U_U_26_norm_mix_normal_base, 
                                                      Rhat_delta_newDPd_U_U_26_norm_mix_normal_base, precDPd_delta_new_U_U_26_norm_mix_normal_base,
                                                      
                                                      tau2_DPd_U_U_26_norm_mix_normal_base, LB_tau2_DPd_U_U_26_norm_mix_normal_base, UB_tau2_DPd_U_U_26_norm_mix_normal_base,
                                                      precDPd_mu_U_U_26_norm_mix_normal_base, precDPd_tau2_U_U_26_norm_mix_normal_base, prec_base_mu_DPd_U_U_26_norm_mix_normal_base, prec_base_tau2_DPd_U_U_26_norm_mix_normal_base,
                                                      prec_alpha_DPd_U_U_26_norm_mix_normal_base,  alpha_DPd_U_U_26_norm_mix_normal_base , LB_alpha_DPd_U_U_26_norm_mix_normal_base, UB_alpha_DPd_U_U_26_norm_mix_normal_base,
                                                      Rhat_muDPd_U_U_26_norm_mix_normal_base, Rhat_tauDPd_U_U_26_norm_mix_normal_base, Rhat_base_mu_DPd_U_U_26_norm_mix_normal_base,
                                                      Rhat_base_tau_DPd_U_U_26_norm_mix_normal_base, Rhat_alpha_DPd_U_U_26_norm_mix_normal_base, 
                                                      numclus_DPd_U_U_26_norm_mix_normal_base, LB_K_DPd_U_U_26_norm_mix_normal_base, UB_K_DPd_U_U_26_norm_mix_normal_base,
                                                      Pr_delta_new_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base,
                                                      Pr_mu_DPd_U_U_26_norm_mix_norm_mix_norm_mix_normal_base)

extra_col = row.names(DPdresults_U_U_26_norm_mix_normal_base)
DPdresults_U_U_26_norm_mix_normal_base$extra_col = extra_col

write.csv(resDPd_U_U_26_norm_mix_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPn_U_U_26_normal_base_58\\res_DPd_26_U_U_norm_mix_normal_base_same_tau_comp_pred.csv",row.names=FALSE )  
write.csv(listaDPd_U_U_26_norm_mix_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPn_U_U_26_normal_base_58\\rel_eff_DPd_26_U_U_norm_mix_normal_base_pred.csv",row.names=FALSE )  
write.csv(clust_mean_DPd_U_U_26_norm_mix_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPn_U_U_26_normal_base_58\\clusters_means_DPd_26_U_U_norm_mix_normal_base_pred.csv",row.names=FALSE )  
write.csv(DPdresults_U_U_26_norm_mix_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPn_U_U_26_normal_base_58\\all_res_DPd_26_U_U_norm_mix_normal_base_pred.csv",row.names=FALSE )  


