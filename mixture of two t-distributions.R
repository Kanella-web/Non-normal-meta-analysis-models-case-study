
pre_term_data_58 = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\pre_term_data58.csv")

############ Mixture of two t-distributions model ############
set.seed(2158)
##normal-2t-mixture(Gamma/dir)
library(R2jags)
writeLines("
  model {
    for (i in 1:ns) { 
      
      m[i] ~ dnorm(0,.0001)
      
      # Mean outcome for group 1
      y1[i] ~ dnorm(mu1[i], prec1[i])
      mu1[i] <- m[i]
      
      # Mean outcome for group 2
      y2[i] ~ dnorm(mu2[i], prec2[i])
      mu2[i] <- m[i] + delta12[i]*pooled_sd[i]
      
      # Between-study standardized effect
      delta12[i] ~  dt(mu[z[i]], tau1[z[i]], df[z[i]])  # non-central-t distribution
      
      z[i] ~ dcat(p[])  # Mixture component indicator
      
      # Precision parameters for group 1 and 2 based on observed standard deviations
      prec1[i] <- 1 / sd1[i]*sd1[i]
      prec2[i] <- 1 / sd2[i]*sd2[i]

      # Calculate pooled standard deviation
      pooled_sd[i] <- sqrt(((n1[i] - 1)*pow(sd1[i], 2) + (n2[i] - 1)*pow(sd2[i], 2)) / (n1[i] + n2[i] - 2))
    
    }

    # Priors for mixing proportions
    for (i in 1:N) {

        alpha[i] <- 1  
        
    }
    ###### Assigning a DP prior for p reasures that sum[p]=1 ###
    
    p ~ ddirch(alpha[])


#### Priors
    for (k in 1:N) {
    
      #mu[k] ~ dnorm(0, t1)
      
      mu[k] ~ dnorm(m1, t1)  ##### variation of mu across components
      
      tau1[k] <- 1/ tau_k_sqr[k]
      
      tau_k_sqr[k] = tau_k[k] * tau_k[k]
      
      tau_k[k] ~ dgamma(3,b) ## variation of tau across components
      
      df[k] ~ dexp(0.10)T(2.5, )   

    }
    
    
    m1 ~ dnorm(0,0.0001)
     
    t1 = 1/ t*t
    t ~ dnorm(0, 1)T(0,)
    
    b~dgamma(0.03, 0.03)
    
    # Calculate mean of the mixture distribution
    for (i in 1:N) {
    
      mean_mu[i] <- p[i] * mu[i]
      
    }
    
    total_mean <- sum(mean_mu[])  # Total mean of the mixture
    
 
   # Random effects distribution variance #
  
  
  for(i in 1:N){
  
    var1[i]<-p[i]*tau_k_sqr[i]
  
    sq[i]<-mu[i]-total_mean
    
    var2[i]<-p[i]*sq[i]*sq[i]   #### E[X2]
  }
  
  total_var <-sum(var1[])+sum(var2[]) ### E[X2] - E[X]2
   
   
  # Summary statistics #
  for(i in 1:ns){
   for (j in 1:N){
    SC[i,j] <- equals(j,z[i])
   }
  }
  # total clusters K#
  for (j in 1:N){
   cl[j] <-step(sum(SC[,j])-1)
   }
  K<-sum(cl[])
  
  for(i in 1:ns){
   rho[i] <- (exp(delta12[i]) / (1+exp(delta12[i])))
  }
  
  for(i in 1:ns){
   for(j in 1:ns){
    equalsmatrix[i,j]<-equals(rho[i],rho[j])
   }
  equalsres[i]<-sum(equalsmatrix[i,])
  }
  
  pred ~ dcat(p[1:N])
   delta_new ~ dnorm(mu[pred],tau1[pred])
  
}",
           con = "model.finite_mixture_t_dir_weights_tau_comp.txt")

modfile3 = 'model.finite_mixture_t_dir_weights_tau_comp.txt'

library(R2jags)

dati1  <-list(ns = nrow(pre_term_data_58),
              y1 = pre_term_data_58$mean_FT,
              y2 = pre_term_data_58$mean_EPT.VPT,
              sd1 = pre_term_data_58$sd_FT,
              sd2 = pre_term_data_58$sd_EPT.VPT,
              n1 = pre_term_data_58$n_FT,
              n2 = pre_term_data_58$n_EPT.VPT,
              N = 2)

run.model_fin_mix_t_dir_weigths_tau_comp = jags(
  data = dati1,
  inits = NULL,
  parameters.to.save = c(
    "mu",
    "tau1",
    "delta12",
    "tau_k_sqr",
    "tau_k",
    "p",
    "total_mean",
    "total_var",
    "SC",
    "df",
    "z",
    "K",
    "pred",
    "delta_new"
  ),
  n.chains = 2,
  n.iter = 50000,
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile3
)

rel_eff_fin_mix_t_dw = as.data.frame(run.model_fin_mix_t_dir_weigths_tau_comp$BUGSoutput$summary)
rel_eff_fin_mix_t_dw1 = round(rel_eff_fin_mix_t_dw, 2)
View(rel_eff_fin_mix_t_dw1)

write.csv(rel_eff_fin_mix_t_dw1, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_t_mix_58\\dir_weights\\gamma_comp_var\\res_all.csv",row.names=FALSE )  

######## cluster assignments ######
fdd <- grepl("z", row.names(rel_eff_fin_mix_t_dw))
mdd <- rel_eff_fin_mix_t_dw[fdd,]
##### only study 49 is identified as an outlier #########
clust <- round(mdd$`50%`, 0)

cluster_mix_t_1 = which(clust == 1)
cluster_mix_t_2 = which(clust == 2)

write.csv(cluster_mix_t_1, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_t_mix_58\\dir_weights\\gamma_comp_var\\cluster1_assign_fin_mix_t_dw.csv",row.names=FALSE )  
write.csv(cluster_mix_t_2, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_t_mix_58\\dir_weights\\gamma_comp_var\\cluster2_assign_fin_mix_t_dw.csv",row.names=FALSE )  

### Pr of mu<0 
Pr_mu_fin_mix_t_dir_weigths = mean(run.model_fin_mix_t_dir_weigths_tau_comp$BUGSoutput$sims.matrix[  ,"total_mean"] < 0 )

### Pr of mu<0 
Pr_delta_new_fin_mix_t_dir_weigths = mean(run.model_fin_mix_t_dir_weigths_tau_comp$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )



########## CLUSTER ASSIGNEMENT ##########
rel_eff_fin_mix_t_dw$ind <- row.names(rel_eff_fin_mix_t_dw)

############# extraction of the parameters of interest ###########
f55sc_t_mix <- grepl("SC", row.names(rel_eff_fin_mix_t_dw))
m55sc_t_mix <- rel_eff_fin_mix_t_dw[f55sc_t_mix,]
SC_t_mix <- m55sc_t_mix$mean

#### BASED ON  THE max(SC[i,j]) ####
##### Distribution of data points to clusters ###

extra_col = rownames(m55sc_t_mix)
m55sc_t_mix$names = extra_col

prob_list <- data.frame(Column10 = m55sc_t_mix[, 10], Column1 = m55sc_t_mix[, 1])

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

write.csv(max_values, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_t_mix_58\\dir_weights\\gamma_comp_var\\cluster_prob_pred.csv",row.names=FALSE )  
write.csv(max_indices, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_t_mix_58\\dir_weights\\gamma_comp_var\\cluster_assignment.csv",row.names=FALSE )  

###### study-specific effects
fdd <- grepl("delta12", row.names(rel_eff_fin_mix_t_dw))
mdd <- rel_eff_fin_mix_t_dw[fdd,]
rel_eff_fin_mix_t_dw_dcv <- mdd$`50%`
LB_rel_eff_fin_mix_t_dw_dcv <- mdd$`2.5%`
UB_rel_eff_fin_mix_t_dw_dcv <- mdd$`97.5%`
sd_rel_eff_fin_mix_t_dw_dcv <- mdd$sd
Rhat_rel_eff_fin_mix_t_dw_dcv <- mdd$Rhat

mean_cl1 = mean(rel_eff_fin_mix_t_dw_dcv[cluster_mix_t_1])
mean_cl2 = mean(rel_eff_fin_mix_t_dw_dcv[cluster_mix_t_2])
mean_cl = cbind.data.frame(mean_cl1 , mean_cl2)

write.csv(mean_cl, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_t_mix_58\\dir_weights\\gamma_comp_var\\mean_cl.csv",row.names=FALSE )  

# Calculate densities for both clusters 
density1 <- density(rel_eff_fin_mix_t_dw_dcv[cluster_mix_t_1])
density2 <- density(rel_eff_fin_mix_t_dw_dcv[cluster_mix_t_2])

# Plot the first density
plot(density1, col = "red", lwd = 2, main = "Density Plot of Two Clusters",
     xlab = "Value", ylab = "Density", xlim= c(-2,0), ylim = c(0,5))

# Add the second density
lines(density2, col = "blue", lwd = 2)



plot(density(rel_eff_fin_mix_t_dw_dcv))
library(metafor)

m_a_fin_mix_t_dw = rma(yi = rel_eff_fin_mix_t_dw_dcv , vi = (sd_rel_eff_fin_mix_t_dw_dcv)^2, 
                       method = "REML",slab = pre_term_data_58$Study)
forest(m_a_fin_mix_t_dw,cex = 0.5)

#### Mean of each component
f33_fin_mix_t_dir_weigths <- grepl("mu", row.names(rel_eff_fin_mix_t_dw))
m33_fin_mix_t_dir_weigths <- rel_eff_fin_mix_t_dw[f33_fin_mix_t_dir_weigths,]
cl_mu_fin_mix_t_dw <- m33_fin_mix_t_dir_weigths$`50%`
LB_cl_mu_fin_mix_t_dw <- m33_fin_mix_t_dir_weigths$`2.5%`
UB_cl_mu_fin_mix_t_dw <- m33_fin_mix_t_dir_weigths$`97.5%`
Rhat_cl_mu_fin_mix_t_dw <- m33_fin_mix_t_dir_weigths$Rhat
prec_cl_mu_fin_mix_t_dw <- UB_cl_mu_fin_mix_t_dw - LB_cl_mu_fin_mix_t_dw

#### Variance of each component
f44_fin_mix_t_dir_weigths <- grepl("tau_k_sqr", row.names(rel_eff_fin_mix_t_dw))
m44_fin_mix_t_dir_weigths <- rel_eff_fin_mix_t_dw[f44_fin_mix_t_dir_weigths,]
m44_fin_mix_t_dir_weigths <- m44_fin_mix_t_dir_weigths[!grepl("tau1", row.names(m44_fin_mix_t_dir_weigths)),]
cl_tau2_fin_mix_t_dw <- m44_fin_mix_t_dir_weigths$`50%`
LB_cl_tau2_fin_mix_t_dw <- m44_fin_mix_t_dir_weigths$`2.5%`
UB_cl_tau2_fin_mix_t_dw <- m44_fin_mix_t_dir_weigths$`97.5%`
Rhat_cl_tau_fin_mix_t_dw <- m44_fin_mix_t_dir_weigths$Rhat
prec_cl_tau2_fin_mix_t_dw <- UB_cl_tau2_fin_mix_t_dw - LB_cl_tau2_fin_mix_t_dw


###### TO VISUALIZE THE ONLY OCCUPIED CLUSTER2 ##########
c1_fin_mix_t_dw = rma(yi = rel_eff_fin_mix_t_dw_dcv[cluster_mix_t_1] , vi = (sd_rel_eff_fin_mix_t_dw_dcv[cluster_mix_t_1])^2, 
                      method = "REML", slab = pre_term_data_58$Study[cluster_mix_t_1])
##### DENSITY PLOT #######
plot(density(rel_eff_fin_mix_t_dw_dcv[cluster_mix_t_1]))
##### FOREST PLOT #######
forest(c1_fin_mix_t_dw, cex = 0.5)

###### TO VISUALIZE THE ONLY OCCUPIED CLUSTER2 ##########
c2_fin_mix_t_dw = rma(yi = rel_eff_fin_mix_t_dw_dcv[cluster_mix_t_2] , vi = (sd_rel_eff_fin_mix_t_dw_dcv[cluster_mix_t_2])^2, 
                      method = "REML", slab = pre_term_data_58$Study[cluster_mix_t_2])
##### DENSITY PLOT #######
plot(density(rel_eff_fin_mix_t_dw_dcv[cluster_mix_t_2]))
##### FOREST PLOT #######
forest(c2_fin_mix_t_dw, cex = 0.5)
########### SAVE THE REST OF THE RESULTS ###########

f11_fin_mix_t_dir_weigths <- grepl("total_mean", row.names(rel_eff_fin_mix_t_dw))
m11_fin_mix_t_dir_weigths <- rel_eff_fin_mix_t_dw[f11_fin_mix_t_dir_weigths,]
mu_fin_mix_t_dw<- m11_fin_mix_t_dir_weigths$`50%`
LB_mu_fin_mix_t_dw <- m11_fin_mix_t_dir_weigths$`2.5%`
UB_mu_fin_mix_t_dw <- m11_fin_mix_t_dir_weigths$`97.5%`
Rhat_mu_fin_mix_t_dw <- m11_fin_mix_t_dir_weigths$Rhat
prec_mu_fin_mix_t_dw <- UB_mu_fin_mix_t_dw - LB_mu_fin_mix_t_dw


f17_fin_mix_t_dir_weigths <- grepl("delta_new", row.names(rel_eff_fin_mix_t_dw))
m17_fin_mix_t_dir_weigths <- rel_eff_fin_mix_t_dw[f17_fin_mix_t_dir_weigths,]
delta_new_fin_mix_t_dw<- m17_fin_mix_t_dir_weigths$`50%`
LB_delta_new_fin_mix_t_dw <- m17_fin_mix_t_dir_weigths$`2.5%`
UB_delta_new_fin_mix_t_dw <- m17_fin_mix_t_dir_weigths$`97.5%`
Rhat_delta_new_fin_mix_t_dw <- m17_fin_mix_t_dir_weigths$Rhat
prec_delta_new_fin_mix_t_dw <- UB_delta_new_fin_mix_t_dw - LB_delta_new_fin_mix_t_dw


f22_fin_mix_t_dir_weigths <- grepl("total_var", row.names(rel_eff_fin_mix_t_dw))
m22_fin_mix_t_dir_weigths <- rel_eff_fin_mix_t_dw[f22_fin_mix_t_dir_weigths,]
tau2_fin_mix_t_dw <- m22_fin_mix_t_dir_weigths$`50%`
LB_tau2_fin_mix_t_dw <- m22_fin_mix_t_dir_weigths$`2.5%`
UB_tau2_fin_mix_t_dw <- m22_fin_mix_t_dir_weigths$`97.5%`
Rhat_tau_fin_mix_t_dw <- m22_fin_mix_t_dir_weigths$Rhat
prec_tau2_fin_mix_t_dw <- UB_tau2_fin_mix_t_dw -  LB_tau2_fin_mix_t_dw

fp_fin_mix_t_dir_weigths <- grepl("p", row.names(rel_eff_fin_mix_t_dw))
mp_fin_mix_t_dir_weigths <- rel_eff_fin_mix_t_dw[fp_fin_mix_t_dir_weigths,]
mp_fin_mix_t_dir_weigths <- mp_fin_mix_t_dir_weigths[!grepl("pop", row.names(mp_fin_mix_t_dir_weigths)),]
mp_fin_mix_t_dir_weigths <- mp_fin_mix_t_dir_weigths[!grepl("pred", row.names(mp_fin_mix_t_dir_weigths)),]
pi_fin_mix_t_dir_weigths <- mp_fin_mix_t_dir_weigths$mean


lista_fin_mix_t_dir_weigths <- cbind.data.frame(rel_eff_fin_mix_t_dw_dcv,sd_rel_eff_fin_mix_t_dw_dcv,
                                                LB_rel_eff_fin_mix_t_dw_dcv,
                                                UB_rel_eff_fin_mix_t_dw_dcv,Rhat_rel_eff_fin_mix_t_dw_dcv)

############ Overall results ##########
res_fin_mix_t_dw <- cbind.data.frame(mu_fin_mix_t_dw, LB_mu_fin_mix_t_dw,UB_mu_fin_mix_t_dw, 
                                     delta_new_fin_mix_t_dw, LB_delta_new_fin_mix_t_dw,UB_delta_new_fin_mix_t_dw,
                                     tau2_fin_mix_t_dw, LB_tau2_fin_mix_t_dw, UB_tau2_fin_mix_t_dw,
                                     prec_mu_fin_mix_t_dw,  prec_delta_new_fin_mix_t_dw,
                                     prec_tau2_fin_mix_t_dw, 
                                     Rhat_mu_fin_mix_t_dw,Rhat_delta_new_fin_mix_t_dw, Rhat_tau_fin_mix_t_dw, 
                                     Pr_mu_fin_mix_t_dir_weigths,Pr_delta_new_fin_mix_t_dir_weigths )

##### Clusters results ##########
clusters_fin_mix_t_dw = cbind.data.frame(cl_mu_fin_mix_t_dw, LB_cl_mu_fin_mix_t_dw, UB_cl_mu_fin_mix_t_dw,
                                         cl_tau2_fin_mix_t_dw, LB_cl_tau2_fin_mix_t_dw, UB_cl_tau2_fin_mix_t_dw,
                                         prec_cl_mu_fin_mix_t_dw, prec_cl_tau2_fin_mix_t_dw,
                                         pi_fin_mix_t_dir_weigths) 

extra_col_cl = row.names(clusters_fin_mix_t_dw)
clusters_fin_mix_t_dw$extra_col = extra_col_cl


write.csv(res_fin_mix_t_dw, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_t_mix_58\\dir_weights\\gamma_comp_var\\res_fin_mix_t_dw.csv",row.names=FALSE )  
write.csv(clusters_fin_mix_t_dw, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_t_mix_58\\dir_weights\\gamma_comp_var\\clusters_fin_mix_t_dw.csv",row.names=FALSE )  
write.csv(lista_fin_mix_t_dir_weigths, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_t_mix_58\\dir_weights\\gamma_comp_var\\rel_eff_fin_mix_t_dw.csv",row.names=FALSE )  


