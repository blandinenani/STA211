
library(rjags)
library(coda)
library(ggmcmc)
library(mvtnorm)
library(tidyverse) # Ensemble de packages pour la manip des donnees
library(latex2exp) # permet d'ajouter des equations latex aux figures
rm(list=ls()) # Nettoyage de l'environnement de travail


######premières approches#########################
#################################################


qweibull <- function(u, kappa=1, alpha=1){
  return( ((-log(1 - u))/alpha)^(1/kappa) )
}

rweibull <- function(n, kappa=1, alpha=1, t0){
  U = runif( n )
  X = qweibull( U, kappa, alpha )
  X = sapply(X, function(x) min(t0, x))
  X=sort(X)
  p=which.max(X)
  X[n+1]=p
  return(X)
}


inverse_wei<-function(q,alpha,kappa){
  return ((-(1/alpha)*log(1-q))^(1/kappa))
}

##pour chaque kappa on cherche le max de alpha!!

##kappa fixe : 



log_vraisemblance_en_k<-function(kappa,t,p){
  return (p*(log(p/sum(t^kappa))+log(kappa))+(kappa-1)*sum(log(t[1:p]))-p)
}


max_vraisemblance_en_k<-function(t,p){
  sol<-optimize(log_vraisemblance_en_k,maximum=TRUE, interval=c(0.01,100),t=t,p=p)
  kappa_hat<-sol$maximum
  alpha_hat<-p/sum(t^kappa_hat)
  return (c(alpha_hat,kappa_hat))
}

test_weibull<-read.csv("donnees_Weibull_censure.csv")
x=test_weibull$x

##prenons t0=4 
t0=4
x = sapply(x, function(x) min(t0, x))
x=sort(x)
p=which.max(x)

max=max_vraisemblance_en_k(x,p)
alpha_hat<-max[1]
kappa_hat<-max[2]
cat("Estimation d'alpha pour les données tests :",alpha_hat)
cat("Estimation de kappa pour les données tests :", kappa_hat)


##QUESTION 7 :

n=20
alpha_re = numeric(10000)
kappa_re = numeric(10000)
q_60 = numeric(10000)
for (i in 1:10000){
  y<-rweibull(n,kappa_hat,alpha_hat,t0=4)
  ##Réestimation des paramètres
  max_k=max_vraisemblance_en_k(y[1:n],y[n+1])
  alpha_re[i]=max_k[1]
  kappa_re[i]=max_k[2]
  q_60[i]=qweibull(0.6,alpha_re[i],kappa_re[i])
}



# Calculer les intervalles de confiance à 95% pour α, κ et q60%
alpha_IC <- quantile(alpha_re, c(0.025, 0.975))
kappa_IC <- quantile(kappa_re, c(0.025, 0.975))
q60_IC <- quantile(q_60, c(0.025, 0.975))

# Afficher les intervalles de confiance obtenus
cat("Intervalle de confiance à 95% pour alpha : [", alpha_IC[1], ", ", alpha_IC[2], "]\n")
cat("Intervalle de confiance à 95% pour kappa : [", kappa_IC[1], ", ", kappa_IC[2], "]\n")
cat("Intervalle de confiance à 95% pour q60% : [", q60_IC[1], ", ", q60_IC[2], "]\n")


#on import à nouveau les données --> nouvelles notations
library(readr)
donnees <- read_csv("donnees_Weibull_censure (1).csv")
p=6
n=nrow(donnees)
t=donnees$x


###question 10: algo MCMC #######################################################
#################################################################################

#Algorithme de Gibbs permettant de simuler une chaîne de Markov pour chaque paramètre (alpha et kappa)

#Tout d'abbord on implemente une fonction qui va nous permettre de calculer
#la densité marginale à postériori de Kappa
#on aura besoin de cette fonction pour calculer le ratio de MH
#on la met sous forme log pour faciliter le calcul

log_kappa_p=function(kappa,t,p){dgamma(kappa,10^-3,10^-3,log=TRUE)+p*log(kappa)+(kappa-1)*sum(log(t[1:p]))-(10^-3+p)*log(10^-3+sum(t^kappa))}


#Ensuite on implément l'algorithme MCMC


MCMC=function(kappa_0,alpha_0,k,niter,p,t)
{
  chain_kappa= c(1:niter)
  chain_alpha=c(1:niter)
  chain_kappa[1]=kappa_0
  chain_alpha[1]=alpha_0
  accept=c(1:niter)
  kappa_curr=kappa_0
  alpha_curr=alpha_0
  
  for(i in 2:niter)
  {
    #ETAPE 1: actualisation de kappa par MH
    kappa_cand = rnorm(1, mean=kappa_curr, sd= k)
    log_ratio=log_kappa_p(kappa_cand,t,p)-log_kappa_p(kappa_curr,t,p)
    u=runif(1)
    #comme on est passé au log on compare à 0
    m=min(0,log_ratio)
    if(isTRUE(log(u)<m)) 
    {
      kappa_curr= kappa_cand
      accept[i]=1
      chain_kappa[i]=kappa_curr
    }
    else
    {
      accept[i]=0
      chain_kappa[i]=kappa_curr
    }
    
    #ETAPE 2: actualisation de alpha en utilisant le kappa actualisé
    alpha_curr=rgamma(1,10^-3+p,10^-3+sum(t^kappa_curr))
    chain_alpha[i]=alpha_curr
  }
  
  return(list(chain1=chain_alpha,chain2=chain_kappa,accept=accept))
}


###question 11: choix du saut k #######################################################
#################################################################################

#on va tester les taux d'acceptation opbtenus pour différentes valeurs de saut k

#on prépare les tests qu'on va mener
G <- 10000
list_saut <- seq(1, 301, 10)
list_accept <- sapply(list_saut, function(saut) {
  #on initialise chaque paramètre à 1
  Res <- suppressWarnings(MCMC(kappa_0 = 1, alpha_0 = 1, k = saut, niter = G, p = 6, t = donnees$x))
  mean(Res$accept)
})

plot(list_saut,list_accept,type="l",xlab = "k",ylab="taux d'acceptation",main="Taux d'acceptation en fonction du paramètre du saut k")
#plus on augmente le saut moins on accepte

#on détermine par lecture graphique la valeur correspondant à 
#un taux d'acceptation de 40%
#on choisie une valeur de saut égale à 4.1
abline(v=4.1,lty=2)

###question 12: execution des chaînes #######################################################
#############################################################################################

#on conserve la valeur de saut obtenue précédement
saut=4.1

#nous cherchons maintenant à déterminer le nombre d'iteration à élmininer qui correspondent 
#aux itération nécessaires avant stabilisation

G=20000
Res1=suppressWarnings(MCMC(kappa_0=2,alpha_0=2,k=saut,niter=G,p=p,t=donnees$x))
Res2=suppressWarnings(MCMC(kappa_0=3,alpha_0=1, k=saut,niter=G,p=p,t=donnees$x))
Res3=suppressWarnings(MCMC(kappa_0=1,alpha_0=5,k=saut,niter=G,p=p,t=donnees$x))
par(mfrow=c(1,3))
plot(Res1$chain1,type="l",ylim=c(-2,5))
lines(Res1$chain2,col="blue")
plot(Res2$chain1,type="l",ylim=c(-2,5))
lines(Res2$chain2,col="blue")
plot(Res3$chain1,type="l",ylim=c(-2,5))
lines(Res3$chain2,col="blue")

#analyse visuelle
#les chaînes semblent converger vers une même valeur

#on supperpose les chaînes
par(mfrow=c(1,2))
plot(Res1$chain1,type="l",main="Les Chaines alpha")
lines(Res2$chain1,col="blue")
lines(Res3$chain1,col="red")
plot(Res1$chain2,type="l",main="Les Chaines kappa")
lines(Res2$chain2,col="blue")
lines(Res3$chain2,col="red")
#elles se confondent après une période de chauffe
#DONC il n'y a aucun indice de non convergence

#statistique de Gelman-Rubin

library(ggmcmc)
library(coda)
sta1=mcmc(cbind(Res1$chain1,Res1$chain2))
sta2=mcmc(cbind(Res2$chain1,Res2$chain2))
sta3=mcmc(cbind(Res3$chain1,Res3$chain2))
sta=mcmc.list(list(sta1,sta2,sta3))
gelman.diag(sta)

# Potential scale reduction factors:
#   
#   Point est. Upper C.I.
# [1,]          1          1
# [2,]          1          1
# 
# Multivariate psrf
# 
# 1

#Règle standard : R < 1.05 + stabilité ⇒ Pas de problème de convergence majeur
#diagnostiqué selon ce critère

gelman.plot(sta)
abline(v=3000,lty=2)
#On voit graphiquement que 20000 itérations dont 3000 itérations de temps de chauffe
#suffisent pour approcher correctement la loi a posteriori 
#on a donc G_0=3000


###question 13: analyser des autocorrélations intra-chaînes #################################
#############################################################################################

#on cherche maintenant à determiner la taille de l'échantillon nécessaire à bien
#representer la loi à postériori

autocorr.diag(sta)
autocorr.plot(sta)

#Plus les chaînes de Markov sont auto-corrélées (problème de mélangeance
#= "slow mixing"), plus il faudra générer de valeurs après le temps de chauffe pour
#bien approcher la loi a posteriori (Calculs plus longs !)

#ici: on remarque qu'il y a des corrélation donc il faudra faire plus de tirage dans 
#la chaîne que de tirage indépendant que l'on souhaite


###question 14: détermination de la taille nécessaire de l'échantillon ######################
#############################################################################################

#on retire le temps de chauffe

G_0=3000
sta1_n=mcmc(cbind(Res1$chain1[G_0:G],Res1$chain2[G_0:G]))
sta2_n=mcmc(cbind(Res2$chain1[G_0:G],Res2$chain2[G_0:G]))
sta3_n=mcmc(cbind(Res3$chain1[G_0:G],Res3$chain2[G_0:G]))
sta_n=mcmc.list(list(sta1_n,sta2_n,sta3_n))

#calcul ESS

effectiveSize(sta_n)
# var1     var2 
# 7985.759 5705.583 
#les echantillons effectifs sont nettement plus faibles que les échantillons originaux 
#ce qui est logique car on a obsrervé des corrélations

#utilisation de l'erreur de MC pour choisir la taille de l'échantillon

summ<- summary(sta_n)
summ$statistics[,4]< 0.05*summ$statistics[,2]
#Vérifier que Time-series SE (= Erreur de Monte-Carlo estimée) est inférieure à 5% de SD
#(= Ecart-type empirique a posteriori) pour chaque paramètre
#ICI: l'écart est approprié on peut concerver cette taille d'échantillon

### question 15: statistques et représentations des lois des paramètres #####################
#############################################################################################

#statistiques

summary(sta_n)
# Iterations = 1:17001
# Thinning interval = 1 
# Number of chains = 3 
# Sample size per chain = 17001 
# 
# Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean      SD  Naive SE Time-series SE
# [1,] 0.08376 0.07868 0.0003484      0.0008814
# [2,] 1.83277 0.63336 0.0028045      0.0083916

#représentation

plot(sta_n,density=TRUE)

### question 15: intervalles de crédibilité #################################################
#############################################################################################

summary(sta_n)

# 2. Quantiles for each variable:
#   
#   2.5%     25%     50%    75%  97.5%
# [1,] 0.005394 0.02893 0.05986 0.1132 0.2956
# [2,] 0.775022 1.36898 1.77032 2.2313 3.2112

