###
##D_j
##gg animates
# check with the ak and how it corresponds to s values or groups
# Tune lambda or degf
# read_csv tidyverse setting
# pipeline to validate 
# compare f vs a 
# generate y_j's automatically so that we can scale things up
# compare PEP prior vs PEP post in a plot
# use ggplot!!!!
# PEP as prior
# color points by Assignment of correct and incorrect sampling


# treat all sigma the same as the truth
# check sigma's values
# change the wijk part of p(y)

# understanding of aks

# rescale with arclength problem

# ak vs x1 and x2

# browser in R etc debugging in R

# print statements for j(\theta) and l(\theta) the log likelihood


rm(list = ls(all.names = TRUE))
library(mvtnorm)
library(numDeriv)
library(Matrix)
library(ggplot2)
library(splines)
library(tidyverse)

PEP<-read_csv("ev_prefiltering.csv")
PEP<-PEP$PEP



n<-100
num_of_data <- n
p <- 2
num_it <- 60
degf <- 5
temp <- c(22,28)
power <- c(2.3,1.5) 


set.seed(40)
pep_sim<-matrix(sample(PEP,n*p,replace = TRUE), nrow = n, ncol = p)
#pep_sim<-matrix(0, nrow = n, ncol = p)
s <- runif(num_of_data,0, 1)
x<-matrix(nrow = n, ncol = p)
y<-matrix(nrow = n, ncol = p)
z<-matrix(nrow = n, ncol = p)
y1<-matrix(nrow = n, ncol = p)

#x is the reference RT
x[,1] <- 3*(temp[1]-10)^2*s^power[1]/100 + 10
x[,2] <- 3*(temp[2]-10)^2*s^power[2]/100 + 10


prior <- array(dim = c(n,p,n))
for (i in 1:n) {
  for (j in 1:p) {
    prob=rep(pep_sim[i,j]/(n-1),n)
    prob[i]=1-pep_sim[i,j]
    prior[i,j,]=prob
    class=rmultinom(1, 1, prob)
    y[i,j]=t(class)%*%x[,j]+rnorm(1,mean = 0,sd=0.25)
    z[i,j]=which(class==1)
  }
}

plot(x=y[,1],y=y[,2])
plot(x=x[,1],y=x[,2])


y1[,1] <- 3*(temp[1]-10)^2*sort(s)^power[1]/100 + 10-mean(y[,1])
y1[,2] <- 3*(temp[2]-10)^2*sort(s)^power[2]/100 + 10-mean(y[,2])
#s<-seq(-pi/2,pi/2,length.out=100)

# x[,1] <- 3*(temp[1]-10)^2*s^2/100 + 10
# x[,2] <- 3*(temp[2]-10)^2*s^2/100 + 10
# x[,1] <- 3*(temp[1]-10)^2*s^power[1]/100 + 10 + rnorm(num_of_data, mean = 0, sd = 0.25)
# x[,2] <- 3*(temp[2]-10)^2*s^power[2]/100 + 10 + rnorm(num_of_data, mean = 0, sd = 0.25)
# plot(x=x[,1],y=x[,2])
y <- scale(y, scale = F)
plot(y[, 1],
     y[, 2],
     bty = 'n',
     xlab = expression(x[1]),
     ylab = expression(x[2]))
legend("topleft", legend = c("Initialize"), bty = "n")




svdx <- svd(y)
d<-svdx$v
abline(a=0, b=d[2,1]/d[1,1])

## plot projections of each point onto line
proj1 <- y%*%d[,1]%*%t(d[,1])
segments(x0=y[,1],y0=y[,2],
         x1=proj1[,1],y1=proj1[,2])


ak_initial<-function(fi,mu,d){
  ak<-t(d)%*%(fi-mu)/(t(d)%*%d)
  return(ak)
}
n=nrow(y)
p=ncol(y)
vk<-rep(1/n,n)
mu<-colMeans(y)
Fun<-proj1#n by p matrix representing fj(ak)
a<-apply(Fun, 1, ak_initial, mu=mu,d=d[,1])
aini<-a
alist<-list()

sigma<-t(replicate(n,diag(cov(y))))#should be n by p (col for each ak)

Q<-function(a,fits,W,X){#f is a list of fits, is empty for first iteration
  n=length(a)
  p<-ncol(X)
  Fun<-matrix(data = NA,nrow = n,ncol = p)
  if(length(fits)==0){#first iteration
    # mu<-colMeans(X)
    # d<-svdx$v[,1]
    # Fun<-t(d%*%t(a)+mu)
    # Sigvec<-diag(cov(X))
    # Qval=0
    # for (k in 1:n) {
    #   S=(t(t(X)-Fun[k,]))^2
    #   second_term<-rep(0.5*sum(log(Sigvec)),n)
    #   Qval<-Qval+sum(S%*%(-1/Sigvec)-second_term)
    # }
  }else{
    for (j in 1:p) {
      Fun[,j]=predict(fits[[j]],a)$y
    }
    Qval<-0
    for (k in 1:n) {
      S=W[i,j,]*(t(t(X)-Fun[k,]))^2
      Sigvec=colSums(W[,,k]*S)/colSums(W[,,k])#p by 1 vector of sigma_j^2(ak)
      second_term<-sum(W[,,k]%*%log(Sigvec)/2)
      Qval<-Qval+sum(S%*%(-1/Sigvec))-second_term
    }
  }
  return(-Qval)
}


#step2 Iteration
#likelihood function
u=0
fits<-list()
SMtx<-list()
#PenMtx<-list()
# Qval<-rep(0,num_it+1)
# Qval[1]<--Q(a,fits,W=NA,x)#Return value of Q is the negative of Q value in EM
# Qval1<-Qval[1]
for (h in 1:num_it) {
  # Qval0<-Qval1
  W2<-array(dim = c(n,p,n))#modify the calculation
  for (i in 1:n) {
    for (j in 1:p) {
      gyij<-0
      for (k1 in 1:n) {
        gyij<-gyij+dnorm(y[i,j],Fun[k1,j],sd = sqrt(sigma[k1,j]))*prior[i,j,k1]
      }
      for (k in 1:n) {
        loggyijak<-log(dnorm(y[i,j],mean = Fun[k,j],sd = sqrt(sigma[k,j])))
        #gy(y_h)
        W2[i,j,k]<-loggyijak-log(gyij)+log(prior[i,j,k])
      }
    }
  }
  W2=exp(W2)
  # for (i in 1:n) {
  #   for (j in 1:p) {
  #     W2[i,j,]=W2[i,j,]/sum(W2[i,j,])
  #   }
  # }
  # if(h!=1){
  #   a<-nlm(f=Q, p=a, fits=fits, W=W2, X=y,iterlim = 1)$estimate
  # }
  #now update f and sigma and vk
   a2<-max(a)
   a1<-min(a)
  weights<-t(apply(W2, 2:3, sum))
  
  for (j in 1:p) {
    Dmtx<-diag(1/weights[,j])#it's the inverse of D
    ybar<-t(W2[,j,])%*%y[,j]
    yvalue<-Dmtx%*%ybar
    fit <- smooth.spline(x=a, y=yvalue, w=weights[,j], df=degf)##fixed df because as iteration goes on the fit will be interplotation.
    Fun[,j]=fitted(fit)
    fits[[j]]=fit
  }
  
################################################################################  
  for (k in 1:n) {
    S=(t(t(y)-Fun[k,]))^2
    sigma[k,]=colSums(W2[,,k]*S)/colSums(W2[,,k])#p by 1 vector of sigma_j^2(ak)
  }
  vk<-weights/(n*p)
  u<-u+1
  # plot data and the principal curve for a sequence of a's
  plot(y[,1], y[,2], bty='n',
       xlab=expression(x[1]),
       ylab=expression(x[2]))
  legend("topleft", legend=c("Proj 1"), bty="n")
  points(x=Fun[,1],y=Fun[,2],col="red")
  
  points(x=y[2,1],y=y[2,2],col="blue")
  # plot(x[,1], x[,3], bty='n',
  #      xlab=expression(x[1]),
  #      ylab=expression(x[3]))
  # legend("topleft", legend=c("Proj 2"), bty="n")
  # points(x=Fun[,1],y=Fun[,3],col="red")
  # 
  # 
  # plot(x[,2], x[,3], bty='n',
  #      xlab=expression(x[2]),
  #      ylab=expression(x[3]))
  # legend("topleft", legend=c("Proj 3"), bty="n")
  # seq_a <- seq(min(a),max(a),length.out=100)
  # #lines(predict(fits[[1]], seq_a)$y, predict(fits[[2]], seq_a)$y)
  # points(x=Fun[,2],y=Fun[,3],col="red")
  #lines(Fun[,1], Fun[,2])
  # 
  # 
  # 
  
  # Qval1<--Q(a,fits,W2,x)
  # Qval[h+1]<-Qval1
   arc_length<-sum(diag(dist(Fun)))
   a<-arc_length*(a-a1)/(a2-a1)
   alist[[h]]<-a
}

apply(W2, 1:2, which.max)
lines(y1[,1],y1[,2])
ranks <- matrix(nrow = n,ncol = p)
diffs <- matrix(nrow = n,ncol = p)
rdiffs <- matrix(nrow = n,ncol = p)
for (i in 1:n) {
  for (j in 1:p) {
    ranks[i,j]=rank(W2[i,j,])[z[i,j]]
    max.value<-max(W2[i,j,])
    diffs[i,j]=abs(W2[i,j,z[i,j]]-max.value)
    rdiffs[i,j]=abs(W2[i,j,z[i,j]]-max.value)/max.value
  }
}

mean(ranks) # 47.755 when n =50  93.5 when n =100
mean(diffs) # 0.160311 n=50   0.211 when n =100 significant because divided by 50
mean(rdiffs) #                0.294 when n= 100

colMeans(sigma) #sigma doesn't converge in this simplified model
lapply(alist, range)
