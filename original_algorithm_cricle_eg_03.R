rm(list = ls(all.names = TRUE))
library(mvtnorm)
library(numDeriv)
library(Matrix)
library(ggplot2)
library(splines)
#how to improve the calculation of Derivative of Q or the Newton Method
#x3 should be related to s
#How do I know the newton method works


#step0 generate some bivariate data
num_of_data <- 100
num_it <- 20
degf <- 10
set.seed(24)
# s <- runif(num_of_data,-pi / 2, pi / 2)
s<-seq(-pi/2,pi/2,length.out=num_of_data)
epsilon1 <- rnorm(num_of_data)
epslion2 <- rnorm(num_of_data)



x1 <- 3 * sin(s) + 0.5 * epsilon1
x2 <- 3 * cos(s) + 0.5 * epslion2
#x3<-runif(num_of_data,-1,1)
z1 <- 3 * sin(sort(s))-mean(x1)
z2 <- 3 * cos(sort(s))-mean(x2)



x <- scale(cbind(x1, x2), scale = F)




plot(x[, 1],
     x[, 2],
     bty = 'n',
     xlab = expression(x[1]),
     ylab = expression(x[2]))
legend("topleft", legend = c("Initialize"), bty = "n")



#step1 initiation (WARNING: This is scaled data, can not use svd when it's not scaled)
svdx <- svd(x)
d<-svdx$v
clip(alim_[1],alim_[2],alim_[1],alim_[2])
abline(a=0, b=d[2,1]/d[1,1])

## plot projections of each point onto line
proj1 <- x%*%d[,1]%*%t(d[,1])
segments(x0=x[,1],y0=x[,2],
         x1=proj1[,1],y1=proj1[,2])


ak_initial<-function(fi,mu,d){
  ak<-t(d)%*%(fi-mu)/(t(d)%*%d)
  return(ak)
}
n=nrow(x)
p=ncol(x)
vk<-rep(1/n,n)
mu<-colMeans(x)
Fun<-proj1#n by p matrix representing fj(ak)
a<-apply(Fun, 1, ak_initial, mu=mu,d=d[,1])
#a<- with(svdx, as.numeric(u[,1]*d[1]))
sigma<-t(replicate(n,diag(cov(x))))#should be n by p (col for each ak)
Q<-function(a,fits,W,X){#f is a list of fits, if empty for first iteration
  n=length(a)
  p<-ncol(X)
  Fun<-matrix(data = NA,nrow = n,ncol = p)
  if(length(fits)==0){#first iteration
    mu<-colMeans(X)
    d<-svdx$v[,1]
    Fun<-t(d%*%t(a)+mu)
    Sigvec<-diag(cov(X))
    Qval=0
    for (k in 1:n) {
      S=(t(t(X)-Fun[k,]))^2
      second_term<-rep(0.5*sum(log(Sigvec)),n)
      Qval<-Qval+sum(S%*%(-1/Sigvec)-second_term)
    }
  }else{
    weights<-colSums(W)
    for (j in 1:p) {
      Fun[,j]=predict(fits[[j]],a)$y
    }
    Qval<-0
    for (k in 1:n) {#didn't include weights Wik.
      S=(t(t(X)-Fun[k,]))^2 #May be calculate once
      Sigvec=t(t(W[,k])%*%S)/sum(W[,k])#p by 1 vector of sigma_j^2(ak)
      second_term<-0.5*sum(log(Sigvec))*sum(W[,k])
      Qval<-Qval+t(W[,k])%*%S%*%(-1/Sigvec)-second_term#add wik here before summation
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
Qval<-rep(0,num_it+1)
Qval[1]<--Q(a,fits,W=NA,x)#Return value of Q is the negative of Q value in EM
Qval1<-Qval[1]

a_matrix <- matrix(nrow = num_it, ncol = n)
a_matrix[1,]<-a


for (h in 1:num_it) {
  Qval0<-Qval1
  W2<-matrix(nrow = n,ncol = n)
  for (i in 1:n) {
    gyi<-0
    for (k1 in 1:n) {
      gyi<-gyi+dmvnorm(x[i,],Fun[k1,],diag(sigma[k1,]))*vk[k1]
    }
    for (k in 1:n) {
      loggyiak<-log(dmvnorm(x[i,],Fun[k,],diag(sigma[k,])))
      #gy(y_h)
      
      W2[i,k]<-loggyiak-log(gyi)+log(vk[k])
    }
  }
  W2=exp(W2)
  rowSums(W2)
  if(h!=1){
    a<-nlm(f=Q, p=a, fits=fits, W=W2, X=x,iterlim = 1)$estimate

  }
  #now update f and sigma and vk
  a2<-max(a)
  a1<-min(a)
  weights<-colSums(W2)
  Dmtx<-diag(1/weights)
  for (j in 1:p) {
    yvalue<-Dmtx%*%t(W2)%*%x[,j]
    fit <- smooth.spline(x=a, y=yvalue, w=weights, df=degf)##fixed df because as iteration goes on the fit will be interplotation.
    Fun[,j]=fitted(fit)
    fits[[j]]=fit
  }
  for (k in 1:n) {
    S=(t(t(x)-Fun[k,]))^2
    sigma[k,]=t(t(W2[,k])%*%S)/sum(W2[,k])#p by 1 vector of sigma_j^2(ak)
  }
  vk<-colSums(W2)/n
  u<-u+1
  # plot data and the principal curve for a sequence of a's
  plot(x[,1], x[,2], bty='n',
       xlab=expression(x[1]),
       ylab=expression(x[2]))
  legend("topleft", legend=c("Proj 1"), bty="n")
  points(x=Fun[,1],y=Fun[,2],col="red")
  lines(z1,z2)
  
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
  Qval1<--Q(a,fits,W2,x)
  Qval[h+1]<-Qval1
  arc_length<-sum(diag(dist(Fun)))
  a<-arc_length*(a-a1)/(a2-a1)
  a_matrix[h+1,]=a
}
colMeans(sigma)
rowMeans(a_matrix)
rank(a_matrix[1,])==rank(a_matrix[2,])
