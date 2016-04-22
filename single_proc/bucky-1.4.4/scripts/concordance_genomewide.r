#------------ Examples of use -----------
#
# first, source the code contained in the file:
# source("concordance_genomewide.r")
#
# then enter the posterior distribution of the samplewide CF
# for the clade of interest. 5 examples below, all of which from
# an analysis on 10 genes:
# sample.pp = c(0,0,0,0,0,0,0,.2,.3,.4,.1)
# sample.pp = c(0,0,0,0,0,0,0,0,0,0,1)
# sample.pp = c(0,0,0,0,0,0,0,0,0,0.000060,0.999940)  # clades:2 in 8
# sample.pp = c(0.999320,0.00059,0.00009,0,0,0,0,0,0,0,1)  # 2 in 8
# sample.pp = c(0.999320,0.00059,0.00009,0,0,0,0,0,0,0,1)  # 2 in 8
# N=100  # size of the grid between 0 and 1
# mat2.6 = dpPosteriorWeights(alpha=1,n=10,prob.clade(8,2),N=N)
# xx = 0:N/N
# genome.pp = mat2.6 %*% sample.pp
# with the line above, values in genome.pp sum up to 1. In order to get
# a density which continuously integrate to 1, one would multiply by N.
# plot(xx, genome.pp,type="h")
# distribution.summary(genome.pp)
# For the mean: use analytic formula
# genomewide(sample.pp,alpha=1,n=10,8,2,conf.level=.95)

# Older method: need to specify a genome size N
# sample.pp = c(0,0,0,0,0,0,0,.2,.3,.4,.1)
# mat4.5 = dpPosteriorWeights(alpha=1,N=6000,n=10,prob.clade(9,4))
# xx = 0:6000/6000
# genome.pp = mat4.5 %*% sample.pp
# distribution.summary(genome.pp)
# plot(xx, genome.pp,type="h")

UB = function(s){
# Number of unrooted binary trees with s taxa
 s = floor(s)
 if (s<4){return(1)}
 else{    return( prod(2*(4:s)-5))}
}

logUB = function(s){
# Logarithm of the number of unrooted binary trees with s taxa
 s = floor(s)
 if (s<4){ return(0)}
 else{ 
  return( sum(log(2*(4:s)-5)))
 }
}

prob.clade = function(Ntax,Nclade){
 # returns the prior probability of a clade.
 # Ntax = total # of taxa, Nclade= # of taxa in one part of the clade.
 exp(logUB(Nclade+1)+logUB(Ntax-Nclade+1)-logUB(Ntax))
}

lA = function(alpha,n,nonly=TRUE){
# calculates log(A_n(alpha)/A_n(1)) = log( alpha*...*(alpha+n-1) / n! )
# if nonly is TRUE. Otherwise, it returns the aforementioned value
# for all integers between 1 and n.
 if (nonly){
  res=sum(log(1+(alpha-1)/1:n))
 } else{
  res=c(0,cumsum(log(1+(alpha-1)/1:n)))
 }
 return(res)
}

dpPosteriorWeights = function(alpha,N,n,pclade=0,Ntax=0,Nclade=0){
# N = total number of genes in the genome
# n = number of genes in the sample.
# pclade = prior probability of the clade (or feature of interest). If
# not provided, then Ntax and Nclade must be provided.
# Ntax = total number of taxa
# Nclade = # of taxa on one side of the bipartition 
# output: (N+1) x (n+1) matrix of weights.
# The vector of genome-wide posterior probabilities
# P{k genes out of N have the clade|Data} (k varying between 0 and N)
# will then be obtained by multiplying the matrix of weights by the
# vector of sample-wide posterior probabilities
# P{j genes out of n have the clade|Data} (j varying between 0 and n)

 if (pclade == 0 && (Ntax == 0 || Nclade == 0)){
  cat("Ntax or Nclade is missing, as well as the clade probability (pclade)\n")
  return(NULL)
 } else {
  if (pclade ==0){ pclade = prob.clade(Ntax,Nclade)}
  wts1 = sapply(alpha*pclade+(0:n) , lA, n=N-n,nonly=F)
  wts2 = sapply(alpha*(1-pclade)+(0:n), lA, n=N-n,nonly=F)
  wts3 = lA(alpha+n,N-n,nonly=T)
  # these are (N-n+1)x(n+1) matrices
  wts4 = matrix(0,N+1,n+1)
  for (j in 0:n){
   wts4[1+j+0:(N-n),j+1] = exp(wts1[,j+1] + wts2[(N-n):0+1,n-j+1] - wts3)
  }
  return(wts4)
 }
}

genomewide = function(samplewide,alpha,N,n,Ntax,Nclade,conf.level=.95){
# samplewide = vector of posterior probabilities for the
#              sample-wide concordance factors. 
# alpha: prior level of discordance
# N = total number of genes in the genome
# n = number of genes in the sample.
# Ntax = total number of taxa
# Nclade = # of taxa in one side of the bipartition
# output = vector of posterior probabilities for the
#          genome-wide concordance factors. 

 mat = dpPosteriorWeights(alpha,N,n,pclade=0,Ntax,Nclade)
 genomePP = mat %*% samplewide
 plot((0:N)/N,genomePP,type="h",xlab="genome-wide concordance factor",
      ylab="posterior probability")
 distribution.summary(genomePP,conf.level=conf.level)
}

genomewide = function(samplewide,alpha,n,Ntax,Nclade,conf.level=.95){
# overwrites previous function: assumes infinite genome
# samplewide = vector of posterior probabilities for the
#              sample-wide concordance factors. 
# alpha: prior level of discordance
# n = number of genes in the sample.
# Ntax = total number of taxa
# Nclade = # of taxa in one side of the bipartition
# output = posterior mean and sd for the
#          genome-wide concordance factor. 

 if (n != length(samplewide)-1) 
  stop(paste("n (",n,") does not correspond to 'samplewide'", sep=""))
 #N=1000 # number of points for the grid
 pclade = prob.clade(Ntax,Nclade)
 #mat = dpPosteriorWeights(alpha,n,pclade=pclade,Ntax,Nclade,N=1000)
 #genomePP = mat %*% samplewide
 #plot((0:N)/N,genomePP,type="h",xlab="genome-wide concordance factor",
 #     ylab="posterior probability")
 #distribution.summary(genomePP,conf.level=conf.level)
 m = alpha/(alpha+n)*pclade + n/(alpha+n)*sum(samplewide*seq(0,1,by=1/n))
 m2j=(alpha*pclade+0:n)*(alpha*pclade+0:n +1)/((alpha+n)*(alpha+n+1))
 m2 = sum(samplewide*m2j)
 s  = sqrt(m2-m^2)
 return(list(mean=m,sd=s))
}


distribution.summary = function(probs, values=0, conf.level=.95){
# returns mean, standard deviation and 95% confidence limits
# of a probability distribution. probs = vector of probabilities.
# values: if 0 (default), the values of the random variable are
# assumed to be 0,1/N,2/N,...,1 where (N+1) is the length of the
# vector "probs" (adapted to concordance factors).
# Otherwise, "values" needs to be a vector of same length as "probs".  
 probs = probs/sum(probs) 
 if (length(values)==1 && values==0) {
  Ngenes = length(probs)-1
  values = 0:Ngenes/Ngenes
 }
 N = length(values);
 a = (1-conf.level)/2;
 u = cumsum(probs);
 Spinf = min(which(u >= a))
 Spinf = values[Spinf]
 Spsup = min(which(u >= 1-a),N)
 Spsup = values[Spsup]
 mu = weighted.mean(values,probs);
 va = weighted.mean(values^2,probs) - mu^2;

 return(list(mean=mu,std.dev=sqrt(va), lower=Spinf,upper=Spsup));
}


dpPosteriorWeights = function(alpha,n,pclade=0,Ntax=0,Nclade=0,N=1000){
# overwrites previous definition of the function. 
# assumes infinite number of genes in the genome
# N = number of points in the grid
# n = number of genes in the sample.
# pclade = prior probability of the clade (or feature of interest). If
# not provided, then Ntax and Nclade must be provided.
# Ntax = total number of taxa
# Nclade = # of taxa on one side of the bipartition 
# output: (N-1) x (n+1) matrix of weights.
# The vector of genome-wide posterior probabilities
# P{proportion p of genes have the clade|Data} (p in the grid)
# will then be obtained by multiplying the matrix of weights by the
# vector of sample-wide posterior probabilities
# P{j genes out of n have the clade|Data} (j varying between 0 and n)

 if (pclade == 0 && (Ntax == 0 || Nclade == 0))
  stop("Ntax or Nclade is missing, as well as the clade probability (pclade)\n")
 if (pclade ==0){ pclade = prob.clade(Ntax,Nclade)}
 qclade=1-pclade
 lq=matrix(NA,N+1,n+1)
 logalphap = log(alpha*pclade+0:n)
 logalphaq = log(alpha*qclade+0:n)
 ii      = seq(1,N-1)+1
 xinside = seq(1,N-1)/N
 u = log(xinside)-log(1-xinside)
 lq[ii,1+1] = alpha*pclade*log(xinside)+(alpha*qclade+n-2)*log(1-xinside)
 lq[ii,0+1] = lq[ii,1+1]-u+logalphap[0+1]-logalphaq[n]
 for (j in 2:n) lq[ii,j+1]=lq[ii,j]+u-logalphap[j]+logalphaq[n-j+1]
 bestj = round(n/2+alpha*(1/2-pclade))
# allz = apply(exp(lq),2,sum,na.rm=T)/N; 
# z = allz[max(min( n , bestj+1),2)]
 z = sum(exp(lq[,bestj+1]),na.rm=T)/N  

 wts = exp(lq-log(z))
 wts[1,  1:n +1]=0; 
 wts[N+1,1:n   ]=0;
 wts[1,1]    =N-sum(wts[2:(N+1),1])
 wts[N+1,n+1]=N-sum(wts[1:N,n+1])
 b = beta(alpha*pclade+1,alpha*qclade+n-1)
 #cat("logalphap=",logalphap,"\nlogalphaq=",logalphaq,"\n")
 #cat("u=",u,"\nlq=",lq,"\n")
 cat("bestj=",bestj,", beta=",b,"\n")
 cat("logz=",log(z),", z=",z,"\n")
 cat("max precision:",(z-b)/b,"\n")
 matplot(wts,type="l",lty=1,xaxt="n",ylim=c(0,max(wts[2:N,2:n])))
 axis(side=1,at=seq(0,N,length.out=11),labels=seq(0,1,by=.1))
 return(wts)
}
#dpPosteriorWeights(alpha=2,n=5,pclade=.01)

