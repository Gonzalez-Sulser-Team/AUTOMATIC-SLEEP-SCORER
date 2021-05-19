
#library(signal)
library(multitaper)

# generates and returns a function to calculate multitapper spectrum
rat.genmtspec <- function(n,k,nw=k/2,fa=1,fb=n/2) {
  tapper <- dpss(n,k,nw,F)$v
  return( function(v) {
    v <- v - mean(v)
    x <- Mod(fft(v*tapper[,1])[fa:fb])^2
    if (k>1) for (i in 2:k) x <- x + Mod(fft(v*tapper[,i])[fa:fb])^2
    x/n
  } )
}

# aplies a function f to segments of s of size ws and overlaping we
rat.apply <- function(s,ws=30,we=0,f=function(chunk) chunk, ...) {
  if(is.null(dim(s))) dim(s) <- c(1,length(s))
  lmax <- ncol(s)
  N <- floor(lmax/ws)
  r <- ((1-we):(ws+we))
  v <- f(s[,pmax(r,1)], ...)
  m <- matrix(v,length(v),N)
  pb <- txtProgressBar(style=3,width=20)
  for(i in 2:N){
    m[,i] <- f(s[,pmin(pmax(r+(i-1)*ws,1),lmax)], ...)
    setTxtProgressBar(pb, i/N)
  }
  close(pb)
  if(nrow(m)==1) as.vector(m) else m
}

# returns a raw multitapper spectrogram
rat.rawmtspec <- function(x,sr=250.4,fr=c(1,20),w=15,e=0,k=10) {
  l <- w+2*e
  rat.apply(x,w*sr,e*sr,rat.genmtspec(l*sr,k,k*.75,1+l*fr[1],1+l*fr[2]))
}

# returns a log10 scaled multitapper amplitud spectrogram
rat.mtspec <- function(x,sr=250.4,fr=c(1,20),w=5,e=0,k=6,a=1e-6,nm=T,nsd=F,nrd=F) {
  l <- w+2*e
  log10(rat.apply(x,w*sr,e*sr,rat.genmtspec(l*sr,k,k*.75,1+l*fr[1],1+l*fr[2]))+a)/2
}

# normalizes a spectrogram relatively to a given quantile
rat.norm <- function(x,q=.05) {
  t(apply(x,1,function(v) v-quantile(v,q)))
}

# returns a matrix of the euclidian distances of each element (column) of x to its z neighbours (z is the column displacement)
rat.ldist <- function(x,z=-5:5) {
  n <- ncol(x)
  t( sapply( z,
             function(i) {
               k <- 1:n + i
               if (i<0) k[1:-i] <- (1-i):2
               else if (i>0) k[n+1+(-1:-i)] <- n + (-i:-1)
               sqrt(colSums((x-x[,k])^2))
             }
  ))
}

# returns a locally smoothed version of x replacing each column by the mean of its k closest neigbours at d range
rat.lsmth <- function(x,d=5,k=d+1,f=1) {
  l <- 2*d+1
  n <- ncol(x)
  i <- sapply( -d:d,
               function(i) {
                 k <- 1:n + i
                 if (i<0) k[1:-i] <- (1-i):2
                 else if (i>0) k[n+1+(-1:-i)] <- n + (-i:-1)
                 k
               }
  )
  m <- rat.ldist(x,-d:d)
  i <- apply(rbind(t(i),m),2,function(v) v[1:l][order(v[(l+1):(2*l)])])
  r <- x[,i[1,]]*f
  if (k>1) for (j in 2:k) r <- r + x[,i[j,]]
  r/(k+f-1)
}

# band activity
rat.bact <- function(b,p=.5) {
  i = round(nrow(b)*p)
  a <- apply(b,2,function(v) cumsum(sort(v,decreasing=T)))
  a[i,]/i
}

# point (column) number of neighbours at d distance (euclidian)
rat.nn <- function(m,d=.1) {
  w <- apply(m,2,function(p) sum(sqrt(colSums(sweep(m,1,p)^2))<d))
}

# dynamic distance
rat.dd <- function(m,l=2) {
  n <- ncol(m)
  sqrt(colSums((m[,c((1+l):n,rep(n,l))] - m[,c(rep(1,l),1:(n-l))])^2))/(2*l+1)
}

# point (column) of maximum density (euclidian)
rat.fmaxd <- function(m,q=1/3) {
  w <- apply(m,2,function(p) quantile(sqrt(colSums(sweep(m,1,p)^2)),q))
  m[,which.min(w)]
}

# covariance matrix of m with respect to centroid v
rat.fcov <- function(m,v) { m <- sweep(m,1,v); m%*%t(m)/ncol(m) }

# mahalanobis distance of points p (columns) to cluster
rat.fmhd <- function(p,ccnt,ccvm) {
  p <- sweep(p,1,ccnt); sqrt(colSums(p*(solve(ccvm)%*%p)))
}

# MAIN processing of eeg and emg signals (250Hz)
rat.proc <- function(eeg,emg) {
  x <<- rat.spectr(eeg,emg)
  p <<- rat.3dprojection(x)
  r <<- rat.clusters(p$r,p$s,p$d)
  s <<- rat.state(r$mhd.r,r$mhd.s)
}

# spectrograms of eeg and emg signals, sampling rate sr, window size w, starting frequency fe1...
rat.spectr <- function(eeg,emg,sr=250.4,w=5) {
  # eeg and emg spectrograms:
  x.eeg <- rat.norm(rat.mtspec(eeg,sr=sr,fr=c(0,sr/2),w=w))
  x.emg <- rat.norm(rat.mtspec(emg,sr=sr,fr=c(0,sr/2),w=w))
  # smoothed versions:
  s.eeg <- rat.lsmth(x.eeg)
  s.emg <- rat.lsmth(x.emg)
  list(eeg=list(r=x.eeg,s=s.eeg),emg=list(r=x.emg,s=s.emg))
}

# frequency ranges
rat.3dvars <- list(muscle=c(60,90),theta=c(5.8,8.2),total=c(1,20))

# 3d projection: (muscle,theta,delta+sigma)
rat.3dprojection <- function(x,w=5) {
  i <- lapply(rat.3dvars, function(p) p*w+1)
  ms <- i$muscle[1]:i$muscle[2]
  th <- i$theta[1]:i$theta[2]
  ds <- c(i$total[1]:(th[1]-1),(th[2]+1):i$total[2]) # total except theta
  # 3d projection variables: r[1,]=muscle, r[2,]=theta, r[3,]=delta+sigma
  r <- rbind(2*colMeans(x$emg$r[ms,]),colMeans(x$eeg$r[th,]),colMeans(x$eeg$r[ds,])) # real 3d projection variables
  s <- rbind(2*colMeans(x$emg$s[ms,]),colMeans(x$eeg$s[th,]),colMeans(x$eeg$s[ds,])) # smoothed ones
  # distance between real and smoothed projections
  d <- sqrt(colSums((r-s)^2))
  # distance between successive smoothed projections
  sd <- sqrt(colSums((s[,-1]-s[,-ncol(s)]))^2)
  sd <- (c(sd[1],sd)+c(sd,sd[length(sd)]))/2
  list(r=r,s=s,d=d,sd=sd)
}

# clusters of the 3d projection, centroids, mahalanobis distance to centroids
rat.clusters <- function(r,s,d,qty=rep(1,ncol(r))) {
  dd <- rat.dd(s)
  d.i <- which(dd<quantile(dd,.8) & qty==1) # stable points
  # delta/sigma threshold
  ds.q <- quantile(s[3,d.i],c(.1,.6))
  ds.h <- hist(pmin(pmax(s[3,d.i],ds.q[1]),ds.q[2]),40,plot=F)
  ds.t <- ds.h$mids[which.min(ds.h$counts)]
  ds.l <- d.i[s[3,d.i]<ds.t] # stable points under delta/sigma threshold
  # muscle threshold
  ms.q <- quantile(s[1,d.i],c(.1,.6))
  ms.h <- hist(pmin(pmax(s[1,d.i],ms.q[1]),ms.q[2]),40,plot=F)
  ms.t <- ms.h$mids[which.min(ms.h$counts)]
  ms.l <- d.i[s[1,d.i]<ms.t] # stable points under muscle threshold
  # delta/sigma threshold for low muscle
  lds.q <- quantile(s[3,ms.l],c(.1,.6))
  lds.h <- hist(pmin(pmax(s[3,ms.l],lds.q[1]),lds.q[2]),40,plot=F)
  lds.t <- lds.h$mids[which.min(lds.h$counts)]
  # muscle threshold for low delta/sigma
  lms.q <- quantile(s[1,ds.l],c(.1,.6))
  lms.h <- hist(pmin(pmax(s[1,ds.l],lms.q[1]),lms.q[2]),40,plot=F)
  lms.t <- lms.h$mids[which.min(lms.h$counts)]
  # tentative states, 0:W, 1:NREM, 2:REM, 3:?, 4:Transition
  e <- 0*d.i+3 # default state is 3:?
  e[s[1,d.i]<ms.t] <- 1
  e[s[1,d.i]<lms.t & s[3,d.i]<lds.t] <- 2
  e[s[1,d.i]>=lms.t & s[3,d.i]<ds.t] <- 0
  bs <- rep(4,ncol(r)) # 4 for discarded epochs (not in d.i)
  bs[d.i] <- e
  bs.c <- bs
  # for (j in 1:3) {
  #   k <- which(bs==j)
  #   v <- rat.nn(s[,k])
  #   bs.c[k[v<quantile(v,.1)]] <- 4 # 4 for excentric epochs
  # }
  # first aproximation to centroids of clusters
  # centroids of clusters (median)
  ccnt <- sapply(0:2,function(x) apply(s[,bs.c==x],1,median))
  # clusters shape (covariance matrix)
  ccvm <- lapply(0:2,function(x) rat.fcov(s[,bs.c==x],ccnt[,x+1]))
  # mahalanobis distances to clusters
  mhd.s <- t(sapply(0:2,function(x) rat.fmhd(s,ccnt[,x+1],ccvm[[x+1]])))
  mhd.r <- t(sapply(0:2,function(x) rat.fmhd(r,ccnt[,x+1],ccvm[[x+1]])))
  list(mhd.s=mhd.s,mhd.r=mhd.r,d.i=d.i,bs=bs,bs.c=bs.c)
}

# behavioural state (0:W, 1:REM, 2:NREM, 7:None, 8:LQ, 9:Undefined)
rat.state <- function(mhd.r,mhd.s,qty=rep(1,ncol(mhd.r))) {
  L <<- ncol(mhd.r) # length
  o.r <- apply(mhd.r,2,order)-1 # ordered states according to closest mahalanobis distance (real)
  o.s <- apply(mhd.s,2,order)-1 # ordered states according to closest mahalanobis distance (smoothed)
  d.r <- apply(mhd.r,2,sort) # ordered mahalanobis distances (real)
  d.s <- apply(mhd.s,2,sort) # ordered mahalanobis distances (smoothed)
  
  #MARK Abnormalities likely to be swd or W
  s_mark <- o.r[1,] # the state is the one with minimum real distance (W, REM, NREM), but it is
  s_mark[d.r[1,]>3.5 & d.s[1,]>1] <- 4 # Undefined for large distances
  s_mark[o.r[1,]!=o.s[1,] | qty==0] <- 0 # Inconsistency or low qty means W
  s_mark[d.r[1,]>4 & d.s[1,]>2] <- 5 # and assigned to None if distance too large
  
  # Interrupted REM which is NOT REM
  s_like <- s_mark # Marked EEG abnormalities likely to be swd to correct surround
  k <- c(1,which(s_like[2:L]!=s_like[1:(L-1)])+1) # starting points of new states
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(s_like[k]==0 & (j-k)<4 & (j-k)>1) #W length =2
  i <- which (s_like[k-1]==1 & s_like[j+1]==2) #W length =2 preceded by NR, followed by R
  for (m in i) s_like[k[m]:k[m]+3] <- 0 #To be considered as W

   # Interrupted REM 1
  # s_like <- s_mark # Marked EEG abnormalities likely to be swd to correct surround
  k <- c(1,which(s_like[2:L]!=s_like[1:(L-1)])+1) # starting points of new states
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(s_like[k]==2 & (j-k)<2) #R states of length 1  
  i <- which(s_like[k-1]==1 & s_like[j+1]==0) #R states of length 1 preceded by NR followed W
  for (m in i) s_like[k[m]+1] <- 2 # to be considered as R
  
  #k <- c(1,which(s_like[2:L]!=s_like[1:(L-1)])+1) # starting points of new states
  #l <- length(k)
  #k <- k[c(k[2:l],L+1)-k==1 & s_like[k]==2] # R states of length 1
  #k <- k[s_like[k-1]==1] # R states of length 1 preceded by NR 
  #k <- k[s_like[k+1]==0] # R states of length 1 preceded by NR followed W
  #s_like[k+1] <- 2 # will be considered as REM
  
  # Interrupted REM 2
  k <- c(1,which(s_like[2:L]!=s_like[1:(L-1)])+1) # starting points of new states
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(s_like[k]==5 & (j-k)<2) #5 length =1  
  i <- which(s_like[k-1]==1) #5 length =1 preceded by NR 
  i <- which(s_like[k+1]==0 & s_like[j+2]==2) #5 length =1 preceded by NR, followed by 0 and then R
  for (m in i) s_like[k[m]:k[m]+1] <- 2 # 5 and W to be considered as R
  
  # Interrupted REM 3
  k <- c(1,which(s_like[2:L]!=s_like[1:(L-1)])+1) # starting points of new states
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(s_like[k]==5 & (j-k)<2) #5 length =1  
  i <- which(s_like[k-1]==1) #5 length =1 preceded by NR 
  i <- which(s_like[k+1]==4 & s_like[j+2]==2) #5 length =1 preceded by NR, followed by 4 and then R
  for (m in i) s_like[k[m]:k[m]+1] <- 2 # 5 and 4 to be considered as R
  
  # Small 4 surrounded by NR
  #s_like <- s_mark # Marked EEG abnormalities likely to be swd to correct surround
  k <- c(1,which(s_like[2:L]!=s_like[1:(L-1)])+1) # starting points of new states 
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(s_like[k]==4 & (j-k)<3 & s_like[k-1]==1 & s_like[j+1]==1) #4 length <3 surrounded by NR 
  for (m in i) s_like[k[m]:j[m]] <- 1 #To be considered a NR
  
  # Correct wrong NREM inside SWD IF SEIZURES ==1
  if (seizures==1) { 
  # s_like <- s_mark # Marked EEG abnormalities likely to be swd to correct surround
  k <- c(1,which(s_like[2:L]!=s_like[1:(L-1)])+1) # starting points
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(s_like[k]==1 & (j-k)<4 & s_like[k-1]>3 & s_like[j+1]>3) #NREM length <5 surrounded by 5 or 4
  for (m in i) s_like[k[m]:j[m]] <- 0 #To be considered as W
  }
  
  # first approximation: remaining abnormalities to W 
  s1 <- s_like
  for (i in 1:length(s_like)) (
  if (s_like[i]>3) s1[i]<-0  # abnormalities changed to W
  else if  (s_like[i]<3) s_like[i]<-s1[i]
  )
  
  # seizure correction:
  s4 <- s1
   for (i in 1:length(swd)) (
    if (swd[i]==6) s4[i]<-0
    else if  (swd[i]==4) s4[i]<-s1[i]
  )
  
  # second approximation:
  s2 <- s4
  k <- c(1,which(s2[2:L]!=s2[1:(L-1)])+1) # starting points of new states
#  if (s2[1]==0) { s2[1:(k[2]-1)] <- s2[k[2]]; k <- k[-2] } # undefined at start replaced by following state
#  l <- length(k)
#  if (s2[L]==0) { s2[k[l]:L] <- s2[k[l]-1]; k <- k[-l] } # undefined at end replaced by previous state
 
  # NR singles surrounded by same state  
  #k <- c(1,which(s2[2:L]!=s2[1:(L-1)])+1) # starting points of new states
  #l <- length(k)
  #k <- k[c(k[2:l],L+1)-k==1 & s2[k]==1] # NR states of length 1
  #k <- k[s2[k-1]==s2[k+1]] # NR states of length 1 surrounded by the same state
  #s2[k] <- s2[k-1] # are replaced by that state
 
  k <- c(1,which(s2[2:L]!=s2[1:(L-1)])+1) # starting points of new states
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(s2[k]==1 & (j-k)<2) # NR states of length 1
  i <- which(s2[k-1]==s2[k+1]) # NR states of length 1 surrounded by the same state 
  for (m in i) s2[k[m]] <- s2[k[m]+1] # are replaced by that state
  
  
  # W singles surrounded by same state 
  #k <- c(1,which(s2[2:L]!=s2[1:(L-1)])+1) # starting points of new states
  #l <- length(k)
  #k <- k[c(k[2:l],L+1)-k==1 & s2[k]==0] # W states of length 1
  #k <- k[s2[k-1]==s2[k+1]] # NR states of length 1 surrounded by the same state
  #s2[k] <- s2[k-1]  # are replaced by that state
  
  k <- c(1,which(s2[2:L]!=s2[1:(L-1)])+1) # starting points of new states
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(s2[k]==0 & (j-k)<2) # W states of length 1
  i <- which(s2[k-1]==s2[k+1]) # W states of length 1 surrounded by the same state 
  for (m in i) s2[k[m]] <- s2[k[m]+1] # are replaced by that state
  
  # R singles surrounded by same state replaced by that state
  #k <- c(1,which(s2[2:L]!=s2[1:(L-1)])+1) # starting points of new states
  #l <- length(k)
  #k <- k[c(k[2:l],L+1)-k==1 & s2[k]==2] # R states of length 1
  #k <- k[s2[k-1]==s2[k+1]] # R states of length 1 surrounded by the same state
  #s2[k] <- s2[k-1]  # are replaced by that state
  
  k <- c(1,which(s2[2:L]!=s2[1:(L-1)])+1) # starting points of new states
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(s2[k]==2 & (j-k)<2) # NR states of length 1
  i <- which(s2[k-1]==s2[k+1]) # NR states of length 1 surrounded by the same state 
  for (m in i) s2[k[m]] <- s2[k[m]+1] # are replaced by that state
  
  # Isolated small NREM inside W
  #k <- c(1,which(s2[2:L]!=s2[1:(L-1)])+1) # starting points
  #j <- c(k[-1]-1,L) # corresponding ending points
  #k <- k[-c(1,length(k))] # first and last state not considered
  #j <- j[-c(1,length(j))]
  #i <- which(s2[k]==1 & (j-k)<4 & s2[k-1]==0 & s2[j+1]==0) #NREM length <4 surrounded by W
  #for (m in i) s2[k[m]:j[m]] <- 0 #To be considered as W
  #k <- c(1,which(s2[2:L]!=s2[1:(L-1)])+1) # starting points of new states
  
  
  k <- which(c(s2[1],s2[-L])!=s2 & c(s2[-1],s2[L])!=s2) # isolated states
  s2[k] <- 9 # are replaced by undefined
  k <- c(1,which(s2[2:L]!=s2[1:(L-1)])+1) # new starting points of new states
  j <- c(k[-1]-1,L) # corresponding ending points
  i <- which(s2[k]==9) # undefined states
  i <- i[s2[k[i]-1]==s2[j[i]+1]] # undefined states of any length surrounded by the same state
  for (m in i) s2[k[m]:j[m]] <- s2[k[m]-1] # are replaced by that state
  
  # third aproximation: 
  s3 <- s2
  k <- which(s2==9) # remaining undefined states
  for (i in k) if (s3[i-1]!=2) s3[i] <- s3[i-1] # states except REM are extended forward over undefined epochs
  k <- which(s2==9) # remaining undefined states
  for (i in rev(k)) s3[i] <- s3[i+1] # states are extended backward over undefined epochs
  
  # Isolated small REM 
  k <- c(1,which(s3[2:L]!=s3[1:(L-1)])+1) # starting points
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(s3[k]==2 & (j-k)<3 & s3[j+1]==0) #Small REM length <3 followed by W
  for (m in i) s3[k[m]:j[m]] <- 0 #To be considered as W
  list(o.r=o.r,o.s=o.s,d.r=d.r,d.s=d.s,s1=s1,s2=s2,s3=s3,s4=s4,s_mark=s_mark,s_like=s_like)
}
