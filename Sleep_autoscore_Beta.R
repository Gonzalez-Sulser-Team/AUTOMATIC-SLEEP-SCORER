# -*- coding: utf-8 -*-
# AUTOMATIC SLEEP SCORER - BETA VERSION
#Original version by Dr. Alejandro Bassi and Dr. Javier Díaz in https://doi.org/10.3389/fncel.2017.00302 - Sleep & Chronobiology Lab - ICBM - Faculty of Medicine - Universidad de Chile
#Ported to R by Dr. Alejandro Bassi (abassi@dcc.uchile.cl)
#Adapted by Dr. Ingrid Buller (ingrid.buller@ed.ac.uk)
#@Gonzalez-Sulser-Team CDBS - SIDB - University of Edinburgh 2021
#
#
# test:
# rat.read()   # reads the signals (required first)
# rat.xtest()  # generates and displays the espectrogram (12 x 2h rows)
# rat.test()   # generates and displays the 3d projection with clusters (blue:NREM, red:REM, green:W)
# rat.tests()  # generates the behavioural states (requires running first rat.test())
#              # the result in variable s$s1 (0:W, 1:NREM, 2:REM, 7:None)
# rat.draw()   # to display the results in 3D space



setwd("/Users/ibull/Documents/POSTDOC/ANALISYS/EEG RESULTS/For Igor") #SET SOURCE DIR FOR CODE AND LIBRARY
getwd()
source("mgg.R")
source("sleep_score_library.R")
library(rgl)

rat.patchsgn <- function(s) {
  l <- length(s)
  k <- match(T,s!=0)
  if (k>1) s[1:(k-1)] <- s[k]
  k <- match(T,rev(s)!=0)
  if (k>1) s[(l-k+1):l] <- s[l-k]
  i <- which(s[1:(l-1)]!=0 & s[2:l]==0)
  j <- which(s[2:l]!=0 & s[1:(l-1)]==0) + 1
  d <- j-i
  n <- length(i)
  if (n>0) for (k in 1:n) s[(i[k]+1):(j[k]-1)] <- s[i[k]] + (s[j[k]]-s[i[k]])*(1:(d[k]-1))/d[k]
  s
}

rat.choose <-function() {
  setwd("/Users/ibull/Documents/POSTDOC/ANALISYS/EEG RESULTS/For Igor")
  getwd()
  folder<<-getwd()
  animal<<-paste(rat_line, "_", animal_id, sep="")
  dn<<-paste(animal, "_", animal_day, sep="")
  d<<-paste(folder, rat_line, animal, sep="/")
}


rat.find <-function (d=paste(folder, rat_line, animal, sep="/")) {
  rats <<- dir(d,"_Channels",full.names=T,recursive=T)
}

#rat.readspec <- function(f.eeg=rats[9],f.emg=rats[13]) {
#rat.readspec <- function() {
#  f.eeg<-paste(dn, eeg_val,".txt", sep="")
#  f.emg=paste(dn, emg_val,".txt", sep="")
#  eeg <- read.table(f.eeg)
#  emg <- read.table(f.emg)
#  rat.data <<- rbind(rat.patchsgn(eeg[[1]]),rat.patchsgn(emg[[1]]),(eeg[[1]]!=0)+0)
#}



#IF CHANNELS ON .csv FORMAT
rat.read <- function() {
  data_ch<-paste(d, "/", dn, "_Channels.csv", sep="")
  f_chan <- read.csv(data_ch)
  eeg <<- as.data.frame(f_chan[,1])
  emg <<- as.data.frame(f_chan[,2])
  #eeg <- read.table(f.eeg)
  #emg <- read.table(f.emg)
  rat.data <<- rbind(rat.patchsgn(eeg[[1]]),rat.patchsgn(emg[[1]]),(eeg[[1]]!=0)+0)
}


#IF CHANNELS ON .TXT FORMAT
#rat.read <- function(f.eeg=rats[9],f.emg=rats[13]) {
#rat.read <- function() {
#  f.eeg<-paste(dn, "eeg_OK.txt", sep="")
#  f.emg=paste(dn, "emg_OK.txt", sep="")
#  eeg <- read.table(f.eeg)
#  emg <- read.table(f.emg)
#  rat.data <<- rbind(rat.patchsgn(eeg[[1]]),rat.patchsgn(emg[[1]]),(eeg[[1]]!=0)+0)
#}

rat.xtest <- function() {
  qty <<- colMeans(matrix(rat.data[3,],250.4*5,ncol(rat.data)/250.4/5))
  x <<- rat.spectr(rat.data[1,],rat.data[2,])
  rat.xdraw(x$eeg$s[3:101,])
}

rat.test <- function(k=3) {
  p <<- rat.3dprojection(x)
  r <<- rat.clusters(p$r,p$s,p$d,qty>2/3)
  s <- numeric(ncol(p$s))+4; for (i in 1:3) s[r$mhd.s[i,]<1] <- i;
  if (k==1) plot3d(t(p$s),xlim=c(-.1,1.1),ylim=c(-.1,1.1),zlim=c(-.1,1.1),col=c(3,4,2,1,5)[s])
  else if(k==2) plot3d(t(p$s[,r$d.i]),xlim=c(-.1,1.1),ylim=c(-.1,1.1),zlim=c(-.1,1.1),col=c(3,4,2,1)[r$bs[r$d.i]+1])
  else if(k==3) plot3d(t(p$s),xlim=c(-.1,1.1),ylim=c(-.1,1.1),zlim=c(-.1,1.1),col=c(3,4,2,1,5)[r$bs.c+1])
}

rat.tests <- function() {
  s <<- rat.state(r$mhd.r,r$mhd.s,qty>2/3)
}

rat.xrform <- function(x,n=2,l=60*12,e=0) {
  if (e>0) x <- rbind(x,matrix(NA,e,ncol(x)))
  l <- n*l
  y <- x[,1:l]
  for (i in 1:(ncol(x)/l-1)) y <- rbind(x[,l*i+1:l],y)
  y
}

rat.xdraw = function(m,q=c(.05,.95)) {
  mgg.img(m=rat.xrform(m,e=5),q=q,xb=30*12,yb=Inf,col=viridis(256)) +
    labs(title="Spectrogram",x="epochs (5s)",y="freq. 0.4-20Hz 0.2Hz res.")
  grDevices::dev.new()
  print(last_plot())
  Sys.sleep(.1)
}

rat.draw = function(b=0,l=1440,i=1:l) {
  g <- mgg.plot()
  g <- mgg.img(g,x$eeg$s[5:101,b*l+i],ylim=c(-1,-.1),col=viridis(256))
  g <- mgg.img(g,x$emg$s[251:301,b*l+i],ylim=c(-2,-1.1),col=viridis(256))
  g <- mgg.lin(g,s$s3[b*l+i]) %>% mgg.grid(yb=1)
  g
  grDevices::dev.new()
  print(last_plot())
}

 
correct.A <- function(){
  sleep.score <<- s$s3
  sleep.score <<- as.data.frame(sleep.score)
  s.score <<- sleep.score
  for (i in 2:nrow(sleep.score)) (
    if ((sleep.score[i-1,1]=="7")&(sleep.score[i,1]=="2")) s.score[i,1]<-"8"
    else if  ((sleep.score[i-1,1]=="0")&(sleep.score[i,1]=="2")) s.score[i,1]<-"8"
    else s.score[i,1]<-sleep.score[i,1]
  )
  step.A<<-sum(s.score=="8")
  s.score <<- as.data.frame(s.score)
  print(sum(s.score == "8" ))
}

correct.B<- function(){
  s.correct <<- s.score
  for (i in 2:nrow(sleep.score)) (
    if ((sleep.score[i,1]=="2")&(s.score[i,1]=="8")) s.correct[i,1]<-"8"
    else if ((s.correct[i-1,1]=="8")&(sleep.score[i,1]=="2")) s.correct[i,1]<-"8"
    else s.correct[i,1]<-s.correct[i,1]
    )
  step.B<<-sum(s.correct=="8")
  s.correct <<- as.data.frame(s.correct)
  print(sum(s.correct == "8" ))
}

correct.C <- function(){
  s.dge <<- s.correct
  for(i in 1:nrow(s.dge)) (
   if((s.dge[i,1] == "7" ) | (s.dge[i,1] =="8")) s.dge[i,1]<- "0"
   else s.dge[i,1]<- s.dge[i,1]
  )
  step.C<<-sum(s.dge >"2")
  
  # Isolated small NREM inside W
  dge.ok<<-data.matrix(sapply(s.dge, as.numeric))
  L <<- 17280 # length
  k <- c(1,which(dge.ok[2:L]!=dge.ok[1:(L-1)])+1) # starting points
  j <- c(k[-1]-1,L) # corresponding ending points
  k <- k[-c(1,length(k))] # first and last state not considered
  j <- j[-c(1,length(j))]
  i <- which(dge.ok[k]==1 & (j-k)<4 & dge.ok[k-1]==0 & dge.ok[j+1]==0) #NREM length <4 surrounded by W
  for (m in i) dge.ok[k[m]:j[m]] <- "0" #To be considered as W
  
  #WRITE .CSV WITH SLEEP SCORE (0=WAKE, 1=NON REM, 2=REM, EPOCHS 5 SECS)
  dge.ok <<- as.data.frame(dge.ok)
  print(sum(dge.ok > "2" ))
  setwd(d)
  getwd() 
  dge_out <<- paste(dn,"-","dge_ok.csv",sep="")
  write.table(dge.ok, file = (dge_out), row.names=F, na="",col.names=T, sep=",")
  
  #WRITE .CSV WITH SLEEP SCORE AND SWD (0=WAKE, 1=NON REM, 2=REM, 4=SWD, EPOCHS 5 SECS)
  # add seizures to score for coherence analysis
  if (seizures=="1") {
  dge_seiz <- dge.ok
  for (i in 1:length(swd)) (
    if (swd[i]==6) dge_seiz[i]<-4
    else if  (swd[i]==4) dge_seiz[i]<-dge.ok[i]
  )
  dge_seiz <<- as.data.frame(dge_seiz)
  print(sum(dge_seiz > "2" ))
  setwd(d)
  getwd() 
  dge_out <<- paste(dn,"-","dge_swd.csv",sep="")
  write.table(dge_seiz, file = (dge_out), row.names=F, na="",col.names=T, sep=",")
  }
}

rat.hypno <- function(){ #Plot Hypnogram per hour
  hour<-rep(c(0:23),each=720)
  hour<-as.data.frame(hour)
  dge.ok<-as.matrix.data.frame(dge.ok)
  hour<-as.matrix.data.frame(hour)
  per_hour<-data.frame(hour,dge.ok)
  factor(per_hour$hour)
  factor(per_hour$sleep.score)
  require(dplyr)
  hr<-per_hour %>% count(per_hour$hour, per_hour$sleep.score)
  states <- c("0" = "WAKE", "1" = "NonREM", "2" = "REM")
  ggplot(hr,aes(x=hr$`per_hour$hour`,y=hr$n,fill=hr$`per_hour$sleep.score`)) + geom_bar(stat="identity") + ggtitle("States by Hour") + xlab("Hours since lights on") + ylab("epochs") + labs(fill="States") + scale_fill_manual(values=c("green","blue","red"), labels=states)
  grDevices::dev.new()
  print(last_plot())
}


###SEIZURE CORRECTION
swd.read <- function() {
  if (seizures=="1") {
    f=paste(d,"/", "seiz", dn,"_DGE_SWDs.csv")
    abs<-read.csv(file=f)
    swd<<-data.matrix(sapply(abs, as.numeric)) 
  } else if (seizures=="0") {
    abs<-as.data.frame(rep(4,17280))
    swd<<-data.matrix(sapply(abs, as.numeric)) 
  }
}


###-----------------------------------------------
###SPECTRAL POWER BY STATE


## Main function rat spec
rat.spctr_by_state <- function(xp=x$eeg$s,st=sdge,w=5) {
  df <- data.frame(hz=((1:nrow(xp))-1)/5)
  for (i in sort(unique(st))) df[[paste0("s_",i)]] <- rowMeans(xp[,st==i])
  return(df)
}

## Get PW spectrum by state
pw.spct_state <- function() {
  df <- rat.spctr_by_state()
  pwspectrum<<-df
 # spect_out <<- paste(animal_day,"-","spec.csv",sep="")
  if (max(sdge==2)){
    d1<<-rep(animal_day,times=4)
  }else if(max(sdge==4)){
    d1<<-rep(animal_day,times=5)
  }  
  spec_rat <<-rbind(d1,pwspectrum)
  #write.csv(df, file = (spect_out), row.names=F)
}

spectral.by.state <- function() {
  if (seizures=="1") {
    ss<<-dge_seiz
  } else {
    ss<<-dge.ok
  }  
  sdge<<-data.matrix(sapply(ss, as.numeric)) 
  pw.spct_state()
  spec_rat<<-as.data.frame(spec_rat)
  setwd(d)
  getwd() 
  spect_ratout <<- paste(dn,"-","pw_spectrum",sep="")
  #assign(spect_ratout,spec_rat, envir = parent.frame())
  spect_go <<- paste(spect_ratout,".csv",sep="")
  write.table(spec_rat, file = (spect_go), row.names=F, na="",col.names=T, sep=",")
  #write.csv(spect_ratout, file = (spect_go), row.names=F)
  print(spect_go)
}


rat.get <- function(){
  rat.choose()
  swd.read()
  #rat.find()
  rat.read()
  rat.xtest()
  rat.test()
  rat.tests()
}

rat.correct <- function(){
  correct.A()
  correct.B()
  correct.C()
  print (step.A)
  print (step.B)
  print (step.C)
  rat.hypno()
}

sleep.autoscore<- function(){
  rat.get()
  rat.correct()
  spectral.by.state()
}

