

rm(list=ls())



 
path.main<-getwd()

path.out<-paste(path.main,"/Latex/",sep="")
path.sweave<-paste(path.main,"/Sweave/",sep="")
#------------------------------

recompute_results<-

paper<-"BIP_predictor"

script <- paste(path.sweave,paper,sep="")
## enforce par(ask=FALSE)
options(device.ask.default=FALSE)
## create a LaTeX file
Sweave(script,output=paste(path.out,paper,".tex",sep=""))

#Stangle(script)

