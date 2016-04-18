rm(list=ls()) 

# loading dependencies

library(Matrix)
library(slam)
library(Rcplex)
library(Rgraphviz)

# parsing files and generating network structure

ilp_options<-read.delim(file="ilp_options.txt",header=FALSE)
mipgap<-as.numeric(as.character((ilp_options$V2[ilp_options$V1=="mipgap"])))
timelimit<-as.numeric(as.character((ilp_options$V2[ilp_options$V1=="timelimit"])))
populatelim<-as.numeric(as.character((ilp_options$V2[ilp_options$V1=="number_solutions"])))
network_file<-as.character(ilp_options$V2[ilp_options$V1=="network_file"])
net<-read.delim(file=network_file,header=FALSE)
nreactions<-nrow(net)
species_file<-as.character(ilp_options$V2[ilp_options$V1=="species_file"])
spc<-read.delim(file=species_file, header=FALSE) 
nspc<-nrow(spc)
inputs_file<-as.character(ilp_options$V2[ilp_options$V1=="inputs_file"])
inputs<-read.delim(file=inputs_file, header=FALSE) 
measurements_file<-as.character(ilp_options$V2[ilp_options$V1=="measurements_file"])
measurements<-read.delim(file=measurements_file, header=FALSE)

nexp<-(nrow(measurements)-1)
ninp<-ncol(inputs)
nsignals<-ncol(measurements)
inputs_hash<-array(NA, dim=c(ninp, 1))

for (i in 1:ninp) {
  inputs_hash[[inputs[1,i]]]<-i
}

upactivations<-list()
dactivations<-list()
upinhibitions<-list()
dinhibitions<-list()
for (i in 1:nspc) {
  upactivations[[i]]<-list()
  dactivations[[i]]<-list()
  upinhibitions[[i]]<-list()
  dinhibitions[[i]]<-list()
}

for (i in 1:nreactions) {
  if ((net$V2[i]=="1") | (net$V2[i]=="2")) {  
    upactivations[[(net$V3[i])]]<-c(upactivations[[(net$V3[i])]], i)
    dactivations[[net$V1[i]]]<-c(dactivations[[net$V1[i]]], i)
  } else {
    upinhibitions[[(net$V3[i])]]<-c(upinhibitions[[(net$V3[i])]], i)
    dinhibitions[[net$V1[i]]]<-c(dinhibitions[[net$V1[i]]],i)
  }
}

source("optimize.r")