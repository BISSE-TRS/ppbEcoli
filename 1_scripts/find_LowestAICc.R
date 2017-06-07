## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## USE AT OWN RISK!!
## This code accompanies the paper
## Predicting pathogenicity behavior in E. coli population through a state-dependent model and TRS profiling
## by Krzysztof Bartoszek, Marta Majchrzak, Sebastian Sakowski, Anna B. Kubiak-Szeligowska, Ingemar Kaj and Pawel Parniewski
## one runs it in R as
## source("run_BiSSEanalysis.R")  ## estimate the parameters of the BiSSE model
## source("run_SummaryBiSSEparams.R") ## summarize the parameters
## in R' working directory the following files are required to be
## strainsData.RData the trait data
## the data has to be a matrix called strainsData
## each row corresponds to a strain, the rows have to be named with rows
## corresponding to the rows on the phylogeny
## the columns correspond to the various binary traits, each column has to be named
## the first columns describe the phylogenetic families to which the strains belong to
## each strain belongs to only one family
## the next column the actual pathogenicity traits
## the files containing the phylogenetic trees:
## CGG+GTG.txt, GTG.txt, CGG.txt ## all strains, different files for different starters
## K_CGG+GTG.txt, K_GTG.txt, K_CGG.txt ## only K strains
## U_CGG+GTG.txt, U_GTG.txt, U_CGG.txt ## only U straints
## =================
## the code takes a long time to run and will create in the working directory the following files:
## outputFromWrapper.txt : real time output from the estimation procedure (run_BiSSEanalysis.R), can be used to monitor the code's progress
## tmpoutput.txt : output from run_BiSSEanalysis.R, reporting what model is considered now, can be used for monitoring the code's progress
## BiSSEestimates.RData :  file created at the end of run_BiSSEanalysis.R containg the results of the estimation, can be huge
## BiSSEestimatesSummary.txt : a nicely formatted summary of the estimation
## BiSSEestimatesSummary.csv : model parameters and summary statistics for the different traits and setups in a csv file

## search through considered models to find the one with the lowest AICc score

fbest.no.death<-function(lresults,varnames,bfixed.starter=TRUE){
    lsearch<-list()
    if (!bfixed.starter){ltosearch<-c(lresults$all.tree.estimates$CGG$all.estimates,lresults$all.tree.estimates$GTG$all.estimates,lresults$all.tree.estimates$CGGGTG$all.estimates)}
    else{ltosearch<-lresults$all.estimates}
    vAICcvals<-sapply(ltosearch,function(lres){
	res<-Inf
	if ((!is.null(lres$model.par.est))&&(!is.na(lres$model.par.est[1]))){
	    names(lres$model.par.est)<-varnames
	    if ((lres$model.par.est["mu0"]==0)&&(lres$model.par.est["mu1"]==0)){res<-lres$AICc}
	}
	res
    },simplify=TRUE)
    minmodel<-which(vAICcvals==min(vAICcvals))
    if (length(minmodel)>1){cat("\n");cat("WARNING: more than one model with identically lowest AICc! Choosing an arbitrary one!");minmodel<-minmodel[1]}
    ltosearch[[minmodel]]
}

