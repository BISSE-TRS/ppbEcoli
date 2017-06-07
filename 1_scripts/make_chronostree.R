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

## make the dendrograms ultramatric, i.e. all tips contemporary

chronos.tree<-function(apetree,method,treename,tips.keep){
    library(ape)
    tips.drop<-setdiff(apetree$tip.label,tips.keep)
    apetree<-drop.tip(apetree,tips.drop)
    chrapetree<-NA
    tryCatch({
	    chrapetree<- chronos(apetree,model=method)
    },error=function(e){sink("outputFromWrapper.txt",append=TRUE);print(paste("Error with chronos: ",e,"for tree",treename,"and model",method));sink()})	
    chrapetree
}


