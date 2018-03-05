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

## the wrapper function around make.bisse()

library(diversitree)

BiSSE.model.analysis<-function(parmodel,phylfam,starter,strain.family,dfdata,binvar,ltrees){
    
    sink("outputFromWrapper.txt",append=TRUE)
    chrapetree<-ltrees[[which(names(ltrees)==paste(phylfam,strain.family,starter,sep="."))]]
    
    
    bnoerror<-FALSE ## the tree chronos is done earlier
    if ((!is.na(chrapetree))&&(!is.null(chrapetree))){bnoerror<-TRUE}
    
    res<-list(estpar=NA,lnLik=-Inf,AIC=Inf,AICc=Inf,BIC=Inf,bisse.mle.res=NA,bisse.obj=NA,chrapetree=chrapetree,vbindata=NA,parmodel=parmodel,starter=starter,strain.family=strain.family)

    if (bnoerror){
	tryCatch({
	    vbindata<-c(dfdata[,binvar])[chrapetree$tip.label]
	    names(vbindata)<-chrapetree$tip.label
	    
	    print(c(starter,strain.family))
	    print(vbindata)
	    print("=====================================")
	    
	    
	    ## setup the object for BiSSE depending on the assumptions on the parameters
	    bisse.obj<-make.bisse(chrapetree,vbindata)
	    bisse.eq<-function(pars,bisse.obj,fpar.trans){bisse.obj(fpar.trans(pars))}
	    ## estimate the parameters
	    bisse.mle.res<-find.mle(bisse.eq,parmodel$par0,method="optim",bisse.obj=bisse.obj,fpar.trans=parmodel$func.par.transform)
    
	    ## extract what is needed
	    names(bisse.mle.res$par)<-names(parmodel$par0)
	    estpar<-bisse.mle.res$par
	    model.par.est<-parmodel$func.par.transform(estpar)
    
	    ## extract likelihood and calculate information criteria
	    lnLik<-bisse.mle.res$lnLik
	    k<-length(parmodel$par0)
	    n<-length(vbindata)

	    AIC<- (-2)*lnLik+2*k
	    AICc<- AIC+ 2*k*(k+1)/(n-k-1)
	    BIC<- (-2)*lnLik + k*log(n)

	    res<-list(model.par.est=model.par.est,lnLik=lnLik,AIC=AIC,AICc=AICc,BIC=BIC,bisse.mle.res=bisse.mle.res,bisse.obj=bisse.obj,chrapetree=chrapetree,vbindata=vbindata,parmodel=parmodel,starter=starter,strain.family=strain.family,estpar=estpar)
	},error=function(e){print(paste("Error with BiSSE calculations: ",e,"for starter",starter,"and",strain.family,"and",binvar))})
    }
    sink()
    res
}

