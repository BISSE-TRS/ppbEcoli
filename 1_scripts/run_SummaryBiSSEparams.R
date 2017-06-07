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

## report the output of the estimation in a nice way and calculate summary statistics

source("find_LowestAICc.R")

cdatafile<-"BiSSEestimates.RData" ## file created at the end of run_BiSSEanalysis.R
coutfile<-"BiSSEestimatesSummary.txt" ## file to write out the summary of the estimation
coutestcsv<-"BiSSEestimatesSummary.csv" ## csv file to store the estimated parameters and summary statistics

if (!is.element("parests",ls())){load(cdatafile)}
vphylfam<-c("A","B1+C","B2+F","D+E")

varnames<-c("lambda0","lambda1","mu0","mu1","q01","q10")
numofstats<-2*19


fsumstats<-function(vpars,noext=FALSE){
## function to calculate the summary statistics
    vsumstatnames<-c("ASlim0","ASlim1","lambda0/lambda1","mu0/mu1","q01/q10","lambda0/mu0","lambda1/mu1","lambda0/mu1","lambda1/mu0","lambda0/q01","lambda0/q10","lambda1/q01","lambda1/q10")
    vstats<-c(vpars,c(rep(NA,length(vsumstatnames))))
    k<-length(vpars);i<-1
    
## These formulae will NOT work if q01 or q10 is 0, in our case this did not happen

    L0<-vpars["lambda0"]-vpars["mu0"]
    L1<-vpars["lambda1"]-vpars["mu1"]
    q01<-vpars["q01"]
    q10<-vpars["q10"]
    if ((q01==0)||(q01==0)){cat("WARNING, one of the transistion rates is 0! Resulting limits are WRONG!");cat("\n")}
    H<- L0-L1-q01+q10
    sqrtDelta<-sqrt(H^2+4*q01*q10)
    vstats[k+i]<-(H+sqrtDelta)/(H+sqrtDelta+2*q01) ;i<-i+1
    vstats[k+i]<-(2*q01)/(H+sqrtDelta+2*q01) ;i<-i+1
    
    vstats[k+i]<-vpars["lambda0"]/vpars["lambda1"];i<-i+1
    vstats[k+i]<-vpars["mu0"]/vpars["mu1"];i<-i+1
    vstats[k+i]<-vpars["q01"]/vpars["q10"];i<-i+1
    vstats[k+i]<-vpars["lambda0"]/vpars["mu0"];i<-i+1
    vstats[k+i]<-vpars["lambda1"]/vpars["mu1"];i<-i+1
    vstats[k+i]<-vpars["lambda0"]/vpars["mu1"];i<-i+1
    vstats[k+i]<-vpars["lambda1"]/vpars["mu0"];i<-i+1
    vstats[k+i]<-vpars["lambda0"]/vpars["q01"];i<-i+1
    vstats[k+i]<-vpars["lambda0"]/vpars["q10"];i<-i+1    
    vstats[k+i]<-vpars["lambda1"]/vpars["q01"];i<-i+1
    vstats[k+i]<-vpars["lambda1"]/vpars["q10"];i<-i+1
        
    names(vstats)<-c(names(vpars),vsumstatnames)
    if (noext){names(vstats)<-paste("NoExtinction.",names(vstats),sep="")}
    vstats
}

## description of various rate setups
model.list<-list(
list(model.desc="All parameters free",func.par.transform=function(pars){pars},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,mu1=0.1,q01=0.1,q10=0.1)),
list(model.desc="Transition rates equal: q01=q10",func.par.transform=function(pars){c(pars,pars[5])},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,mu1=0.1,q01=0.1)),
list(model.desc="Death and transition rates equal: mu0=mu1, q01=q10",func.par.transform=function(pars){c(pars[1:3],pars[3],pars[4],pars[4])},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,q01=0.1)),
list(model.desc="Birth and transition rates equal: lambda0=lambda1, q01=q10",func.par.transform=function(pars){c(pars[1],pars[1],pars[2:4],pars[4])},par0=c(lambda0=0.1,mu0=0.1,mu1=0.1,q01=0.1)),
list(model.desc="Birth, death and transition rates equal: lambda0=lambda1, mu0=mu1, q01=q10",func.par.transform=function(pars){c(rbind(pars,pars))},par0=c(lambda0=0.1,mu0=0.1,q01=0.1)),
list(model.desc="Death rates equal: mu0=mu1",func.par.transform=function(pars){c(pars[1:3],pars[3],pars[4:5])},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,q01=0.1,q10=0.1)),
list(model.desc="Birth rates equal: lambda0=lambda1",func.par.transform=function(pars){c(pars[1],pars[1],pars[2:5])},par0=c(lambda0=0.1,mu0=0.1,mu1=0.1,q01=0.1,q10=0.1)),
list(model.desc="Birth and death rates equal: lambda0=lambda1, mu0=mu1 ",func.par.transform=function(pars){c(pars[1],pars[1],pars[2],pars[2],pars[3:4])},par0=c(lambda0=0.1,mu0=0.1,q01=0.1,q10=0.1)),

list(model.desc="mu0=0",func.par.transform=function(pars){c(pars[1:2],0,pars[3:5])},par0=c(lambda0=0.1,lambda1=0.1,mu1=0.1,q01=0.1,q10=0.1)),
list(model.desc="mu0=0 and transition rates equal: q01=q10",func.par.transform=function(pars){c(pars[1:2],0,pars[3:4],pars[4])},par0=c(lambda0=0.1,lambda1=0.1,mu1=0.1,q01=0.1,q10=0.1)),
list(model.desc="mu0=0 and birth, transition rates equal: lambda0=lambda1, q01=q10",func.par.transform=function(pars){c(pars[1],pars[1],0,pars[2:3],pars[3])},par0=c(lambda0=0.1,mu1=0.1,q01=0.1)),
list(model.desc="mu0=0 and birth rates equal: lambda0=lambda1",func.par.transform=function(pars){c(pars[1],pars[1],0,pars[2:4])},par0=c(lambda0=0.1,mu1=0.1,q01=0.1,q10=0.1)),

list(model.desc="mu1=0",func.par.transform=function(pars){c(pars[1:3],0,pars[4:5])},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,q01=0.1,q10=0.1)),
list(model.desc="mu1=0 and transition rates equal: q01=q10",func.par.transform=function(pars){c(pars[1:3],0,pars[4],pars[4])},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,q01=0.1,q10=0.1)),
list(model.desc="mu1=0 and birth, transition rates equal: lambda0=lambda1, q01=q10",func.par.transform=function(pars){c(pars[1],pars[1:2],0,pars[3],pars[3])},par0=c(lambda0=0.1,mu0=0.1,q01=0.1)),
list(model.desc="mu1=0 and birth rates equal: lambda0=lambda1",func.par.transform=function(pars){c(pars[1],pars[1:2],0,pars[3:4])},par0=c(lambda0=0.1,mu0=0.1,q01=0.1,q10=0.1)),

list(model.desc="mu0=0, mu1=0",func.par.transform=function(pars){c(pars[1:2],0,0,pars[3:4])},par0=c(lambda0=0.1,lambda1=0.1,q01=0.1,q10=0.1)),
list(model.desc="mu0=0, mu1=0 and transition rates equal: q01=q10",func.par.transform=function(pars){c(pars[1:2],0,0,pars[3],pars[3])},par0=c(lambda0=0.1,lambda1=0.1,q01=0.1)),
list(model.desc="mu0=0, mu1=0 and birth, transition rates equal: lambda0=lambda1, q01=q10",func.par.transform=function(pars){c(pars[1],pars[1],0,0,pars[2],pars[2])},par0=c(lambda0=0.1,q01=0.1)),
list(model.desc="mu0=0, mu1=0 and birth rates equal: lambda0=lambda1",func.par.transform=function(pars){c(pars[1],pars[1],0,0,pars[3:4])},par0=c(lambda0=0.1,q01=0.1,q10=0.1))
)


sink(coutfile)

vvirfacts<-names(parests)
for (i in 1:length(vvirfacts)){
    parests[[i]]$virfactname<-vvirfacts[i]
}

## create a nice output by extracting the properties of the considered
## setup and call the function to calculate the summary statistics
res.all<-sapply(parests,function(lvirfact.data,varnames,model.list,fsumstats,numofstats,vphylfam){
    virfact<-lvirfact.data$virfactname
    lvirfact.data<-lvirfact.data[-which(names(lvirfact.data)=="virfactname")]
    cat("Virulence factor: ");cat(virfact);cat("\n")
    allphylfam.res<-vector("list",length(vphylfam))
    phylfam.id<-1
    for (Phylfam in vphylfam){
	cat("Phylogenetic family: ");cat(Phylfam);cat("\n")
	lvirfact<-lvirfact.data[[which(names(lvirfact.data)==Phylfam)]]
	if (!is.element("strain.family",names(lvirfact$UK$best.model))){lvirfact$UK$best.model$strain.family<-"UK"}
	if (!is.element("strain.family",names(lvirfact$U$best.model))){lvirfact$U$best.model$strain.family<-"U"}
	if (!is.element("strain.family",names(lvirfact$K$best.model))){lvirfact$K$best.model$strain.family<-"K"}	
	res.fam<-sapply(lvirfact,function(ltreetyperes.allstarters,varnames,virfact,model.list,fsumstats,numofstats,Phylfam){	
	    vstarters<-c("CGG","GTG","CGGGTG")
    	    res.str<-vector("list",length(vstarters))
    	    starter.id<-1
    	    for (starter in vstarters){
        	res<-c(paste(virfact,".",Phylfam,".",ltreetyperes.allstarters$best.model$strain.family,".",starter,sep=""),rep(NA,numofstats))
        	names(res)<-c("id","lambda0","lambda1","mu0","mu1","q01","q10","ASlim0","ASlim1","lambda0/lambda1","mu0/mu1","q01/q10","lambda0/mu0","lambda1/mu1","lambda0/mu1","lambda1/mu0","lambda0/q01","lambda0/q10","lambda1/q01","lambda1/q10"   ,"NoExtinction.lambda0","NoExtinction.lambda1","NoExtinction.mu0","NoExtinction.mu1","NoExtinction.q01","NoExtinction.q10","NoExtinction.ASlim0","NoExtinction.ASlim1","NoExtinction.lambda0/lambda1","NoExtinction.mu0/mu1","NoExtinction.q01/q10","NoExtinction.lambda0/mu0","NoExtinction.lambda1/mu1","NoExtinction.lambda0/mu1","NoExtinction.lambda1/mu0","NoExtinction.lambda0/q01","NoExtinction.lambda0/q10","NoExtinction.lambda1/q01","NoExtinction.lambda1/q10")
        	ltreetyperes<-ltreetyperes.allstarters$all.tree.estimates[[starter]]
        	cat("Phylogenetic family: ");cat(Phylfam);
        	cat(", Virulence factor: ");cat(virfact);
        	cat(",  Strain family: ");cat(ltreetyperes.allstarters$best.model$strain.family);
        	cat(", Starter: ");cat(starter);
        	if (is.infinite(ltreetyperes$best.model$AICc)){
            	    cat("\n");cat("Tree ultrametrization or estimation procedure failed for this virulence factor starter");
        	}else{  
		    cat("  best model");cat("\n")
		    model.index<-which(sapply(model.list,function(lmod,estmodfunc){isTRUE(all.equal(lmod$func.par.transform,estmodfunc))},estmodfunc=ltreetyperes$best.model$parmodel$func.par.transform,simplify=TRUE))
		    cat("Model description : ");
		    cat(model.list[[model.index]]$model.desc);cat("\n")
	
		    cat("Best estimated parameters :");cat("\n")
		    vpars<-ltreetyperes$best.model$model.par.est
		    names(vpars)<-varnames
		    cat(names(vpars));cat("\n")
		    cat(vpars);cat("\n")
	
		    cat("Summary statistics of parameters :");cat("\n")
		    vsumstats<-fsumstats(vpars)
		    cat(names(vsumstats));cat("\n")
		    cat(vsumstats);cat("\n");cat("\n")
		    res<-c(paste(virfact,".",Phylfam,".",ltreetyperes$best.model$strain.family,".",starter,sep=""),vsumstats)	
				
	    	    noextinction.best.model<-fbest.no.death(ltreetyperes,varnames,TRUE)
            	    cat("\n");cat("\n")
            	    cat("Best model when no extinction assumed, virulence factor: ");cat(virfact);cat(",  strain family: ");cat(noextinction.best.model$strain.family);
            	    cat(", starter: ");cat(noextinction.best.model$starter);cat(" :\n")
            	    if (!is.null(noextinction.best.model$model.par.est)){				
			cat("\n");cat("\n")
			names(noextinction.best.model$model.par.est)<-varnames
		        model.index<-which(sapply(model.list,function(lmod,estmodfunc){isTRUE(all.equal(lmod$func.par.transform,estmodfunc))},estmodfunc=noextinction.best.model$parmodel$func.par.transform,simplify=TRUE))
			cat("Model description : ");
			cat(model.list[[model.index]]$model.desc);cat("\n")
	
			cat("Best estimated parameters :");cat("\n")
			vpars<-noextinction.best.model$model.par.est
			names(vpars)<-varnames
			cat(names(vpars));cat("\n")
			cat(vpars);cat("\n")
	
			cat("Summary statistics of parameters :");cat("\n")
			vsumstats<-fsumstats(vpars,TRUE)
			cat(names(vsumstats));cat("\n")
			cat(vsumstats);
		     }else{
                	cat("Model with no extinction failed")
                	names(vsumstats)<-paste("NoExtinction.",names(vsumstats),sep="")
                	vsumstats<-rep(NA,length(vsumstats))
                    }
		    res<-c(res,vsumstats)		
		}
		res.str[[starter.id]]<-res
        	starter.id<-starter.id+1
        	cat("\n");cat("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");cat("\n")
	    }
	    cat("\n");cat("*****************************************");cat("\n");cat("\n")
	    res.str
	},varnames=varnames,virfact=virfact,model.list=model.list,fsumstats=fsumstats,numofstats=numofstats,Phylfam=Phylfam,simplify=FALSE,USE.NAMES = FALSE)
	cat("========================================================================");cat("\n");cat("\n");cat("\n");cat("\n")
	allphylfam.res[[phylfam.id]]<-res.fam
	phylfam.id<-phylfam.id+1
    }
    allphylfam.res
},varnames=varnames,model.list=model.list,fsumstats=fsumstats,numofstats=numofstats,vphylfam=vphylfam,simplify=FALSE)
sink()


## create a nice table with the summary statistics to be written out to a csv file

res.all.org<-res.all
mEcoliBiSSEres<-rep(NA,numofstats)
n<-1
res.all<-unlist(res.all,recursive=FALSE)
res.all<-unlist(res.all,recursive=FALSE)
res.all<-unlist(res.all,recursive=FALSE)
names(res.all)<-NULL



for (i in 1:length(res.all)){
    k<-1 
    mests<-as.numeric(res.all[[i]][-1]) 
    mEcoliBiSSEres<-rbind(mEcoliBiSSEres,mests)
    rownames(mEcoliBiSSEres)[(n+1):(n+k)]<-res.all[[i]][1]
    n<-n+k
}

colnames(mEcoliBiSSEres)<-names(res.all[[1]])[-1]
mEcoliBiSSEres<-mEcoliBiSSEres[-1,]

## save the parameter summaries to a csv file
write.table(mEcoliBiSSEres,file=coutestcsv,sep=";", col.names=NA,quote=FALSE)
