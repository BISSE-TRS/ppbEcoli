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

## the code to run the BiSSE analysis using the wrapper around make.bisse()

library(diversitree)
library(parallel)
## load the file with the strains data
## the data has to be a matrix called strainsData
## each row corresponds to a strain, the rows have to be named with rows
## corresponding to the rows on the phylogeny
## the columns correspond to the various binary traits, each column has to be named
## the first columns describe the phylogenetic families to which the strains belong to
## each strain belongs to only one family
## the next column the actual pathogenicity traits
load("strainsData.RData") 

## The phylogenetic trees will be read in later
## the files containing the trees have to be called
## CGG+GTG.txt, GTG.txt, CGG.txt ## all strains, different files for different starters
## K_CGG+GTG.txt, K_GTG.txt, K_CGG.txt ## only K strains
## U_CGG+GTG.txt, U_GTG.txt, U_CGG.txt ## only U straints

## here we can lump together some phylogenetic families if we analyze
## them separately
vphylfam<-c("A","B1+C","B2+F","D+E")

## choose the traits to analyze
vbinvars<-c("hly1","iroN","iutA","fyuA","FimG","papC","sat","sfa","usp","tsh","escV","cnf","astA")
## removed aggR,pic,stx1,stx2 as they caused errors

## strain families
vstrain.family<-c("UK","U","K")
## different starters for constructing the dendrograms
vstarter<-c("CGG","GTG","CGGGTG")


## setup the various rate combinations we are interested in
## i.e. set some to be 0, others to be equal to each other
model.list<-list(
list(func.par.transform=function(pars){pars},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,mu1=0.1,q01=0.1,q10=0.1)),
list(func.par.transform=function(pars){c(pars,pars[5])},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,mu1=0.1,q01=0.1)),
list(func.par.transform=function(pars){c(pars[1:3],pars[3],pars[4],pars[4])},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,q01=0.1)),
list(func.par.transform=function(pars){c(pars[1],pars[1],pars[2:4],pars[4])},par0=c(lambda0=0.1,mu0=0.1,mu1=0.1,q01=0.1)),
list(func.par.transform=function(pars){c(rbind(pars,pars))},par0=c(lambda0=0.1,mu0=0.1,q01=0.1)),
list(func.par.transform=function(pars){c(pars[1:3],pars[3],pars[4:5])},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,q01=0.1,q10=0.1)),
list(func.par.transform=function(pars){c(pars[1],pars[1],pars[2:5])},par0=c(lambda0=0.1,mu0=0.1,mu1=0.1,q01=0.1,q10=0.1)),
list(func.par.transform=function(pars){c(pars[1],pars[1],pars[2],pars[2],pars[3:4])},par0=c(lambda0=0.1,mu0=0.1,q01=0.1,q10=0.1)),

list(func.par.transform=function(pars){c(pars[1:2],0,pars[3:5])},par0=c(lambda0=0.1,lambda1=0.1,mu1=0.1,q01=0.1,q10=0.1)),
list(func.par.transform=function(pars){c(pars[1:2],0,pars[3:4],pars[4])},par0=c(lambda0=0.1,lambda1=0.1,mu1=0.1,q01=0.1,q10=0.1)),
list(func.par.transform=function(pars){c(pars[1],pars[1],0,pars[2:3],pars[3])},par0=c(lambda0=0.1,mu1=0.1,q01=0.1)),
list(func.par.transform=function(pars){c(pars[1],pars[1],0,pars[2:4])},par0=c(lambda0=0.1,mu1=0.1,q01=0.1,q10=0.1)),

list(func.par.transform=function(pars){c(pars[1:3],0,pars[4:5])},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,q01=0.1,q10=0.1)),
list(func.par.transform=function(pars){c(pars[1:3],0,pars[4],pars[4])},par0=c(lambda0=0.1,lambda1=0.1,mu0=0.1,q01=0.1,q10=0.1)),
list(func.par.transform=function(pars){c(pars[1],pars[1:2],0,pars[3],pars[3])},par0=c(lambda0=0.1,mu0=0.1,q01=0.1)),
list(func.par.transform=function(pars){c(pars[1],pars[1:2],0,pars[3:4])},par0=c(lambda0=0.1,mu0=0.1,q01=0.1,q10=0.1)),

list(func.par.transform=function(pars){c(pars[1:2],0,0,pars[3:4])},par0=c(lambda0=0.1,lambda1=0.1,q01=0.1,q10=0.1)),
list(func.par.transform=function(pars){c(pars[1:2],0,0,pars[3],pars[3])},par0=c(lambda0=0.1,lambda1=0.1,q01=0.1)),
list(func.par.transform=function(pars){c(pars[1],pars[1],0,0,pars[2],pars[2])},par0=c(lambda0=0.1,q01=0.1)),
list(func.par.transform=function(pars){c(pars[1],pars[1],0,0,pars[2:3])},par0=c(lambda0=0.1,q01=0.1,q10=0.1))
)

## read in the phylogenetic trees
apeCGGGTG.UK<-read.tree("CGG+GTG.txt")  
tree.heigth<-max(node.depth.edgelength(apeCGGGTG.UK),na.rm=TRUE)
apeCGGGTG.UK$edge.length<-apeCGGGTG.UK$edge.length/tree.heigth

apeCGG.UK<-read.tree("CGG.txt")  
tree.heigth<-max(node.depth.edgelength(apeCGG.UK),na.rm=TRUE)
apeCGG.UK$edge.length<-apeCGG.UK$edge.length/tree.heigth

apeGTG.UK<-read.tree("GTG.txt")  
tree.heigth<-max(node.depth.edgelength(apeGTG.UK),na.rm=TRUE)
apeGTG.UK$edge.length<-apeGTG.UK$edge.length/tree.heigth

apeCGGGTG.K<-read.tree("K_CGG+GTG.txt")  
tree.heigth<-max(node.depth.edgelength(apeCGGGTG.K),na.rm=TRUE)
apeCGGGTG.K$edge.length<-apeCGGGTG.K$edge.length/tree.heigth

apeCGG.K<-read.tree("K_CGG.txt")  
tree.heigth<-max(node.depth.edgelength(apeCGG.K),na.rm=TRUE)
apeCGG.K$edge.length<-apeCGG.K$edge.length/tree.heigth

apeGTG.K<-read.tree("K_GTG.txt")  
tree.heigth<-max(node.depth.edgelength(apeGTG.K),na.rm=TRUE)
apeGTG.K$edge.length<-apeGTG.K$edge.length/tree.heigth

apeCGGGTG.U<-read.tree("U_CGG+GTG.txt")  
tree.heigth<-max(node.depth.edgelength(apeCGGGTG.U),na.rm=TRUE)
apeCGGGTG.U$edge.length<-apeCGGGTG.U$edge.length/tree.heigth

apeCGG.U<-read.tree("U_CGG.txt")  
tree.heigth<-max(node.depth.edgelength(apeCGG.U),na.rm=TRUE)
apeCGG.U$edge.length<-apeCGG.U$edge.length/tree.heigth

apeGTG.U<-read.tree("U_GTG.txt")
tree.heigth<-max(node.depth.edgelength(apeGTG.U),na.rm=TRUE)
apeGTG.U$edge.length<-apeGTG.U$edge.length/tree.heigth



ltrees.org=list(
UK.CGG=apeCGG.UK,
UK.GTG=apeGTG.UK,
UK.CGGGTG=apeCGGGTG.UK,
U.CGG=apeCGG.U,
U.GTG=apeGTG.U,
U.CGGGTG=apeCGGGTG.U,
K.CGG=apeCGG.K,
K.GTG=apeGTG.K,
K.CGGGTG=apeCGGGTG.K
)

## run in parallel
if (!exists("cl")){cl <- makeCluster(getOption("cl.cores", 8),outfile="")}

## make the different trees for the different phylogenetic families
## and make them ultrametric by different methods
sink("outputFromWrapper.txt",append=TRUE)
ltrees.chr.relaxed<-list()
ltrees.chr.correlated<-list()
ltrees.chr.discrete<-list()
vnames<-c()
for (PhylFam in vphylfam){
    vstrainskeep<-c()
    ## extract the phylogenetic families
    vfamgrs<-strsplit(PhylFam,"+",fixed=TRUE)[[1]] 
    ## and only keep those strains from the appropriate phylogenetic family
    for(famgr in vfamgrs){vstrainskeep<-c(vstrainskeep,rownames(strainsData[which(strainsData[,famgr]==1),]))}
    vnames<-c(vnames,paste(PhylFam,names(ltrees.org),sep="."))
    ltrees.chr.relaxed<-c(ltrees.chr.relaxed,parSapply(cl=cl,1:length(ltrees.org),function(i,ltrees,method,tips.keep){source("make_chronostree.R");chronos.tree(ltrees[[i]],method,names(ltrees)[i],tips.keep=tips.keep)},ltrees=ltrees.org,method="relaxed",tips.keep=vstrainskeep,simplify=FALSE))
    ltrees.chr.correlated<-c(ltrees.chr.correlated,parSapply(cl=cl,1:length(ltrees.org),function(i,ltrees,method,tips.keep){source("make_chronostree.R");chronos.tree(ltrees[[i]],method,names(ltrees)[i],tips.keep=tips.keep)},ltrees=ltrees.org,method="correlated",tips.keep=vstrainskeep,simplify=FALSE))
    ltrees.chr.discrete<-c(ltrees.chr.discrete,parSapply(cl=cl,1:length(ltrees.org),function(i,ltrees,method,tips.keep){source("make_chronostree.R");chronos.tree(ltrees[[i]],method,names(ltrees)[i],tips.keep=tips.keep)},ltrees=ltrees.org,method="discrete",tips.keep=vstrainskeep,simplify=FALSE))
}
names(ltrees.chr.relaxed)<-vnames 
names(ltrees.chr.correlated)<-vnames 
names(ltrees.chr.discrete)<-vnames 
sink()
ltrees.ultras<-list("rel"=ltrees.chr.relaxed,"corr"=ltrees.chr.correlated,"disc"=ltrees.chr.discrete)

numtrees<-length(vnames)

## check if there were any errors during making trees ultrametric
## and choose all error-free trees

mNAtrees<-matrix(NA,ncol=4,nrow=numtrees)
colnames(mNAtrees)<-c("rel","corr","disc","numgood")
rownames(mNAtrees)<- vnames 
mNAtrees[,"rel"]<-as.numeric(sapply(ltrees.chr.relaxed,function(x){!is.na(x[1])},simplify=TRUE))
mNAtrees[,"corr"]<-as.numeric(sapply(ltrees.chr.correlated,function(x){!is.na(x[1])},simplify=TRUE))
mNAtrees[,"disc"]<-as.numeric(sapply(ltrees.chr.discrete,function(x){!is.na(x[1])},simplify=TRUE))
mNAtrees[,"numgood"]<-apply(mNAtrees,1,sum,na.rm=TRUE)
sink("outputFromWrapper.txt",append=TRUE);print(mNAtrees);sink()

vnumOKtrees<-which(mNAtrees[,"numgood"]!=0)

vnumnNAs<-nrow(mNAtrees)-c(sum(sapply(ltrees.chr.relaxed,function(x){is.na(x[1])},simplify=TRUE)),sum(sapply(ltrees.chr.correlated,function(x){is.na(x[1])},simplify=TRUE)),sum(sapply(ltrees.chr.discrete,function(x){is.na(x[1])},simplify=TRUE)))
names(vnumnNAs)<-c("rel","corr","disc")

ltrees<-vector("list",length(vnames))
names(ltrees)<-vnames
ltrees[which(mNAtrees[,"numgood"]==0)]<-NA
maxmodel<-order(vnumnNAs,decreasing=TRUE)

vultramethod<-rep(0,length(ltrees))
names(vultramethod)<-names(ltrees)

if(maxmodel[1]==length(vnumOKtrees)){
    ltrees[vnumOKtrees]<-ltrees.ultras[[maxmodel[1]]][vnumOKtrees]
    vultramethod[vnumOKtrees]<-maxmodel[1]
}else{
    vmaxnNAtrees<-which(mNAtrees[,maxmodel[1]]!=0)
    ltrees[vmaxnNAtrees]<-ltrees.ultras[[maxmodel[1]]][vmaxnNAtrees]
    vultramethod[vmaxnNAtrees]<-maxmodel[1]
    vnumOKtrees<-setdiff(vnumOKtrees,vmaxnNAtrees)
    if(maxmodel[2]==length(vnumOKtrees)){
	ltrees[vnumOKtrees]<-ltrees.ultras[[maxmodel[2]]][vnumOKtrees]
	vultramethod[vnumOKtrees]<-maxmodel[2]
    }else{
	vmaxnNAtrees<-which(mNAtrees[,maxmodel[2]]!=0)
	ltrees[vmaxnNAtrees]<-ltrees.ultras[[maxmodel[2]]][vmaxnNAtrees]
	vultramethod[vmaxnNAtrees]<-maxmodel[2]
	vnumOKtrees<-setdiff(vnumOKtrees,vmaxnNAtrees)
	if(maxmodel[3]==length(vnumOKtrees)){
	    ltrees[vnumOKtrees]<-ltrees.ultras[[maxmodel[3]]][vnumOKtrees]
	    vultramethod[vnumOKtrees]<-maxmodel[3]
	}else{
	    vmaxnNAtrees<-which(mNAtrees[,maxmodel[3]]!=0)
	    ltrees[vmaxnNAtrees]<-ltrees.ultras[[maxmodel[3]]][vmaxnNAtrees]
	    vultramethod[vmaxnNAtrees]<-maxmodel[3]
	}
    }
}


vultramethod<-as.character(vultramethod)
vultramethod[which(vultramethod=="1")]<-"rel"
vultramethod[which(vultramethod=="2")]<-"corr"
vultramethod[which(vultramethod=="3")]<-"disc"
vultramethod[which(vultramethod=="0")]<- NA


sink("outputFromWrapper.txt",append=TRUE);print(vultramethod);sink()



dfdata<-strainsData


## estimate the BiSSE model parameters for the different rate setups
## it is parallalized over the different traits, i.e. binary variables
parests<-parSapply(cl=cl,vbinvars,function(binvar,model.list,vstarter,vstrain.family,vphylfam,dfdata,ltrees){
	    source("BiSSE_wrapperFunction.R")	    
	    modelests.phylfam<-vector("list",length=length(vphylfam))
	    names(modelests.phylfam)<-vphylfam
	    u<-1
	    for (phylfam in vphylfam){ ## for each phylogenetic family
	    	modelests.strainfamily<-vector("list",length=length(vstrain.family))
		names(modelests.strainfamily)<-vstrain.family
		j<-1	
		for(strain.family in vstrain.family){## for each strain family
		    modelests.tree<-vector("list",length=length(vstarter))
		    names(modelests.tree)<-vstarter
		    k<-1
		    best.model<-list(AICc=Inf)
		    for(starter in vstarter){## for each starter, i.e. different tree
    			modelests<-vector("list",length=length(model.list))
			i<-1
			best.model.tree<-list(AICc=Inf)
			sink("tmpoutput.txt",append=TRUE);print(paste("Doing ",binvar," ",phylfam," ",starter," ",strain.family,sep=""));sink()
			for(parmodel in model.list){		
			## do the actual estimation and find the one with lowest AICc
			    modelests[[i]]<-BiSSE.model.analysis(parmodel,phylfam,starter,strain.family,dfdata,binvar,ltrees)
			    if (modelests[[i]]$AICc< best.model$AICc){best.model<-modelests[[i]]}
			    if (modelests[[i]]$AICc< best.model.tree$AICc){best.model.tree<-modelests[[i]]}
			    i<-i+1
			}
			modelests.tree[[k]]<-list(best.model=best.model,best.model.tree=best.model.tree,all.estimates= modelests)
			k<-k+1
		    }
		    modelests.strainfamily[[j]]<-list(best.model=best.model,all.tree.estimates= modelests.tree)
		    j<-j+1
		}		
		modelests.phylfam[[u]]<-modelests.strainfamily    
		u<-u+1
	    }
	    modelests.phylfam
},model.list=model.list,vstarter=vstarter,vstrain.family=vstrain.family,vphylfam=vphylfam,dfdata=dfdata,ltrees=ltrees,simplify=FALSE)
names(parests)<-vbinvars

## WARNING: the below file can be VERY large
save(parests,model.list,vstarter,vstrain.family,dfdata,ltrees,ltrees.org,mNAtrees,vultramethod,ltrees.ultras,strainsData,file="BiSSEestimates.RData")


