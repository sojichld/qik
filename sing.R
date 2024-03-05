getwd()
R.Version()
library("MCMC.qpcr")
library("dplyr")
library("ggplot2")


HPDsummary <-
function(model,data,xgroup=NULL,genes=NA,relative=FALSE,
log.base=2,summ.plot=TRUE,ptype="z",...) {
	
#model=msr3;data=qs;xgroup=NULL;genes=NA;relative=FALSE;log.base=2;summ.plot=TRUE;ptype="z"

	mm=model;base=log.base
	gene.results=list()
	trms=attr(terms(mm$Fixed$formula),"term.labels")[grep("gene:",attr(terms(mm$Fixed$formula),"term.labels"))]
	trms=sub("gene:","",trms)
	sols=colnames(mm$Sol)
	if (is.na(genes[1])) { 
		genes=sub("gene","",sols[grep("gene.*$",sols)])
		genes=unique(sub(":.*","",genes))	
	}
	facts =list()
	for (t in trms) {
		if (!grepl(":",t)) { facts =append(facts,list(levels(data[,t])))}
	}
	names(facts)=trms[1:length(facts)]
	nfactors=length(facts)
	if (nfactors>2) { 
		stop("not implemented for more than 2 crossed factors")
	}
	gsols=c();fac1=c();fac2=c();samps=c();skips=c()
	
	globalInt=rep(0,length(mm$Sol[,1]))	
	if(sols[1]=="(Intercept)") { 
		globalInt=rep(mean(mm$Sol[,"(Intercept)"]),length(mm$Sol[,1]))
	} 
	
	interaction=0
	if (nfactors==2) {
		sol=paste("gene",genes[2],":",names(facts)[1],facts[[1]][2],":",names(facts)[2],facts[[2]][2],sep="")
		if (sol %in% sols) {
			interaction=1
		}
	}
	
	for (gene in genes) {
		sol=paste("gene",gene,sep="")
		if (!(sol %in% sols)) colnames(mm$Sol)[1]=sol
	}
	sols=colnames(mm$Sol)	

	for (gene in genes) {
		fac1g=c();fac2g=c();sampsg=c();skip=FALSE
		for (lev1 in 1:length(facts[[1]])) {
			if (nfactors==2) {
				for (lev2 in 1:length(facts[[2]])) {
					if (lev1==1 & lev2==1) { 
						sol=paste("gene",gene,sep="")
						if(sum(grep(sol,sols))==0) {
							skip=TRUE
							next
						}
						samp=(globalInt+mm$Sol[,sol])*as.numeric(!relative)
						int0=samp
					} else {
						if (lev2==1) { 
							sol=paste("gene",gene,":",names(facts)[1],facts[[1]][lev1],sep="")
							if(sum(grep(sol,sols))==0) {
								skip=TRUE
								next
							}
							
							samp=int0+mm$Sol[,sol]
							int2=mm$Sol[,sol]
						} else {
							if (lev1==1) { 
								sol=paste("gene",gene,":",names(facts)[2],facts[[2]][lev2],sep="") 
								if(sum(grep(sol,sols))==0) {
									skip=TRUE
									next
								}
								samp=int0+mm$Sol[,sol]
								int1=mm$Sol[,sol]
							} else {
								sol=paste("gene",gene,":",
								names(facts)[1],facts[[1]][lev1],
								":",names(facts)[2],facts[[2]][lev2],sep="") 
								if(sum(grep(sol,sols))==0 & interaction==1) {
									skip=TRUE
									next
								}
								if (interaction==1) {
									samp=int0+int1+int2+mm$Sol[,sol]
								} else { 
									samp=int0+int1+int2
								}
							}
						}
					}
		#			print(paste(lev1,lev2,sol))
					if (skip) { 
#						genes=genes[!(genes %in% gene)]						
						next 
					}
					gsols=append(gsols,gene)
					fac1g=append(fac1g,facts[[1]][lev1])
					fac2g=append(fac2g,facts[[2]][lev2])
					sampsg=cbind(sampsg,samp)
				}
			} else {
				if (lev1==1) { 
						sol=paste("gene",gene,sep="") 
						if(sum(grep(sol,sols))==0) {
							skip=TRUE
							next
						}
						samp=(globalInt+mm$Sol[,sol])*as.numeric(!relative)
						int0=samp
				} else {
					sol=paste("gene",gene,":",
					names(facts)[1],facts[[1]][lev1],sep="")
					if(sum(grep(sol,sols))==0) {
						skip=TRUE
						next
					}
					samp=int0+mm$Sol[,sol]
				}
				if (skip) { 
#					genes=genes[!(genes %in% gene)]
					next 
				}
				gsols=append(gsols,gene)
				fac1g=append(fac1g,facts[[1]][lev1])
				sampsg=cbind(sampsg,samp)
			}
		}
		if (skip) { 
			skips=append(skips,gene)
			next 
		}
		fac1=append(fac1,fac1g)
		samps=cbind(samps,sampsg)
		if (nfactors==2) { fac2=append(fac2,fac2g) }
		sampsg=data.frame(sampsg)
		if (nfactors==2) {
			gres=matrix(nrow=length(fac1g),ncol=length(fac1g),
			dimnames=list("pvalue"=paste(names(facts)[1],fac1g,":",
			names(facts)[2],fac2g,sep=""),
			"difference"=paste(names(facts)[1],fac1g,":",
			names(facts)[2],fac2g,sep="")))
		} else { 
			gres=matrix(nrow=length(fac1g),ncol=length(fac1g),
			dimnames=list("pvalue"=fac1g,"difference"=fac1g))
		}
		for (i in 1:(length(fac1g)-1)) {
			for (j in (i+1):length(fac1g)) {
				diff=sampsg[,j]-sampsg[,i]
				gres[j,i]=mcmc.pval(diff,ptype=ptype)
				gres[i,j]=mean(diff)/log(base)
			}	
		}		
		gene.results=append(gene.results,list(gres))
		names(gene.results)[length(gene.results)]=gene		
	}

	genes=genes[!(genes %in% skips)]
	for (ge in skips) { 
		gsols=gsols[-c(grep(ge,gsols))]		
	}
	big.summary=data.frame(cbind("gene"=gsols,"f1"=fac1))
	names(big.summary)[2]=names(facts[1])
	if(nfactors==2) {
		big.summary$f2=fac2
		names(big.summary)[3]=names(facts[2])
	}
	samps=apply(samps,2, function(x){ return(x/log(base)) })
	mns=apply(samps,2,mean)
	sds=apply(samps,2,sd)
	lower=apply(samps,2,function(x) { return(quantile(x,0.05)) })
	upper=apply(samps,2,function(x) { return(quantile(x,0.95)) })
	big.summary=cbind(big.summary,"mean"=mns,
	"sd"=sds,"lower"=lower,"upper"=upper)
	ploo=NULL

	if(relative) {
		if (nfactors==2) { 
			remov=which(
				big.summary[,2]==facts[[1]][1] & big.summary[,3]==facts[[2]][1]
			)
			big.summary=big.summary[-remov,] 
		} else { 	
			big.summary=big.summary[big.summary[,2]!=facts[[1]][1],] 
			}
	}
	if (is.null(xgroup)) {
		xgroup=names(facts[1])
		if(nfactors==2) { facet=names(facts[2]) }
	} else {
		if (xgroup==names(facts[1])) {
			if(nfactors==2) { facet=names(facts[2]) }
		} else {
			facet=names(facts[1])
		}
	}

	if(summ.plot) { 
		if (!relative) { 
			if (nfactors==2) {
				ploo=summaryPlot(big.summary,xgroup=xgroup,
				facet=facet,type="line",log.base=log.base,...) 
			} else {
				ploo=summaryPlot(big.summary,xgroup=names(facts[1]),type="line",log.base=log.base,...) 
			}
		} else {
			if (nfactors==2) {
				ploo=summaryPlot(big.summary,xgroup=xgroup,
				facet=facet,type="bar",log.base=log.base,...) 
			} else {
				ploo=summaryPlot(big.summary,xgroup=names(facts[1]),type="bar",log.base=log.base,...) 
			}
		}
	}

	return(list("summary"=big.summary,"geneWise"=gene.results,"ggPlot"=ploo))
}



getwd()
myc2 <- read.csv("Cq/myc2_ts.csv")
ubc1 <- read.csv("Cq/ubc2_ts.csv")
tuba2<- read.csv("Cq/tuba2_ts.csv")
tubb1<- read.csv("Cq/tubb1_ts.csv")
myc2 <- within(myc2, Biological.Set.Name[Content == 'NTC'] <- "Control")
ubc1 <- within(ubc1, Biological.Set.Name[Content == 'NTC'] <- "Control")
tuba2 <- within(tuba2, Biological.Set.Name[Content == 'NTC'] <- "Control")
tubb1 <- within(tubb1, Biological.Set.Name[Content == 'NTC'] <- "Control")
mergeda <- merge(x=myc2, y=ubc1, by = "Well")
merged1 <- data.frame(mergeda["Well"],mergeda["Biological.Set.Name.x"],mergeda["Cq.x"],mergeda["Cq.y"])
colnames(merged1) <- c('sample', 'Time', "Myc2", "Ubc2")
mergedb <- merge(x=tubb1, y=tuba2, by = "Well")
merged2 <- data.frame(mergedb["Well"],mergedb["Biological.Set.Name.x"],mergedb["Cq.x"],mergedb["Cq.y"])
colnames(merged2) <- c('sample', 'Time', "Tubb1", "Tuba2")
mergedc <- merge(x=merged1, y=merged2, by = "sample")
df <- mergedc[c("sample", "Time.x", "Myc2", "Ubc2", "Tuba2","Tubb1")]
df <- df[order(df$Time.x),]
df <- df %>% mutate_all(~ifelse(is.nan(.), NA, .))
get_avg <- function(a_vec){
    vec <- c()
    for (i in seq(1,length(a_vec), by=2)){
        item <- (a_vec[i] + a_vec[i+1])/2
        vec <-c(vec,item)
    }
    return(vec)
}
#TUBA2
dilutTB2 <- c(0,-1,-2,-3,-5)
stdTB2 <- c(7,8,9,10,12)
valTB2 <- c(21.76,21.76,25.56,25.53, 29.89, 28.78, 34.35, 32.70,38.44,38.44)
#TUBB1
dilutTA1 <- c(-1,-2,-3,-4)
stdTA1 <- c(1,2,3,4)
valTA1 <- c(26.59,26.80,30.10,29.93,33.36,33.18,37.35,37.45)
#UBC2
dilutBC2 <- c(-1,-2,-3,-4)
stdBC2 <- c(1,2,3,4)
valBC2 <- c(22.83,23.48,26.24,26.01,29.36,29.41,33.54,33.41)

dilutions <- list(dilutBC2,dilutTB2,dilutTA1)
std <- list(stdBC2,stdTB2,stdTA1)
val <- list(valBC2,valTB2,valTA1)
names <-c("Tuba2","Tubb1","Ubc2")


for (i in seq(1,3)){
    v <- get_avg(val[[i]])
    typeof(v)
    temp <- data.frame(dilutions[[i]],std[[i]],v)
    colnames(temp) <- c("dilut", "std", "v")
    ax <- unlist(temp["dilut"]) 
    ay <- unlist(temp["v"])
    bf <- lm(ay ~ ax)
    cf <- coef(bf)
    Intercept <- cf[1]
    Slope <- cf[2]
    E=10^(-1/Slope)
    temp
    plot(ax, ay, pch = 16, cex = 2.3, col = "green", main = names[i], xlab="Log10Concentration",ylab="Cq") 
    abline(Intercept,Slope, col="green", lwd=3)
    print(paste0("Slope:",Slope))
    print(paste0("E:",E))

}

#    E = 10^(-1/slope)
gene <- c("Myc2", "Ubc2", "Tubb1", "Tuba2")
efficiency <- c(1.95, 1.96, 1.91, 1.97)
eff <- data.frame(gene, efficiency)
colnames(eff) <- c("gene", "efficiency")
show(eff)


qs=cq2counts(data=df, effic=eff, genecols=c(3:6),condcols=c(1:2))

dl=cq2log(data=df,genecols=(3:6),condcols=c(1:2),effic=eff)

classic=mcmc.qpcr.classic(fixed="Time.x",
data=dl,
controls=c("Ubc2","Tubb1","Tuba2"),
pr=T,
pl=T
)
diagnostic.mcmc(model = classic,
col="grey50",
cex=0.8)

summary(classic)


plot(classic)

s1=HPDsummary(model=classic,data=dl)


s0=HPDsummary(model=classic,data=dl,relative=TRUE)

HPDsummary <-
function(model,data,xgroup=NULL,genes=NA,relative=FALSE,
log.base=2,summ.plot=TRUE,ptype="z",...) {
	
#model=msr3;data=qs;xgroup=NULL;genes=NA;relative=FALSE;log.base=2;summ.plot=TRUE;ptype="z"

	mm=model;base=log.base
	gene.results=list()
	trms=attr(terms(mm$Fixed$formula),"term.labels")[grep("gene:",attr(terms(mm$Fixed$formula),"term.labels"))]
	trms=sub("gene:","",trms)
	sols=colnames(mm$Sol)
	if (is.na(genes[1])) { 
		genes=sub("gene","",sols[grep("gene.*$",sols)])
		genes=unique(sub(":.*","",genes))	
	}
	facts =list()
	for (t in trms) {
		if (!grepl(":",t)) { facts =append(facts,list(levels(data[,t])))}
	}
	names(facts)=trms[1:length(facts)]
	nfactors=length(facts)
	if (nfactors>2) { 
		stop("not implemented for more than 2 crossed factors")
	}
	gsols=c();fac1=c();fac2=c();samps=c();skips=c()
	
	globalInt=rep(0,length(mm$Sol[,1]))	
	if(sols[1]=="(Intercept)") { 
		globalInt=rep(mean(mm$Sol[,"(Intercept)"]),length(mm$Sol[,1]))
	} 
	
	interaction=0
	if (nfactors==2) {
		sol=paste("gene",genes[2],":",names(facts)[1],facts[[1]][2],":",names(facts)[2],facts[[2]][2],sep="")
		if (sol %in% sols) {
			interaction=1
		}
	}
	
	for (gene in genes) {
		sol=paste("gene",gene,sep="")
		if (!(sol %in% sols)) colnames(mm$Sol)[1]=sol
	}
	sols=colnames(mm$Sol)	

	for (gene in genes) {
		fac1g=c();fac2g=c();sampsg=c();skip=FALSE
		for (lev1 in 1:length(facts[[1]])) {
			if (nfactors==2) {
				for (lev2 in 1:length(facts[[2]])) {
					if (lev1==1 & lev2==1) { 
						sol=paste("gene",gene,sep="")
						if(sum(grep(sol,sols))==0) {
							skip=TRUE
							next
						}
						samp=(globalInt+mm$Sol[,sol])*as.numeric(!relative)
						print("check0")
						print(samp)
						int0=samp
					} else {
						if (lev2==1) { 
							sol=paste("gene",gene,":",names(facts)[1],facts[[1]][lev1],sep="")
							if(sum(grep(sol,sols))==0) {
								skip=TRUE
								next
							}
							
							samp=int0+mm$Sol[,sol]
							print(samp)
							int2=mm$Sol[,sol]
						} else {
							if (lev1==1) { 
								sol=paste("gene",gene,":",names(facts)[2],facts[[2]][lev2],sep="") 
								if(sum(grep(sol,sols))==0) {
									skip=TRUE
									next
								}
								samp=int0+mm$Sol[,sol]
								int1=mm$Sol[,sol]
							} else {
								sol=paste("gene",gene,":",
								names(facts)[1],facts[[1]][lev1],
								":",names(facts)[2],facts[[2]][lev2],sep="") 
								if(sum(grep(sol,sols))==0 & interaction==1) {
									skip=TRUE
									next
								}
								if (interaction==1) {
									samp=int0+int1+int2+mm$Sol[,sol]
								} else { 
									samp=int0+int1+int2
								}
							}
						}
					}
		#			print(paste(lev1,lev2,sol))
					if (skip) { 
#						genes=genes[!(genes %in% gene)]						
						next 
					}
					gsols=append(gsols,gene)
					fac1g=append(fac1g,facts[[1]][lev1])
					fac2g=append(fac2g,facts[[2]][lev2])
					sampsg=cbind(sampsg,samp)
				}
			} else {
				if (lev1==1) { 
						sol=paste("gene",gene,sep="") 
						if(sum(grep(sol,sols))==0) {
							skip=TRUE
							next
						}
						samp=(globalInt+mm$Sol[,sol])*as.numeric(!relative)
						print("check1")
						print(samp)
						int0=samp
				} else {
					sol=paste("gene",gene,":",
					names(facts)[1],facts[[1]][lev1],sep="")
					if(sum(grep(sol,sols))==0) {
						skip=TRUE
						next
					}
					samp=int0+mm$Sol[,sol]
					print("check2")
					print(samp)
				}
				if (skip) { 
#					genes=genes[!(genes %in% gene)]
					next 
				}
				gsols=append(gsols,gene)
				fac1g=append(fac1g,facts[[1]][lev1])
				sampsg=cbind(sampsg,samp)
			}
		}
		if (skip) { 
			skips=append(skips,gene)
			next 
		}
		fac1=append(fac1,fac1g)
		samps=cbind(samps,sampsg)
		if (nfactors==2) { fac2=append(fac2,fac2g) }
		sampsg=data.frame(sampsg)
		if (nfactors==2) {
			gres=matrix(nrow=length(fac1g),ncol=length(fac1g),
			dimnames=list("pvalue"=paste(names(facts)[1],fac1g,":",
			names(facts)[2],fac2g,sep=""),
			"difference"=paste(names(facts)[1],fac1g,":",
			names(facts)[2],fac2g,sep="")))
		} else { 
			gres=matrix(nrow=length(fac1g),ncol=length(fac1g),
			dimnames=list("pvalue"=fac1g,"difference"=fac1g))
		}
		for (i in 1:(length(fac1g)-1)) {
			for (j in (i+1):length(fac1g)) {
				diff=sampsg[,j]-sampsg[,i]
				gres[j,i]=mcmc.pval(diff,ptype=ptype)
				gres[i,j]=mean(diff)/log(base)
			}	
		}		
		gene.results=append(gene.results,list(gres))
		names(gene.results)[length(gene.results)]=gene		
	}

	genes=genes[!(genes %in% skips)]
	for (ge in skips) { 
		gsols=gsols[-c(grep(ge,gsols))]		
	}
	big.summary=data.frame(cbind("gene"=gsols,"f1"=fac1))
	names(big.summary)[2]=names(facts[1])
	if(nfactors==2) {
		big.summary$f2=fac2
		names(big.summary)[3]=names(facts[2])
	}
		print("check3")
		print(samps)
	samps=apply(samps,2, function(x){ return(x/log(base)) })
	mns=apply(samps,2,mean)
	sds=apply(samps,2,sd)
	lower=apply(samps,2,function(x) { return(quantile(x,0.05)) })
	upper=apply(samps,2,function(x) { return(quantile(x,0.95)) })
	big.summary=cbind(big.summary,"mean"=mns,
	"sd"=sds,"lower"=lower,"upper"=upper)
	ploo=NULL

	if(relative) {
		if (nfactors==2) { 
			remov=which(
				big.summary[,2]==facts[[1]][1] & big.summary[,3]==facts[[2]][1]
			)
			big.summary=big.summary[-remov,] 
		} else { 	
			big.summary=big.summary[big.summary[,2]!=facts[[1]][1],] 
			}
	}
	if (is.null(xgroup)) {
		xgroup=names(facts[1])
		if(nfactors==2) { facet=names(facts[2]) }
	} else {
		if (xgroup==names(facts[1])) {
			if(nfactors==2) { facet=names(facts[2]) }
		} else {
			facet=names(facts[1])
		}
	}

	if(summ.plot) { 
		if (!relative) { 
			if (nfactors==2) {
				ploo=summaryPlot(big.summary,xgroup=xgroup,
				facet=facet,type="line",log.base=log.base,...) 
			} else {
				ploo=summaryPlot(big.summary,xgroup=names(facts[1]),type="line",log.base=log.base,...) 
			}
		} else {
			if (nfactors==2) {
				ploo=summaryPlot(big.summary,xgroup=xgroup,
				facet=facet,type="bar",log.base=log.base,...) 
			} else {
				ploo=summaryPlot(big.summary,xgroup=names(facts[1]),type="bar",log.base=log.base,...) 
			}
		}
	}

	return(list("summary"=big.summary,"geneWise"=gene.results,"ggPlot"=ploo))
}


data(beckham.data)
data(beckham.eff)

qs$treatment.time=as.factor(paste(qs$tr,qs$time,sep="."))


qs$treatment.time=relevel(qs$treatment.time,ref="control.0h")

naive=mcmc.qpcr(data=qs, fixed="treatment.time")


s1=HPDsummary(model=naive,data=qs)

#mm=naive
#data=qs
mm=classic
data=dl
genes=NA
relative=TRUE
summ.plot=TRUE
ptype="z"
base=2
	gene.results=list()
	trms=attr(terms(mm$Fixed$formula),"term.labels")[grep("gene:",attr(terms(mm$Fixed$formula),"term.labels"))]
	trms=sub("gene:","",trms)
	sols=colnames(mm$Sol)
	if (is.na(genes[1])) { 
		genes=sub("gene","",sols[grep("gene.*$",sols)])
		genes=unique(sub(":.*","",genes))	
	}
	facts =list()
	for (t in trms) {
		if (!grepl(":",t)) { facts =append(facts,list(levels(data[,t])))}
	}
	names(facts)=trms[1:length(facts)]
	nfactors=length(facts)
	if (nfactors>2) { 
		stop("not implemented for more than 2 crossed factors")
	}
	gsols=c();fac1=c();fac2=c();samps=c();skips=c()
	
	globalInt=rep(0,length(mm$Sol[,1]))	
	if(sols[1]=="(Intercept)") { 
		globalInt=rep(mean(mm$Sol[,"(Intercept)"]),length(mm$Sol[,1]))
	} 
	
	interaction=0
	if (nfactors==2) {
		sol=paste("gene",genes[2],":",names(facts)[1],facts[[1]][2],":",names(facts)[2],facts[[2]][2],sep="")
		if (sol %in% sols) {
			interaction=1
		}
	}


    for (gene in genes) {
		sol=paste("gene",gene,sep="")
		if (!(sol %in% sols)) colnames(mm$Sol)[1]=sol
	}
	sols=colnames(mm$Sol)	


    install.packages("vscDebugger")


    for (gene in genes) {
		fac1g=c();fac2g=c();sampsg=c();skip=FALSE
		for (lev1 in 1:length(facts[[1]])) {
			if (nfactors==2) {
				for (lev2 in 1:length(facts[[2]])) {
					if (lev1==1 & lev2==1) { 
						sol=paste("gene",gene,sep="")
						if(sum(grep(sol,sols))==0) {
							skip=TRUE
							next
						}
						samp=(globalInt+mm$Sol[,sol])*as.numeric(!relative)
						typeof(samp)
						int0=samp
					} else {
						if (lev2==1) { 
							sol=paste("gene",gene,":",names(facts)[1],facts[[1]][lev1],sep="")
							if(sum(grep(sol,sols))==0) {
								skip=TRUE
								next
							}
							
							samp=int0+mm$Sol[,sol]
							typeof(samps)
							int2=mm$Sol[,sol]
						} else {
							if (lev1==1) { 
								sol=paste("gene",gene,":",names(facts)[2],facts[[2]][lev2],sep="") 
								if(sum(grep(sol,sols))==0) {
									skip=TRUE
									next
								}
								samp=int0+mm$Sol[,sol]
								int1=mm$Sol[,sol]
							} else {
								sol=paste("gene",gene,":",
								names(facts)[1],facts[[1]][lev1],
								":",names(facts)[2],facts[[2]][lev2],sep="") 
								if(sum(grep(sol,sols))==0 & interaction==1) {
									skip=TRUE
									next
								}
								if (interaction==1) {
									samp=int0+int1+int2+mm$Sol[,sol]
									print(typeof(samps))
								} else { 
									samp=int0+int1+int2
									print(typeof(samps))
								}
							}
						}
					}
		#			print(paste(lev1,lev2,sol))
					if (skip) { 
#						genes=genes[!(genes %in% gene)]						
						next 
					}
					gsols=append(gsols,gene)
					fac1g=append(fac1g,facts[[1]][lev1])
					fac2g=append(fac2g,facts[[2]][lev2])
					sampsg=cbind(sampsg,samp)
					print(typeof(samps))
				}
			} else {
				if (lev1==1) { 
						sol=paste("gene",gene,sep="") 
						if(sum(grep(sol,sols))==0) {
							skip=TRUE
							next
						}
						samp=(globalInt+mm$Sol[,sol])*as.numeric(!relative)
						typeof(samps)
						int0=samp
						print(typeof(samps))
				} else {
					sol=paste("gene",gene,":",
					names(facts)[1],facts[[1]][lev1],sep="")
					if(sum(grep(sol,sols))==0) {
						skip=TRUE
						next
					print(typeof(samps))
					}
					samp=int0+mm$Sol[,sol]
					print(typeof(samps))
				}
				if (skip) { 
#					genes=genes[!(genes %in% gene)]
					next 
				}
				gsols=append(gsols,gene)
				fac1g=append(fac1g,facts[[1]][lev1])
				sampsg=cbind(sampsg,samp)
				print(typeof(samps))
			}
		}
		if (skip) { 
			skips=append(skips,gene)
			next 
		}
		fac1=append(fac1,fac1g)
		samps=cbind(samps,sampsg)
		print("oh")
		print(typeof(samps))
		if (nfactors==2) { fac2=append(fac2,fac2g) }
		sampsg=data.frame(sampsg)
		if (nfactors==2) {
			gres=matrix(nrow=length(fac1g),ncol=length(fac1g),
			dimnames=list("pvalue"=paste(names(facts)[1],fac1g,":",
			names(facts)[2],fac2g,sep=""),
			"difference"=paste(names(facts)[1],fac1g,":",
			names(facts)[2],fac2g,sep="")))
			print(typeof(samps))
		} else { 
			gres=matrix(nrow=length(fac1g),ncol=length(fac1g),
			dimnames=list("pvalue"=fac1g,"difference"=fac1g))
		}
		for (i in 1:(length(fac1g)-1)) {
			for (j in (i+1):length(fac1g)) {
				diff=sampsg[,j]-sampsg[,i]
				gres[j,i]=mcmc.pval(diff,ptype=ptype)
				gres[i,j]=mean(diff)/log(base)
				print(typeof(samps))
			}	
		}		
		print(typeof(samps))
		print("ohn")
		gene.results=append(gene.results,list(gres))
		names(gene.results)[length(gene.results)]=gene		
		print(typeof(samps))
	}
			print(typeof(samps))


            print(typeof(samps))


            genes=genes[!(genes %in% skips)]
	for (ge in skips) { 
		gsols=gsols[-c(grep(ge,gsols))]		
	}
	big.summary=data.frame(cbind("gene"=gsols,"f1"=fac1))
	names(big.summary)[2]=names(facts[1])
	if(nfactors==2) {
		big.summary$f2=fac2
		names(big.summary)[3]=names(facts[2])
	}
		print("check3")
		print(samps)
	samps=apply(samps,2, function(x){ return(x/log(base)) })
	mns=apply(samps,2,mean)
	sds=apply(samps,2,sd)
	lower=apply(samps,2,function(x) { return(quantile(x,0.05)) })
	upper=apply(samps,2,function(x) { return(quantile(x,0.95)) })
	big.summary=cbind(big.summary,"mean"=mns,
	"sd"=sds,"lower"=lower,"upper"=upper)
	ploo=NULL

    if(relative) {
		if (nfactors==2) { 
			remov=which(
				big.summary[,2]==facts[[1]][1] & big.summary[,3]==facts[[2]][1]
			)
			big.summary=big.summary[-remov,] 
		} else { 	
			big.summary=big.summary[big.summary[,2]!=facts[[1]][1],] 
			}
	}
	if (is.null(xgroup)) {
		xgroup=names(facts[1])
		if(nfactors==2) { facet=names(facts[2]) }
	} else {
		if (xgroup==names(facts[1])) {
			if(nfactors==2) { facet=names(facts[2]) }
		} else {
			facet=names(facts[1])
		}
	}

    if(summ.plot) { 
		if (!relative) { 
			if (nfactors==2) {
				ploo=summaryPlot(big.summary,xgroup=xgroup,
				facet=facet,type="line",log.base=log.base,...) 
			} else {
				ploo=summaryPlot(big.summary,xgroup=names(facts[1]),type="line",log.base=log.base,...) 
			}
		} else {
			if (nfactors==2) {
				ploo=summaryPlot(big.summary,xgroup=xgroup,
				facet=facet,type="bar",log.base=log.base,...) 
			} else {
				ploo=summaryPlot(big.summary,xgroup=names(facts[1]),type="bar",log.base=log.base,...) 
			}
		}
	}

    return(list("summary"=big.summary,"geneWise"=gene.results,"ggPlot"=ploo))