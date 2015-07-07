library(magrittr)
library(dplyr)

#infer haplotypes

###input 'altervcf' output. 
###estimate genetic distances between markers

#Lcm_vcf<-function(cmma,vcf1,vcf2,l1,l2,inva=FALSE,core="cc"){

#output posterior probability of distance given lnorm prior on mb/cm relationship
#bayesian version of Lcm_vcf that uses a lognormal prior on the relationship between physical and genetic distance. 

bay_Lcm_vcf<-function(cmma,vcf1,vcf2,l1,l2,inva=FALSE,core="cc",ml=log(.5),sdl=.4,bay=TRUE){
	dist<-c(vcf1[,2],vcf2[,2])
	dist2<-(max(dist)-min(dist))/1000000/(cmma*100)
	prior<-dlnorm(dist2,meanlog=ml,sdlog=sdl,log=TRUE)
	lik<-Lcm_vcf(cmma=cmma,vcf1=vcf1,vcf2=vcf2,l1=l1,l2=l2,inva=inva,core=core)
	post<-lik
	if(bay==TRUE){post<-lik+prior}
	return(post)
	}

#output posterior probabilities across a range of distances
	
#given a VCF file for a scaffold and the cross types in onemap format (all must be one cross type, but as of now it expects a vector)
	#calculate map distances between adjacent loci
	#ml is meanlog parameter from dlnorm. gives mean of prior distribution. sdl is sdlog from dlnorm, gives sd of prior distribution. 

map.dist<-function(vcf,crosstype,ml=log(.5),sdl=0.4,bay=TRUE){
	
	outdist<-c()
	for(i in 2:length(crosstype)){
		lik1<-optimize(f=bay_Lcm_vcf,interval=c(0,.3),maximum=TRUE,
			vcf1=vcf[i-1,],
			vcf2=vcf[i,],
			l1=crosstype[i-1],
			l2=crosstype[i],
			inva=TRUE,core="cc",ml=ml,sdl=sdl,bay=bay)
		lik2<-optimize(f=bay_Lcm_vcf,interval=c(0,.3),maximum=TRUE,
			vcf1=vcf[i-1,],
			vcf2=vcf[i,],
			l1=crosstype[i-1],
			l2=crosstype[i],
			inva=FALSE,core="cc",ml=ml,sdl=sdl,bay=bay)
		ml2<-rbind(c(lik1$maximum,lik1$objective),c(lik2$maximum,lik2$objective))	
		ml2w<-which.max(c(lik1$objective,lik2$objective))
		outdist<-rbind(outdist,ml2[ml2w,])	
		}

	outdist
	
	}
	
###following test lines functions show a histogram of random draws from the prior, ML distances vs physical distances, and max a posteriori distances vs physical distances. 

sub<-grepl("Scaffold9860_",rownames(vcf))&onemapcross=="D2.15"
altervcf(vcf[sub,],final.parent[sub])->test	
par(mfrow=c(1,3))
xxx<-log(.8);yyy<-.4
hist(rlnorm(10000,meanlog=xxx,sdlog=yyy),breaks=100)
plot(diff(test[[4]]),map.dist(test[[5]],test[[6]],ml=xxx,sdl=yyy,bay=FALSE)[,1])
lines(x=c(0,1500000),y=c(0,.03))
plot(diff(test[[4]]),map.dist(test[[5]],test[[6]],ml=xxx,sdl=yyy,bay=TRUE)[,1])
lines(x=c(0,1500000),y=c(0,.03))
exp(xxx-(yyy^2))
exp(xxx+(yyy^2)/2)

map.dist(test[[5]],test[[6]],ml=xxx,sdl=yyy,bay=TRUE)->test2

#given map distances, calculate probability of all possible (0,1)-crossover gametes

gamprob<-function(dtab){
	hald<-function(d){
		0.5*(1-exp(-2*d))
		}
	pr<-c()
	for(i in 1:length(dtab[,1])){
		r<-hald(dtab[i,1])
		nr<-1-hald(dtab[-i,1])
		prtemp<-sum(log(r,base=10),log(nr,base=10))
		pr<-c(pr,prtemp)
		}
	pr<-c(sum(log(1-hald(dtab[,1]),base=10)),pr)
	pr
	}

gamprob(test2)->test3



#given altervcf output, calculate probability of read data given gametes

gamprob_readprob<-function(avcf,gamprobs){
	
	lgen<-avcf[[1]]
	lchr<-avcf[[2]]
	lrat<-avcf[[3]]
	
	
	lengam<-length(lgen[,1])
	nind<-length(lgen[1,])
	ochr<-matrix(NA,nrow=lengam,ncol=nind)
	
	gps<-matrix(-Inf,nrow=90,ncol=5)
	
			
	for(i in 1:(lengam)){
		for(j in 1:nind){
			gam1<-c(rep(1,i-1),rep(0,lengam-(i-1)))
			gam2<-c(rep(0,i-1),rep(1,lengam-(i-1)))
			
			lik1a<-(-sum(lrat[lchr[,j]!=gam1,j],na.rm=TRUE))
			lik2a<-(-sum(lrat[lchr[,j]!=gam2,j],na.rm=TRUE))
			
			lik1<-lik1a+gamprobs[i]
			lik2<-lik2a+gamprobs[i]
			
			if(lik1>gps[j,1]){
				ochr[,j]<-gam1
				gps[j,2]<-gps[j,1]
				gps[j,1]<-lik1
				gps[j,3]<-lik1a
				gps[j,4]<-gamprobs[i]
				gps[j,5]<-i
				}else{if(lik1>gps[j,2]){
							gps[j,2]<-lik1
							}
						}
			if(lik2>gps[j,1]){
				ochr[,j]<-gam2
				gps[j,2]<-gps[j,1]
				gps[j,1]<-lik2
				gps[j,3]<-lik2a
				gps[j,4]<-gamprobs[i]
				gps[j,5]<-i
				}else{if(lik2>gps[j,2]){
							gps[j,2]<-lik2
							}
						}
			
			}
		}
	avcf[[7]]<-ochr
	avcf[[8]]<-gps
	return(avcf)
	}

gamprob_readprob(test,test3)->test4

##use output from gamprob_readprob to correct genotypes, add a new matrix to the end. 
	#use crosstype to put into "final.offspring" format. 
	#then write a script to actually replace lines in final.offspring
	
phase_correct<-function(gprpout){
	
	lengam<-length(gprpout[[1]][,1])
	outmat<-matrix(NA,nrow=lengam,ncol=length(gprpout[[1]][1,]))
	
	raw<-as.matrix(gprpout[[1]])
	phased<-as.matrix(gprpout[[7]])
	crosstype<-gprpout[[6]]
	
	refmaj<-crosstype=="abaa"|crosstype=="aaab"
	altmaj<-crosstype=="abbb"|crosstype=="bbba"
	
	if(any((!refmaj)&(!altmaj))){stop("big problems here")}
	
	for(i in 1:lengam){
		inv<-cor(raw[i,],phased[i,],use="complete")<0
		if(inv){
			phased[i,]<-as.numeric(!phased[i,])
			}
		if(refmaj[i]){
			outmat[i,]<-phased[i,]
			}
		if(altmaj[i]){
			ch<-which(phased[i,]==0)
			phased[i,ch]<-2
			outmat[i,]<-phased[i,]
			}
				
		}
	
	rownames(outmat)<-rownames(gprpout[[5]])
	
	gprpout[[9]]<-outmat
	
	
	return(gprpout)
	
	}


#enveloping function for above subfunctions. 
	#requires a subset of tightly linked, ordered markers in vcf format with final.parent style cross types (e.g. "abaa","abbb")
	#all markers must be either male or female informative, no double hets, no mixture
	
impute_genotypes<-function(vcf_sub, crosstypes,prior_meanlog=log(0.5),prior_sd=0.4){
	
	step1<-altervcf(vcf_sub,crosstypes)
	step2<-map.dist(step1[[5]],step1[[6]],ml=prior_meanlog,sdl=prior_sd,bay=TRUE)
	step3<-gamprob(step2)
	step4<-gamprob_readprob(step1,step3)
	step5<-phase_correct(step4)
	names(step5)<-c("raw_genos","chromo_assignment","likelihood_ratios","scaf_position","vcf_info","crosstype","corrected_genos","likelihoods_per_gamete","final.offspring_format")
	return(step5)
	
	}

plot.mapdist<-function(ig_out,prior_meanlog=log(0.5),prior_sd=0.4){
	
	par(mfrow=c(1,3))
	xxx<-prior_meanlog
	yyy<-prior_sd
	mdl<-map.dist(test[[5]],test[[6]],ml=xxx,sdl=yyy,bay=FALSE)
	mdb<-map.dist(test[[5]],test[[6]],ml=xxx,sdl=yyy,bay=TRUE)
	hist(rlnorm(10000,meanlog=xxx,sdlog=yyy),breaks=100,xlab="prior expectation of megabases per centimorgan",main="random samples from prior distribution")
	plot(diff(ig_out[[4]]),mdl[,1],xlab="physical distance",ylab="genetic distance",main="maximum likelihood")
	lines(x=c(0,1500000),y=c(0,.03))
	plot(diff(test[[4]]),mdb[,1],xlab="physical distance",ylab="genetic distance",main="bayesian")
	lines(x=c(0,1500000),y=c(0,1500000/1000000/exp(xxx)/100),lwd=2)
	quan<-qlnorm(p=c(.025,.975),meanlog=xxx,sdlog=yyy)
	lines(x=c(0,1500000),y=c(0,1500000/1000000/quan[1]/100),lty=2,lwd=2,col="red")
	lines(x=c(0,1500000),y=c(0,1500000/1000000/quan[2]/100),lty=2,lwd=2,col="red")	
	cat("mode of prior is ",exp(xxx-(yyy^2)),"\n")
	cat("mean of prior is ",exp(xxx+(yyy^2)/2),"\n")
	colnames(mdb)<-c("genetic distance","posterior")
	print(mdb)
	
	}

plot.mapdata<-function(ig_out,ind=c(1,20),psize=5){
	
	plot(NULL,ylim=ind,xlim=c(1,length(ig_out[[4]])))
	
	for(i in 1:length(ig_out[[3]][,1])){
		al<-ig_out[[3]]
		al[is.na(al)]<-0
		al[al<0]<-0
		al[al>10]<-10
		al2<-ig_out[[2]]
		al2[is.na(al2)]<-0
		points(x=rep(i,90),y=1:90,col=rgb(al2[i,],0,0,al[i,]/max(al,na.rm=TRUE)),pch=ig_out[[1]][i,]*2+18,cex=psize)
		}

	
	}



####full test

sub<-grepl("Scaffold9994_",rownames(vcf))&onemapcross=="D1.10"
test<-impute_genotypes(vcf[sub,],final.parent[sub])
plot.mapdist(test)
plot.mapdata(test,c(61,90),psize=4)


###backtranslate results to same format as final.offspring given crosstype

final.offspring.corrected_0.5_0.4<-final.offspring

nrecM_0.5_0.4<-c()
bpdM_0.5_0.4<-c()
nrecF_0.5_0.4<-c()
bpdF_0.5_0.4<-c()
for(i in 1:length(junk[,1])){
	print(i)	
	scaf<-paste(junk[i,1],"_",sep="")
	if(junk[i,2]>2){
		sub<-grepl(scaf,rownames(vcf))&onemapcross=="D1.10"
		test<-impute_genotypes(vcf[sub,],final.parent[sub],prior_meanlog=log(0.5),prior_sd=0.4)
		tab<-table(colSums(test[[7]]))
		print(tab)
		tab<-tab[-1]
		tab<-tab[-length(tab)]
		nrecM_0.5_0.4<-c(nrecM_0.5_0.4,sum(tab))
		bpdM_0.5_0.4<-c(bpdM_0.5_0.4,diff(range(test[[4]])))
		for(j in 1:length(test[[9]][,1])){
			namevec<-rownames(test[[9]])
			final.offspring.corrected_0.5_0.4[namevec[j],]<-test[[9]][namevec[j],]
			}		
		}

	if(junk[i,3]>2){
		sub<-grepl(scaf,rownames(vcf))&onemapcross=="D2.15"
		test<-impute_genotypes(vcf[sub,],final.parent[sub],prior_meanlog=log(0.5),prior_sd=0.4)
		tab<-table(colSums(test[[7]]))
		print(tab)
		tab<-tab[-1]
		tab<-tab[-length(tab)]
		nrecF_0.5_0.4<-c(nrecF_0.5_0.4,sum(tab))
		bpdF_0.5_0.4<-c(bpdF_0.5_0.4,diff(range(test[[4]])))
		for(j in 1:length(test[[9]][,1])){
			namevec<-rownames(test[[9]])
			final.offspring.corrected_0.5_0.4[namevec[j],]<-test[[9]][namevec[j],]		
			}
		}


	}


final.offspring.corrected_1.0_0.4<-final.offspring

nrecM_1.0_0.4<-c()
bpdM_1.0_0.4<-c()
nrecF_1.0_0.4<-c()
bpdF_1.0_0.4<-c()
for(i in 1:length(junk[,1])){
	print(i)	
	scaf<-paste(junk[i,1],"_",sep="")
	if(junk[i,2]>2){
		sub<-grepl(scaf,rownames(vcf))&onemapcross=="D1.10"
		test<-impute_genotypes(vcf[sub,],final.parent[sub],prior_meanlog=log(1),prior_sd=0.4)
		tab<-table(colSums(test[[7]]))
		print(tab)
		tab<-tab[-1]
		tab<-tab[-length(tab)]
		nrecM_1.0_0.4<-c(nrecM_1.0_0.4,sum(tab))
		bpdM_1.0_0.4<-c(bpdM_1.0_0.4,diff(range(test[[4]])))
		for(j in 1:length(test[[9]][,1])){
			namevec<-rownames(test[[9]])
			final.offspring.corrected_1.0_0.4[namevec[j],]<-test[[9]][namevec[j],]
			}		
		}

	if(junk[i,3]>2){
		sub<-grepl(scaf,rownames(vcf))&onemapcross=="D2.15"
		test<-impute_genotypes(vcf[sub,],final.parent[sub],prior_meanlog=log(1),prior_sd=0.4)
		tab<-table(colSums(test[[7]]))
		print(tab)
		tab<-tab[-1]
		tab<-tab[-length(tab)]
		nrecF_1.0_0.4<-c(nrecF_1.0_0.4,sum(tab))
		bpdF_1.0_0.4<-c(bpdF_1.0_0.4,diff(range(test[[4]])))
		for(j in 1:length(test[[9]][,1])){
			namevec<-rownames(test[[9]])
			final.offspring.corrected_1.0_0.4[namevec[j],]<-test[[9]][namevec[j],]		
			}
		}


	}


recrateM<-cbind(nrecM_1.0_0.4,bpdM_1.0_0.4,nrecM_0.5_0.4,bpdM_0.5_0.4)
recrateF<-cbind(nrecF_1.0_0.4,bpdF_1.0_0.4,nrecF_0.5_0.4,bpdF_0.5_0.4)


final.offspring.corrected_1.5_0.4<-final.offspring

nrecM_1.5_0.4<-c()
bpdM_1.5_0.4<-c()
nrecF_1.5_0.4<-c()
bpdF_1.5_0.4<-c()
for(i in 1:length(junk[,1])){
	print(i)	
	scaf<-paste(junk[i,1],"_",sep="")
	if(junk[i,2]>2){
		sub<-grepl(scaf,rownames(vcf))&onemapcross=="D1.10"
		test<-impute_genotypes(vcf[sub,],final.parent[sub],prior_meanlog=log(1.5),prior_sd=0.4)
		tab<-table(colSums(test[[7]]))
		print(tab)
		tab<-tab[-1]
		tab<-tab[-length(tab)]
		nrecM_1.5_0.4<-c(nrecM_1.5_0.4,sum(tab))
		bpdM_1.5_0.4<-c(bpdM_1.5_0.4,diff(range(test[[4]])))
		for(j in 1:length(test[[9]][,1])){
			namevec<-rownames(test[[9]])
			final.offspring.corrected_1.5_0.4[namevec[j],]<-test[[9]][namevec[j],]
			}		
		}

	if(junk[i,3]>2){
		sub<-grepl(scaf,rownames(vcf))&onemapcross=="D2.15"
		test<-impute_genotypes(vcf[sub,],final.parent[sub],prior_meanlog=log(1.5),prior_sd=0.4)
		tab<-table(colSums(test[[7]]))
		print(tab)
		tab<-tab[-1]
		tab<-tab[-length(tab)]
		nrecF_1.5_0.4<-c(nrecF_1.5_0.4,sum(tab))
		bpdF_1.5_0.4<-c(bpdF_1.5_0.4,diff(range(test[[4]])))
		for(j in 1:length(test[[9]][,1])){
			namevec<-rownames(test[[9]])
			final.offspring.corrected_1.5_0.4[namevec[j],]<-test[[9]][namevec[j],]		
			}
		}


	}

recrateM<-cbind(recrateM,nrecM_1.5_0.4,bpdM_1.5_0.4)
recrateF<-cbind(recrateF,nrecF_1.5_0.4,bpdF_1.5_0.4)
