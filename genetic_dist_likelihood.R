##Likelihood of parent phases and genetic distance for two locus crosses, haldane mapping function

#two het-hom loci
	#not going to return the actual phase for these functions, just the likelihoods and distances for the two phases
#(1/2) (1- exp(-2d)). 
	
#gls are two, 2-value vectors
Lcm_dp_ind<-function(cm,gl1,gl2,inv=FALSE){
	
	if(any(is.na(c(gl1,gl2)))){
		return(0)
		}
	gl1<-10^gl1
	gl2<-10^gl2
	prR<-0.5*(1-exp(-2*cm))
	
	if(inv){prR<-(1-prR)}
	
	p1<-(.5*gl1[1]*gl2[1]*(1-prR))+
	(.5*gl1[2]*gl2[2]*(1-prR))+
	(.5*gl1[1]*gl2[2]*(prR))+
	(.5*gl1[2]*gl2[1]*(prR))

	return(log(p1))
	
	}

#requires two 2-column GL inputs, one for each locus, containing the two relevant genotypes 

Lcm_dp_sites<-function(cmm,gl1m,gl2m,invm=FALSE){

	gls<-cbind(gl1m,gl2m)
	
	lik<-apply(gls,MAR=1,FUN=function(x){Lcm_dp_ind(cm=cmm,gl1=x[1:2],gl2=x[3:4],inv=invm)})
	
	return(sum(lik))
	
	}

#given two VCF lines and parent genotypes for each locus

Lcm_vcf<-function(cmma,vcf1,vcf2,l1,l2,inva=FALSE){
	
	subc<-rbind(c(1,2),c(1,2),c(2,3),c(2,3))
	rownames(subc)<-c("aaab","abaa","abbb","bbba")
	
	gl1a<-vcf1[16:105] %>% gsub(".*:","",.) %>% strsplit(split=",",x=.) %>% do.call(rbind,.)
	class(gl1a)<-"numeric"
	gl1a<-gl1a[,subc[l1,]]
	gl2a<-vcf2[16:105] %>% gsub(".*:","",.) %>% strsplit(split=",",x=.) %>% do.call(rbind,.)
	class(gl2a)<-"numeric"	
	gl2a<-gl2a[,subc[l2,]]	
	
	lik<-Lcm_dp_sites(cmm=cmma,gl1m=gl1a,gl2m=gl2a,invm=inva)
	
	return(lik)
	
	}

#out<-matrix(nrow=14,ncol=4)
#for(i in 2:14){
#	out[i-1,1]<-tempvcf[i-1,2]
#	out[i-1,2]<-tempvcf[i,2]
#	t1<-optimize(f=Lcm_vcf, interval=c(0.0001,.5),vcf1=tempvcf[i-1,],vcf2=tempvcf[i,],l1=tempp[i-1],l2=tempp[i],inva=FALSE,maximum=TRUE)
#	t2<-optimize(f=Lcm_vcf, interval=c(0.0001,.5),vcf1=tempvcf[i-1,],vcf2=tempvcf[i,],l1=tempp[i-1],l2=tempp[i],inva=TRUE,maximum=TRUE)
#	if(t1$objective>t2$objective){
#		out[i-1,3]<-t1$maximum
#		out[i-1,4]<-t1$objective
#		}else{
#			out[i-1,3]<-t2$maximum
#			out[i-1,4]<-t2$objective	
#			
#			}
#	
#	
#	}


L_pairwise<-function(vcftab,scafname,cross1,cross2){
	tempp<-final.parent[vcf[,1]==scafname&(final.parent==cross1|final.parent==cross2)]
	tempvcf<-vcf[vcf[,1]==scafname&(final.parent==cross1|final.parent==cross2),]
	t(combn(1:length(tempp),2))->pjunk
	out<-matrix(nrow=length(pjunk[,1]),ncol=4)
	for(i in 1:length(pjunk[,1])){
	
		out[i,1]<-tempvcf[pjunk[i,1],2]
		out[i,2]<-tempvcf[pjunk[i,2],2]
		t1<-optimize(f=Lcm_vcf, interval=c(0.0001,.5),vcf1=tempvcf[pjunk[i,1],],vcf2=tempvcf[pjunk[i,2],],l1=tempp[pjunk[i,1]],l2=tempp[pjunk[i,2]],inva=FALSE,maximum=TRUE)
		t2<-optimize(f=Lcm_vcf, interval=c(0.0001,.5),vcf1=tempvcf[pjunk[i,1],],vcf2=tempvcf[pjunk[i,2],],l1=tempp[pjunk[i,1]],l2=tempp[pjunk[i,2]],inva=TRUE,maximum=TRUE)
		if(t1$objective>t2$objective){
			out[i,3]<-t1$maximum
			out[i,4]<-t1$objective
			}else{
				out[i,3]<-t2$maximum
				out[i,4]<-t2$objective	
			
				}	
		}
	out<-as.data.frame(out)	
	out<-cbind(scafname,out)
	return(out)		
	}

L_pairwise(vcf,"Scaffold10000",cross1="abaa",cross2="abbb")->out
sum(out[,4])/(sum(out[,3]-out[,2]))*1000000
plot(NULL,ylim=c(0,.5),xlim=c(0,6100000))
	for(i in 1:length(out[,2])){
	lines(out[i,2:3],rep(out[i,4],2),pch=20,type="b")
	}

names(table(vcf[final.parent=="abaa"|final.parent=="abbb",1])[table(vcf[final.parent=="abaa"|final.parent=="abbb",1])>1])->subscaf1
names(table(vcf[final.parent=="aaab"|final.parent=="bbba",1])[table(vcf[final.parent=="aaab"|final.parent=="bbba",1])>1])->subscaf2

library(magrittr)
allpwise1<-c()	
for(i in 1:length(subscaf1)){
	allpwise1<-rbind(allpwise1,L_pairwise(vcf,subscaf1[i],cross1="abaa",cross2="abbb"))
	print(i)
	}

allpwise2<-c()	
for(i in 1:length(subscaf2)){
	allpwise2<-rbind(allpwise2,L_pairwise(vcf,subscaf2[i],cross1="aaab",cross2="bbba"))
	print(i)
	}

summarize_dist_byloc<-function(pwise){
	
	pwise1<-cbind(pwise[,c(1,2,4,5)],pwise[,3]-pwise[,2])
	pwise2<-cbind(pwise[,c(1,3,4,5)],pwise[,3]-pwise[,2])

	names(pwise1)<-letters[1:5]
	names(pwise2)<-letters[1:5]
	pwise<-rbind(pwise1,pwise2)
	
	pwise2<-group_by(pwise,a,b) %>% summarize(., avdist=sum(c)/sum(e),dist=mean(c),lik=mean(d),count=n())
	
	as.data.frame(pwise2,stringsAsFactors=FALSE)
	
	}

summarize_dist_byloc<-function(pwise){
	
	pwise1<-cbind(pwise[,c(1,2,4,5)],pwise[,3]-pwise[,2])
	pwise2<-cbind(pwise[,c(1,3,4,5)],pwise[,3]-pwise[,2])

	names(pwise1)<-letters[1:5]
	names(pwise2)<-letters[1:5]
	pwise<-rbind(pwise1,pwise2)
	
	pwise2<-group_by(pwise,a,b) %>% summarize(., avdist=sum(c)/sum(e),dist=mean(c),lik=mean(d),count=n())
	
	as.data.frame(pwise2,stringsAsFactors=FALSE)
	
	}

summarize_dist_byscaf<-function(pwise){
	
	pwise2<-group_by(pwise,scafname) %>% summarize(., avdist=sum(V3)/sum(V2-V1),sumdist=sum(V2-V1),lik=mean(V4),count=n())
	
	as.data.frame(pwise2,stringsAsFactors=FALSE)
	
	
	}

