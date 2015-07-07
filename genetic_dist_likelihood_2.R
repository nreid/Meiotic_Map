#gls are 2 3-value vectors.
#for double coupling or double repulsion. separate function for coupling/repulsion. 
Lcm_dp_ind_4het_cc<-function(cm,gl1,gl2,inv=FALSE){
	
	if(any(is.na(c(gl1,gl2)))){
		return(0)
		}
	gl1<-10^gl1
	gl2<-10^gl2
	prR<-0.5*(1-exp(-2*cm))
	
	if(inv){prR<-(1-prR)}
	
	p1<-
	gl1[1]*gl2[1]*.25*(1-prR)^2+
	gl1[2]*gl2[2]*.5*(1-prR)^2+
	gl1[3]*gl2[3]*.25*(1-prR)^2+

	gl1[1]*gl2[2]*.25*2*(1-prR)*prR +
	gl1[2]*gl2[3]*.25*2*(1-prR)*prR +
	gl1[2]*gl2[1]*.25*2*(1-prR)*prR +
	gl1[3]*gl2[2]*.25*2*(1-prR)*prR +

	gl1[1]*gl2[3]*.25*prR^2 +
	gl1[3]*gl2[1]*.25*prR^2 +
	gl1[2]*gl2[2]*.5*prR^2


	
	return(log(p1))
	
	}


Lcm_dp_ind_4het_cr<-function(cm,gl1,gl2){
	
	if(any(is.na(c(gl1,gl2)))){
		return(0)
		}
	gl1<-10^gl1
	gl2<-10^gl2
	prR<-0.5*(1-exp(-2*cm))
	
	
	p1<-
	gl1[1]*gl2[2]*.25*((1-prR)^2+prR^2)+
	gl1[2]*gl2[3]*.25*((1-prR)^2+prR^2)+
	gl1[2]*gl2[1]*.25*((1-prR)^2+prR^2)+
	gl1[3]*gl2[2]*.25*((1-prR)^2+prR^2)+

	gl1[1]*gl2[1]*.25*(1-prR)*prR+
	gl1[3]*gl2[3]*.25*(1-prR)*prR+
	gl1[1]*gl2[3]*.25*(1-prR)*prR+
	gl1[3]*gl2[1]*.25*(1-prR)*prR+
	gl1[2]*gl2[2]*.5*2*(1-prR)*prR
	
		
	return(log(p1))
	
	}


#het1 is a het/hom cross (2-value vector). het2 is a hethet cross (3-value vector. 

Lcm_dp_ind_3het<-function(cm,het1,het2,inv=FALSE){
	
	if(any(is.na(c(het1,het2)))){
		return(0)
		}
	het1<-10^het1
	het2<-10^het2
	prR<-0.5*(1-exp(-2*cm))
	
	if(inv){prR<-(1-prR)}
	
	p1<-
	het1[1]*het2[1]*.25*(1-prR)+
	het1[2]*het2[3]*.25*(1-prR)+
	het1[1]*het2[3]*.25*prR+
	het1[2]*het2[1]*.25*prR+
	het1[1]*het2[2]*.25+
	het1[2]*het2[2]*.25

	return(log(p1))
	
	}


#gls are two, 2-value vectors

Lcm_dp_ind_2het<-function(cm,gl1,gl2,inv=FALSE){
	
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


Lcm_dp_sites_multi<-function(cmm,gl1m,gl2m,invm=FALSE, ctype){

	if(ctype=="2het"){
		gls<-cbind(gl1m,gl2m)
		if(length(gls[1,])!=4){stop("ya done messed up a-aron!")}
		lik<-apply(gls,MAR=1,FUN=function(x){Lcm_dp_ind_2het(cm=cmm,gl1=x[1:2],gl2=x[3:4],inv=invm)})
		}
		
	if(ctype=="3het"){
		if(length(gl1m[1,])==2){gls<-cbind(gl1m,gl2m)}
		else{gls<-cbind(gl2m,gl1m)}
		if(length(gls[1,])!=5){stop("ya done messed up a-aron!")}
		lik<-apply(gls,MAR=1,FUN=function(x){Lcm_dp_ind_3het(cm=cmm,het1=x[1:2],het2=x[3:5],inv=invm)})
		}

	if(ctype=="4het_cc"){
		gls<-cbind(gl1m,gl2m)
		if(length(gls[1,])!=6){stop("ya done messed up a-aron!")}
		lik<-apply(gls,MAR=1,FUN=function(x){Lcm_dp_ind_4het_cc(cm=cmm,gl1=x[1:3],gl2=x[4:6],inv=invm)})
		}

	if(ctype=="4het_cr"){
		gls<-cbind(gl1m,gl2m)
		if(length(gls[1,])!=6){stop("ya done messed up a-aron!")}
		lik<-apply(gls,MAR=1,FUN=function(x){Lcm_dp_ind_4het_cr(cm=cmm,gl1=x[1:3],gl2=x[4:6])})
		}
		
		
	return(sum(lik))
	
	}

###optimization needs to be done inside of this. repeatedly processing text is CRAZY. 
Lcm_vcf<-function(cmma,vcf1,vcf2,l1,l2,inva=FALSE,core="cc"){

	subc<-list(c(1,2),c(1,2),c(3,2),c(3,2),c(1:3))
	names(subc)<-c("aaab","abaa","abbb","bbba","abab")

	gl1a<-vcf1[16:105] %>% gsub(".*:","",.) %>% strsplit(split=",",x=.) %>% do.call(rbind,.)
	class(gl1a)<-"numeric"
	gl2a<-vcf2[16:105] %>% gsub(".*:","",.) %>% strsplit(split=",",x=.) %>% do.call(rbind,.)
	class(gl2a)<-"numeric"	
	gl1a<-gl1a[,subc[[l1]]]
	gl2a<-gl2a[,subc[[l2]]]	

	if(sum(l1=="abab",l2=="abab")==2){
		if(core=="cc"){
			ct<-"4het_cc"
			lik<-Lcm_dp_sites_multi(cmm=cmma,gl1m=gl1a,gl2m=gl2a,invm=inva,ctype=ct)
			}
		if(core=="cr"){
			ct<-"4het_cr"
			lik<-Lcm_dp_sites_multi(cmm=cmma,gl1m=gl1a,gl2m=gl2a,invm=inva,ctype=ct)
			}
		}
	if(sum(l1=="abab",l2=="abab")==1){
		ct<-"3het"
		lik<-Lcm_dp_sites_multi(cmm=cmma,gl1m=gl1a,gl2m=gl2a,invm=inva,ctype=ct)
		}	
	if(l1==l2&sum(l1=="abab",l2=="abab")==0){
		ct<-"2het"
		lik<-Lcm_dp_sites_multi(cmm=cmma,gl1m=gl1a,gl2m=gl2a,invm=inva,ctype=ct)
		}
	if(sum(c("abaa","abbb")%in%c(l1,l2))==2){
		ct<-"2het"
		lik<-Lcm_dp_sites_multi(cmm=cmma,gl1m=gl1a,gl2m=gl2a,invm=inva,ctype=ct)
		}
	if(sum(c("aaab","bbba")%in%c(l1,l2))==2){
		ct<-"2het"
		lik<-Lcm_dp_sites_multi(cmm=cmma,gl1m=gl1a,gl2m=gl2a,invm=inva,ctype=ct)
		}
	return(lik)	
	}

Lcm_vcf_surface<-function(vcf1,vcf2,l1,l2,inve=FALSE,core="cc"){
	x<-seq(from=.0001,to=.5,by=.001)
	y<-unlist(lapply(x,FUN=Lcm_vcf,vcf1=vcf1,vcf2=vcf2,l1=l1,l2=l2,inva=inve,core=core))
	cbind(x,y)
	}

plot(Lcm_vcf_surface(vcf1=vcf["Scaffold9860_156093",],vcf2=vcf["Scaffold9860_373981",],l1=final.parent["Scaffold9860_156093"],l2=final.parent["Scaffold9860_373981"],inve=TRUE,core="cc"))
Lcm_vcf(cmma=.03,vcf1=vcf["Scaffold9860_156093",],vcf2=vcf["Scaffold9860_373981",],l1=final.parent["Scaffold9860_156093"],l2=final.parent["Scaffold9860_373981"],inva=TRUE,core="cc")
Lcm_vcf_optim(vcf1=vcf["Scaffold9860_156093",],vcf2=vcf["Scaffold9860_373981",],l1=final.parent["Scaffold9860_156093"],l2=final.parent["Scaffold9860_373981"])

Lcm_vcf_optim<-function(vcf1,vcf2,l1,l2){

	subc<-list(c(1,2),c(1,2),c(3,2),c(3,2),c(1:3))
	names(subc)<-c("aaab","abaa","abbb","bbba","abab")

	gl1a<-vcf1[16:105] %>% gsub(".*:","",.) %>% strsplit(split=",",x=.) %>% do.call(rbind,.)
	class(gl1a)<-"numeric"
	gl2a<-vcf2[16:105] %>% gsub(".*:","",.) %>% strsplit(split=",",x=.) %>% do.call(rbind,.)
	class(gl2a)<-"numeric"	
	gl1a<-gl1a[,subc[[l1]]]
	gl2a<-gl2a[,subc[[l2]]]	

	if(sum(l1=="abab",l2=="abab")==2){
		ct<-"4het_cc"
		lik1<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=FALSE,ctype=ct,maximum=TRUE)
		lik2<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=TRUE,ctype=ct,maximum=TRUE)
		ct<-"4het_cr"
		lik3<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=TRUE,ctype=ct,maximum=TRUE)
		
		lik<-rbind(c(lik1$maximum,lik1$objective),c(lik2$maximum,lik2$objective),c(lik3$maximum,lik3$objective))
		lik<-lik[which.max(lik[,2]),]
		}
	if(sum(l1=="abab",l2=="abab")==1){
		ct<-"3het"
		lik1<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=FALSE,ctype=ct,maximum=TRUE)
		lik2<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=TRUE,ctype=ct,maximum=TRUE)
		lik<-rbind(c(lik1$maximum,lik1$objective),c(lik2$maximum,lik2$objective))
		lik<-lik[which.max(lik[,2]),]
		}	
	if(l1==l2&sum(l1=="abab",l2=="abab")==0){
		ct<-"2het"
		lik1<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=FALSE,ctype=ct,maximum=TRUE)
		lik2<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=TRUE,ctype=ct,maximum=TRUE)
		lik<-rbind(c(lik1$maximum,lik1$objective),c(lik2$maximum,lik2$objective))
		lik<-lik[which.max(lik[,2]),]
		}
	if(sum(c("abaa","abbb")%in%c(l1,l2))==2){
		ct<-"2het"
		lik1<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=FALSE,ctype=ct,maximum=TRUE)
		lik2<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=TRUE,ctype=ct,maximum=TRUE)
		lik<-rbind(c(lik1$maximum,lik1$objective),c(lik2$maximum,lik2$objective))
		lik<-lik[which.max(lik[,2]),]
		}
	if(sum(c("aaab","bbba")%in%c(l1,l2))==2){
		ct<-"2het"
		lik1<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=FALSE,ctype=ct,maximum=TRUE)
		lik2<-optimize(f=Lcm_dp_sites_multi,interval=c(.0001,.5),gl1m=gl1a,gl2m=gl2a,invm=TRUE,ctype=ct,maximum=TRUE)
		lik<-rbind(c(lik1$maximum,lik1$objective),c(lik2$maximum,lik2$objective))
		lik<-lik[which.max(lik[,2]),]
		
		}
	if(any(c("abaa","abbb")%in%c(l1,l2))&any(c("aaab","bbba")%in%c(l1,l2))){return(NULL)}
	lik<-as.vector(c(vcf1[,1],vcf1[,2],vcf2[,2],lik,l1,l2))
	return(lik)
	}



L_pairwise<-function(vcftab,crosses,scafname){
	lvec<-vcftab[,1]==scafname
	tempp<-crosses[lvec]
	tempvcf<-vcf[lvec,]
	t(combn(1:length(tempp),2))->pjunk
	#out<-matrix(nrow=length(pjunk[,1]),ncol=7)
	out<-c()
	for(i in 1:length(pjunk[,1])){
		lik<-Lcm_vcf_optim(tempvcf[pjunk[i,1],],tempvcf[pjunk[i,2],],tempp[pjunk[i,1]],tempp[pjunk[i,2]])
		out<-rbind(out,lik)
		}
	rownames(out)<-NULL
	out<-data.frame(out,stringsAsFactors=FALSE)	
	out[,2]<-as.numeric(out[,2])
	out[,3]<-as.numeric(out[,3])
	out[,4]<-as.numeric(out[,4])
	out[,5]<-as.numeric(out[,5])
	return(out)
	}


Lcm_vcf_optim(vcf["Scaffold1_1976295",],vcf["Scaffold1_2700034",],final.parent["Scaffold1_1976295"],final.parent["Scaffold1_2700034"])

L_pairwise(vcf,final.parent,"Scaffold0")->test

names(table(vcf[,1])[table(vcf[,1])>=8])->subscaf
allpwise3<-c()	
for(i in 1:length(subscaf)){
	allpwise3<-rbind(allpwise3,L_pairwise(vcf,final.parent,subscaf[i]))
	print(i)
	}
quit(save="yes")


###look through allpwise3 and find pairs (excluding ones with abab loci) and find loci with very low likelihoods. 
	###check to see if they are legit, or if there are problems. 
	#list for exclusion: "Scaffold9924_1576440", #wrong cross type. others seem fine. 

which(names(final.parent)=="Scaffold9924_1576440")
final.parent[-4615]->final.parent
final.offspring[-4615,]->final.offspring
vcf[-4615,]->vcf
allpwise3<-allpwise3[!(allpwise3[,1]=="Scaffold9924"&(allpwise3[,2]==1576440|allpwise3[,3]==1576440)),]
	
###look through allpwise3 and find markers unlinked to their scaffolds. assign them new scaffold ids. 

	#list for unlinking: 
	
unlink<-c("Scaffold9940_56468","Scaffold9926_520606","Scaffold9926_56851",
	names(final.parent[grep("Scaffold9916",names(final.parent))]), "Scaffold9911_5017096","Scaffold9859_34668",
	"Scaffold9859_862773","Scaffold173_171962","Scaffold173_2810","Scaffold173_1039607","Scaffold1_1740178")

names(final.parent)[which(names(final.parent)=="Scaffold9940_56468")]<-"Scaffold9940b_56468"
rownames(final.offspring)[which(rownames(final.offspring)=="Scaffold9940_56468")]<-"Scaffold9940b_56468"

x<-28
re<-"b_"
names(final.parent)[which(names(final.parent)==unlink[x])]<-gsub("_",re,unlink[x])
rownames(final.offspring)[which(rownames(final.offspring)==unlink[x])]<-gsub("_",re,unlink[x])


system.time(
	optimize(f=Lcm_vcf,interval=c(0,.5),vcf1=vcf[5,],vcf2=vcf[10,],l1=final.parent[5],l2=final.parent[10],inv=FALSE,maximum=TRUE)
	)
system.time(	
	Lcm_vcf_optim(vcf1=vcf[5,],vcf2=vcf[10,],l1=final.parent[5],l2=final.parent[10])
	)



#ratio is centimorgans/megabase for each pairwise comparison, i.e. 1 for 1 centimorgan per megabase
#distvec is the pairwise distance between all marker pairs
#logicvec is goes in "inva"

Lcm_vcf_gdist<-function(ratio,pairwise,distvec,logicvec,vcf,cross){
	ratio<-ratio*1e-08
	outlik<-0
	for(i in 1:length(distvec)){
		cm<-distvec[i]*ratio
		ind1<-pairwise[i,1]
		ind2<-pairwise[i,2]
		tmp<-Lcm_vcf(cm,vcf[ind1,],vcf[ind2,],cross[ind1],cross[ind2],inva=logicvec[i])
		outlik<-outlik+tmp
		
		}
	
	outlik
	
	}

vcf[((final.parent=="abaa"|final.parent=="abbb")&vcf[,1]=="Scaffold0"),]->p01
final.parent[((final.parent=="abaa"|final.parent=="abbb")&vcf[,1]=="Scaffold0")]->p01f
vcf[((final.parent=="aaab"|final.parent=="bbba")&vcf[,1]=="Scaffold0"),]->p02
final.parent[((final.parent=="aaab"|final.parent=="bbba")&vcf[,1]=="Scaffold0")]->p02f

pairwise_test<-t(combn(dim(p01)[1],2))

log_test<-c()
test_dist<-c()
for(i in 1:(length(p01[,1])-1)){
	for(j in (i+1):length(p01[,1])){
		a<-optimize(f=Lcm_vcf,interval=c(0,.5),vcf1=p01[i,],vcf2=p01[j,],l1=p01f[i],l2=p01f[j],inv=FALSE,maximum=TRUE)$objective
		b<-optimize(f=Lcm_vcf,interval=c(0,.5),vcf1=p01[i,],vcf2=p01[j,],l1=p01f[i],l2=p01f[j],inv=TRUE,maximum=TRUE)$objective
		log_test<-c(log_test,b>a)
		test_dist<-c(test_dist,p01[j,2]-p01[i,2])
		}
	}

optimize(f=Lcm_vcf_gdist,interval=c(.1,4),pairwise_test,test_dist,log_test,p01,p01f,maximum=TRUE)

sum((allpwise3[1:541,4])/((allpwise3[1:541,3]-allpwise3[1:541,2]))*1000000*((allpwise3[1:541,3]-allpwise3[1:541,2])))/sum(((allpwise3[1:541,3]-allpwise3[1:541,2])))
