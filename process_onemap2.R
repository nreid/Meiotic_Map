###In this script we are testing whether or not we can leverage linkage information from markers already clustered into scaffolds
###to generate linkage groups. our radseq data are too sparse and error prone to effectively use joinmap or onemap naively

library(onemap)
library(magrittr)
library(stringr)

####read in ALL onemap data:
fname<-"onemap.out.replace.2.csv"

onemapdata<-read.table(fname,skip=1,stringsAsFactors=FALSE,sep=" ")

#create a vector of scaffold names corresponding to each locus
scafnames<-onemapdata[,1]
gsub("\\*","",scafnames) %>% gsub("_[0-9]*", "",.) -> scafnames

#create an object containing the names of the n scaffolds with the largest numbers of markers
first60<-names(head(sort(table(scafnames),decreasing=TRUE),n=168))
nscaf<-length(first60)


onemapsub <- read.outcross(".",fname)
onemapsub$n.mar->nmarkers

#calculate 2 point recombination frequencies. 
onemapsub.rec <- rf.2pts(onemapsub)

	#####the output of rf.2pts is a large object. 
	#####$analysis contains a 3 dimensional array
	#####the dimensions are [nmarkers,4,c("Theta","LODs")]
	#####theta is the ML recombination fraction, LODs are the LOD scores 
	
#The next step is to pull out all the LOD scores for the recombination frequencies

#table containing numbers for all marker pairs. pcom = "pairwise comparisons of markers"
pcom<-(t(combn(1:nmarkers,2)))

##internal onemap function used to retrieve data for marker pairs
acum<-function (w){
    if (w < 0) 
        stop("'w' should be equal to or higher than zero")
    w * (w + 1)/2
	}	


#each pair of markers has 4 LOD scores for the four possible phases. get the highest of these four for each marker pair
lodout.1<-((apply(onemapsub.rec$analysis[,,"LODs"],MAR=1,FUN=max)))

#reorder and truncate LOD scores
lodout.2<-lodout.1[unlist(lapply(pcom[,2]-2,FUN=acum))+pcom[,1]]
lodout.3<-lodout.2
lodout.3[lodout.2>10]<-10
lodout.3[lodout.2<3]<-0

#plot LOD scores
plot(NULL,xlim=c(1,nmarkers),ylim=c(1,nmarkers))
points(pcom[,1],pcom[,2],pch=15,col=rgb(1,0,0,((lodout.3)/10)),cex=.2)
abline(v=(1:nmarkers)[!duplicated(gsub("(?<=[0-9]_).*","",onemapsub.rec$marnames,perl=TRUE))])

#create a square matrix containing LOD scores. dunno if this one works. didn't run it yet. 
LODmat<-matrix(nrow=nmarkers,ncol=nmarkers)
for(i in 1:length(lodout.3)){
	LODmat[pcom[i,1],pcom[i,2]]<-lodout.3[i]
	LODmat[pcom[i,2],pcom[i,1]]<-lodout.3[i]
	}
dimnames(LODmat)<-list(onemapsub.rec$marnames,onemapsub.rec$marnames)

colfunc<- colorRampPalette(c("white", "red"))
image(LODmat,zlim=c(3,10),col=colfunc(7))

#####now we calculate the mean LOD score between pairs of scaffolds

######pcomscaf enumerates pairwise comparisons of markers, but contains only scaffold names, not marker positions
pcomscaf<-cbind(onemapsub.rec$marnames[pcom[,1]],onemapsub.rec$marnames[pcom[,2]])
pcomscaf<-gsub("_[0-9]*$","",pcomscaf)

pcomcross<-cbind(onemapsub$segr.type[pcom[,1]],onemapsub$segr.type[pcom[,2]])
pcomcross<-(pcomcross[,1]=="D2.15"&pcomcross[,2]=="D1.10")|(pcomcross[,2]=="D2.15"&pcomcross[,1]=="D1.10")
pcomscaf<-pcomscaf[!pcomcross,]
pcom<-pcom[!pcomcross,]
lodout.4<-lodout.3[!pcomcross]
###pairscafs contains all pairwise comparisons of scaffolds in the data
gsub("^","_",first60)->first60.2
pairscafs<-t(combn(first60.2,2))

####iterate over all pairwise comparisons of scaffolds, calculating mean lod scores. 
### the number of comparisons should be (nscafs*nscafs-1)/2

#meanlod<-c()
#for(i in 1:length(pairscafs[,1])){
	
#	lvec1<-pcomscaf[,1]==pairscafs[i,1]&pcomscaf[,2]==pairscafs[i,2]
#	lvec2<-pcomscaf[,1]==pairscafs[i,2]&pcomscaf[,2]==pairscafs[i,1]
#	meanlod<-c(meanlod, mean(lodout.3[lvec1|lvec2]))
#	if((i%%10)==0){print(i)}
#	}

meanlod<-c()
len<-length(first60.2)
it<-0
for(i in 1:(len-1)){
	lvec1<-pcomscaf[,1]==first60.2[i]|pcomscaf[,2]==first60.2[i]
	#print(sum(lvec1))
	templod<-lodout.4[lvec1]
	tempscaf<-pcomscaf[lvec1,]
	for(j in (i+1):len){
		lvec2<-tempscaf[,1]==first60.2[j]|tempscaf[,2]==first60.2[j]
		#print(sum(lvec2))
		meanlod<-c(meanlod,mean(templod[lvec2]))
		it<-it+1
		if((it%%10)==0){print(it)}
		}
		
#	lvec1<-pcomscaf[,1]==pairscafs[i,1]&pcomscaf[,2]==pairscafs[i,2]
#	lvec2<-pcomscaf[,1]==pairscafs[i,2]&pcomscaf[,2]==pairscafs[i,1]
#	meanlod<-c(meanlod, mean(lodout.3[lvec1|lvec2]))
#	if((i%%10)==0){print(i)}
	}



#### to visualize scaffold linkages... need to play with the rgb settings and/or scale meanlod to get the best results. 
sn<-(t(combn(1:nscaf,2)))
plot(sn[,1],sn[,2],pch=15,col=rgb(1,0,0,meanlod/9))

#####Now create clusters of scaffolds based on the mean LOD scores connecting them. 

###given a scaffold, or set of scaffolds, mean lod scores, and a Xx2 matrix of scaffold pairs, 
	###this recursive function generates a cluster of all scaffolds connected to at greater than a threshold

growclust<-function(clust, pairwise, lods,thresh=.1){
	
	startl<-length(clust)
	scafs<-rbind(pairwise[pairwise[,1]%in%clust&lods>thresh,], pairwise[pairwise[,2]%in%clust&lods>thresh,])
	scafs<-unique(as.vector(scafs))
	if(length(scafs)==0){return(clust)}	
	clust<-scafs
	endl<-length(clust)
	
	if(startl<endl){growclust(clust,pairwise,lods,thresh)}
	
	else{return(clust)}
	
	}

lgs<-cbind(gsub("^","_",first60),NA)
tout<-c()
lgroup<-1

while(sum(is.na(lgs[,2]))>0){
	
	tout<-lgs[is.na(lgs[,2]),1][1]
	tout<-growclust(tout,pairscafs,meanlod,thresh=.5)
	lgs[lgs[,1]%in%tout,2]<-lgroup
	lgroup<-lgroup+1
	
	}

lgs[order(lgs[,2]),]->lgs

#####Now we want to reorder the square matrix of LOD scores to see if our scaffold clustering is working well. 

#####this function reorders the square LODmat based on the linkage groups assigned in scaf_links (lgs above)
reorder.LODmat<-function(marker_LODs, scaf_links,addall=FALSE){
	
	loci<-rownames(marker_LODs)
	lgs<-unique(scaf_links[,2])
	ngroups<-length(lgs)
	
	neworder<-c()
	newloci<-c()
	
	for(i in lgs){
		
		scafs<-scaf_links[scaf_links[,2]==i,1]
		
		for(j in scafs){
			
			so<-grep(paste(j,"_",sep=""),loci)
			neworder<-c(neworder,so)
			newloci<-c(newloci,gsub("^",paste("L",i,sep=""),loci[so]))
			}
		
		}
	if(addall){
		leftout<-!(1:length(loci))%in%neworder
		neworder<-c(neworder,which(leftout))
		newloci<-c(newloci,loci[leftout])
		}
	marker_LODs<-marker_LODs[neworder,neworder]
	rownames(marker_LODs)<-newloci
	colnames(marker_LODs)<-newloci
	return(marker_LODs)
	}

LODmat2<-reorder.LODmat(LODmat,lgs)
nmarkers2<-dim(LODmat2)[1]
colfunc<- colorRampPalette(c("white", "red"))
image(LODmat2,zlim=c(1,10),col=colfunc(10),x=1:nmarkers2,y=1:nmarkers2)
abline(v=(1:nmarkers2)[!duplicated(gsub("_.*","",rownames(LODmat2),perl=TRUE))]-.5,col="blue",lwd=2)
abline(v=(1:nmarkers2)[!duplicated(str_extract(rownames(LODmat2),"SCFLD[0-9]+"))]-.5)

LODmat<-matrix(nrow=nmarkers,ncol=nmarkers)
for(i in 1:length(lodout.4)){
	LODmat[pcom[i,1],pcom[i,2]]<-lodout.4[i]
	LODmat[pcom[i,2],pcom[i,1]]<-lodout.4[i]
	}
dimnames(LODmat)<-list(onemapsub.rec$marnames,onemapsub.rec$marnames)

LODmat3<-reorder.LODmat(LODmat,lgs,addall=TRUE)
image(LODmat3,zlim=c(1,10),col=colfunc(10),x=1:5777,y=1:5777)
image(LODmat3[1:5777,3271:5777],zlim=c(1,10),col=colfunc(10),x=1:5777,y=3271:5777)



bounds<-c(0,(1:nmarkers2)[!duplicated(gsub("_.*","",rownames(LODmat2),perl=TRUE))][-1]-1,3270)
amat<-c()
for(i in 1:(length(bounds)-1)){
	amat<-cbind(amat,apply(LODmat3[,(bounds[i]+1):bounds[i+1]],MAR=1,FUN=mean,na.rm=TRUE))
	}
amatmax<-apply(amat,MAR=1,FUN=max)
amatmaxind<-apply(amat,MAR=1,FUN=which.max)
amatdiff<-amatmax-apply(amat,MAR=1,FUN=function(x){sort(x,decreasing=TRUE)[2]})

plot(hclust(dist(t(amat[which((amatdiff/amatmax)<.8),]))))

image(LODmat3[grep("L20_|L26_",colnames(LODmat3)),grep("L20_|L26_",colnames(LODmat3))],
	zlim=c(1,10),col=colfunc(10))
	
image(LODmat3[grep("L24_|L6_",colnames(LODmat3)),grep("L24_|L6_",colnames(LODmat3))],
	zlim=c(1,10),col=colfunc(10))
	
#lodmat contains lod scores per locus, sorted by LG with unassigned markers tacked on (from reorder.LODmat)
#lodthresh is the min average lod to assign a marker to an LG. diffthresh is the min difference from the highest average lod to the second highest average lod
assignsome<-function(lodmat,lodthresh,diffthresh){
	

	ind<-grepl("L[0-9]+_",colnames(lodmat))
	nmara<-sum(ind)
	marnames<-colnames(lodmat)
	nmart<-length(marnames)
	lgs<-unique(gsub("_.*","",rownames(lodmat[ind,])))
	
	cat(sum(ind), " markers already assigned.\n", sum(!ind), " markers to be assigned to ", length(lgs), " linkage groups.\n",sep="")

	bounds<-c(0,(1:nmara)[!duplicated(gsub("_.*","",rownames(lodmat[ind,]),perl=TRUE))][-1]-1,nmara)
	amat<-c()
	for(i in 1:(length(bounds)-1)){
		amat<-cbind(amat,apply(lodmat[,(bounds[i]+1):bounds[i+1]],MAR=1,FUN=mean,na.rm=TRUE))
		}
	amatmax<-apply(amat,MAR=1,FUN=max)
	amatmaxind<-apply(amat,MAR=1,FUN=which.max)
	amatdiff<-amatmax-apply(amat,MAR=1,FUN=function(x){sort(x,decreasing=TRUE)[2]})
		
	assignable<-!ind&amatmax>=lodthresh&amatdiff/amatmax>=diffthresh
	cat(sum(assignable,na.rm=TRUE), " markers can be assigned.\n",
		sum(!assignable&!ind,na.rm=TRUE), " markers have ambiguous assignments.\n",
		sum(amatmax==0), " markers are hopeless.\n",sep="")
	
	for(i in 1:length(lgs)){
		cat("assigning ", sum(assignable&amatmaxind==i,na.rm=TRUE), " markers to ", lgs[i],".\n",sep="")
		marnames[which(assignable&amatmaxind==i)]<-paste(lgs[i],"_",marnames[which(assignable&amatmaxind==i)],sep="")
		
		}
	
	colnames(lodmat)<-marnames
	rownames(lodmat)<-marnames
	
	ord<-sort(marnames[grep("^L",marnames)])
	ord<-c(ord,sort(marnames[grep("^_",marnames)]))
	lodmat[ord,ord]
	
	}
	
assignsome(LODmat3,1,.95)->LODmat4
image(LODmat4,zlim=c(1,10),col=colfunc(10),x=1:5777,y=1:5777)
abline(v=(1:5320)[!duplicated(gsub("_.*","",rownames(LODmat4),perl=TRUE))]-.5,col="blue",lwd=2)

#join LGs L20/L26 and L6/L24
LODmat5<-LODmat4
newmarnames<-gsub("L26_","L20_",colnames(LODmat5))
newmarnames<-gsub("L24_","L6_",newmarnames)
colnames(LODmat5)<-newmarnames
rownames(LODmat5)<-newmarnames
ord<-sort(newmarnames[grep("^L",newmarnames)])
ord<-c(ord,sort(newmarnames[grep("^_",newmarnames)]))
LODmat5<-LODmat5[ord,ord]

assignsome(LODmat5,.7,.9)->LODmat5
nass<-5679
bounds<-c(0,(1:nass)[!duplicated(gsub("_.*","",rownames(LODmat5[1:nass,]),perl=TRUE))][-1]-1,nass)
amat<-c()
for(i in 1:(length(bounds)-1)){
	amat<-cbind(amat,apply(LODmat5[,(bounds[i]+1):bounds[i+1]],MAR=1,FUN=mean,na.rm=TRUE))
	}
amatmax<-apply(amat,MAR=1,FUN=max)
amatmaxind<-apply(amat,MAR=1,FUN=which.max)
amatdiff<-amatmax-apply(amat,MAR=1,FUN=function(x){sort(x,decreasing=TRUE)[2]})

image(LODmat5,zlim=c(1,10),col=colfunc(10),x=1:5777,y=1:5777)
abline(v=(1:nass)[!duplicated(gsub("_.*","",rownames(LODmat5),perl=TRUE))]-.5,col="blue",lwd=2)



image(LODmat5,zlim=c(1,10),col=colfunc(10),x=1:5777,y=1:5777)
abline(v=(1:5600)[!duplicated(gsub("_.*","",rownames(LODmat5),perl=TRUE))]-.5,col="blue",lwd=2)

subl<-grep("L6_|L8_",colnames(LODmat5))
image(LODmat5[subl,subl],zlim=c(1,10),col=colfunc(10))


# method: the distance measure to be used. This must be one of
#         ‘"euclidean"’, ‘"maximum"’, ‘"manhattan"’, ‘"canberra"’,
#         ‘"binary"’ or ‘"minkowski"’.  Any unambiguous substring can
#         be given.
#USE CANBERRA OR BINARY

subl<-grep("L25_",colnames(LODmat5))
hclust(dist(LODmat5[subl,subl],method="canberra"))->subcl
image(LODmat5[subl,subl][subcl$order,subcl$order],zlim=c(1,10),col=colfunc(10))

gsub("__","_",colnames(LODmat5)[1:5679]) %>% strsplit(split="_",x=.) %>% do.call(rbind,.) ->final.lgs

