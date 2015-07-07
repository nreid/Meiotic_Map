#input vcf lines (single parent cross markers only. e.g. D1.10)
#input inferred cross type (e.g. "aaab")
	#infer parent gametes
#output three tables
	#1 parent gamete identity
	#2 offspring genotype (0=homo,1=het)
	#3 likelihood ratio to next most likely genotype

altervcf<-function(vcf,crosstype){
	
	vcf2<-vcf
	vcf<-as.matrix(vcf)
	#out1 will be genotypes
	out1<-vcf[,16:105]
	out1<-gsub(":.*","",out1)
	#out2 will be likelihood ratios
	out2<-vcf[,16:105]
	out2<-gsub(".*:","",out2)
	
	for(i in 1:length(crosstype)){
		
		tmplik<-do.call(rbind,strsplit(split=",",x=out2[i,]))
		class(tmplik)<-"numeric"
		
		if(crosstype[i]=="abaa"|crosstype[i]=="aaab"){
			
			tmplik<-tmplik[,1:2]
			
			out1[i,]<-gsub("0/0","0",out1[i,])
			out1[i,]<-gsub("0/1","1",out1[i,])
			out1[i,]<-gsub("1/1","1",out1[i,])
			
			for(j in 1:90){
				if(out1[i,j]==0){
					out2[i,j]<-tmplik[j,1]-tmplik[j,2]
					}else{out2[i,j]<-tmplik[j,2]-tmplik[j,1]}
				}
			
			}

		if(crosstype[i]=="abbb"|crosstype[i]=="bbba"){
			
			tmplik<-tmplik[,3:2]
			
			out1[i,]<-gsub("0/0","1",out1[i,])
			out1[i,]<-gsub("0/1","1",out1[i,])
			out1[i,]<-gsub("1/1","0",out1[i,])
							
			for(j in 1:90){
				if(out1[i,j]==0){
					out2[i,j]<-tmplik[j,1]-tmplik[j,2]
					}else{out2[i,j]<-tmplik[j,2]-tmplik[j,1]}
				}
						
			}


		
		}
	
	class(out1)<-"numeric"
	class(out2)<-"numeric"

	#out3 will be parental gamete identity
	out3<-out1
	for(i in 2:length(crosstype)){
		
		ind<-cor(out1[1,],out1[i,],use="complete")
		if(ind<0){
			switch<-out1[i,]==0
			out3[i,switch]<-1
			out3[i,!switch]<-0
			}
		
		
		}
	
	list(out1,out3,out2,as.numeric(vcf[,2]),vcf2,crosstype)
	
	}
	
sub<-grepl("Scaffold9860_",rownames(vcf))&onemapcross=="D2.15"
altervcf(vcf[sub,],final.parent[sub])->test	

plot(NULL,ylim=c(71,90),xlim=c(1,length(test[[4]])))
#plot(NULL,ylim=c(1,20),xlim=c(0,3))

for(i in 1:length(test[[3]][,1])){
	al<-test[[3]]
	al[is.na(al)]<-0
	al[al<0]<-0
	al[al>10]<-10
	al2<-test[[2]]
	al2[is.na(al2)]<-0
	points(x=rep(i,90),y=1:90,col=rgb(al2[i,],0,0,al[i,]/max(al,na.rm=TRUE)),pch=test[[1]][i,]*2+18,cex=5)
	}


#scaffold9860 
#M -14 -49 -70 -74 -83, 5443516 bases, 9.567958e-09, 1.04mb/cm
#F -18, 156093 bases, 1.805211e-09, 5539519mb/cm


