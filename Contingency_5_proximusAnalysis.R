#Takes the binary strings (rdat_I.csv) and performs analysis using the Proximus algorithm. Exports a file (dataFromR2.txt) to the directory of the  “Contingency_6_nameAndFilter” sketch (see below) that contains the results of the analysis.
library("scales")
library("cba")
setwd(getSrcDirectory(function(){})[1])

lousedatI<-read.csv("Contingency_4_prepDataForR/rdat_I.csv")
tk<-ncol(lousedatI)
namsI<-colnames(lousedatI)
lousedatI$patDis<-rescale(lousedatI$patDis,to=c(-1,1))
tl<-nrow(lousedatI)
lousedat2<-matrix(nrow=(tk-2),ncol=tl)
colnames(lousedat2)<-lousedatI$sym
rownames(lousedat2)<-namsI[3:tk]
for(k in 3:tk){
	for(l in 1:tl){
		if(lousedatI[,namsI[k]][l]==1){
			lousedat2[(k-2),l]<-TRUE
		}else{
			lousedat2[(k-2),l]<-FALSE
		}
	}
}
save(lousedat2,file="lousedat2.rdata")
load("lousedat2.rdata")
proxI<-proximus(lousedat2,max.radius=7,min.size=1)
pSumI<-summary(proxI)
patIDs<-which(pSumI$pattern[,1]>1 & (pSumI$pattern[,4]<0.05 | pSumI$pattern[,6]>0.95))
tm<-length(patIDs)
cat("\tid\tsize\trelative error\tFrobenius norm\tJaccard similarity\tmembers",file="Contingency_6_nameAndFilter/dataFromR2.txt",sep="\n",append=FALSE)
m<-0
for(m in 1:tm){
	#/Users/Ian/documents/processing/NewPipe2/
	cat(m,file="Contingency_6_nameAndFilter/dataFromR2.txt",sep=" ",append=TRUE)
	cat("\t",file="Contingency_6_nameAndFilter/dataFromR2.txt",sep=" ",append=TRUE)
	cat(patIDs[m],file="Contingency_6_nameAndFilter/dataFromR2.txt",sep="\t",append=TRUE)
	cat("\t",file="Contingency_6_nameAndFilter/dataFromR2.txt",sep=" ",append=TRUE)
	cat(pSumI$pattern[patIDs[m],1],file="Contingency_6_nameAndFilter/dataFromR2.txt",sep="\t",append=TRUE)
	cat("\t",file="Contingency_6_nameAndFilter/dataFromR2.txt",sep=" ",append=TRUE)
	cat(pSumI$pattern[patIDs[m],4],file="Contingency_6_nameAndFilter/dataFromR2.txt",sep="\t",append=TRUE)
	cat("\t",file="Contingency_6_nameAndFilter/dataFromR2.txt",sep=" ",append=TRUE)
	cat(pSumI$pattern[patIDs[m],5],file="Contingency_6_nameAndFilter/dataFromR2.txt",sep="\t",append=TRUE)
	cat("\t",file="Contingency_6_nameAndFilter/dataFromR2.txt",sep=" ",append=TRUE)
	cat(pSumI$pattern[patIDs[m],6],file="Contingency_6_nameAndFilter/dataFromR2.txt",sep="\t",append=TRUE)
	cat("\t",file="Contingency_6_nameAndFilter/dataFromR2.txt",sep=" ",append=TRUE)
	cat(proxI$rownames[proxI$a[[patIDs[m]]]$x],file="Contingency_6_nameAndFilter/dataFromR2.txt",sep=" ",append=TRUE)
		cat("\n",file="Contingency_6_nameAndFilter/dataFromR2.txt",sep=" ",append=TRUE)
}