# source("http://bioconductor.org/biocLite.R") #To download DESeq package (you can comment these lines out, they only need to be run once ever)
 # biocLite( "MCMCglmm")
install.packages("MCMC.qpcr") # for trellisByGene function
library(MCMC.qpcr)
library(MCMC.OTU)
library(dplyr)
library(tidyr)
setwd("~/Dropbox/State Year 1/Symbiodinium/State Year 1 Symbiont info")
dat<-read.csv("041719EMits-OTU.csv",header=T)


####organizing original seqs table, where many columns of unique seq variants were BLASTn as the same hit
#Using dplyr to select columns that start with the sequence string since there are over 150 C15's, then summing the numbers in each row and creating a total column
datns<- dplyr::select(dat,starts_with("ns"))
datns <-cbind(datns, nssum = rowSums(datns))

datsp<- dplyr::select(dat,starts_with("sp"))
datsp<-datsp[,-c(1)]
datsp <- cbind(datsp, spsum = rowSums(datsp))

datD2<- dplyr::select(dat,starts_with("D2"))
datD2 <- cbind(datD2, D2sum = rowSums(datD2))

datD<- dplyr::select(dat,starts_with("D"))
datD<-datD[,-c(1:140)]#remove D2 columns
datD <- cbind(datD, Dsum = rowSums(datD))

datC1b<- dplyr::select(dat,starts_with("C1b"))
datC1b <- cbind(datC1b, C1bsum = rowSums(datC1b))

datC15<- dplyr::select(dat,starts_with("C15"))
datC15 <- cbind(datC15, C15sum = rowSums(datC15))

datB14<- dplyr::select(dat,starts_with("B14"))
datB14 <- cbind(datB14, B14sum = rowSums(datB14))

datB<- dplyr::select(dat,starts_with("B"))
datB<-datB[,-c(2743:2773)]#remove B14 columns
which(colnames(datB)=="B14")
datB <- cbind(datB, Bsum = rowSums(datB))

#column bind totalled seq count columns and sample column into new data frame
itseqs<-data.frame(cbind(dat[,c(1:4)], datsp$spsum, datD$Dsum, datD2$D2sum, datB$Bsum, datB14$B14sum, datC15$C15sum, datC1b$C1bsum))
head(itseqs)
write.csv(itseqs,file="2018_TotalITSseqs.csv",row.names=FALSE)

###MCMC.otu pipeline
itseqs$Genet<-factor(itseqs$Genet)
itseqs$Rep<-factor(itseqs$Rep)
colnames(itseqs)[5:11] <- c('Symbiodiniaceae spp.','D','D2','B','B14','C15','C1b')
colnames(itseqs)[1] <- c('sample')
pclidat<-itseqs[itseqs$Species=='PCLI',]
ofavdat<-itseqs[itseqs$Species=='OFAV',]
pclidat$Genet=factor(pclidat$Genet)
ofavdat$Genet=factor(ofavdat$Genet)
pclidat<-pclidat[,-c(2)]#remove Species columns
ofavdat<-ofavdat[,-c(2)]#remove Species columns

# numdat<-data.frame(apply(dat,2,function(x){as.numeric(x)}))
# row.names(numdat)<-row.names(dat)
# numdat<-numdat[,1:(ncol(numdat)-2)]
# colnames(numdat)<-c(1:ncol(numdat))
# goods=purgeOutliers(numdat, count.columns=c(1:(ncol(numdat))), otu.cut=0.001) 
goods=purgeOutliers(pclidat, count.columns=4:10, otu.cut=0.001) 
head(goods)

# what is the proportion of samples with data for these otus?
withData=apply(goods[,4:length(goods[1,])],2,function(x){sum(x>0)/length(x)})
hist(withData)
#apply(goods[,6:length(goods[1,])],2,function(x){sum(x>0)/length(x)})

# what percentage of global total counts each otu represents?
#apply(goods[,6:length(goods[1,])],2,function(x){sum(x)/sum(goods[,6:length(goods[1,])])})
props=apply(goods[,4:length(goods[1,])],2,function(x){sum(x)/sum(goods[,4:length(goods[1,])])})
barplot(sort(props,decreasing=T),xaxt="n",log="y")

# stacking the data; adjust otu.columns and condition.columns values for your data
gss=otuStack(data=goods,count.columns=c(4:length(goods[1,])),condition.columns=c(1:3))
#gss$Genet=factor(gss$Genet, levels=c("1","2","3","5","8","11","12","27","32","61","125","126","132"))

#write.csv(gss,file="202007)porII_ITSstacked.csv",row.names=FALSE)

# fitting the model. Replace the formula specified in 'fixed' with yours, add random effects if present. 
# See ?mcmc.otu for these and other options. 
#island by inner outer
###gon###
mm=mcmc.otu(
	fixed="Genet",
	random="Rep",
	data=gss,
	nitt=35000,thin=50,burnin=5000,singular.ok=T # a long MCMC chain to improve modeling of rare otus
	)
	
	# selecting the otus that were modeled reliably
# (otus that are too rare for confident parameter estimates are discarded) 
acpass=otuByAutocorr(mm,gss)

# calculating differences and p-values between all pairs of factor combinations
smm0=OTUsummary(mm,gss,summ.plot=FALSE) 

# adjusting p-values for multiple comparisons:
smmA=padjustOTU(smm0)

# significant otus at FDR<0.05:
sigs=signifOTU(smmA)
sigs

# plotting the significant ones
pdf(file="202208_2018st_pcli_sigOTUs-oneplot.pdf",width=6,height=6)
smm1=OTUsummary(mm,gss,otus=sigs,xgroup="Genet")
dev.off()
#smm1=otusummary(mm,gss,otus=sigs,xgroup="species")

# trellis by otu (some massaging to the summary object first to make trellisByGene eat it):
smg=smm1
names(smg$summary)[1]="gene"
pdf(file="202208_2018st_plci_sigOTUs.pdf",width=6,height=6)
trellisByGene(smg,xFactor="Genet",groupFactor="Genet",nrow=1,legendPos = 'none')
dev.off()

# table of log10-fold changes and p-values: this one goes into supplementary info in the paper
smmA$otuWise[sigs]
