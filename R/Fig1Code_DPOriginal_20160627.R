setwd("C:/Users/David Paez/Documents/Bighorn")
# read in data - this is the first worksheet in 'LostineMoviPastLung_150630.xlsx' with shortened column names.
bdat <- read.csv("LostineSample2016Final.csv")

####Data formatting
for(i in 1:ncol(bdat)) if(class(bdat[,i])=="character") bdat[bdat[,i]=="",i]<-NA
for(i in 1:ncol(bdat)) if(class(bdat[,i])=="character") {
	cat("\n");cat(names(bdat)[i],"\n");print(table(bdat[,i]))}
for(i in grep("date",names(bdat))) bdat[,i] <- as.Date(bdat[,i],"%m/%d/%Y") # encode date
######


################Sum Stats for paper
length(unique(bdat$id)[table(bdat$id)>1])###individuals with more than 1 sample
length(unique(bdat$id)[table(bdat$id)==1])###individuals with only 1 sample
length(unique(bdat$id))#total number of unique ids

##number of sheep that were less than 4 at capture
head(bdat)
sheep <- unique(bdat$id)

MinAgesOut <-matrix(nrow=0, ncol=3)
colnames(MinAgesOut) <- c("id", "entry_date", "age_entry")
for (i in 1:length(sheep)){
tmp <- bdat[bdat$id == sheep[i],c(1,20,19)]
MinAgesOut<- rbind(MinAgesOut, tmp[which(tmp$entry_date==min(tmp$entry_date)),][1,])
}

length(which(MinAgesOut$age_entry<4))
length(which(MinAgesOut$age_entry>4 & MinAgesOut$age_entry<5))
length(which(MinAgesOut$age_entry>5))

AgeOut<-MinAgesOut[order(MinAgesOut$age_entry),]
row.names(AgeOut)<-1:nrow(AgeOut)	
#AgeOut
#############################


##############removing indeterminate samples
bdat <- bdat[!is.na(bdat$movi_qpcr),]

bdat$init_cap_date <- tapply(bdat$capture_date,bdat$id,min)[bdat$id]
bdat$days_since_init_cap <- as.numeric(bdat$capture_date)-bdat$init_cap_date
bdat$yrs_since_1st_cap <- bdat$days_since_init_cap/365
bdat$movi <- as.numeric(bdat$movi_qpcr=="Detected")
bdat$movi[bdat$movi_qpcr=="Indeterminate"]<-NA
bdat$nsamp <- tapply(bdat$movi,bdat$id,length)[bdat$id]

bdat$class <- tapply(bdat$movi,bdat$id,function(x){
  if(sum(!is.na(x))<=1) return(NA)
  x <- x[!is.na(x)]
  if(all(x>0)) return(1)
  if(all(x<1)) return(0)
  if(max(x)==1 &  min(x)==0) return(.5)
  return(2)
})[bdat$id]

bdat$class_color <- "blue"
bdat$class_color[bdat$class==1] <- "red"
bdat$class_color[bdat$class==.5] <- "dark gray"
bdat$class_color[is.na(bdat$class)] <- NA


######################################################################

# figure showing captures and movi status by month ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bdat$month<-format(bdat$capture_date,"%b-%Y")
bdat$month<-factor(bdat$month,levels=unique(format(sort(bdat$capture_date),"%b-%Y"))) # this is needed for plotting
table(bdat$movi,bdat$month)


###############################################################
###############################################################

# figure showing movi status by individual over time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

indivs <- unique(bdat$id[order(bdat$age_entry)])
bdat$ID <- factor(bdat$id,levels=indivs)
bdat$ID2 <- as.character(bdat$ID)
bdat$id_no <- as.numeric(bdat$ID)###
bdat$day <- bdat$days_since_init_cap

id_labs <- levels(bdat$ID)
id_pts <- 1:max(bdat$id_no)
summary(bdat$day)
bdat$pch<-16
bdat$pch[is.na(bdat$movi)]<-NA

bdat$pch2<-1
bdat$pch2[is.na(bdat$movi)]<-NA

bdat$col<-"beige"
bdat$col[is.na(bdat$movi)]<-"grey88"
bdat$col[!is.na(bdat$movi) & bdat$movi==1]<-"firebrick"
bdat <- bdat[order(bdat$capture_date),]
bdat$Movi <- c("neg","pos")[bdat$movi+1]
bdat$Movi[is.na(bdat$movi)]<-"ind"
bdat$Movi <- factor(bdat$Movi,levels=c("neg","pos","ind"))
n<-table(bdat$ID2[!is.na(bdat$movi)])

brks <- seq(0,20,by=2)
mids <- seq(1,19,by=2)
bdat$abin<-cut(bdat$age_capture,brks)


bdat$class <- tapply(bdat$movi,bdat$id,function(x){
  if(length(x)<=1) return(NA)
  x <- x[!is.na(x)]
  if(all(x>0)) return(1)
  if(all(x<1)) return(0)
  if(max(x)==1 &  min(x)==0) return(.5)
  return(2)
})[bdat$id]

bdat2 <- bdat[bdat$nsamp>1,];sum(is.na(bdat2$class[!duplicated(bdat2$id)]));
	p2 <- table(bdat2$class[!duplicated(bdat2$id)]);p2<-c(0,cumsum(p2/sum(p2)))

bdat3 <- bdat[bdat$nsamp>2,];sum(is.na(bdat3$class[!duplicated(bdat2$id)]));
	table(bdat3$class);p3 <- table(bdat3$class[!duplicated(bdat3$id)]);
	p3<-c(0,cumsum(p3)/sum(p3))

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

##########Main figure
#png("Fig1SLW.png",10,8,units="in",res=300)

NI<-1:11
INT <-12:(35)
PI <- (36):(48)

par(oma=c(3.5,4,.5,.5))
for(i in 3:1){
  tmp <- bdat[bdat$class==c(0,.5,1)[i] & !is.na(bdat$class) & bdat$nsamp>1 & !is.na(bdat$movi_qpcr),] # data by class
  if(i==3){
	tmp <- tmp[tmp$id!="99L03",]
	}
  indivs<-unique(as.character(tmp$id)[order(tmp$age_entry)])
  tmp$ID <- factor(tmp$id,levels=indivs)
  tmp$ID2 <- as.character(tmp$ID)
  tmp$id_no<-as.numeric(tmp$ID)
  id_labs <- levels(tmp$ID)
  id_pts <- 1:max(tmp$id_no)

	#tmp$cex <- ifelse(is.na(tmp$Lung), .6, ifelse(tmp$Lung == 0, 1, ifelse(tmp$Lung >=1 & tmp$Lung <= 5, 1.4,
	#		ifelse(tmp$Lung >= 6 & tmp$Lung <= 10, 1.8, ifelse(tmp$Lung > 11, 2.3, 10)))))##if using cex for lunworms
  
  
  par(fig=c(0,1,p2[i],p2[i+1]),mar=c(.25,.5,.25,.5),new=i<3)
  plot(1~1,type="n",ylim=c(1,max(tmp$id_no)),xlim=c(-.25,max(bdat$yrs_since_1st_cap)),yaxt="n",xlab="",ylab="",xaxt="n")
  points(id_no~yrs_since_1st_cap,col=col,pch=pch,data=tmp, cex=1.5)
  points(id_no~yrs_since_1st_cap,data=tmp, pch=pch2, cex=1.5)
  #text(y=tmp$id_no, x=tmp$yrs_since_1st_cap, labels= tmp$Lung, cex=.55, font=2)##if using counts for lungworms


if(0){###if using individual id, change to 1
	axis(side=2,at=id_pts,labels=id_labs,las=1,cex.axis=.8,xpd=NA)
}
else{
if(i==1){
	axis(side=2,at=id_pts,labels=NI,las=1,cex.axis=.8,xpd=NA)
	legend("bottomright", legend="Never detected", bty="n")
  }
if(i==2){
	axis(side=2,at=id_pts,labels=INT,las=1,cex.axis=.8,xpd=NA)
	legend("bottomright", legend="Intermittently detected", bty="n")
  }

if(i==3){
	axis(side=2,at=id_pts,labels=PI,las=1,cex.axis=.8,xpd=NA)
	legend("bottomright", legend="Always detected", bty="n")

  }
}

  if(i==3){ legend("topright",legend=c("Detected","Not detected"),
			pch=c(21,21),inset=c(.02,.02),xpd=NA,cex=.8,
	title=expression(italic('M. ovipneumoniae')),
	pt.bg=c("firebrick","beige"))
	
	}
  if(i==1) {axis(side=1,xpd=NA);title(xlab="Years since initial capture",xpd=NA,line=2.7, cex.lab=2)}
}
#mtext("Individual",side=2,line=2.25,outer=T, cex=2)



#dev.off()
####


#############Main figure ends here



if(0){
#####################3Lung thing that raina wanted
#########################################removing intermittend obs and making lung larvae a numeric variabe
if(0){
bdatx<- bdat[bdat$movi_qpcr!="Indeterminate",]
row.names(bdatx)<-1:nrow(bdatx)
bdatx <- bdatx[!is.na(bdaxt$id),]

Lung<-as.numeric(as.character(bdat$lung_larvae))
bdat$Lung<-ifelse(bdat$lung_larvae=="negative" | bdat$lung_larvae=="Negative", 0,
	 Lung)
}
######################################################################


Lungvsclass<-bdat[bdat$class==c(0,.5,1) & !is.na(bdat$class) & bdat$nsamp>2,] # data by class
par(mfrow=c(1,2))
boxplot(Lungvsclass$Lung~Lungvsclass$class, ylab="lung worm count")
boxplot(log(Lungvsclass$Lung+0.01)~Lungvsclass$class, ylab="log lung worm count")

library(lme4)

bdat$class<-as.factor(bdat$class)
mod<-glmer(Lung~ yrs_since_1st_cap + (1|class/id), "poisson", data=bdat[!is.na(bdat$class),])
summary(mod)

mod1<-glmer(Lung~ yrs_since_1st_cap + (1|class), "poisson", data=bdat[!is.na(bdat$class),])
summary(mod1)

mod2<-glmer(Lung~ yrs_since_1st_cap + (1|id), "poisson", data=bdat[!is.na(bdat$class),])
summary(mod2)

anova(mod, mod1, test="Chi")
anova(mod, mod2, test="Chi")

####Lung vs age

modx <- glmer(Lung~age_capture + (1|id), data=bdat[!is.na(bdat$class),], family="poisson")
summary(modx)

pchisq(2*(as.numeric(logLik(mod))-as.numeric(logLik(mod2))), df=1, lower=F)

###############################################
##########Excluding less than 2 obs#################
###############################################
###############################################
png("fig2_capture_history_by_indiv_grouped_n_gt_2.png",5,6,units="in",res=150)
par(oma=c(3.5,4,.5,.5))
for(i in 3:1){
  tmp <- bdat[bdat$class==c(0,.5,1)[i] & !is.na(bdat$class) & bdat$nsamp>2,] # data by class
  indivs<-unique(as.character(tmp$id)[order(tmp$age_entry)])
  tmp$ID <- factor(tmp$id,levels=indivs)
  tmp$ID2 <- as.character(tmp$ID)
  tmp$id_no<-as.numeric(tmp$ID)
  id_labs <- levels(tmp$ID)
  id_pts <- 1:max(tmp$id_no)
  
  
  par(fig=c(0,1,p3[i],p3[i+1]),mar=c(.25,.5,.25,.5),new=i<3)
  plot(1~1,type="n",ylim=c(1,max(tmp$id_no)),xlim=c(-.25,max(bdat$yrs_since_1st_cap)),yaxt="n",xlab="",ylab="",xaxt="n")
  points(id_no~yrs_since_1st_cap,col=col,pch=16,data=tmp)
  points(id_no~yrs_since_1st_cap,data=tmp)
  axis(side=2,at=id_pts,labels=id_labs,las=1,cex.axis=.8,xpd=NA)
  
  if(i==3) legend("topright",legend=c("detected","not detected","indeterminate"),col=c(2,4,1),pch=16,inset=c(.02,.02),xpd=NA,cex=.8,title="MOVI status")
  if(i==1) {axis(side=1,xpd=NA);title(xlab="years since intial capture",xpd=NA,line=2.25)}
}
mtext("individual",side=2,line=3.75,outer=T)
dev.off()

par(mar=c(5,5,1,1))
plot(1~1,type="n",ylim=c(1,max(bdat$id_no)),xlim=c(-.25,max(bdat$yrs_since_1st_cap)),yaxt="n",xlab="years since initial capture",ylab="")
# for(i in id_pts){
#   tmp <- as.data.frame(rbind(bdat[bdat$id_no==i,c("day","movi","col")],c(max(bdat$day),NA,1)))
#   for(j in 1:(nrow(tmp)-1)) segments(x0=tmp$day[j],y0=i,x1=tmp$day[j+1],y1=i,col=tmp$col[j],lwd=3)
# }
points(id_no~yrs_since_1st_cap,col=col,pch=16,data=bdat)
points(id_no~yrs_since_1st_cap,data=bdat)
axis(side=2,at=id_pts,labels=id_labs,las=1,cex.axis=.8)
mtext("individual",side=2,line=3.75)
# axis(side=1,at=axis_days,labels=axis_months,las=3)
# mtext("date",side=1,line=5.5)
legend("topright",legend=c("detected","not detected","indeterminate"),col=c(2,4,1),pch=16,inset=c(.02,.01),xpd=NA,cex=.8,title="MOVI status")
dev.off()

} ##if 0 loop