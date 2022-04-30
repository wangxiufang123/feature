a = list.files("C:/Users/wxf/Desktop/GSE5680_family.xml")
dir = paste("C:/Users/wxf/Desktop/GSE5680_family.xml/",a,sep="")            
n = length(dir)     
merge.data = read.table(file = dir[1],header=F,dec = ".")
for (i in 2:n){
  new.data = read.table(file = dir[i], header=F, dec = ".")[,2]
  merge.data = cbind(merge.data,new.data)
}
write.csv(merge.data,file = "C:/Users/wxf/Desktop/GSE5680_family.xml/merge.csv",row.names=F)  
dim(merge.data)

data=read.csv("C:/Users/wxf/Desktop/GSE5680_family.xml/merge.csv",header = F)
data0=t(data[,-1])
data0=as.data.frame(data0)
colnames(data0)=data[,1]
row.names(data0)=NULL
coln=colnames(data0)
y_index=which(coln=="1384040_at")
x_index=which(coln!="1384040_at")
p=dim(data0)[2]-1
n=dim(data0)[1]
d0=round(n/log(n))
d1=n
datax=data0[,-y_index]
datay=data0[,y_index]
#######################
c=12
GMC_gene=rep(0,p)
for (l in 1:p) {
  GMC_gene[l]=cal_GMC(datax[,l],datay,c)
}
GMC=abs(GMC_gene)
GMC_sel=which(rank(-GMC)<=d0)
colnames(datax)[GMC_sel]

GMC_sel1=which(rank(-GMC)<=d1)
colnames(datax)[GMC_sel1]

#################
CD=rep(0,p)
for (l in 1:p) {
  CD[l]=cal_hat_CD(datax[,l],datay)
}

CD_sel=which(rank(-as.numeric(CD))<=d0)
colnames(datax)[CD_sel]
CD_sel1=which(rank(-as.numeric(CD))<=d1)
colnames(datax)[CD_sel1]
##################
#SIRS
SIRS=rep(0,p)
for (l in 1:p) {
  SIRS[l]=cal_hat_SIRS(datax[,l],datay)
}
SIRS_sel=which(rank(-SIRS)<=d0)
colnames(datax)[SIRS_sel]
SIRS_sel1=which(rank(-SIRS)<=d1)
colnames(datax)[SIRS_sel1]
#################
#tSIS
tau=abs(cor(datax,datay,method = "kendall"))
tau=as.numeric(tau)
tSIS_sel=which(rank(-tau)<=d0)
colnames(datax)[tSIS_sel]
tSIS_sel1=which(rank(-tau)<=d1)
colnames(datax)[tSIS_sel1]
#################
#SIS
corr=abs(cor(datax,datay))
corr=as.numeric(corr)
SIS_sel=which(rank(-corr)<=d0)
colnames(datax)[SIS_sel]
SIS_sel1=which(rank(-corr)<=d1)
colnames(datax)[SIS_sel1]

##DC

dc=rep(0,p)
for (l in 1:p) {
 dc[l]=cal_hat_dcorr(datax[,l],datay)
}
DC_sel=which(rank(-as.numeric(dc))<=d0)
DC_sel1=which(rank(-as.numeric(dc))<=d1)
colnames(datax)[DC_sel]
colnames(datax)[DC_sel1]

##MDC
MDC=rep(0,p)
for (l in 200:p) {
  MDC[l]=cal_hat_MDC(datax[,l],datay)
}
MDC_sel=which(rank(-abs(MDC))<=d0)
MDC_sel1=which(rank(-abs(MDC))<=d1)
colnames(datax)[MDC_sel]
colnames(datax)[MDC_sel1]

######################################
##############################
xindex_bbs14=which(coln!="1384040_at")
bbs14=which(coln=="1384040_at")
d2=500
d3=25
M=100
datax14=data0[,-bbs14]
datay14=data0[,bbs14]
library(glmnet)
#############
c=20
d2=500
size_gmc=rep(0,M)
mspe_gmc=rep(0,M)
for (m in 1:M){
  GMC_bbs14=rep(0,p)
  for (l in 1:p) {
    GMC_bbs14[l]=cal_GMC(datax14[,l],datay14,c)
  }
  GMC=abs(GMC_bbs14)
  GMC_sel=which(rank(-GMC)<=d2)
  colnames(datax14)[GMC_sel]
  datax_screened=datax14[,GMC_sel]
  
  #split data 100:20
  set.seed(m)
  train_index=sample(1:120,size = 100,replace = F)
  trainx=datax_screened[train_index,]
  trainy=datay14[train_index]
  testx=datax_screened[-train_index,]
  testy=datay14[-train_index]

  model1=cv.glmnet(as.matrix(trainx),trainy,standardize = T ,nfolds = 5,keep = T)
  lam_opt=which(model1$lambda==model1$lambda.min)
  size_gmc[m]=model1$nzero[lam_opt]#model size
  pre_y=predict.glmnet(model1$glmnet.fit,as.matrix(testx),type="response",s=model1$lambda.min)
  mspe_gmc[m]=sum((pre_y-testy)^2)#prediction mean squared error (PE).
}
mean(size_gmc)
mean(mspe_gmc)

#################
size_cd=rep(0,M)
mspe_cd=rep(0,M)
for (m in 1:M){
  CD_bbs14=rep(0,p)
  for (l in 1:p) {
    CD_bbs14[l]=cal_hat_CD(datax14[,l],datay14)
  }
  
  CD_sel=which(rank(-as.numeric(CD_bbs14))<=d2)
  datax_screened1=datax14[,CD_sel]
  set.seed(m)
  train_index=sample(1:120,size = 100,replace = F)
  trainx1=datax_screened1[train_index,]
  trainy=datay14[train_index]
  testx1=datax_screened1[-train_index,]
  testy=datay14[-train_index]
  
  model2=cv.glmnet(as.matrix(trainx1),trainy,standardize = T ,nfolds = 5,keep = T)
  lam_opt2=which(model2$lambda==model2$lambda.min)
  size_cd[m]=model2$nzero[lam_opt2]#model size
  pre_y2=predict.glmnet(model2$glmnet.fit,as.matrix(testx1),type="response",s=model2$lambda.min)
  mspe_cd[m]=sum((pre_y2-testy)^2)
  
}
mean(size_cd)
mean(mspe_cd)


#SIS
size_sis=rep(0,M)
mspe_sis=rep(0,M)
for (m in 1:M){
  corr=abs(cor(datax14,datay14))
  corr=as.numeric(corr)
  SIS_sel=which(rank(-corr)<=d2)
  
  datax_screened_SIS=datax14[,SIS_sel]
  set.seed(m)
  train_index=sample(1:120,size = 100,replace = F)
  trainx_SIS=datax_screened_SIS[train_index,]
  trainy=datay14[train_index]
  testx_SIS=datax_screened_SIS[-train_index,]
  testy=datay14[-train_index]
  
  model3=cv.glmnet(as.matrix(trainx_SIS),trainy,standardize = T ,nfolds = 5,keep = T)
  lam_opt3=which(model3$lambda==model3$lambda.min)
  size_sis[m]=model3$nzero[lam_opt3]#model size
  pre_y3=predict.glmnet(model3$glmnet.fit,as.matrix(testx_SIS),type="response",s=model3$lambda.min)
  mspe_sis[m]=sum((pre_y3-testy)^2)#
  
}

mean(size_sis)
mean(mspe_sis)


#tSIS
size_tsis=rep(0,M)
mspe_tsis=rep(0,M)
for (m in 1:M){
  tau=abs(cor(datax14,datay14,method = "kendall"))
  tau=as.numeric(tau)
  tSIS_sel=which(rank(-tau)<=d2)
  
  datax_screened_tSIS=datax14[,tSIS_sel]
  set.seed(m)
  train_index=sample(1:120,size = 100,replace = F)
  trainx_tSIS=datax_screened_tSIS[train_index,]
  trainy=datay14[train_index]
  testx_tSIS=datax_screened_tSIS[-train_index,]
  testy=datay14[-train_index]
  
  model4=cv.glmnet(as.matrix(trainx_tSIS),trainy,standardize = T ,nfolds = 5,keep = T)
  lam_opt4=which(model4$lambda==model4$lambda.min)
  size_tsis[m]=model4$nzero[lam_opt4]#model size
  pre_y4=predict.glmnet(model4$glmnet.fit,as.matrix(testx_tSIS),type="response",s=model4$lambda.min)
  mspe_tsis[m]=sum((pre_y4-testy)^2)#
}
mean(size_tsis)
mean(mspe_tsis)


#SIRS
size_sirs=rep(0,M)
mspe_sirs=rep(0,M)
for (m in 1:M){
  SIRS=rep(0,p)
  for (l in 1:p) {
    SIRS[l]=cal_hat_SIRS(datax14[,l],datay14)
  }
  SIRS_sel=which(rank(-SIRS)<=d2)
  
  datax_screened_SIRS=datax14[,SIRS_sel]
  set.seed(m)
  train_index=sample(1:120,size = 100,replace = F)
  trainx_SIRS=datax_screened_SIRS[train_index,]
  trainy=datay14[train_index]
  testx_SIRS=datax_screened_SIRS[-train_index,]
  testy=datay14[-train_index]
  
  model5=cv.glmnet(as.matrix(trainx_SIRS),trainy,standardize = T ,nfolds = 5,keep = T)
  lam_opt5=which(model5$lambda==model5$lambda.min)
  size_sirs[m]=model5$nzero[lam_opt5]#model size
  pre_y5=predict.glmnet(model5$glmnet.fit,as.matrix(testx_SIRS),type="response",s=model5$lambda.min)
  mspe_sirs[m]=sum((pre_y5-testy)^2)#
}
mean(size_sirs)
mean(mspe_sirs)


##DC
size_dc=rep(0,M)
mspe_dc=rep(0,M)
for (m in 1:M){
  DC_bbs14=rep(0,p)
  for (l in 1:p) {
    DC_bbs14[l]=cal_hat_dcorr(datax14[,l],datay14)
  }
  
  DC_sel=which(rank(-as.numeric(DC_bbs14))<=d2)
  
  datax_screened_DC=datax14[,DC_sel]
  set.seed(m)
  train_index=sample(1:120,size = 100,replace = F)
  trainx_DC=datax_screened_DC[train_index,]
  trainy=datay14[train_index]
  testx_DC=datax_screened_DC[-train_index,]
  testy=datay14[-train_index]
  
  model6=cv.glmnet(as.matrix(trainx_DC),trainy,standardize = T ,nfolds = 5,keep = T)
  lam_opt6=which(model6$lambda==model6$lambda.min)
  size_dc[m]=model6$nzero[lam_opt6]#model size
  pre_y6=predict.glmnet(model6$glmnet.fit,as.matrix(testx_DC),type="response",s=model6$lambda.min)
  mspe_dc[m]=sum((pre_y6-testy)^2)#
  
}
mean(size_dc)
mean(mspe_dc)


##MDC
MDC_bbs14=rep(0,p)
for (l in 26114:p) {
  MDC_bbs14[l]=cal_hat_MDC(datax14[,l],datay14)
}
MDC_sel=which(rank(-abs(MDC_bbs14))<=d2)

datax_screened_MDC=datax14[,MDC_sel]
size_mdc=rep(0,200)
mspe_mdc=rep(0,200)
for (m in 1:200) {
  set.seed(m)
  train_index=sample(1:120,size = 100,replace = F)
  trainx_MDC=datax_screened_MDC[train_index,]
  testx_MDC=datax_screened_MDC[-train_index,]
  trainy=datay14[train_index]
  testy=datay14[-train_index]
  model7=cv.glmnet(as.matrix(trainx_MDC),trainy,standardize = T ,nfolds = 5,keep = T)
 
  lam_opt7=which(model7$lambda==model7$lambda.min)
  size_mdc[m]=model7$nzero[lam_opt7]#model size
  pre_y7=predict.glmnet(model7$glmnet.fit,as.matrix(testx_MDC),type="response",s=model7$lambda.min)
  mspe_mdc[m]=sum((pre_y7-testy)^2)#

}
mean(size_mdc)
mean(mspe_mdc)





###RF
library(randomForest)
train_index=sample(1:120,size = 100,replace = F)
trainx_rf=datax14[train_index,]
testx_rf=datax14[-train_index,]
trainy=datay14[train_index]
testy=datay14[-train_index]
train=cbind(trainy,trainx_rf)
modelrf=randomForest(trainy~.,data = train,importance=T,xtest=testx_rf,ytest=testy)

#a new corr
library(dplyr)
cal_ncorr=function(x,y){
  n=length(x)
  dataf=data.frame(x=x,y=y)
  dataf_slice=arrange(dataf,x)
  dataf_slice[,2]=as.numeric(rank(dataf_slice[,2]))
  y1=dataf_slice[,2]
  s=sum(abs(dataf_slice[-1,2]-y1[-n]))
  ncorr=1-3*s/(n^2-1)
  return(ncorr)
}

ncorr_gene=rep(0,p)
for (l in 1:p) {
  ncorr_gene[l]=cal_ncorr(datax[,l],datay)
}

d0=25
d1=120
ncorr=abs(ncorr_gene)
ncorr_sel=which(rank(-ncorr)<=d0)
colnames(datax)[ncorr_sel]

ncorr_sel1=which(rank(-ncorr)<=d1)
colnames(datax)[ncorr_sel1]


d2=500
size_nc=rep(0,M)
mspe_nc=rep(0,M)

  nc_bbs14=rep(0,p)
  for (l in 1:p) {
    nc_bbs14[l]=cal_ncorr(datax14[,l],datay14)
  }
  
  nc_sel=which(rank(-as.numeric(DC_bbs14))<=25)
  
  datax_screened_nc=datax14[,nc_sel]
for (m in 1:M){
  set.seed(m)
  train_index=sample(1:120,size = 100,replace = F)
  trainx_nc=datax_screened_nc[train_index,]
  trainy=datay14[train_index]
  testx_nc=datax_screened_nc[-train_index,]
  testy=datay14[-train_index]
  
  model=cv.glmnet(as.matrix(trainx_nc),trainy,standardize = T ,nfolds = 5,keep = T)
  lam_opt=which(model$lambda==model$lambda.min)
  size_nc[m]=model$nzero[lam_opt]#model size
  pre_y=predict.glmnet(model$glmnet.fit,as.matrix(testx_nc),type="response",s=model$lambda.min)
  mspe_nc[m]=sum((pre_y-testy)^2)#
  
}
mean(size_nc)
mean(mspe_nc)
