library(mlr3verse)
library(parallel)
library(caret)
#(RE_m_all)
data=readRDS('your_data.rds')

head(data[,1:5])
dim(data)

set.seed(666)
final=list()
#len=15000
p<-"
 pr=sum(data$label==1)>sum(data$label==0)
 train=data[unlist(Resampling[-xxx]),]
 test=data[unlist(Resampling[xxx]),]
 COR=cor(train$label,train)
 COR=abs(COR)
 #Feature-selection
 n=c(500,1000,2000,5000,10000,ifelse(15000>ncol(data),ncol(data)-1,15000))
 n=n+1
 n2=length(a)*length(n)
 gg=data.frame(cor_numbers=1:n2,learner_name=1:n2,AUC=1:n2,ACC=1:n2,PCCs=1:n2,BAC=1:n2,PRAUC=1:n2)
 for(i in 1:length(a)){
  for(ii in 1:length(n)){
   id=order(COR,decreasing=T)[1:(n[ii])]
   #id=which(COR>n[ii])
   otrain=train[,c(id)]
   otest=test[,c(id)]
   otrain$label=as.factor(otrain$label)
   otest$label=as.factor(otest$label)
   t=as_task_classif(otrain,target='label')
   tt=as_task_classif(otest,target = 'label')
   l=lrn(a[i],predict_type = 'prob')
   l$train(t)
   result=l$predict(tt)$prob[,2]
   pre=ifelse(result>0.5,1,0)
   testDataY=as.numeric(l$predict(tt)$truth)-1
   accuracy_class1 <- sum(pre[testDataY == 1] == 1) / sum(testDataY == 1)
   accuracy_class0 <- sum(pre[testDataY == 0] == 0) / sum(testDataY == 0)
   x=(i-1)*length(n)+ii
   gg[x,1]=n[ii]
   gg[x,6]=(accuracy_class1 + accuracy_class0) / 2
   gg[x,2]=a[i]
   gg[x,3]=ROCR::performance(ROCR::prediction(result,testDataY),'auc')@y.values[[1]]
   gg[x,5]=cor(testDataY,result,method='pearson')
   testDataY=as.factor(testDataY)
   gg[x,7]=ifelse(pr,mlr3measures::prauc(testDataY, 1-result, '0'),mlr3measures::prauc(testDataY, result, '1'))
   pre=as.factor(pre)
   gg[x,4]=caret::confusionMatrix(pre, testDataY)$overall['Accuracy'][[1]]
   
  }
 }
 saveRDS(gg,outname[xxx])
"
set.seed(666)

#10repeat
for(zsj in 1:10){
  t1=Sys.time()
  
  #5k-folds
  nfolds=5
  
  Resampling=createFolds(data$label,k=nfolds)
  outname=sapply(1:5,function(x){paste0('ree',x,'.rds')})
  a=c('classif.C50','classif.gausspr','classif.gbm','classif.kknn','classif.liblinear',
      'classif.cv_glmnet','classif.ranger','classif.svm','classif.log_reg','classif.naive_bayes')
  y=sapply(1:nfolds,function(x) gsub('xxx',x,p))
  yy=mclapply(1:nfolds,function(x) eval(parse(text =y[x])),mc.cores=nfolds)
  
  f=list.files(pattern = 'ree')
  myfiles = do.call(rbind, lapply(f, function(x) readRDS(x)))
  myfiles2=aggregate(myfiles[,3:7],by=list(learner_name=myfiles$learner_name,cor_numbers=myfiles$cor_numbers),mean)
  myfiles2=myfiles2[order(myfiles2$AUC),]
  print(myfiles2)
  saveRDS(myfiles2,paste0('No',zsj,'.rds'))
  final[[zsj]]=myfiles2 
  t2=Sys.time()
  print(paste0('No',zsj,':',t2-t1))
}


final=do.call(rbind,final)
SD=aggregate(final[,3:7],by=list(learner_name=final$learner_name,cor_numbers=final$cor_numbers),sd)
Mean=aggregate(final[,3:7],by=list(learner_name=final$learner_name,cor_numbers=final$cor_numbers),mean)
result=list(mean=Mean,sd=SD)
