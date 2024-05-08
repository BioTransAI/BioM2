#This is the preliminary version of the BioM2 parameter selection, which will be revised subsequently.


library(mlr3verse)
library(parallel)
library(caret)


TrainData=readRDS('your_data.rds')
# First column name must be 'label', and the rest are the features (e.g., CpGs).
#label cg21870274 cg09499020 cg16535257 cg00168193
#0     0.0057     0.0002    -0.0313     0.0002
#0    -0.0317    -0.0444    -0.0578    -0.0160
#1    -0.0341    -0.0541    -0.0056    -0.0230
#1     0.0811    -0.0029     0.0049     0.0274
#1    -0.0187     0.0475     0.1168     0.0169
#0    -0.0158     0.0032    -0.0173     0.0133



FeatureAnno=readRDS('your_Feature_annotation_data.rds')
# The data frame must contain the two column names 'ID' and 'entrezID' .

#ID entrezID symbol
#cg00000029     5934   RBL2
#cg00000109    64778 FNDC3B
#cg00000155   221927  BRAT1
#cg00000221   162282 ANKFN1
#cg00000236     7419  VDAC3
#cg00000289       87  ACTN1




pathlistDB=readRDS('your_pathway_list.rds')
# The name of each subset of the list is the ID of the pathway, and each subset contains a vector of gene entrezIDs.

#List of 15719
#$ GO:0000002: chr [1:31] "142" "291" "1763" "1890" ...
#$ GO:0000003: chr [1:1513] "2" "18" "49" "51" ...
#$ GO:0000012: chr [1:12] "142" "1161" "2074" "3981" ...
#$ GO:0000017: chr [1:2] "6523" "6524"
#$ GO:0000018: chr [1:131] "60" "86" "142" "604" ...
#$ GO:0000019: chr [1:7] "2068" "4292" "4361" "7014" ...
#$ GO:0000022: chr [1:9] "4926" "6795" "9055" "9212" ...
#$ GO:0000023: chr [1:3] "2548" "2595" "8972"



baseModel=function ( trainData, testData, predMode = c( "classification","regression","probability" ),
                     classifier,paramlist=NULL, inner_folds=10){
  predMode=match.arg(predMode)
  
  if(!is.null(testData)){
    if (colnames(trainData)[1] != "label") {
      stop("The first column of the 'trainData' must be the 'label'!")
    }
    if (colnames(testData)[1] != "label") {
      stop("The first column of the 'testData' must be the 'label'!")
    }
    
    if( predMode == 'probability'){
      classifier=paste0('classif.',classifier,'')
      trainData[,1]=as.factor(trainData[,1])
      testData[,1]=as.factor(testData[,1])
      trainData=as_task_classif(trainData,target='label')
      testData=as_task_classif(testData,target='label')
      model=lrn(classifier,predict_type = "prob")
      
      if(!is.null(paramlist)){
        at = auto_tuner(
          tuner = tnr("grid_search", resolution = 10, batch_size = 5),
          learner =  model,
          search_space = paramlist,
          resampling =rsmp("cv", folds =5),
          measure = msr("classif.acc")
        )
        at$train(trainData)
        model$param_set$values = at$tuning_result$learner_param_vals[[1]]
        model$train(trainData)
        predict=model$predict(testData)$prob[,2]
        return(predict)
      }else{
        sink(nullfile())
        model$train(trainData)
        sink()
        predict=model$predict(testData)$prob[,2]
        return(predict)
      }
    }
    else if( predMode == 'regression'){
      classifier=paste0('regr.',classifier,'')
      trainData[,1]=as.numeric(trainData[,1])
      testData[,1]=as.numeric(testData[,1])
      trainData=as_task_regr(trainData,target='label')
      testData=as_task_regr(testData,target='label')
      model=lrn(classifier)
      if(!is.null(paramlist)){
        at = auto_tuner(
          tuner = tnr("grid_search", resolution = 5, batch_size = 5),
          learner =  model,
          search_space = paramlist,
          resampling = rsmp("cv", folds =5),
          measure = msr("regr.mae")
        )
        at$train(trainData)
        model$param_set$values = at$tuning_result$learner_param_vals[[1]]
        model$train(trainData)
        predict=model$predict(testData)$response
        return(predict)
      }else{
        model$train(trainData)
        predict=model$predict(testData)$response
        return(predict)
      }
    }
  }
  else if(is.null(testData)){
    if (colnames(trainData)[1] != "label") {
      stop("The first column of the 'trainData' must be the 'label'!")
    }
    if( predMode == 'probability'){
      classifier=paste0('classif.',classifier,'')
      trainData[,1]=as.factor(trainData[,1])
      trainData=as_task_classif(trainData,target='label')
      model=lrn(classifier,predict_type = "prob")
      sink(nullfile())
      #set.seed(seed)
      rr=resample(trainData, model, rsmp("cv", folds = inner_folds))$prediction()
      sink()
      re=as.data.frame(as.data.table(rr))[,c(1,5)]
      re=re[order(re$row_ids),][,2]
      return(re)
    }
  } 
}

Step2_FeartureSelection=function(Step2_FeartureSelection_Method=NULL,data=NULL,label=NULL,cutoff=NULL,preMode=NULL,classifier=NULL){
  if(preMode=='probability' | preMode=='classification'){
    
    if(Step2_FeartureSelection_Method=='cor'){
      print(paste0('     Using  <<  correlation  >> , you choose cutoff ==>> ',cutoff))
      up=ifelse(classifier=='lda',0.999,100)
      corr=sapply(1:length(data),function(x) cor(data[[x]],label,method='pearson'))
      index=order(corr,decreasing = T)[1:cutoff]
      index2=which(corr>0)
      index=intersect(index,index2)
      print(paste0('     |> Final number of pathways >>>',length(index),'......Min correlation of pathways>>>',round(min(corr[index]),digits = 3)))
      #index=which(corr>cutoff & corr < up )
      return(index)
    }else if(Step2_FeartureSelection_Method=='wilcox.test'){
      
      data=do.call(cbind,data)
      data=cbind(label=label,data)
      data=as.data.frame(data)  
      
      data_0=data[which(data$label==unique(data$label)[1]),]
      data_1=data[which(data$label==unique(data$label)[2]),]
      pvalue=unlist(mclapply(2:ncol(data),function(x) wilcox.test(data_0[,x],data_1[,x])$p.value,mc.cores=60))
      if(cutoff < length(pvalue)){
        index=order(pvalue)[1:cutoff]
      }else{
        index=order(pvalue)
      }
      
      #index=which(pvalue < cutoff)
      print(paste0('     |> Final number of pathways >>>',length(index),'......Max p-value of pathways>>>',round(max(pvalue[index]),digits = 3)))
      return(index)
      
    }else if(Step2_FeartureSelection_Method=='RemoveHighcor'){
      if(!is.null(label)){
        corr=sapply(1:length(data),function(x) cor(data[[x]],label,method='pearson'))
        data=do.call(cbind,data)
        corm=cor(data)
        unindex=findCorrelation(corm,cutoff =cutoff)
        index=which(corr > 0)
        index=setdiff(index,unindex)
        print(paste0('     |> Final number of pathways >>>',length(index),'......Min correlation of pathways>>>',round(min(corr[index]),digits = 3)))
        return(index)
      }else{
        data=do.call(rbind,data)
        label=data[,1]
        corr=sapply(2:ncol(data),function(x) cor(data[,x],label,method='pearson'))
        index=which(corr > 0)
        data=data[,-1]
        corm=cor(data)
        unindex=findCorrelation(corm,cutoff =cutoff)
        index=setdiff(index,unindex)
        index=index+1
        return(index)
      }
    }else if(Step2_FeartureSelection_Method=='RemoveLinear'){
      if(!is.null(label)){
        corr=sapply(1:length(data),function(x) cor(data[[x]],label,method='pearson'))
        data=do.call(cbind,data)
        unindex=findLinearCombos(data)$remove
        index=which(corr > 0)
        index=setdiff(index,unindex)
        print(paste0('     |> Final number of pathways >>>',length(index),'......Min correlation of pathways>>>',round(min(corr[index]),digits = 3)))
        return(index)
      }else{
        data=do.call(rbind,data)
        label=data[,1]
        corr=sapply(2:ncol(data),function(x) cor(data[,x],label,method='pearson'))
        index=which(corr > 0)
        data=data[,-1]
        unindex=findLinearCombos(data)$remove
        index=setdiff(index,unindex)
        index=index+1
        return(index)
      }
      
    }else{
      if(!is.null(label)){
        up=ifelse(classifier=='lda',0.99,100)
        corr=sapply(1:length(data),function(x) cor(data[[x]],label,method='pearson'))
        
        print(paste0('     |> Final number of pathways >>>',length(order(corr,decreasing=T)[which(corr > 0 & corr < up)]),'......Min correlation of pathways>>>',round(min(corr[which(corr > 0 & corr < up)]),digits = 3)))
        index=which(corr > 0 & corr < up )
        return(index)
      }else{
        data=do.call(rbind,data)
        label=data[,1]
        corr=sapply(2:ncol(data),function(x) cor(data[,x],label,method='pearson'))
        index=which(corr > 0)
        index=index+1
        return(index)
      }
      
    }
  }
  if(preMode=='regression'){
    
  }
}


#Setting the basic parameters
nfolds=5
classifier='liblinear'
predMode = "probability"
PathwaySizeUp=200
PathwaySizeDown=20
MinfeatureNum_pathways=10
Add_UnMapped='Yes'
Unmapped_num=c(0,100,300,500,1000)
Add_FeartureSelection_Method='wilcox.test'
Inner_CV='None'
inner_folds=10
Step1_FeartureSelection_Method='cor'
Step2_FeartureSelection_Method='RemoveHighcor'
classifier2=NULL


# Window doesn't work, you need to use it on the server.
cores=30


#stage-1 feature_selection
off=c(0.3,0.1)

#stage-2 feature-selection
cutoff2=c(0,seq(0.7,0.9,0.1))

#add_unmapped_num
Unmapped_num=c(0,100,300,500,1000)
Unmap='w'

#Setting up the Learner(liblinear ->>  L2-log)
aa='liblinear'
classifier=aa
classifier2=aa

# Whether nested resampling was used 
# It is the most standard practice, but it consumes more time and the results are not always better
# If you want to use this set the parameter to 'Yes'

Inner_CV='None'
inner_folds=10


HOPE=list()
for(luck in 1:length(off)){
  set.seed(666)
  cutoff=off[luck]
  Resampling=createFolds(TrainData$label,k=nfolds)
  final=list()
  t1=Sys.time()
  pr=sum(TrainData$label==1)>sum(TrainData$label==0)
  for(xxx in 1:nfolds){
    print('Step1: ReadData')
    trainData=TrainData[unlist(Resampling[-xxx]),]
    testData=TrainData[unlist(Resampling[xxx]),]
    geneNum_pathways=sapply(1:length(pathlistDB),function(i) length(pathlistDB[[i]]))
    pathlistDB_sub=pathlistDB[which(geneNum_pathways > PathwaySizeDown & geneNum_pathways < PathwaySizeUp )]
    print(paste0('     |>Total number of pathways==>>',length(pathlistDB_sub)))
    
    
    print('Step2: FeartureSelection-features')
    print(paste0('      Using <<  correlation  >>',' ,and you choose cutoff:',cutoff))
    Cor=cor(trainData$label,trainData)
    Cor=ifelse(Cor>0,Cor,-Cor)
    names(Cor)=colnames(trainData)
    Cor_names=names(Cor)
    Cor_cutoff=Cor[which(Cor>cutoff)]
    Cor_cutoff_names=names(Cor_cutoff)
    
    featureAnno=FeatureAnno[FeatureAnno$ID %in% colnames(trainData),]
    MinfeatureNum_pathways2=MinfeatureNum_pathways+1
    featureNum_pathways=mclapply(1:length(pathlistDB_sub),function(x){
      id=c('label',featureAnno$ID[which(featureAnno$entrezID %in% pathlistDB_sub[[x]])])
      if(length(id)>MinfeatureNum_pathways2){
        id2=id[which(id %in% Cor_cutoff_names)]
        if(length(id2)<MinfeatureNum_pathways2){
          a=Cor[id]
          id2=names(a)[order(a,decreasing = T)[1:MinfeatureNum_pathways2]]
          return(id2)
        }else{
          return(id2)
        }
      }else{
        return(id)
      }
    } ,mc.cores=60)
    
    lens=sapply(1:length(featureNum_pathways),function(x) length(featureNum_pathways[[x]]))
    
    print('Step3: MergeData')
    trainDataList=mclapply(1:length(featureNum_pathways),function(x) trainData[,featureNum_pathways[[x]]] ,mc.cores=30)
    testDataList=mclapply(1:length(featureNum_pathways),function(x) testData[,featureNum_pathways[[x]]] ,mc.cores=30)
    names(trainDataList)=names(pathlistDB_sub)
    names(testDataList)=names(pathlistDB_sub)
    trainDataList=trainDataList[which(lens>MinfeatureNum_pathways)]
    testDataList=testDataList[which(lens>MinfeatureNum_pathways)]
    featureNum_pathways2=sapply(1:length(trainDataList),function(i2) length(trainDataList[[i2]]))
    print(paste0('     |> Total number of selected pathways==>>',length(trainDataList)))
    print(paste0('     |> Min features number of pathways==>>',min(featureNum_pathways2)-1,'.......','Max features number of pathways==>>',max(featureNum_pathways2)-1))
    
    
    #(FeartureSelection-pathways)
    print('Step4: Reconstruction')
    if(Inner_CV =='Yes'){
      print('     |> Using Inner CV ~ ~ ~')
      train=mclapply(1:length(trainDataList),function(i4) baseModel(trainData =trainDataList[[i4]],testData =NULL,predMode =predMode,classifier = classifier,inner_folds=inner_folds),mc.cores=cores)
      test=mclapply(1:length(testDataList),function(i5) baseModel(trainData =trainDataList[[i5]],testData =testDataList[[i5]],predMode =predMode,classifier = classifier),mc.cores=cores)
    }else{
      
      train=mclapply(1:length(trainDataList),function(i4) baseModel(trainData =trainDataList[[i4]],testData =trainDataList[[i4]],predMode = predMode ,classifier = classifier),mc.cores=cores)
      test=mclapply(1:length(testDataList),function(i5) baseModel(trainData =trainDataList[[i5]],testData =testDataList[[i5]],predMode = predMode ,classifier = classifier),mc.cores=cores)
      
    }
    
    print('     <<< Reconstruction Done! >>>     ')
    
    #(FeartureSelection-pathways)
    print('Step5: FeartureSelection-pathways')
    n=length(Unmapped_num)*length(cutoff2)
    Record=data.frame(Unmapped_num=1:n,learner=1:n,AUC=1:n,
                      ACC=1:n,PCCs=1:n,BAC=1:n,PRAUC=1:n,cutoff=1:n,MCCS=1:n)
    
    dbmap=unique(unlist(pathlistDB_sub))
    annomap=unique(featureAnno$entrezID)
    mapgene=intersect(annomap,dbmap)
    map=featureAnno$ID[which(featureAnno$entrezID %in% mapgene)]
    Unmapped_Train=trainData[,setdiff(colnames(trainData),map)]
    Unmapped_Test=testData[,setdiff(colnames(trainData),map)]
    print(paste0('Unmapped_num: ',ncol(Unmapped_Train)))
    #Unmapped_Train=trainData[,setdiff(colnames(trainData),featureAnno$ID)]
    #Unmapped_Test=testData[,setdiff(colnames(trainData),featureAnno$ID)]
    #print(paste0('Unmapped_num: ',ncol(Unmapped_Train)))
    Unmapped_0=Unmapped_Train[which(Unmapped_Train$label==unique(Unmapped_Train$label)[1]),]
    Unmapped_1=Unmapped_Train[which(Unmapped_Train$label==unique(Unmapped_Train$label)[2]),]
    if(Unmap=='cor'){
      Unmapped_pvalue= abs(cor(Unmapped_Train$label,Unmapped_Train[,-1]))
      Unmapped_pvalue=1/Unmapped_pvalue
    }else{
      Unmapped_0=Unmapped_Train[which(Unmapped_Train$label==unique(Unmapped_Train$label)[1]),]
      Unmapped_1=Unmapped_Train[which(Unmapped_Train$label==unique(Unmapped_Train$label)[2]),]
      Unmapped_pvalue=unlist(mclapply(2:ncol(Unmapped_Train),function(x) wilcox.test(Unmapped_0[,x],Unmapped_1[,x])$p.value,mc.cores=cores))
      
    }
    
    
    for(ii in  1:length(cutoff2)){
      if(cutoff2[ii]==0){
        index=Step2_FeartureSelection(Step2_FeartureSelection_Method='None',data=train,
                                      label=trainDataList[[1]]$label,cutoff=cutoff2[ii],preMode='probability',classifier =classifier )
      }else{
        index=Step2_FeartureSelection(Step2_FeartureSelection_Method=Step2_FeartureSelection_Method,data=train,
                                      label=trainDataList[[1]]$label,cutoff=cutoff2[ii],preMode='probability',classifier =classifier )
        
      }
      
      for(i in 1:length(Unmapped_num)){
        newtrain=do.call(cbind,train[index])
        colnames(newtrain)=names(trainDataList)[index]
        newtrain=cbind(label=trainDataList[[1]]$label,newtrain)
        
        newtest=do.call(cbind,test[index])
        colnames(newtest)=names(trainDataList)[index]
        newtest=cbind(label=testDataList[[1]]$label,newtest)
        colnames(newtest)=gsub(':','',colnames(newtest))
        colnames(newtrain)=gsub(':','',colnames(newtrain))
        if(Add_UnMapped=='Yes' & Unmapped_num[i] > 0 ){
          #Unmapped_Data=AddUnmapped(train=trainData,test=testData,Add_FeartureSelection_Method=Add_FeartureSelection_Method,
          #     Unmapped_num=Unmapped_num[i],len=ncol(newtrain),anno=featureAnno)
          if(is.null(Unmapped_num)){
            Unmapped_num=len-1
          }
          if(length(Unmapped_pvalue)<Unmapped_num[i]){
            Unmapped_Num=length(Unmapped_pvalue)
          }else{
            Unmapped_Num=Unmapped_num[i]
          }
          Unmapped_id=order(Unmapped_pvalue)[1:Unmapped_Num]
          Unmapped_id=Unmapped_id+1
          Unmapped_Train2=Unmapped_Train[,Unmapped_id]
          Unmapped_Test2=Unmapped_Test[,Unmapped_id]
          newtrain=cbind(newtrain,Unmapped_Train2)
          newtest=cbind(newtest,Unmapped_Test2)
          #newtrain=cbind(newtrain,Unmapped_Data$train)
          #newtest=cbind(newtest,Unmapped_Data$test)
          #print(paste0('     |> Merge PathwayFeature and AddFeature ==>>',ncol(newtrain)))
        }
        #print('Step6: Predict and Metric')
        
        if(is.null(classifier2)){
          classifier2=classifier
        }
        
        a=(ii-1)*length(Unmapped_num)+i
        result=baseModel(trainData=newtrain,testData=newtest,predMode ='probability',classifier = classifier2)
        testDataY=testDataList[[1]]$label
        pre=ifelse(result>0.5,1,0)
        accuracy_class1 <- sum(pre[testDataY == 1] == 1) / sum(testDataY == 1)
        accuracy_class0 <- sum(pre[testDataY == 0] == 0) / sum(testDataY == 0)
        Record[a,6]=(accuracy_class1 + accuracy_class0) / 2
        Record[a,1]=Unmapped_num[i]
        Record[a,2]=classifier2
        Record[a,5]=cor(testDataY,result,method='pearson')
        Record[a,9]=ModelMetrics::mcc(predicted = result, actual = testDataY,cutoff = 0.5)
        Record[a,3]=ROCR::performance(ROCR::prediction(result,testDataY),'auc')@y.values[[1]]
        testDataY=as.factor(testDataY)
        Record[a,7]=ifelse(pr,mlr3measures::prauc(testDataY, 1-result, '0'),mlr3measures::prauc(testDataY, result, '1'))
        pre=as.factor(pre)
        Record[a,4]=confusionMatrix(pre, testDataY)$overall['Accuracy'][[1]]
        Record[a,8]=cutoff2[ii]
        #print(i)
      }
      
    }
    
    final[[xxx]]=Record
  }
  hope=do.call(rbind,final)
  t2=Sys.time()
  print(t2-t1)
  hope=aggregate(hope[,c(3:7,9)],by=list(Unmapped_num=hope$Unmapped_num,cutoff2=hope$cutoff),mean)
  print(hope[,c(1:3,6:8)])
  print(paste0('cutoff : ',cutoff))
  print(luck)
  
  HOPE[[luck]]=hope
}


