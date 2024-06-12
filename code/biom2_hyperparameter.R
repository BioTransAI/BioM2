#This is the preliminary version of the BioM2 parameter selection, which will be revised subsequently.


library(mlr3verse)
library(parallel)
library(caret)
library(BioM2)

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


#' BioM2 Hyperparametric Combination
#'
#' @param TrainData The input training dataset. The first column
#' is the label or the output. For binary classes,
#' 0 and 1 are used to indicate the class member.
#' @param pathlistDB A list of pathways with pathway IDs and their
#' corresponding genes ('entrezID' is used).
#' For details, please refer to ( data("GO2ALLEGS_BP") )
#' @param FeatureAnno The annotation data stored in a data.frame for probe
#' mapping. It must have at least two columns named 'ID' and 'entrezID'.
#' (For details, please refer to data( data("MethylAnno") )
#' @param resampling Resampling in mlr3verse.
#' @param nfolds k-fold cross validation ( Only supported when TestData = NULL )
#' @param classifier Learners in mlr3
#' @param predMode The prediction mode. Available options are
#' c('probability', 'classification').
#' @param PathwaySizeUp The upper-bound of the number of genes in each
#' biological pathways.
#' @param PathwaySizeDown The lower-bound of the number of genes in each
#' biological pathways.
#' @param MinfeatureNum_pathways The minimal defined pathway size after mapping your
#' own data to pathlistDB(KEGG database/GO database).
#' @param Add_UnMapped Whether to add unmapped probes for prediction
#' @param Unmapped_num The number of unmapped probes
#' @param Add_FeartureSelection_Method Feature selection methods.
#' @param Inner_CV Whether to perform a k-fold verification on the training set.
#' @param inner_folds k-fold verification on the training set.
#' @param Stage1_FeartureSelection_Method Feature selection methods.
#' @param stage1_cutoff The cutoff used for feature selection threshold. It can be any value
#' between 0 and 1.
#' @param Stage2_FeartureSelection_Method Feature selection methods.
#' @param stage2_cutoff The cutoff used for feature selection threshold. It can be any value
#' between 0 and 1.
#' @param classifier2 Learner for stage 2 prediction(if classifier2==NULL,then it is the same as the learner in stage 1.)
#' @param cores The number of cores used for computation.
#' @param verbose Whether to print running process information to the console
#'
#'
#' @return A data frame contains hyperparameter results
#' @export
#' @import ROCR
#' @import caret
#' @importFrom utils head 
#' @importFrom stats  wilcox.test 

#'
HyBioM2=function(TrainData=NULL,pathlistDB=NULL,FeatureAnno=NULL,resampling=NULL,nfolds=5,classifier='liblinear', predMode = "probability",
                 PathwaySizeUp=200,PathwaySizeDown=20,MinfeatureNum_pathways=10,
                 Add_UnMapped=TRUE,Add_FeartureSelection_Method='wilcox.test',Unmapped_num=300,
                 Inner_CV=TRUE,inner_folds=10,
                 Stage1_FeartureSelection_Method='cor',stage1_cutoff=0.3,
                 Stage2_FeartureSelection_Method='RemoveHighcor',stage2_cutoff=0.8,
                 classifier2=NULL,cores=1,verbose=TRUE){
  re=list()
  if(verbose)print('===================HyBioM2==================')
  for(c1 in 1:length(classifier)){
    stage1_learner=classifier[c1]
    HOPE=list()
    t1=Sys.time()
    for(luck in 1:length(stage1_cutoff)){
      set.seed(666)
      cutoff=stage1_cutoff[luck]
      Resampling=createFolds(TrainData$label,k=nfolds)
      final=list()
      
      pr=sum(TrainData$label==1)>sum(TrainData$label==0)
      for(xxx in 1:nfolds){
        #print('Step1: ReadData')
        trainData=TrainData[unlist(Resampling[-xxx]),]
        testData=TrainData[unlist(Resampling[xxx]),]
        geneNum_pathways=sapply(1:length(pathlistDB),function(i) length(pathlistDB[[i]]))
        pathlistDB_sub=pathlistDB[which(geneNum_pathways > PathwaySizeDown & geneNum_pathways < PathwaySizeUp )]
        #print(paste0('     |>Total number of pathways==>>',length(pathlistDB_sub)))
        
        
        #print('Step2: FeartureSelection-features')
        if(Stage1_FeartureSelection_Method=='cor'){
          #print(paste0('      Using <<  correlation  >>',' ,and you choose cutoff:',cutoff))
          Cor=stats::cor(trainData$label,trainData)
          Cor=ifelse(Cor>0,Cor,-Cor)
          names(Cor)=colnames(trainData)
          Cor_names=names(Cor)
          Cor_cutoff=Cor[which(Cor>cutoff)]
          Cor_cutoff_names=names(Cor_cutoff)
        }else{
          #print(paste0('      Using <<  wilcox.test  >>',' ,and you choose cutoff:',cutoff))
          train_0=trainData[which(trainData$label==0),]
          train_1=trainData[which(trainData$label==1),]
          Cor=unlist(mclapply(1:ncol(trainData),function(x) wilcox.test(train_0[,x],train_1[,x])$p.value,mc.cores=10))
          
          names(Cor)=colnames(trainData)
          Cor_names=names(Cor)
          Cor_cutoff=Cor[which(Cor<cutoff)]
          Cor_cutoff_names=names(Cor_cutoff)
        }
        
        featureAnno=FeatureAnno[FeatureAnno$ID %in% colnames(trainData),]
        MinfeatureNum_pathways2=MinfeatureNum_pathways+1
        featureNum_pathways=mclapply(1:length(pathlistDB_sub),function(x){
          id=c('label',featureAnno$ID[which(featureAnno$entrezID %in% pathlistDB_sub[[x]])])
          if(length(id)>MinfeatureNum_pathways2){
            id2=id[which(id %in% Cor_cutoff_names)]
            if(length(id2)<MinfeatureNum_pathways2){
              a=Cor[id]
              if(Stage1_FeartureSelection_Method=='cor'){
                id2=names(a)[order(a,decreasing = T)[1:MinfeatureNum_pathways2]]
              }else{
                id2=names(a)[order(a,decreasing = F)[1:MinfeatureNum_pathways2]]
              }
              return(id2)
            }else{
              return(id2)
            }
          }else{
            return(id)
          }
        } ,mc.cores=cores)
        lens=sapply(1:length(featureNum_pathways),function(x) length(featureNum_pathways[[x]]))
        
        #print('Step3: MergeData')
        trainDataList=mclapply(1:length(featureNum_pathways),function(x) trainData[,featureNum_pathways[[x]]] ,mc.cores=cores)
        testDataList=mclapply(1:length(featureNum_pathways),function(x) testData[,featureNum_pathways[[x]]] ,mc.cores=cores)
        names(trainDataList)=names(pathlistDB_sub)
        names(testDataList)=names(pathlistDB_sub)
        trainDataList=trainDataList[which(lens>MinfeatureNum_pathways)]
        testDataList=testDataList[which(lens>MinfeatureNum_pathways)]
        featureNum_pathways2=sapply(1:length(trainDataList),function(i2) length(trainDataList[[i2]]))
        #print(paste0('     |> Total number of selected pathways==>>',length(trainDataList)))
        #print(paste0('     |> Min features number of pathways==>>',min(featureNum_pathways2)-1,'.......','Max features number of pathways==>>',max(featureNum_pathways2)-1))
        
        
        #(FeartureSelection-pathways)
        #print('Step4: Reconstruction')
        if(Inner_CV){
          #print('     |> Using Inner CV ~ ~ ~')
          train=mclapply(1:length(trainDataList),function(i4) baseModel(trainData =trainDataList[[i4]],testData =NULL,predMode =predMode,classifier = classifier,inner_folds=inner_folds),mc.cores=cores)
          test=mclapply(1:length(testDataList),function(i5) baseModel(trainData =trainDataList[[i5]],testData =testDataList[[i5]],predMode =predMode,classifier = classifier),mc.cores=cores)
          #pred=mclapply(1:length(trainDataList),function(i4) baseModel2(trainData =trainDataList[[i4]],testData =testDataList[[i4]],classifier = stage1_learner,inner_folds=inner_folds),mc.cores=cores)
          #train=lapply(1:length(pred),function(x) pred[[x]]$predtrain)
          #test=lapply(1:length(pred),function(x) pred[[x]]$predtest)
        }else{
          
          train=mclapply(1:length(trainDataList),function(i4) baseModel(trainData =trainDataList[[i4]],testData =trainDataList[[i4]],predMode = predMode ,classifier = stage1_learner),mc.cores=cores)
          test=mclapply(1:length(testDataList),function(i5) baseModel(trainData =trainDataList[[i5]],testData =testDataList[[i5]],predMode = predMode ,classifier = stage1_learner),mc.cores=cores)
          
        }
        
        #print('     <<< Reconstruction Done! >>>     ')
        
        #(FeartureSelection-pathways)
        #print('Step5: FeartureSelection-pathways')
        
        
        dbmap=unique(unlist(pathlistDB_sub))
        annomap=unique(featureAnno$entrezID)
        mapgene=intersect(annomap,dbmap)
        map=featureAnno$ID[which(featureAnno$entrezID %in% mapgene)]
        Unmapped_Train=trainData[,setdiff(colnames(trainData),map)]
        Unmapped_Test=testData[,setdiff(colnames(trainData),map)]
        #print(paste0('Unmapped_num: ',ncol(Unmapped_Train)))
        Unmapped_0=Unmapped_Train[which(Unmapped_Train$label==unique(Unmapped_Train$label)[1]),]
        Unmapped_1=Unmapped_Train[which(Unmapped_Train$label==unique(Unmapped_Train$label)[2]),]
        if(Add_FeartureSelection_Method=='cor'){
          Unmapped_pvalue= abs(stats::cor(Unmapped_Train$label,Unmapped_Train[,-1]))
          Unmapped_pvalue=1/Unmapped_pvalue
        }else{
          Unmapped_0=Unmapped_Train[which(Unmapped_Train$label==unique(Unmapped_Train$label)[1]),]
          Unmapped_1=Unmapped_Train[which(Unmapped_Train$label==unique(Unmapped_Train$label)[2]),]
          Unmapped_pvalue=unlist(mclapply(2:ncol(Unmapped_Train),function(x) wilcox.test(Unmapped_0[,x],Unmapped_1[,x])$p.value,mc.cores=cores))
          
        }
        
        Record2=list()
        sink(nullfile())
        for(ii in  1:length(stage2_cutoff)){
          
          if(stage2_cutoff[ii]==0){
            index=Stage2_FeartureSelection(Stage2_FeartureSelection_Method='None',data=train,
                                           label=trainDataList[[1]]$label,cutoff=stage2_cutoff[ii],preMode='probability',classifier =classifier,cores=cores)
          }else{
            index=Stage2_FeartureSelection(Stage2_FeartureSelection_Method=Stage2_FeartureSelection_Method,data=train,
                                           label=trainDataList[[1]]$label,cutoff=stage2_cutoff[ii],preMode='probability',classifier =classifier,cores=cores)
            
          }
          
          if(is.null(classifier2)){
            stage2_learners=stage1_learner
          }else{
            stage2_learners=classifier2
          }
          
          n=length(Unmapped_num)*length(stage2_learners)
          Record=data.frame(Unmapped_num=1:n,stage2_learner=1:n,AUC=1:n,
                            PCC=1:n,BAC=1:n,cutoff=1:n)
          for(i in 1:length(Unmapped_num)){
            newtrain=do.call(cbind,train[index])
            colnames(newtrain)=names(trainDataList)[index]
            newtrain=cbind(label=trainDataList[[1]]$label,newtrain)
            
            newtest=do.call(cbind,test[index])
            colnames(newtest)=names(trainDataList)[index]
            newtest=cbind(label=testDataList[[1]]$label,newtest)
            colnames(newtest)=gsub(':','',colnames(newtest))
            colnames(newtrain)=gsub(':','',colnames(newtrain))
            if(Add_UnMapped==TRUE & Unmapped_num[i] > 0 ){
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
            }
            
            
            for(z in 1:length(stage2_learners)){
              classifier_2=stage2_learners[z]
              a=(i-1)*length(stage2_learners)+z
              result=baseModel(trainData=newtrain,testData=newtest,predMode ='probability',classifier = classifier_2)
              testDataY=testDataList[[1]]$label
              pre=ifelse(result>0.5,1,0)
              accuracy_class1 <- sum(pre[testDataY == 1] == 1) / sum(testDataY == 1)
              accuracy_class0 <- sum(pre[testDataY == 0] == 0) / sum(testDataY == 0)
              Record[a,5]=(accuracy_class1 + accuracy_class0) / 2
              Record[a,1]=Unmapped_num[i]
              Record[a,2]=classifier_2
              Record[a,4]=stats::cor(testDataY,result,method='pearson')
              Record[a,3]=ROCR::performance(ROCR::prediction(result,testDataY),'auc')@y.values[[1]]
              Record[a,6]=stage2_cutoff[ii]
            }
          }
          Record2[[ii]]=Record
          
          
          
        }
        sink()
        Record2=do.call(rbind,Record2) 
        final[[xxx]]=Record2
      }
      
      
      hope=do.call(rbind,final)
      
      hope=aggregate(hope[,c(3:5)],by=list(Unmapped_num=hope$Unmapped_num,stage2_cutoff=hope$cutoff,stage2_learner=hope$stage2_learner),mean)
      hope$stage1_cutoff=stage1_cutoff[luck]
      #if(verbose)print(hope)
      #print(paste0('stage1_cutoff : ',cutoff))
      #print(luck)
      
      HOPE[[luck]]=hope
    }
    HOPE=do.call(rbind,HOPE)
    HOPE$stage1_learner=stage1_learner
    HOPE=HOPE[,c('stage1_learner','stage2_learner','stage1_cutoff','stage2_cutoff','Unmapped_num','AUC','BAC','PCC')]
    if(verbose)print(HOPE)
    re[[c1]]=HOPE
    t2=Sys.time()
    if(verbose)print(t2-t1)
    if(verbose)print(' ')
  }
  re=do.call(rbind,re)
  re=re[order(re$AUC,decreasing = T),]
  return(re)
}



#Selection of stage-1 basemodels
classifier1=c('liblinear','svm')


#stage-1 feature_selection
stage1_cutoff=c(0.3,0.5)


#Number of unmapped features
Unmapped_num=c(5,10)

Step2_FeartureSelection_Method='wilcox.test'

#stage-2 feature_selection
stage2_cutoff=c(0.9,0.8)

#Selection of stage-2 basemodels(The default is the same as the stage-1)
classifier2=NULL

pathlistDB=readRDS('/home/zhangshunjie/BioMM/GO2ALLEGS_BP.rds')
FeatureAnno=readRDS('/home/zhangshunjie/BioMM/extra_3/m_anno.rds')
TrainData=readRDS('/home/zhangshunjie/BioMM/extra_3/Final_m_data.rds')[1:40,1:50000]

#A data frame contains hyperparameter results
result=HyBioM2(TrainData=TrainData,pathlistDB=pathlistDB,FeatureAnno=FeatureAnno,resampling=NULL,nfolds=2,classifier=classifier1,
           PathwaySizeUp=200,PathwaySizeDown=150,MinfeatureNum_pathways=10,
           Add_FeartureSelection_Method='wilcox.test',Unmapped_num=Unmapped_num,
           Inner_CV=F,inner_folds=10,
           Stage1_FeartureSelection_Method='cor',stage1_cutoff=stage1_cutoff,
           Stage2_FeartureSelection_Method='RemoveHighcor',stage2_cutoff=stage2_cutoff,
           classifier2=NULL,cores=20,verbose=TRUE)


#View the optimal hyperparameter combination
head(result[order(result$AUC,decreasing = T),c(1,3:6,8)])

#   stage1_learner stage1_cutoff stage2_cutoff Unmapped_num       AUC       PCC
#15            svm         0.001            10            5 0.7373737 0.4020263
#2       liblinear         0.010             5           10 0.7171717 0.3806561
#4       liblinear         0.010            10           10 0.7171717 0.3815476
#6       liblinear         0.001             5           10 0.7171717 0.3798767
#8       liblinear         0.001            10           10 0.7171717 0.3809311
#16            svm         0.001            10           10 0.7121212 0.4094221






