library(GO.db)
library(AnnotationDbi)

max_ancestor=names(which.max(lengths(as.list(GOBPOFFSPRING))))
second_ancestor=as.character(as.list(GOBPCHILDREN)[[max_ancestor]])
third_ancestor=as.character(unique(unlist(as.list(GOBPCHILDREN)[second_ancestor])))
xx=as.list(GOBPANCESTOR)

ancestor_v1=sapply(1:length(xx),function(x){
  a1=second_ancestor[which(second_ancestor %in% unlist(xx[[x]]))]
  names(a1)=rep(names(xx)[x],length(a1))
  a1
})
ancestor_v1=unlist(ancestor_v1)
df=data.frame(GO_id=1:length(ancestor_v1),GO_Ancestor=1:length(ancestor_v1))
df$GO_id=names(ancestor_v1)
df$GO_Ancestor=ancestor_v1
df$GO_id_term=AnnotationDbi::select(GO.db, keys=df$GO_id, columns=c("GOID", "TERM"), keytype="GOID")[,2]
df$GO_Ancestor_term=AnnotationDbi::select(GO.db, keys=df$GO_Ancestor, columns=c("GOID", "TERM"), keytype="GOID")[,2]
df=df[,c(1,3,2,4)]
write.csv(df,'GOBP_ancestor_v1.csv')


ancestor_v2=sapply(1:length(xx),function(x){
  a1=third_ancestor[which(third_ancestor %in% unlist(xx[[x]]))]
  names(a1)=rep(names(xx)[x],length(a1))
  a1
})
ancestor_v2=unlist(ancestor_v2)
df=data.frame(GO_id=1:length(ancestor_v2),GO_Ancestor=1:length(ancestor_v2))
df$GO_id=names(ancestor_v2)
df$GO_Ancestor=ancestor_v2
df$GO_id_term=AnnotationDbi::select(GO.db, keys=df$GO_id, columns=c("GOID", "TERM"), keytype="GOID")[,2]
df$GO_Ancestor_term=AnnotationDbi::select(GO.db, keys=df$GO_Ancestor, columns=c("GOID", "TERM"), keytype="GOID")[,2]
df=df[,c(1,3,2,4)]
write.csv(df,'GOBP_ancestor_v2.csv')