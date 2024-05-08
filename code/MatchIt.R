library(MatchIt)

meta=readrds('your_meta.rds')
#Meta include 'Label' and other covariates of interest
head(meta[,1:5])
#Label Age Female       Cd8t       Cd4t
#201172200060_R08C01     1  26      0 0.06150491 0.19702276
#201172200055_R06C01     1  43      1 0.08574886 0.02965455
#201810150072_R03C01     1  55      0 0.16006435 0.16938068
#201800470034_R07C01     0  27      1 0.06554267 0.14906774
#201172200012_R08C01     0  62      1 0.07471174 0.09312758
#201800470020_R08C01     0  40      0 0.08413080 0.11795726

#Sample imbalance
table(meta$Label)

#0   1 
#68 361


#Set minority class samples to 'TRUE'
meta$Label=as.logical(meta$Label == 0)


#Paired samples based on all covariates(setting method and pairing ratio)
meta_match <- matchit(Label~., data =meta, method="nearest", ratio=1)
summary(meta_match)

#Summary of Balance for All Data:
#  Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#distance            0.2994        0.1320          0.8504     2.8508    0.2727
#Age                47.0588       39.5623          0.5322     1.3622    0.1535
#Female              0.6618        0.6122          0.1048          .    0.0496
#Cd8t                0.0819        0.0740          0.1888     0.9874    0.0656
#Cd4t                0.1595        0.1623         -0.0512     0.6441    0.0541
#Nk                  0.0423        0.0385          0.1069     0.7498    0.0576
#...
#population_6       -0.1631       -0.0255         -0.0372     1.0245    0.0327
#population_7        0.0889        0.0064          0.0289     0.8067    0.0279
#population_8       -0.3241        0.0612         -0.0784     3.7793    0.0344
#population_9        0.3923       -0.0258          0.2495     0.7135    0.0404
#population_10        0.9117       -0.1822          0.2190     5.7112    0.0476
#smoking             1.0109        2.5165         -0.3554     0.7408    0.0910



#Summary of Balance for Matched Data:
# Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
#distance            0.2994        0.2642          0.1789     1.6232    0.0141
#Age                47.0588       44.2206          0.2015     1.5343    0.0817
#Female              0.6618        0.7059         -0.0933          .    0.0441
#Cd8t                0.0819        0.0929         -0.2656     0.7571    0.0579
#Cd4t                0.1595        0.1592          0.0056     0.4130    0.0748
#Nk                  0.0423        0.0444         -0.0577     0.7739    0.0451
#...
#population_6       -0.1631       -0.1822          0.0052     1.0568    0.0278
#population_7        0.0889        0.7119         -0.2185     0.7062    0.0484
#population_8       -0.3241        0.5084         -0.1694     4.9049    0.0457
#population_9        0.3923        0.7959         -0.2409     0.8792    0.0734
#population_10        0.9117        0.3284          0.1168     6.4489    0.0413
#smoking             1.0109        0.8027          0.0492     1.0853    0.0275




#Sample Sizes:
#  Control Treated
#All           361      68
#Matched        68      68
#Unmatched     293       0
#Discarded       0       0

#Extraction of paired samples
match <- match.data(meta_match)
