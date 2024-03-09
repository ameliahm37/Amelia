#Factor Analysis - Amelia Mills Project

#Calculating the correlation matrix for EFA, section 4.2
x <- cor(divorce_data)

install.packages('ggplot2')
install.packages('psych')
library(ggplot2)
library(psych)

#Visual representation of the correlation matrix, section 4.2
heatmap(x, scale='none')

#Bartlett's Test of Sphericity, section 4.2
cortest.bartlett(x, 54)

#Scree plot to calculate the number of factors, section 4.3
scree(x[,1:54], factors=TRUE, pc=FALSE)

#Cumulative variance explained by the factors, section 4.3
results <- prcomp(divorce_data, scale = TRUE)
var_explained = results$sdev^2 / sum(results$sdev^2)
print(var_explained)
sum(var_explained[1:3])

#PAF factor loadings with no rotation, section 4.5
fa1.1 <- fa(r = x, nfactors = 3, rotate = 'none', fm = 'pa') 
summary(fa1.1) 
fa1.1$loadings

#PAF factor loadings with varimax rotation, section 4.5
fa2.1 <- fa(r = x, nfactors = 3, rotate = 'varimax', fm = 'pa') 
summary(fa2.1) 
fa2.1$loadings

#PAF factor loadings with promax rotation, section 4.5
fa3.1 <- fa(r = x, nfactors = 3, rotate = 'promax', fm = 'pa') 
summary(fa3.1) 
fa3.1$loadings

#Calculating the communalities for CFA, section 5.2
fa3.1$communality

#MLE factor loadings with no rotation, section 4.8
fa1.2 <- fa(r = x, nfactors = 3, rotate = 'none', fm = 'ml') 
summary(fa1.2) 
fa1.2$loadings

#MLE factor loadings with varimax rotation, section 4.8
fa2.2 <- fa(r = x, nfactors = 3, rotate = 'varimax', fm = 'ml') 
summary(fa2.2) 
fa2.2$loadings

#MLE factor loadings with promax rotation, section 4.8
fa3.2 <- fa(r = x, nfactors = 3, rotate = 'promax', fm = 'ml') 
summary(fa3.2) 
fa3.2$loadings

if (!requireNamespace("psych", quietly = TRUE)) {
  install.packages("psych")
}
if (!requireNamespace("knitr", quietly = TRUE)) {
  install.packages("knitr")
}

library(psych)
library(knitr)

#Extracting factor loadings
loadings_table1.1 <- as.data.frame(fa1.1$loadings[,1:3])
loadings_table2.1 <- as.data.frame(fa2.1$loadings[,1:3])
loadings_table3.1 <- as.data.frame(fa3.1$loadings[,1:3])
loadings_table1.2 <- as.data.frame(fa1.2$loadings[,1:3])
loadings_table2.2 <- as.data.frame(fa2.2$loadings[,1:3])
loadings_table3.2 <- as.data.frame(fa3.2$loadings[,1:3])

#Print factor loadings table in LaTeX format using kable
kable(loadings_table1.1, format = "latex", digits=c(4,4,4))
kable(loadings_table2.1, format = "latex", digits=c(4,4,4))
kable(loadings_table3.1, format = "latex", digits=c(4,4,4))
kable(loadings_table1.2, format = "latex", digits=c(4,4,4))
kable(loadings_table2.2, format = "latex", digits=c(4,4,4))
kable(loadings_table3.2, format = "latex", digits=c(4,4,4))

#CFA Procedure, Chapter 5
install.packages("lavaan")
library(lavaan)

#Defining the measurement model, section 5.1.1
cfa_model <- '
  #Define latent factors
  f1 =~ V1 + V2 + V3 + V4 + V5 + V7 + V8 + V9 + V10 + V11 + V12
        + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 
        + V22 + V24 + V25 + V26 + V27 + V28 + V29 + V30
  f2 =~ V42 + V43 + V45 + V46 + V47
  f3 =~ V31 + V32 + V33 + V34 + V35 + V36 + V37 + V38 + V39
        + V40 + V41 + V48 + V50 + V51 + V52 + V53 + V54
'
cfa_result <- cfa(cfa_model, data = divorce_data)
summary(cfa_result, fit.measures = TRUE)

#Calculating the modification indices (MIs), section 5.2
mi_results <- modificationindices(cfa_result)
ordered_mi_results <- mi_results[order(mi_results$mi, 
                                       decreasing = TRUE), ]
print(ordered_mi_results)

#Defining the structural model, section 5.2.1
cfa_model2 <- '
  # Define latent factors
  f1 =~ V1 + V2 + V3 + V4 + V5 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22 + V24 + V25 + V26 + V27 + V28 + V29 + V30
  f2 =~ V42 + V43 + V45 + V46 + V47
  f3 =~ V31 + V32 + V33 + V34 + V35 + V36 + V37 + V38 + V39 + V40 + V41 + V48 + V50 + V51 + V52 + V53 + V54
  V22 ~~ V28
  V50 ~~ V51
  V10 ~~ V16
  V17 ~~ V19
  V9 ~~ V15
  V52 ~~ V53
  V25 ~~ V29
  V24 ~~ V26
  V33 ~~ V40
  V18 ~~ V20
  V12 ~~ V18
  V18 ~~ V29
  V21 ~~ V27
  V4 ~~ V10
  V12 ~~ V20
  V14 ~~ V20
  V35 ~~ V36
  V5 ~~ V17
  V38 ~~ V40
  V4 ~~ V16
  V21 ~~ V29
  V11 ~~ V19
  V17 ~~ V30
  V52 ~~ V54
  V14 ~~ V17
  V11 ~~ V26
  V8 ~~ V18
  V20 ~~ V30
  V36 ~~ V39
  V21 ~~ V22
  f1 =~ V34
  V5 ~~ V26
  V5 ~~ V13
  V21 ~~ V40
  V14 ~~ V18
  V5 ~~ V19
  f3 =~ V42
  V26 ~~ V30
  V18 ~~ V21
  V5 ~~ V11
  V11 ~~ V17
  V39 ~~ V41
  V13 ~~ V17
  V40 ~~ V54
  V12 ~~ V29
  V5 ~~ V30
  f1 =~ V42
  f2 =~ V29
  V26 ~~ V27
  V19 ~~ V29
  V22 ~~ V40
  V36 ~~ V38
  V43 ~~ V45
  V18 ~~ V25
  V39 ~~ V40
  f3 =~ V18
  V48 ~~ V51
  V20 ~~ V29
  V29 ~~ V40 
  V1 ~~  V3
  V14 ~~ V26
  V36 ~~ V37
  V32 ~~ V35
  V13 ~~ V26
  V33 ~~ V54
  V11 ~~ V13
  V35 ~~ V40
  V13 ~~ V30
  V9 ~~ V26
  V16 ~~ V22
  V41 ~~ V51
  V2 ~~ V14
  V14 ~~ V30
  V10 ~~ V29
  V16 ~~ V29
  V2 ~~ V29
  V31 ~~ V37
  V8 ~~ V19
  V3 ~~  V4
  V1 ~~ V19
  V18 ~~ V24
  V5 ~~ V24
  V17 ~~ V26
  V18 ~~ V22
  V10 ~~ V22
  V20 ~~ V26
  V18 ~~ V26
  V16 ~~ V25
  f3 =~ V29
  V5 ~~ V29
  V29 ~~ V35
  V12 ~~ V54
  V17 ~~ V28
  V18 ~~ V40
  V11 ~~ V24
  V21 ~~ V35
  V17 ~~ V24
  V5 ~~ V1
  V2 ~~ V24
  V12 ~~ V14
  V8 ~~ V25
  V46 ~~ V36
  V20 ~~ V25
  V16 ~~ V28
  V25 ~~ V46
  V37 ~~ V39
  f2 =~ V37
  V33 ~~ V38
  V1 ~~ V11
  V4 ~~ V22
  V11 ~~ V38
  V11 ~~ V22
  V47 ~~ V40
  V1 ~~ V17
  V19 ~~ V26
  V24 ~~ V30
  V21 ~~ V33
  V53 ~~ V54
  V46 ~~ V39
  V13 ~~ V24
  V14 ~~ V24
  f3 =~  V5
  V11 ~~ V20
  V25 ~~ V27
  V10 ~~ V28
  V9 ~~ V30
  f1 =~ V39
  V35 ~~ V51
  V2 ~~ V20
  V19 ~~ V30
  V13 ~~ V18
  V19 ~~ V24
  V29 ~~ V33
  V42 ~~ V40 
  V5 ~~ V18
  V15 ~~ V20
  V27 ~~ V29
  V4 ~~ V28
  V2 ~~  V3
  V1 ~~  V8
  V8 ~~ V30
  V46 ~~ V52
  V8 ~~ V24
  V36 ~~ V41
  V1 ~~ V30
  V40 ~~ V41
  V38 ~~ V52
  V17 ~~ V29
  f3 =~ V12
  V22 ~~ V39
  f1 =~ V40
  V26 ~~ V41
  V47 ~~ V53
  V1 ~~ V12
  V2 ~~ V48
  V25 ~~ V26
  V5 ~~ V28
  f3 =~ V43
  V8 ~~ V39
  V34 ~~ V50
  V5 ~~ V14
  V42 ~~ V38
  V14 ~~ V29
  V8 ~~ V20
  V3 ~~ V48
  V13 ~~ V19
  f3 =~ V45
  V32 ~~ V34
  V20 ~~ V24
  V47 ~~ V38
  V3 ~~ V26
  V15 ~~ V30
  V46 ~~ V41
  V2 ~~ V18
  V9 ~~ V28
  # Add any additional specifications (e.g., factor variances, covariances)
'

cfa_result2 <- cfa(cfa_model2, data = divorce_data)
summary(cfa_result2, fit.measures = TRUE)

#Fitting the SEMs
fit <- sem(cfa_model, data = divorce_data)
fit2 <- sem(cfa_model2, data = divorce_data)

install.packages("semPlot")
library(semPlot)

#Plot the models, sections 5.1.1 and 5.2.1
semPaths(fit, layout = "tree", whatLabels = "std", 
         edge.label.cex = 0.75)
semPaths(fit2, layout = "tree", whatLabels = "std",
         edge.label.cex = 0.75)
