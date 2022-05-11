This is a brief introduction about how to use the proposed method for detection of significantly differentially abundant features for longitudinal metagenomic sequencing data between different conditions. 1 Data input format for comparison of two different metagenomic conditions:Suppose the data set contains K features, N samples, M time points under two conditions. The elements in a count matrix Y(k, i, j, h),  corresponding to the number of reads (or relative abundance) of feature k at time j in sample i for condition h, where h=1, 2. Data set is tab-delimited format and name of feature should be clarified at the first column. The example format is listed below:##############################################################################                 |      condition1          |             condition2    ##
## -----------------------------------------------------------------------####                 |  time 1  É     time M    |    time 1  É    time M    ##
## -----------------------------------------------------------------------#### feature1.sample1    c1111   É    c11M1     |    c1112   É    c11M2     #### feature1.sample2    c1211   É    c12M1     |    c1212   É    c12M2     ####     .                        .             |            .		  ##
##     .                        .             |            .		  ##
##     .                        .             |            .		  ##
##                 	                                        	  #### feature1.sampleN    c1N11   É    c1NM1     |    c1N12   É    c2NM2     ##
## feature2.sample1    c2111   É    c21M1     |    c2112   É    c21M2     #### feature2.sample2    c2211   É    c22M1     |    c2212   É    c22M2     ####     .                        .             |            .		  ##
##     .                        .             |            .		  ##
##     .                        .             |            .		  ##
##                 	                                        	  #### feature2.sampleN    c2N11   É    c2NM1     |    c2N12   É    c2NM2     ##
##     .                        .             |            .		  ##
##     .                        .             |            .		  ##
##     .                        .             |            .		  ##
############################################################################An example of data set sample for 1000 features count matrix is included on the website: example.csv.2. R commandsOpen R:2.1 Input the source file CorrZIDF.Rsource("CorrZIDF.R")2.2 Load a feature count matrix # 1000 features, 2 conditions, 25 samples across 10 sampling/time pointsnewdat<- read.csv("example.csv", row.names=1, header=TRUE)newdat<- as.matrix(newdat)N<- 25 # number of samples; 
grp<- 2 # numbers of groups/conditions; 
M<- 10 # number of sampling/time points
2.3 Set initial parameterscond <- c(rep(1, N*M),rep(0, N*M))  # condition, 1-disease, 0-control
time <- rep(1:M, N*grp)
id <- rep(1:(N*grp),each=M)
x.count <- cbind(rep(1, N*M*grp), cond)
x.zero <- cbind(rep(1, N*M*grp))
2.4 Analyze the example data

K=dim(newdat)[1]/N  # number of features
pvalue=rep(1, K)

# run the function over each featurefor(k in 1:K){
  print(k)
  # transfer count to long format
  F1 <- rbind(as.matrix(newdat[(N*k-(N-1)):(N*k),1:M]),
            as.matrix(newdat[(N*k-(N-1)):(N*k),(M+1):(2*M)])) 
  rownames(F1) <- colnames(F1) <- NULL
  F1 <- t(F1) 
  y <- c(F1)

  ## get initial est using obs from the 10th time point
  y1 <- y[time==M]; x1.count <- cond[time==M] 

  ## cond for main model and intercept only for zero model
  m1 <- zeroinfl(y1 ~ x1.count | 1) 
  theta.int <- c(m1$coef$zero,m1$coef$count)

  # estimated data under AR1 correlation 
  fit <- zip.corrgee(y,x.zero,x.count,theta.int,id,time, corr= 'AR1')
  pvalue[k]=1-pchisq((fit$theta[3]/fit$Sigma[3])^2,df=1)
}

plot(1:K, pvalue)  ## plot the pvalues

