This is a brief introduction about how to use the proposed method for detection of significantly differentially abundant features for longitudinal metagenomic sequencing data between different conditions. 
## -----------------------------------------------------------------------##
## -----------------------------------------------------------------------##
##     .                        .             |            .		  ##
##     .                        .             |            .		  ##
##                 	                                        	  ##
## feature2.sample1    c2111   �    c21M1     |    c2112   �    c21M2     ##
##     .                        .             |            .		  ##
##     .                        .             |            .		  ##
##                 	                                        	  ##
##     .                        .             |            .		  ##
##     .                        .             |            .		  ##
##     .                        .             |            .		  ##
############################################################################
grp<- 2 # numbers of groups/conditions; 
M<- 10 # number of sampling/time points

time <- rep(1:M, N*grp)
id <- rep(1:(N*grp),each=M)
x.count <- cbind(rep(1, N*M*grp), cond)
x.zero <- cbind(rep(1, N*M*grp))


K=dim(newdat)[1]/N  # number of features
pvalue=rep(1, K)


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

