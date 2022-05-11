#################################################################################
#                                                                               #
#  A distribution-free model for longitudinal zero-inflated count data          #                                    #
#                                                                               #
#  input:                                                                       #
#    tab-delimited count data matrix that contains n features(f), k samples(s), #
#    m time points(t) under two conditions(c) with the following format         #
#-------------------------------------------------------------------------------#
#####                  c_1            |             c_2                     #####
#####          s_1     ...    s_k     |     s_1     ...     s_k             #####
#####      t_1 ... t_m ... t1 ... t_m | t_1 ... t_m ... t_1 ... t__m        #####
##### f_1                             |                                     #####
##### f_2                             |                                     #####
##### ...                             |                                     #####
##### f_n                             |                                     #####
#-------------------------------------------------------------------------------#
#                                                                               #
#  option for main/daf funtcion:                                                #
#  y: long vector, responses for each subject at each time point                #
#  x.zero: design matrix for zero part,                                         #
#          add column of 1's for intercept                                      #
#  x.count: design matrix for count part,                                       #
#           add column of 1's for intercept                                     #
#  theta.int: initial estimators,(zero part first and main part)                #
#  id: subject id                                                               #
#  time: time point                                                             #
#  corr: correlation structure                                                  #
#  eps: converge criterion                                                      #
#                                                                               #
#  output for main funtcion:                                                    #
#  conv: converge indicator, 1 for converged otherwise 0                        #
#  theta: estimates for zero model and main model                               #
#  iter: iterations used to achieve convergenece                                #
#  Sigma: sd of the estimators                                                  #
#  alpha: estimated correlation paramter                                        #
#################################################################################


# load required libraries
require(pscl)     # initial estimator
require(mvtnorm)
require(MASS)


################################################################
#################### self-defined function #####################
################################################################

inv.logit <- function(x.vec,beta) {
    t=x.vec%*%as.matrix(beta)
    c( exp(t)/( 1+exp(t) ) )
}

mu.pois <- function(x.vec,beta) {
    t=x.vec%*%as.matrix(beta)
    c(exp(t))
}

## h1
func1 <- function(b.u,b.v,x.u,x.v){
    mu.zero <- inv.logit(x.u,b.u); mu.count <- mu.pois(x.v,b.v)
    c(mu.zero + (1-mu.zero)*dpois(0,mu.count))
}

## h2
func2 <- function(b.u,b.v,x.u,x.v){
    mu.zero <- inv.logit(x.u,b.u); mu.count <- mu.pois(x.v,b.v)
    c( (1-mu.zero)*mu.count )
}


##
dh <- function(b.u,b.v,x.zero,x.count){ # input must be matrix
    l.u <- ncol(x.zero);l.v <- ncol(x.count); k <- nrow(x.zero)
    Di <- matrix(0,2*k,(l.u+l.v))
    for(i in 1:k){
        mu  <- c( x.zero[i,,drop=F]%*%as.matrix(b.u) )
        mv <- c( x.count[i,,drop=F]%*%as.matrix(b.v) )
        dh1u <-     x.zero[i,,drop=F]*exp(mu)*(1-exp(-exp(mv)))/(1+exp(mu))^2
        dh1v <-     x.count[i,,drop=F]*( -exp(-exp(mv)+mv)/(1+exp(mu)) )
        dh2u <-     x.zero[i,,drop=F]*( -exp(mu+mv)/(1+exp(mu))^2 )
        dh2v <-     x.count[i,,drop=F]*( exp(mv)/(1+exp(mu)) )
        Di[(2*i-1),] <- c(dh1u,dh1v)   # instead of cbind ...
        Di[(2*i),]   <- c(dh2u,dh2v)
    }
    return(Di=Di)
}

################################################################
######                   main function                   #######
######################## model fitting #########################
################################################################
zip.corrgee=function(y,x.zero,x.count,theta.int,id,time,corr,eps=1e-4) {
    
    x.zero <- as.matrix(x.zero); x.count <- as.matrix(x.count)
    l.u <- ncol(x.zero);l.v <- ncol(x.count)
    n <- length(unique(id)); k <- length(unique(time))
    f1 <- as.numeric(y==0)
    f2 <- y*(1-f1)
    
    theta <- as.vector(theta.int)
    iter <- 0
    loop <- 1
    while(loop==1){
        
        b.u <- theta[1:l.u]
        b.v <- theta[(l.u+1):(l.u+l.v)]
        mu <- c(x.zero%*%as.matrix(b.u))
        mv <- c(x.count%*%as.matrix(b.v))
        h1 <- func1(b.u,b.v,x.zero,x.count)
        h2 <- func2(b.u,b.v,x.zero,x.count)
        if(h1[1]=='NaN' | h2[1]=='NaN') {
            theta=rep(NaN,length(as.vector(theta.int)))
            break}
        
        b.u <- theta[1:l.u]
        b.v <- theta[(l.u+1):(l.u+l.v)]
        mu <- c(x.zero%*%as.matrix(b.u))
        mv <- c(x.count%*%as.matrix(b.v))
        h1 <- func1(b.u,b.v,x.zero,x.count)
        h2 <- func2(b.u,b.v,x.zero,x.count)
        if(h1[1]=='NaN' | h2[1]=='NaN') {
            theta=rep(NaN,length(as.vector(theta.int)))
            break}
        
        var.f1 <- h1*(1-h1)
        var.f2 <- exp(mv)*( 1 + exp(mu)*exp(mv)/(1+exp(mu)) )/(1+exp(mu))
        A <- cbind(var.f1,var.f2)
        #defin S and U
        
        # U and Wr for sigma
        S <- cbind(f1-h1,f2-h2)
        U <- Wr <- 0
        H<- cbind(h1, h2)
        
        ############################################################
        ## alpha estimator:
        if (corr=='AR1') {
            tt<- 1:k
            alpha<-prod<- 0
            for(i in 1:n){
                Ai <- diag(c(t(A[id==i,])))
                Hi <- H[id==i,]
                temp <- rep(0,2*k-1)
                temp[2*(1:k)-1]<- (-h1*h2)[id==i]
                if (k==1) {Ai[1,2] <- temp; Ai[2,1] <- temp
                } else { diag( Ai[-2*k,-1] ) <- temp;diag( Ai[-1,-2*k] ) <- temp }
                Si<-S[id==i,]
                for(j in 1:(k-1)){
                    Si1 <- c(t(Si[tt==j,]))
                    Hi1 <- c(t(Hi[tt==j,]))
                    if(Hi1[2]==0) {Hi1[2] = 1e-100}
                    ei1<- Si1/sqrt(Hi1)
                    l<-j+1
                    Si2 <- c(t(Si[tt==l,]))
                    Hi2 <- c(t(Hi[tt==l,]))
                    if(Hi2[2]==0) {Hi2[2] = 1e-100}
                    ei2<- Si2/sqrt(Hi2)
                    prod<- prod+ei1*ei2
                }
            }
            
            #calculate K
            KK<-0
            for(i in 1:n){
                KK=KK+(k-1)
            }
            
            p<-dim(x.count)[2]
            temp_a<-abs(as.vector((prod[1]+prod[2])/(2*(KK-p))))
            if (temp_a=='NaN' | temp_a>100) {break}
            else if(temp_a>1) {alpha = 1}
            else {alpha=temp_a}
            
            ############################################################
            # define working correltaion matrix R
            times <- 1:k
            M1 <- abs(outer(times, times, "-"))
            M2 <- alpha^M1
            J <- matrix(1,2,2)
            M3 <-outer(J,M2, "*")
            R<-matrix(0,2*k,2*k)
            for(i in 1:k){
                for(j in 1:k){
                    R[(i*2-1):(2*i),(j*2-1):(2*j)]=M3[,,i,j]
                }
            }
            for(i in 1:(2*k)){
                if((i/2)%%1!=0){
                    R[i,(i+1)]=0
                }
                else{R[i,(i-1)]=0}
            }
        }
        else if (corr=='exchangeable'){
            tt<- 1:k
            alpha<-prod<- 0
            for(i in 1:n){
                Ai <- diag(c(t(A[id==i,])))
                Hi <- H[id==i,]
                temp <- rep(0,2*k-1)
                temp[2*(1:k)-1]<- (-h1*h2)[id==i]
                if (k==1) {Ai[1,2] <- temp; Ai[2,1] <- temp
                } else { diag( Ai[-2*k,-1] ) <- temp;diag( Ai[-1,-2*k] ) <- temp }
                Si<-S[id==i,]
                for(j in 1:(k-1)){
                    Si1 <- c(t(Si[tt==j,]))
                    Hi1 <- c(t(Hi[tt==j,]))
                    if(Hi1[2]==0) {Hi1[2] = 1e-100}
                    ei1<- Si1/sqrt(Hi1)
                    for(l in 2:k){
                        Si2 <- c(t(Si[tt==l,]))
                        Hi2 <- c(t(Hi[tt==l,]))
                        if(Hi2[2]==0) {Hi2[2] = 1e-100}
                        ei2<- Si2/sqrt(Hi2)
                        prod<- prod+ei1*ei2
                    }
                }
            }
            
            #calculate K
            KK<-0
            for(i in 1:n){
                KK=KK+k*(k-1)
            }
            
            p<-dim(x.count)[2]
            temp_a<-abs(as.vector((prod[1]+prod[2])/(KK-p)^2))
            if (temp_a=='NaN' | temp_a>100) {break}
            else if(temp_a>1) {alpha = 1}
            else {alpha=temp_a}
            
            ############################################################
            # define working correltaion matrix R
            times <- 1:k
            M1 <- matrix(1,k,k); diag(M1)<- 0
            M2 <- alpha^M1
            J <- matrix(1,2,2)
            M3 <-outer(J,M2, "*")
            R<-matrix(0,2*k,2*k)
            for(i in 1:k){
                for(j in 1:k){
                    R[(i*2-1):(2*i),(j*2-1):(2*j)]=M3[,,i,j]
                }
            }
            for(i in 1:(2*k)){
                if((i/2)%%1!=0){
                    R[i,(i+1)]=0
                }
                else{R[i,(i-1)]=0}
            }
        }
        
        ############################################################
        
        U<- Wr<- B <- sand<- 0 ## sw estimator
        
        for(i in 1:n){
            Ai <- diag(c(t(A[id==i,])))
            temp <- rep(0,2*k-1)
            temp[2*(1:k)-1]<- (-h1*h2)[id==i]
            if (k==1) {Ai[1,2] <- temp; Ai[2,1] <- temp
            } else { diag( Ai[-2*k,-1] ) <- temp;diag( Ai[-1,-2*k] ) <- temp }
            
            if(Ai[1,1]=='NaN') { break }
            else{
                e <- eigen(Ai) #A must be symmetirc positive definite
            }
            
            tempV <- e$vectors
            # squart root of Ai
            SA <- tempV %*% diag(sqrt(round(e$values))) %*% t(tempV)
            
            Vi<- SA%*%R%*%SA
            Di <- dh(b.u,b.v,x.zero[id==i,,drop=F],x.count[id==i,,drop=F])
            Si <- c(t(S[id==i,]))
            # theta estimator:
            out <- tryCatch(solve(Vi) %*% Vi, error = function(e) e)
            if(any(class(out) == "error")){ # if Vi is singular, use generlized inverse instead
                U <- U + t(Di)%*%ginv(Vi)%*%Si
                Wr <- Wr + t(Di)%*%ginv(Vi)%*%Di
                
                ## sw estimator:
                B <- B + t(Di)%*%ginv(Vi)%*%Di/n
                sand <- sand + t(Di)%*%ginv(Vi)%*%Si%*%t(Si)%*%ginv(Vi)%*%Di/n
            }
            else{
                U <- U + t(Di)%*%solve(Vi)%*%Si
                Wr <- Wr + t(Di)%*%solve(Vi)%*%Di
                
                ## sw estimator:
                B <- B + t(Di)%*%solve(Vi)%*%Di/n
                sand <- sand + t(Di)%*%solve(Vi)%*%Si%*%t(Si)%*%solve(Vi)%*%Di/n
            }
        }
        
        out2 <- tryCatch(solve(Wr) %*% Wr, error = function(e) e)
        if(any(class(out2) == "error")){ # if Wr is singular, use generlized inverse instead
            theta.new <- theta + ginv(Wr)%*% U
        }
        else{
            theta.new <- theta + solve(Wr)%*% U
        }
        
        if(max(abs(theta.new-theta))<eps){
            loop <- 0}
        
        theta <- theta.new
        
        sigma <- ginv(B)%*%(sand)%*%t(ginv(B))/n
        Sigma=sqrt(diag(sigma))
        
        iter <- iter+1
        if(iter>50){
            cat("Not Converge!\n")
            loop <- 0;converge=0
        } else {converge=1}
    }
    
    return(list(conv=converge,theta=c(theta),iter=iter,alpha=alpha,
    Sigma=sqrt(diag(sigma))))
}
