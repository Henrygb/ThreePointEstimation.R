######################################################
#
# PERT three point estimation distribution
# Henry Bottomley - September 2016
#
# Given 
#    L low case estimate
#    M most likely estimate
#    H high case estimate
#
# use a location-scale Beta distribution 
# with a mean of (L + w*M + H) / (w + 2)
# and a mode of M 
# typically with w = 4 
#
# standard deviation usually not exactly range/6 but instead divide 
# range by something between sqrt(28) and sqrt(50.4) or 
# more generally between sqrt(12+4*w) and sqrt((w+2)^2*(w+3)/(w+1))
#
# See:
# https://en.wikipedia.org/wiki/Three-point_estimation
# https://en.wikipedia.org/wiki/Program_evaluation_and_review_technique
#
########################################################

tpeDetails <- function(L, M, H, w=4){
   location <- L
   scale <- H-L+M*0
   if (min(scale) <= 0){ warning("High not greater than Low", call. = FALSE) }
   if (sum(L>M)>0){ warning("Most likely less than Low", call. = FALSE) }
   if (sum(H<M)>0){ warning("Most likely greater than High", call. = FALSE) }
   betamode <- ifelse( scale==0, 1/2, (M-location) / scale )
   alpha <- 1 + w * betamode
   beta  <- 1 + w * (1-betamode)
   details <- data.frame(location,scale,betamode,alpha,beta)
   details
   }

tpeSummary <- function(L, M, H, w=4){
   details <- tpeDetails(L=L, M=M, H=H, w=w)
   summtpe <- details$location + details$scale * 
              cbind(0,qbeta(0.25,details$alpha,details$beta),
                   qbeta(0.5,details$alpha,details$beta),
                         (w*details$betamode+1)/(w+2),
                         qbeta(0.75,details$alpha,details$beta), 1)
   colnames(summtpe) <- c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
   summtpe
   }

tpeStatistics <- function(L, M, H, w=4){
   details <- tpeDetails(L=L, M=M, H=H, w=w)
   var <- details$scale^2 * details$alpha * details$beta / 
         (details$alpha + details$beta)^2 / 
         (details$alpha + details$beta + 1) 
   statistics <- cbind(L, M, details$location + details$scale * 
                             qbeta(0.5,details$alpha,details$beta),
                       (L+w*M+H)/(w+2), H, var, sqrt(var), details$scale)
   colnames(statistics) <- c("Min.","Mode","Median","Mean","Max.",
                             "Var.","St.Dev.","Range")
   statistics
   }

dtpe <- function(x, L, M, H, w=4, log=FALSE){
    details <- tpeDetails(L=L, M=M, H=H, w=w)
    if(details$scale==0){return(ifelse(x == details$location, NaN, 0))} 
    dbeta(x = (x-details$location) / details$scale, 
            shape1 = details$alpha, shape2 = details$beta, 
            ncp=0, log=log) / details$scale 
    } 

ptpe <- function(q, L, M, H, w=4, lower.tail=TRUE, log.p=FALSE){
    details <- tpeDetails(L=L, M=M, H=H, w=w)
    if(details$scale==0){return(ifelse(x < details$location, 0, 1))} 
    pbeta(q = (q-details$location) / details$scale, 
            shape1 = details$alpha, shape2 = details$beta, 
            ncp=0, lower.tail=lower.tail, log.p=log.p) 
    } 

qtpe <- function(p, L, M, H, w=4, lower.tail=TRUE, log.p=FALSE){
    details <- tpeDetails(L=L, M=M, H=H, w=w)
    details$location + details$scale * qbeta(p=p,
          shape1 = details$alpha, shape2 = details$beta, 
          ncp=0, lower.tail=lower.tail, log.p=log.p)
    }

rtpe <- function(n, L, M, H, w=4){
    details <- tpeDetails(L=L, M=M, H=H, w=w)
    details$location + details$scale * rbeta(n=n,
          shape1 = details$alpha, shape2 = details$beta,
          ncp=0)
    }


#######################################################
# Examples 
#
# see if sample statistics close to predicted
# top row predicted, second row simulated
set.seed(1) # only if want same results every time
dat <- rtpe(1000000, L=10, M=15, H=30) 
datamode <- function(dat){den <- density(dat, bw="nrd"); 
    den$x[min(which(den$y == max(den$y)))]
    } 
rbind( tpeStatistics(L=10, M=15, H=30),
       c( min(dat), datamode(dat), median(dat), mean(dat), max(dat),
          var(dat), sd(dat), max(dat)-min(dat) ) )
rbind( tpeSummary(L=10, M=15, H=30), 
       summary(dat) )
#
# reverse ptpe and qtpe
max(abs(ptpe(qtpe((0:100)/100, L=10, M=15, H=30), L=10, M=15, H=30) - 
    (0:100)/100 ) )
max(abs(qtpe(ptpe(10:30, L=10, M=15, H=30), L=10, M=15, H=30) - 
    (10:30) ) ) # 29 is largest
#
# draw density, CDF and inverse CDF
fund <- function(x){dtpe(x, L=10, M=15, H=30)}
plot(fund, xlim=c(0,40))
funp <- function(x){ptpe(x, L=10, M=15, H=30)}
plot(funp, xlim=c(0,40))
funq <- function(x){qtpe(x, L=10, M=15, H=30)}
plot(funq, xlim=c(0,1))



