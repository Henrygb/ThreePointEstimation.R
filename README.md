# ThreePointEstimation.R
Provides R probability distribution code aiming at PERT three point estimation

PERT three point estimation distribution
Henry Bottomley - initial attempt September 2016

Given 
   L low case estimate
   M most likely estimate
   H high case estimate

use a location-scale Beta distribution 
with a mean of (L + w*M + H) / (w + 2)
and a mode of M 
typically with w = 4 

Functions are 

    dtpe(x, L, M, H, w=4, log=FALSE)
    ptpe(q, L, M, H, w=4, lower.tail=TRUE, log.p=FALSE)
	qtpe(p, L, M, H, w=4, lower.tail=TRUE, log.p=FALSE)
	rtpe(n, L, M, H, w=4)
	
and also 

    tpeDetails(L, M, H, w=4)
    tpeSummary(L, M, H, w=4)	
	tpeStatistics <- function(L, M, H, w=4)
	
standard deviation usually not exactly range/6 but instead divide 
range by something between sqrt(28) and sqrt(50.4) or 
more generally between sqrt(12+4*w) and sqrt((w+2)^2*(w+3)/(w+1))

See:
https://en.wikipedia.org/wiki/Three-point_estimation
https://en.wikipedia.org/wiki/Program_evaluation_and_review_technique

