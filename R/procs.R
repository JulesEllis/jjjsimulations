CBP <- function(pvalues, lambda) {
  pCBP <- p.adjust(pvalues[pvalues<=lambda]/lambda, method = "bonferroni" )
  return(c(length(pCBP), lambda*0.05/length(pCBP),pCBP))
}

Study_FWER <- function(model, sizes, nsim, nsim2){
#This is for studying one of the models P1 - N3.
	if (model == "P1") {
    u <- Study_Randomcorrelations_positive(sizes, nsim, nsim2)}
	if (model == "P2") {
    u <- Study_autocorrelations_positive(sizes, nsim, nsim2)}
	if (model == "P3") {
    u <- Study_PositiveEquiCorrelated(sizes, nsim, nsim2)}
	if (model == "N1") {
    u <- Study_Randomcorrelations(sizes, nsim, nsim2)}
	if (model == "N2") {
    u <- Study_autocorrelations_negative(sizes, nsim, nsim2)}
	if (model == "N3") {
    u <- Study_minimumcorrelations(sizes, nsim, nsim2)}
  u <- data.frame(u)
  colnames(u)<-c("m","minr","meanr","maxr","rej")
	if (model == "P2") {
    colnames(u)<-c("m","a","b","maxr","rej")}
	if (model == "N2") {
    colnames(u)<-c("m","a","b","maxr","rej")}

  u <- data.frame(u)
  um<-aggregate(u, by=list(u$m),FUN=mean)
  us<-aggregate(u, by=list(u$m),FUN=sd)
  uag <-cbind(um$m, um$rej/nsim2, us$rej/nsim2)
  uag<-round(uag,3)
  return(u)
}


Study_Randomcorrelations <-function(sizes, nsim, nsim2) {
#This is for studying a range of correlation matrices of model N1.
  nlambda=10
  alpha = 0.05
  A <-array(0,5)
  for (size in sizes) {
	  for (isim in 1:nsim) {
		  A <- c(A, size, Simulate_Randomcorrelations(nsim2,nlambda,size,alpha))
	  }
	print(size)
  }
  nsizes = length(sizes)
  dim(A)=c(5,nsim*nsizes+1)
  A<-t(A)
  A<-A[-1,]
  return(A)
}


Study_Randomcorrelations_positive <-function(sizes, nsim, nsim2) {
#This is for studying a range of correlation matrices of model P1.
  nlambda=10
  alpha = 0.05
  A <-array(0,5)
  for (size in sizes) {
	  for (isim in 1:nsim) {
		  A <- c(A, size, Simulate_Randomcorrelations_positive(nsim2,nlambda,size,alpha))
	  }
	  print(size) ;
	  save(A,file=paste("P1",size, ".Rda", sep=""))
  }
  nsizes = length(sizes)
  dim(A)=c(5,nsim*nsizes+1)
  A<-t(A)
  A<-A[-1,]
  return(A)
}




Simulate_Randomcorrelations <-function(nsim, nlambda, size, alpha) {
#This is for simulating data with one correlation matrix of model N1.
  nrej <- array(0,nlambda)
  critZ <- CBP_CriticalValues(nlambda, size, alpha)
  qlambda <- 1:nlambda
  qlambda <- qnorm(qlambda/nlambda)

  B <-  rnorm(size^2,0,1)
  dim(B) <- c(size, size)
  B <- scale(B, scale=apply(B, 2, sd))
  B<-B/sqrt(size -1)
  B<-t(B)
  H2<-runif(size,0,1)
  H2<-diag(H2,size)
  H<-sqrt(H2)
  A<-H%*%B
  R<-A %*% t(A)
  L<-diag(R)
  L = sqrt(1-L)
	R<- R+diag(L*L)
  mincor=min(R)
  meancor=mean(R)
  R<-R-diag(2,size)
  maxcor=max(R)
  dim(L) <- c(size,1)
  for (i in 1:nsim) {
    Z <-rnorm(size, 0, 1)
    dim(Z) <- c(size,1)
    E<- rnorm(size,0,1)
    dim(E) <- c(size,1)
    X<- A%*%Z + L*E
      minx = min(X)
	    for (ilambda in 1:nlambda) {
        S = (X <= qlambda[ilambda])
        s = sum(S)
        if (s > 0) {
      	  if (minx < critZ[s, ilambda]){
      		  nrej[ilambda] = nrej[ilambda]+1
      	  }
        }
      }
  }
  return(c(mincor,meancor, maxcor, max(nrej)))
}



Simulate_Randomcorrelations_positive <-function(nsim, nlambda, size, alpha) {
#This is for simulating data with one correlation matrix of model P1.
  nrej <- array(0,nlambda)
  critZ <- CBP_CriticalValues(nlambda, size, alpha)
  qlambda <- 1:nlambda
  qlambda <- qnorm(qlambda/nlambda)

  B <-  rnorm(size^2,0,1)
  B<-abs(B)
  dim(B) <- c(size, size)
  B <- scale(B, center = FALSE, scale=TRUE)
  B<-B/sqrt(size -1)
  B<-t(B)
  H2<-runif(size,0,1)
  H2<-diag(H2,size)
  H<-sqrt(H2)
  A<-H%*%B
  R<-A %*% t(A)
  L<-diag(R)
  L = sqrt(1-L)
  R<- R+diag(L*L)
  mincor=min(R)
  meancor=mean(R)
  R<-R-diag(2,size)
  maxcor=max(R)
  dim(L) <- c(size,1)
  for (i in 1:nsim) {
    Z <-rnorm(size, 0, 1)
    dim(Z) <- c(size,1)
    E<- rnorm(size,0,1)
    dim(E) <- c(size,1)
    X <- A%*%Z + L*E
    minx = min(X)
	  for (ilambda in 1:nlambda) {
      S = (X <= qlambda[ilambda])
      s = sum(S)
      if (s > 0) {
      	if (minx < critZ[s, ilambda]){
      		nrej[ilambda] = nrej[ilambda]+1
      	}
      }
    }
  }
  return(c(mincor,meancor, maxcor, max(nrej)))
}



Study_autocorrelations_positive <-function(sizes, nsim, nsim2) {
#This is for studying a range of correlation matrices of model P2:
  nlambda=10
  minsize = 10
  alpha = 0.05
  direction = 1
  A <-array(0,5)
  for (size in sizes) {
	  for (isim in 1:nsim) {
		  A <- c(A, size, Simulate_autocorrelations(nsim2,nlambda,size,alpha, direction))
	  }
	  print(size);
	  save(A,file=paste("P2",size, ".Rda", sep=""))
  }
  nsizes = length(sizes)
  dim(A)=c(5, nsim*nsizes+1)
  A<-t(A)
  A<-A[-1,]
  return(A)
}



Study_autocorrelations_negative <-function(sizes, nsim, nsim2) {
#This is for studying a range of correlation matrices of model N2:
  nlambda=10
  minsize = 10
  alpha = 0.05
  direction = -1

  A <-array(0,5)
  for (size in sizes) {
	  for (isim in 1:nsim) {
		  A <- c(A, size, Simulate_autocorrelations(nsim2,nlambda,size,alpha, direction))
	  }
	  print(size)
  }
  nsizes = length(sizes)
  dim(A)=c(5, nsim*nsizes+1)
  A<-t(A)
  A<-A[-1,]
  return(A)
}



CBP_CriticalValues <- function(nlambda, size, alpha) {
#This is to compute a fixed array with critical values for the CBP, which makes some procedures faster.
  critZ <-array(0, size * nlambda)
  dim(critZ) <- c(size, nlambda)
  for (i in 1:size){
	  for (ilambda in 1:nlambda){
		  critZ[i, ilambda] <- alpha * (ilambda / nlambda) / i
	  }
  }
  critZ <- qnorm(critZ)
  return(critZ)
}



Simulate_autocorrelations <-function(nsim, nlambda, size, alpha, direction) {
#This is for simulating data with one correlation matrix of model P2 or N2:
  nrej <- array(0,nlambda)

  critZ <- CBP_CriticalValues(nlambda, size, alpha)

  qlambda <- 1:nlambda
  qlambda <- qnorm(qlambda/nlambda)


	apar = runif(1,0,1)
	bpar = runif(1,0,1)

  A <-  rbeta(size-1,apar,bpar)
  for (i in 1:nsim) {
    Z <-rnorm(size, 0, 1)
    X <-array(0,size)
    X[1] = Z[1]
    for (j in 2:size) {
      X[j] = (direction*A[j-1]*X[j-1] + sqrt(1-A[j-1]^2)*Z[j])
	  }
	  minx = min(X)
	  for (ilambda in 1:nlambda) {
      S = (X <= qlambda[ilambda])
      s = sum(S)
      if (s > 0) {
    	  if (minx < critZ[s, ilambda]){
          nrej[ilambda] = nrej[ilambda]+1
        }
      }
    }
  }
  maxcor = max(A)
  return(c(apar, bpar, maxcor, max(nrej)))
}



Study_PositiveEquiCorrelated <-function(sizes, nsim, nsim2) {
#This is for studying a range of correlation matrices of model P3:
  nlambda=10
  alpha = 0.05
  A <-array(0,5)
  for (size in sizes) {
	  for (isim in 1:nsim) {
		  A <- c(A, size, Simulate_PositiveEquiCorrelated(nsim2,nlambda,size,alpha))
	  }
	  print(size)
  }
  nsizes = length(sizes)
  dim(A)=c(5, nsim*nsizes+1)
  A<-t(A)
  A<-A[-1,]
  return(A)
}



Simulate_PositiveEquiCorrelated <-function(nsim, nlambda, size, alpha) {
#This is for simulating data with one correlation matrix of model P3:
  nrej <- array(0,nlambda)
  critZ <- CBP_CriticalValues(nlambda, size, alpha)

  qlambda <- 1:nlambda
  qlambda <- qnorm(qlambda/nlambda)

  cor = runif(1,0,1)
  for (i in 1:nsim) {
    Z <- rnorm(1, 0,1)
    E <-rnorm(size, 0, 1)
    X <- sqrt(cor) * Z + sqrt(1-cor) * E
	  minx = min(X)
	  for (ilambda in 1:nlambda) {
      S = (X <= qlambda[ilambda])
      s = sum(S)
      if (s > 0) {
      	if (minx < critZ[s, ilambda]){
      		nrej[ilambda] = nrej[ilambda]+1
      	}
      }
    }
  }
  return(c(cor, cor, cor, max(nrej)))
}



Study_minimumcorrelations <-function(sizes, nsim, nsim2) {
#This is for studying a range of correlation matrices of model N3:
  nlambda=10
  alpha = 0.05
  A <-array(0,5)
  for (size in sizes) {
	  for (isim in 1:nsim) {
		  A <- c(A, size, Simulate_minimumcorrelations(nsim2,nlambda,size,alpha))
	  }
	  print(size)
  }
  nsizes = length(sizes)
  dim(A)=c(5, nsim*nsizes+1)
  A<-t(A)
  A<-A[-1,]
  return(A)
}



Simulate_minimumcorrelations <-function(nsim, nlambda, size, alpha) {
#This is for simulating data with one correlation matrix of model P3:
  nrej <- array(0,nlambda)
  critZ <- CBP_CriticalValues(nlambda, size, alpha)

  qlambda <- 1:nlambda
  qlambda <- qnorm(qlambda/nlambda)

  for (i in 1:nsim) {
    Z <-rnorm(size, 0, 1)
    X <- Z-mean(Z)
    X<- X/sqrt(1-1/size)
	  minx = min(X)
	  for (ilambda in 1:nlambda) {
      S = (X <= qlambda[ilambda])
      s = sum(S)
      if (s > 0) {
    	  if (minx < critZ[s, ilambda]){
    		  nrej[ilambda] = nrej[ilambda]+1
    	  }
      }
    }
  }
  mincor = -1/(size-1)
  return(c(mincor, mincor, mincor, max(nrej)))
}

# In the procedures below, pvec is the vector with p-values


Bonf <- function(m, pvec, alpha) {
# This is to apply the Bonferroni multiple testing procedure:
# Applies the Bonferroni MTP
  minp = min(pvec)*m
  if (minp < alpha) {
    reject = 1
  } else {
    reject = 0
	}
	return(reject)
}



FGS <- function(m, pvec, alpha, lambda) {
# This is to apply the Finner-Gontscharuk modification of the Storey MTP:
  r= sum(pvec <= lambda)
  kappa=min(1, alpha * (1 - lambda) / lambda)
  n0 = (m - r + kappa) / (1 - lambda)
  if (n0 > 0) {
    adj_alpha = alpha / n0
  } else {
      adj_alpha = 0
  }
  rejectcount = sum(pvec < adj_alpha)

  if (rejectcount > 0){
	  reject=1
  } else {
    reject = 0
  }
  return(reject)
}




Fisher <- function(m, pvec, alpha) {
# Applies the Fisher MTP
  chi2 = -2*sum(log(pvec))
  df = 2 * m
  if (chi2 == 0) {
	  p = 1
  } else {
    p = pchisq(chi2,df,lower.tail = FALSE)
  }
  if (p < alpha) {
	  reject=1
	} else {
	  reject = 0
  }
  return(reject)
}



Tippet <- function(m, pvec, alpha){
# Applies the Tippet (Sidak) MTP
  minp = min(pvec)
  minp_adj = pexp(-2 * log((1 - minp)), m / 2)
  if (minp_adj < alpha) {
	  reject = 1
	} else {
	  reject = 0
  }
  return(reject)
}



LR <- function(m, pvec, alpha) {
# This is to apply the Likelihood Ratio multiple testing procedure for order restricted inference:
  z = qnorm(pmin(0.5,pvec))
  chi2 = sum(z*z)
  df = sum(pvec <0.5)

	if (df ==0) {
		p = 1
	} else {
    p = 0
    for (i in 1:m) {
      a = dbinom(i, m, 0.5)
    	p = p + a * (1 - pchisq(chi2, df))
    }
	}
  if (p < alpha) {
    reject = 1
  } else {
	  reject = 0
  }
  return(reject)
}



integrateIplus <- function(m, pvec){
# AuX function to compute Iplus.
	area = 0
  psort <- sort(pvec)
  edf = array(0,m)
  for (i in 1:m) {
  	for (j in 1:m){
  		if (psort[j] <= psort[i]){
  			edf[i] = edf[i]+1
  		}
  	}
  }
  edf <-edf/m
  psort <- c(0, psort)
  edf <- c(0, edf)

  for (i in 2:(m+1)) {
  	area <- area+integrateIplus2(edf[i - 1], edf[i], psort[i-1], psort[i])
	}
  area <- area+integrateIplus2(edf[m+1], 1, psort[m+1], 1)
  return(area)
}


integrateIplus2 <- function(f1, f2, x1, x2){
# AuX function to compute Iplus.
	if (x1 >= x2) {
		area <-0
	} else {
		if (f1 >= x2) {
			area <- integrateIplus3(f1, f2, x1, x2)
		}
		if (f1 < x1) {
			area <-0
		}
		if (f1 >= x1) {
			if (f1 < x2) {
				x3 <- f1
				area <- integrateIplus3(f1, f2, x1, x3)
			}
		}
	}
	return(area)
}


integrateIplus3 <- function(f1, f2, x1, x2){
# AuX function to compute Iplus.
  if(x2 == 1) {
    y2 <- -1
  } else {
    y2 <- PrimIplus(f1, x2)
  }
  y1 = PrimIplus (f1, x1)
  return(y2 - y1)
}




PrimIplus <- function(f, x){
# Aux function to compute Iplus.
  if (x == 0) {
    y <- 0
  } else {
    y <- f * f * log(x / (1 - x)) + 2 * f * log(1 - x) - x - log(1 - x)
  }
  return(y)
}



Iplus <- function(m, pvec) {
# This is to apply the Iplus multiple testing procedure.
  area = integrateIplus(m, pvec)
  stat = m * area

  rejectcount = 0
  crit <-array(0,7)
  n <-array(0,7)

  n[1] = 2
  n[2] = 3
  n[3] = 5
  n[4] = 10
  n[5] = 15
  n[6] = 20
  n[7] = 30

  crit[1] = 1.954
  crit[2] = 1.912
  crit[3] = 1.896
  crit[4] = 1.876
  crit[5] = 1.87
  crit[6] = 1.87
  crit[7] = 1.859

  critval = approx(n, crit, m,rule=2)$y
  if (stat > critval) {
    reject = 1
  } else {
	  reject = 0
  }
}



Simulate_Power <- function(nsim, ncp0, ncp1,sigma0, sigma1,n0,n1,lambda,alpha){
#This is to simulate the power in a single setting of figures 2 and 3.
  nrej<-array(0,12)
  for (isim in 1:nsim) {
	  z0  = rnorm(n0,ncp0,sigma0)
	  z1  = rnorm(n1,-ncp1,sigma1)
	  z = c(z0,z1)
	  p = pnorm(z)
	  m = n0+n1
	  nrej[1] <- nrej[1] +Bonf(m,p,alpha)
	  nrej[2] <- nrej[2] +FGS(m,p,alpha,lambda)
	  nrej[3] <- nrej[3] + Fisher(m,p,alpha)
	  nrej[4] <- nrej[4] + Tippet(m,p,alpha)
	  nrej[5] <- nrej[5] + Iplus(m,p)
	  nrej[6] <- nrej[6] + LR(m,p, alpha)
	  pc<-p[p<= lambda]/lambda
	  mc <- length(pc)
	  if (mc > 0){
	    nrej[7] <- nrej[7] +Bonf(mc,pc,alpha)
	    nrej[8] <- nrej[8] +FGS(mc,pc,alpha,lambda)
	    nrej[9] <- nrej[9] + Fisher(mc,pc,alpha)
	    nrej[10] <- nrej[10] + Tippet(mc,pc,alpha)
	    nrej[11] <- nrej[11] + Iplus(mc,pc)
	    nrej[12] <- nrej[12] + LR(mc,pc, alpha)
	  }
  }
  return(nrej)
}




Study_Power_Fig2<-function(nsim) {
#This is to study the power in a range of situations, as displayed in figure 2.
  ndpar = 20
  ntest = 12
  lambda = 0.5
  alpha = 0.05
  ncp0vec <- array(2,ndpar)
  ncp1vec <- array(2,ndpar)
  sigma0vec <- array(1,ndpar)
  sigma1vec <- array(1,ndpar)
  n1vec <- array(5,ndpar)
  n0vec <- 1:20
  n0vec<-5*n0vec
  rej = array(0, 0)
  for (idpar in 1:ndpar) {
	  ncp0 = ncp0vec[idpar]
	  ncp1 = ncp1vec[idpar]
	  sigma0 = sigma0vec[idpar]
	  sigma1 = sigma1vec[idpar]
	  n0 = n0vec[idpar]
	  n1 = n1vec[idpar]
    u <- Simulate_Power(nsim, ncp0, ncp1,sigma0, sigma1,n0,n1,lambda,alpha)
    rej<-c(rej,ncp0,ncp1,sigma0,sigma1,n0,n1,u/nsim)
  }
  dim(rej)<-c(6 + ntest,ndpar)
  return(rej)
}




Study_Power_Fig3<-function(nsim) {
#This is to study the power in a range of situations, as displayed in figure 2.
  ndpar = 10
  ntest = 12
  lambda = 0.5
  alpha = 0.05
  ncp0vec <- array(1.5,ndpar)
  ncp1vec <- array(1.5,ndpar)
  sigma0vec <- array(1,ndpar)
  sigma1vec <- array(1,ndpar)
  n0vec <- 1:ndpar
  n0vec<-10*n0vec
  n1vec <- -n0vec+ 10*ndpar
  rej = array(0, 0)
  for (idpar in 1:ndpar) {
	  ncp0 = ncp0vec[idpar]
	  ncp1 = ncp1vec[idpar]
	  sigma0 = sigma0vec[idpar]
	  sigma1 = sigma1vec[idpar]
	  n0 = n0vec[idpar]
	  n1 = n1vec[idpar]
    u <- Simulate_Power(nsim, ncp0, ncp1,sigma0, sigma1,n0,n1,lambda,alpha)
    rej<-c(rej,ncp0,ncp1,sigma0,sigma1,n0,n1,u/nsim)
  }
  dim(rej)<-c(6 + ntest,ndpar)
  return(rej)
}



Simulate_FDR <- function(nsim, ncp0, ncp1,sigma0, sigma1,n0,n1,lambda,alpha){
#This is to simulate the power of BH in a single setting of figures 2 and 3.
  efdru = 0
  etpru = 0
  efdrc = 0
  etprc = 0
  for (isim in 1:nsim) {
	  z0  = rnorm(n0,ncp0,sigma0)
	  z1  = rnorm(n1,-ncp1,sigma1)
	  z = c(z0,z1)
	  p = pnorm(z)
	  m = n0+n1
	  num <- 1:m
	  pnum <- rbind(num, p)
	  cnum <- num[p < lambda]
	  pu <- p.adjust( p, method="BH" )
      punum <-rbind(num,pu)
      rejected <- num[pu < alpha]
      ru = length(rejected)
      fu = length(rejected[rejected <= n0])
      tu = ru - fu
      if (ru> 0) {
        fdru = fu / ru
      } else {
      	fdru = 0
      }
      if (n1> 0) {
      	tpru = tu/n1
      } else {
      	tpru = 0
      }

      if (length(cnum) > 0) {
      	pc <- p.adjust( p[ p < lambda ] /lambda, method="BH" )
      	pcnum <-rbind(cnum,pc)
      	rejected <- cnum[pc < alpha]
      	rc = length(rejected)
      	fc = length(rejected[rejected <= n0])
      	tc = rc - fc
      	if (rc > 0) {
      		fdrc = fc / rc
      	} else {
      		fdrc = 0
      	}
      if (n1> 0) {
      	tprc = tc/n1
      } else {
      	tprc = 0
      }
    } else {
      fdrc = 0
      tprc = 0
    }
	  efdru <- efdru + fdru
	  etpru <- etpru + tpru
	  efdrc <- efdrc + fdrc
	  etprc <- etprc + tprc
  }
  efdru <- efdru / nsim
  etpru <- etpru/ nsim
  efdrc <- efdrc / nsim
  etprc <- etprc / nsim
  return(c(efdru, etpru, efdrc, etprc))
}



Study_FDR<- function(nsim) {
#This is to study the power of the BH in a range of situations, as displayed in figure 2.
  ncp0 = 2
  ncp1 = 2
  sigma0 = 1
  sigma1 = 1
  n0vec <- 1:20
  n0vec <- 5 * n0vec
  n1 = 5
  lambda = 0.5
  alpha = 0.05
  cdr = array(0,0)
  for (n0 in n0vec) {
    dr <- Simulate_FDR (nsim, ncp0, ncp1,sigma0, sigma1,n0,n1,lambda,alpha)
		cdr <- c(cdr, dr)
	}
	return(cdr)
}




Relative_Power <- function(lambda, pi0, delta, epsilon){
#This is a function to compute rho, the ratio of CBP-corrected and uncorrected z-statistics.
  z = qnorm(lambda,0,1)
  ksi = pi0 * pnorm(z-delta) + (1-pi0) * pnorm (z+epsilon)
  return(ksi / lambda)
}



Study_RelativePower <- function(){
#This is to generate Table 2.
  deltavec = c(0.1, 0.2, 0.4, 0.8, 1.6, 3.2)
  pivec = seq(0.1, 0.9, by = 0.1)
  rho<-array(0,9*6)
  dim(rho)<-c(9, 6)

  for (i in 1:9) {
	  for (j in 1:6) {
		  rho[i,j] <- Relative_Power(0.5, pivec[i], deltavec[j], 50)
	  }
  }
  return(rho)
}



StudyMiniMax <- function(npi,minpi,maxpi,ndelta,maxdelta,nepsilon,maxepsilon,nlambda,minlambda) {
#This is to compute the minimax value of lambda.
  ratio <- array(0, dim= c(npi, ndelta, nepsilon, nlambda))
#array with ratios ksi / lambda
  minratios <-array(0, c(npi, ndelta, nepsilon))
#array with minimum ratio per (pi, delta, epsilon) cell, minimized over lambda
  maxim <-array (0, nlambda)
#array with maximum of ratio / minratios per value of lambda
  maxpar <-array (0, c(nlambda, 3))
#array with for each lambda the three parameters (pi, delta, epsilon) at the maximum
  for (ipi in 1:npi) {
    pi0 = minpi + (maxpi - minpi) * (ipi-1) / npi
    pi1 = 1 - pi0
    for (idelta in 1 : ndelta) {
      delta = maxdelta * idelta / ndelta
      for (iepsilon in 1 : nepsilon) {
        epsilon = maxepsilon * iepsilon / nepsilon
        minratio = 1
        for (ilambda in 1:nlambda) {
          lambda = minlambda + (1 - minlambda)* (ilambda-1) / nlambda
          z = qnorm(lambda,0,1)
          ksi = pi0 * pnorm(z-delta) + pi1 * pnorm (z+epsilon)
          r = ksi / lambda
          ratio[ipi, idelta, iepsilon, ilambda] = r
          if (r < minratio) {
            minratio = r
          }
        }
        minratios[ipi, idelta, iepsilon] = minratio
      }
    }
  }

  for (ilambda in 1 : nlambda) {
    maxim[ilambda] = 0
    for (ipi in 1 : npi)
      pi0 = minpi + (maxpi - minpi) * (ipi-1) / npi
      pi1 = 1 - pi0
      lambda = minlambda + (1 - minlambda)* (ilambda-1) / nlambda
      for (idelta in 1 : ndelta) {
        for (iepsilon in 1 : nepsilon) {
          r = ratio[ipi, idelta, iepsilon, ilambda] / minratios[ipi, idelta, iepsilon]
          if (r > maxim[ilambda] ) {
            maxim[ilambda] = r
            maxpar[ilambda, 1] = ipi
            maxpar[ilambda, 2] = idelta
            maxpar[ilambda, 3] = iepsilon
          }
        }
      }
    }
  minimax = 0
  for (ilambda in 1 : nlambda){
    if ((maxim[ilambda] < minimax)||(minimax==0)) {
      minimax = maxim[ilambda]
      minimaxi = ilambda
    }
  }
  minimaxlambda = minlambda + (1 - minlambda)* (minimaxi-1) / nlambda
  minimaxlambda
  # this is the lambda that yields the minimax risk as given by ratio / minratios
  minimax
  # this the minimax risk
  pi = minpi + (maxpi - minpi) * (maxpar[minimaxi, 1]-1) / npi
  delta = maxdelta * maxpar[minimaxi, 2] / ndelta
  epsilon = maxepsilon* maxpar[minimaxi, 3] / nepsilon
  # these are the parameters at which the minimax risk is obtained
  return(c(minimaxlambda, minimax, pi, delta, epsilon))
}



Study_Power_Fig4 <- function(nsim, effects, sizes, lambda, alpha, fraction) {
  #This is to study the power in a range of situations, as displayed in figure 4 of the Supplementary Material of the paper.
  nncp = length(effects)
  nn = length(sizes)
  ndpar = nncp * nncp * nn
  print(ndpar)
  ntest = 4
  sigma0 = 1
  sigma1 = 1
  rej = array(0, 0)
  i=1
  for (ncp0 in effects){
    for (ncp1 in effects) {
      for (n in sizes) {
        n0 = fraction*n
        n1 = (1-fraction)*n
        u <- Simulate_Power2(nsim, ncp0, ncp1,sigma0, sigma1,n0,n1,lambda,alpha)
        rej<-c(rej,ncp0,ncp1,sigma0,sigma1,n0,n1,u/nsim)
        i=i+1
        print(i)
      }}}
  dim(rej)<-c(6 + ntest,ndpar)
  rej<-t(rej)
  rej <-data.frame(rej)
  colnames(rej)<-c("ncp0","ncp1", "sigma0","sigma1","n0", "n1", "Bonferroni","FGS","Conditional Bonferroni","Conditional FGS")
  return(rej)
}

Simulate_Power2 <- function(nsim, ncp0, ncp1,sigma0, sigma1,n0,n1,lambda,alpha){
  #This is to simulate the power in a single setting of figure 4 of the Supplementary Material, but
  #computes only Bonferroni and FGS to save time.
  nrej<-array(0,4)
  for (isim in 1:nsim) {
    z0  = rnorm(n0,ncp0,sigma0)
    z1  = rnorm(n1,-ncp1,sigma1)
    z = c(z0,z1)
    p = pnorm(z)
    m = n0+n1
    nrej[1] <- nrej[1] +Bonf(m,p,alpha)
    nrej[2] <- nrej[2] +FGS(m,p,alpha,lambda)
    pc<-p[p<= lambda]/lambda
    mc <- length(pc)
    if (mc > 0){
      nrej[3] <- nrej[3] +Bonf(mc,pc,alpha)
      nrej[4] <- nrej[4] +FGS(mc,pc,alpha,lambda)
    }
  }
  return(nrej)
}
