#########################################
##### Spatial Multifractal Analysis #####
#########################################

# Griding the data: free boundaries

randtoGrid <- function(X,Y,Z,MX,MY,fun)
  # X: 1st coord, Y: 2nd coord, Z: values, 
  # MX = c(min,max,nbr of X intervals), MY=c(min,max,nbr of Y intervals)
  # fun = "mean" or "sum"
{
  SizeX <- MX[3]
  SizeY <- MY[3]
  result <- matrix(0,ncol=SizeX,nrow=SizeY)
  coordX <- MX[1]+(0:(SizeX-1))*(MX[2]-MX[1])/(SizeX-1)
  coordY <- MY[1]+(0:(SizeY-1))*(MY[2]-MY[1])/(SizeY-1)
  etiquette <- data.frame(Z=Z,coords=rep(0,length(Z)))
  for(i in 1:length(Z)){
    X0 <- max(coordX[which(coordX <= X[i])])
    Y0 <- max(coordY[which(coordY <= Y[i])])
    etiquette$coords[i] <- paste(X0,",",Y0,sep="")
  }
  values <- aggregate(etiquette$Z,by=list(etiquette$coords),FUN=get(fun))
  for(i in 1:SizeX){
    for(j in 1:SizeY){
      text0 <- paste(coordX[i],",",coordY[j],sep="")
      if(any(values$Group.1 == text0)){
        result[j,i] <- values$x[which(values$Group.1==text0)]
      }
    }
  }
  return(result)
}

# Griding the data: fixed boundaries

newrandtoGrid <- function(X,Y,Z,MX,MY,fun){
  # X: 1st coord, Y: 2nd coord, Z: values, 
  # MX = nbr of X intervals, MY = nbr of Y intervals
  # fun = "mean" or "sum" f.e.
  mx <- min(X)
  Mx <- max(X)
  my <- min(Y)
  My <- max(Y)
  intx <- seq(mx,Mx-(Mx-mx)/MX,by=(Mx-mx)/MX)
  inty <- seq(my,My-(My-my)/MY,by=(My-my)/MY)
  coords <- data.frame(Z=Z,etiquette=rep(0,length(Z)))
  result <- matrix(0,nrow=MY,ncol=MX)
  for(i in 1:length(Z)){
    X0 <- max(which(intx <= X[i]))
    Y0 <- max(which(inty <= Y[i]))
    coords$etiquette[i] <- paste(X0,".",Y0,sep="")
  }
  values <- aggregate(coords$Z,by=list(coords$etiquette),FUN=get(fun))
  for(i in 1:MX){
    for(j in 1:MY){
      text0 <- paste(i,".",j,sep="")
      if(any(values$Group.1 == text0)){
        result[j,i] <- values$x[which(values$Group.1==text0)]
      }
    }
  }
  return(result)
}

# Aggregation: Gliding box

GlidingBox <- function(InitGrid,sizes,neigh=1)
{
  InitRows <- nrow(InitGrid)
  InitCols <- ncol(InitGrid)
  Raster <- raster(InitGrid)
  msize <- max(sizes)
  ResultGrid <- rep(list(matrix(nrow=InitRows-2*msize,ncol=InitCols-2*msize,0)),length(sizes))
  for(s in sizes){
    if(s==0){
      ResultGrid[[which(sizes == 0)]] <- InitGrid[(msize+1):(InitRows-msize),(msize+1):(InitCols-msize)]
    } else {
      weightmatrix <- matrix(nrow=2*s+1,ncol=2*s+1,1)
      if(neigh==2){
        for(Index in 0:((2*s+1)*(2*s+1)-2)){
          locrow <- floor(Index/(2*s+1)) + 1
          loccol <- Index %% (2*s+1) + 1
          if((abs(s+1-locrow)+abs(s+1-loccol) > s)){
            weightmatrix[locrow,loccol] <- 0
          }
        }
      }
      result <- focal(x=Raster,w=weightmatrix,fun=sum,pad=T,na.rm=F)
      ResultGrid[[which(sizes == s)]] <- as.matrix(result)[(msize+1):(InitRows-msize),(msize+1):(InitCols-msize)]
    }
  }
  return(ResultGrid)
}

# Aggregation: standard box-counting

GridingAgg <- function(InitGrid,multp,startpos=1)
  # 'multp' is a vector containing the new grid unit sizes as multiples of the initial unit size
  # 'startpos' indicates from which corner to start (counted clockwise)
{
  lm <- length(multp)
  InitRows <- nrow(InitGrid)
  InitCols <- ncol(InitGrid)
  ResultGrid <- rep(list(NULL),lm)
  for(i in 1:lm){
    locmultp <- multp[i]
    locnrow <- floor(InitRows/locmultp)
    locncol <- floor(InitCols/locmultp)
    StepMatrix <- matrix(nrow=locnrow,ncol=locncol,0)
    for(Index in 0:(locncol*locnrow-1)){
      locrow <- floor(Index/locncol) + 1
      loccol <- Index %% locncol + 1
      if(startpos == 1){
        StepMatrix[locrow,loccol] <- sum(InitGrid[((locrow-1)*locmultp+1):(locrow*locmultp),((loccol-1)*locmultp+1):(loccol*locmultp)])
      } else if(startpos == 2){
        StepMatrix[locrow,locncol-loccol+1] <- sum(InitGrid[((locrow-1)*locmultp+1):(locrow*locmultp),(InitCols-(loccol-1)*locmultp):(InitCols-loccol*locmultp+1)]) 
      } else if(startpos == 3){
        StepMatrix[locrow,loccol] <- sum(InitGrid[(InitRows-(locrow-1)*locmultp):(InitRows-locrow*locmultp+1),(InitCols-(loccol-1)*locmultp):(InitCols-loccol*locmultp+1)])
      } else{
        StepMatrix[locrow,loccol] <- sum(InitGrid[(InitRows-(locrow-1)*locmultp):(InitRows-locrow*locmultp+1),((loccol-1)*locmultp+1):(loccol*locmultp)])  
      }
    }     
    ResultGrid[[i]] <- StepMatrix
  }
  return(ResultGrid)
}

# Main multifractal methods

Moment_StdGrid <- function(data,sizes,q,dd)
{
  data <- GridingAgg(data,sizes)
  Zed <- matrix(rep(0,length(q)*length(sizes)),ncol=length(sizes),byrow=TRUE)
  Tau <- rep(0,length(q))
  D <- rep(0,length(q))
  alpha <- rep(0,length(q))
  falpha <- rep(0,length(q))
  for(i in 1:length(q)){
    for(j in 1:length(sizes)){
      Zed[i,j] <- sum(data[[j]][which(data[[j]] != 0)]^q[i])
    }
    fit <- lm(log(Zed[i,]) ~ log(sizes))
    Tau[i] <- fit$coefficients[[2]]
    D[i] <- Tau[i]/(q[i]-1)
  }
  for(i in 2:(length(q)-1)){
    alpha[i] <- (Tau[i+1]-Tau[i-1])/(q[i+1]-q[i-1])
    falpha[i] <- q[i]*alpha[i]-Tau[i]
  }
  alpha[1] <- (Tau[2]-Tau[1])/(q[2]-q[1])
  alpha[length(q)] <- (Tau[length(q)]-Tau[length(q)-1])/(q[length(q)]-q[length(q)-1])
  falpha[1] <- q[1]*alpha[1]-Tau[1]
  falpha[length(q)] <- q[length(q)]*alpha[length(q)]-Tau[length(q)]
  D[match(1,q)] <- (D[match(1,q)+1]+D[match(1,q)-1])/2
  Result <- data.frame(q=q,Dq=D,alpha=alpha,falpha=falpha)
  return(Result)
}

MMoment_MGBox <- function(grid,sizes,q,dd,neigh=1)
{
  data <- GlidingBox(grid,sizes,neigh)
  sizes <- sizes+0.5
  emin <- min(sizes)
  N <- nnzero(data[[1]])
  Zed <- matrix(rep(0,length(q)*(length(sizes)-1)),ncol=length(sizes)-1,byrow=TRUE)
  prealpha <- matrix(rep(0,length(q)*(length(sizes)-1)),ncol=length(sizes)-1,byrow=TRUE)
  Tau <- rep(0,length(q))
  D <- rep(0,length(q))
  alpha <- rep(0,length(q))
  falpha <- rep(0,length(q))
  for(i in 1:length(q)){
    for(j in 2:length(sizes)){
      Zed[i,j-1] <- -log(sum((data[[1]][which(data[[1]] != 0)]/data[[j]][which(data[[1]] != 0)])^q[i])/N)/log(sizes[j]/emin)
      prealpha[i,j-1] <- -sum(((data[[1]][which(data[[1]] != 0)]/data[[j]][which(data[[1]] != 0)])^q[i])*log(data[[1]][which(data[[1]] != 0)]/data[[j]][which(data[[1]] != 0)]))/(sum((data[[1]][which(data[[1]] != 0)]/data[[j]][which(data[[1]] != 0)])^q[i])*log(sizes[j]/emin))
    }
    Tau[i] <- mean(Zed[i,])
  }
  Tau1 <- Tau[which(q==1)]
  for(i in 1:length(q)){
    Tau[i] <- Tau[i]*dd/Tau1-dd
    D[i] <- Tau[i]/(q[i]-1)
    alpha[i] <- dd/Tau1*mean(prealpha[i,])
    falpha[i] <- q[i]*alpha[i]-Tau[i]
  }
  Result <- data.frame(q=q,Dq=D,tau=Tau,alpha=alpha,falpha=falpha)
  return(Result)
}

MHistogramU_2 <- function(data,sizes,n,dd=2,startpos=1)
{
  l <- length(sizes)
  InitRows <- nrow(data)
  InitCols <- ncol(data)
  ResultGrid <- rep(list(NULL),l)
  for(i in 1:l){
    locmultp <- sizes[i]
    locnrow <- floor(InitRows/locmultp)
    locncol <- floor(InitCols/locmultp)
    StepMatrix <- matrix(nrow=locnrow,ncol=locncol,0)
    for(Index in 0:(locncol*locnrow-1)){
      locrow <- floor(Index/locncol) + 1
      loccol <- Index %% locncol + 1
      if(startpos == 1){
        StepMatrix[locrow,loccol] <- sum(data[((locrow-1)*locmultp+1):(locrow*locmultp),((loccol-1)*locmultp+1):(loccol*locmultp)])
      } else if(startpos == 2){
        StepMatrix[locrow,locncol-loccol+1] <- sum(data[((locrow-1)*locmultp+1):(locrow*locmultp),(InitCols-(loccol-1)*locmultp):(InitCols-loccol*locmultp+1)]) 
      } else if(startpos == 3){
        StepMatrix[locrow,loccol] <- sum(data[(InitRows-(locrow-1)*locmultp):(InitRows-locrow*locmultp+1),(InitCols-(loccol-1)*locmultp):(InitCols-loccol*locmultp+1)])
      } else{
        StepMatrix[locrow,loccol] <- sum(data[(InitRows-(locrow-1)*locmultp):(InitRows-locrow*locmultp+1),((loccol-1)*locmultp+1):(loccol*locmultp)])  
      }
    }
    ResultGrid[[i]] <- StepMatrix
  }
  Xarray <- matrix(rep(0,l*n),ncol=n)
  N <- matrix(rep(0,l*n),ncol=n)
  Histo <- matrix(rep(0,l*n),ncol=n)
  GrandX <- lapply(ResultGrid,function(x) log(x))
  minX <- rep(0,l)
  maxX <- rep(0,l)
  DeltaX <- rep(0,l)
  for(i in 1:l){
    minX[i] <- min(GrandX[[i]][which(GrandX[[i]]>=0)])
    maxX[i] <- max(GrandX[[i]])
    DeltaX[i] <- maxX[i]/n-minX[i]/n
  }
  nzero <- numeric()
  for(j in 1:n){
    for(i in 1:l){
      Xarray[i,j] <- minX[i]+(j-1/2)*DeltaX[i]
      N[i,j] <- length(which((GrandX[[i]] >= minX[i]+(j-1/2)*DeltaX[i]-sqrt(DeltaX[i])) & (GrandX[[i]] <= minX[i]+(j-1/2)*DeltaX[i]+sqrt(DeltaX[i]))))
      Histo[i,j] <- log(N[i,j])
    }
    if(min(Histo[,j]) > - 10000){
      nzero <- c(nzero,j)
    }
  }
  alpha <- rep(0,length(nzero))
  falpha <- rep(0,length(nzero))
  for(j in 1:length(nzero)){
    fit <- lm(Xarray[,nzero[j]] ~ log(sizes/dd))
    alpha[j] <- fit$coefficients[[2]]
    fit2 <- lm(Histo[,nzero[j]] ~ log(sizes/dd))
    falpha[j] <- -fit2$coefficients[[2]]
  }
  Result <- data.frame(alpha=alpha,falpha=falpha)
  return(Result)
}

MHistogramC_2 <- function(data,sizes,n,dd=2,startpos=1)
{
  l <- length(sizes)
  InitRows <- nrow(data)
  InitCols <- ncol(data)
  ResultGrid <- rep(list(NULL),l)
  for(i in 1:l){
    locmultp <- sizes[i]
    locnrow <- floor(InitRows/locmultp)
    locncol <- floor(InitCols/locmultp)
    StepMatrix <- matrix(nrow=locnrow,ncol=locncol,0)
    for(Index in 0:(locncol*locnrow-1)){
      locrow <- floor(Index/locncol) + 1
      loccol <- Index %% locncol + 1
      if(startpos == 1){
        StepMatrix[locrow,loccol] <- sum(data[((locrow-1)*locmultp+1):(locrow*locmultp),((loccol-1)*locmultp+1):(loccol*locmultp)])
      } else if(startpos == 2){
        StepMatrix[locrow,locncol-loccol+1] <- sum(data[((locrow-1)*locmultp+1):(locrow*locmultp),(InitCols-(loccol-1)*locmultp):(InitCols-loccol*locmultp+1)]) 
      } else if(startpos == 3){
        StepMatrix[locrow,loccol] <- sum(data[(InitRows-(locrow-1)*locmultp):(InitRows-locrow*locmultp+1),(InitCols-(loccol-1)*locmultp):(InitCols-loccol*locmultp+1)])
      } else{
        StepMatrix[locrow,loccol] <- sum(data[(InitRows-(locrow-1)*locmultp):(InitRows-locrow*locmultp+1),((loccol-1)*locmultp+1):(loccol*locmultp)])  
      }
    }
    ResultGrid[[i]] <- StepMatrix
  }
  Xarray <- matrix(rep(0,l*n),ncol=n)
  N <- matrix(rep(0,l*n),ncol=n)
  Histo <- matrix(rep(0,l*n),ncol=n)
  GrandX <- lapply(ResultGrid,function(x) log(x))
  minX <- rep(0,l)
  maxX <- rep(0,l)
  DeltaX <- rep(0,l)
  for(i in 1:l){
    minX[i] <- min(GrandX[[i]][which(GrandX[[i]]>=0)])
    maxX[i] <- max(GrandX[[i]])
    DeltaX[i] <- maxX[i]/n-minX[i]/n
  }
  nzero <- numeric()
  for(j in 1:n){
    for(i in 1:l){
      Xarray[i,j] <- minX[i]+(j-1/2)*DeltaX[i]
      N[i,j] <- length(which((GrandX[[i]] >= minX[i]+(j-1/2)*DeltaX[i]-DeltaX[i]) & (GrandX[[i]] <= minX[i]+(j-1/2)*DeltaX[i]+DeltaX[i])))
      Histo[i,j] <- log(N[i,j])
    }
    if(min(Histo[,j]) > - 10000){
      nzero <- c(nzero,j)
    }
  }
  alpha <- rep(0,length(nzero))
  falpha <- rep(0,length(nzero))
  for(j in 1:length(nzero)){
    fit <- lm(Xarray[,nzero[j]] ~ log(sizes/dd))
    alpha[j] <- fit$coefficients[[2]]
    fit2 <- lm(Histo[,nzero[j]] ~ log(sizes/dd))
    falpha[j] <- -fit2$coefficients[[2]]
  }
  Result <- data.frame(alpha=alpha,falpha=falpha)
  return(Result)
}

#############################################
##### Independent computation of D0D1D2 #####
#############################################

boxCounting <- function(M,S)
{
  nr <- nrow(M)
  nc <- ncol(M)
  N <- rep(0,length(S))
  for(s in 1:length(S)){
    ss <- S[s]
    for(i in 0:(nr/ss-1)){
      for(j in 0:(nc/ss-1)){
        if(sum(M[i*ss+1:ss,j*ss+1:ss]) > 0){N[s] <- N[s]+1}
      }
    }
  }
  fit <- lm(log(N) ~ log(S))
  return(-fit$coefficients[2])
}

D0D1D2 <- function(Data,sizes)
{
  data <- GridingAgg(Data,sizes)
  vectD1 <- rep(0,length(sizes))
  vectD2 <- rep(0,length(sizes))
  for(i in 1:length(sizes)){
    data[[i]] <- data[[i]]/sum(data[[i]])
    vectD1[i] <- sum(data[[i]][which(data[[i]]>0)]*log(data[[i]][which(data[[i]]>0)])/log(10))
    vectD2[i] <- log(sum(data[[i]]^2))/log(10)
  }
  S <- log(sizes/2)/log(10)
  fit1 <- lm(vectD1 ~ S)
  fit2 <- lm(vectD2 ~ S)
  result <- c(boxCounting(Data,sizes),fit1$coefficients[2],fit2$coefficients[2])
  return(result)
}

#########################
##### Multiplotting #####
#########################

multiplot <- function(...,plotlist=NULL,file,cols=1,layout=NULL) {
  library(grid)
  plots <- c(list(...),plotlist)
  numPlots = length(plots)
  if (is.null(layout)){
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


###################
##### Entropy #####
###################

entropy <- function(physSpace,radii,nbin=3,logb=1){
  # physSpace : list of matrices at finest resolution representing data for each category (normalized to be between 0 and 1)
  # radii : vector of radii representing scales to be tested by increasing order
  # nbin : number of binning intervals for each category's characteristics
  # logb : base for the logarithm
  size <- c(nrow(physSpace[[1]]),ncol(physSpace[[1]]))
  microstates <- matrix("s",nrow=size[1],ncol=size[2])
  maxr <- max(radii)
  for(i in 1:length(physSpace)){
    physSpace[[i]] <- rbind(physSpace[[i]][(maxr+1):2,],physSpace[[i]],physSpace[[i]][(size[1]-1):(size[1]-maxr),])
    physSpace[[i]] <- cbind(physSpace[[i]][,(maxr+1):2],physSpace[[i]],physSpace[[i]][,(size[2]-1):(size[2]-maxr)])
  }
  for(i in 1:length(physSpace)){
    for(j in radii){
      for(k in 1:size[1]){
        for(l in 1:size[2]){
          microstates[k,l] <- paste(microstates[k,l],
                                    as.character(max(0,floor(sum(physSpace[[i]][(maxr+k-j):(maxr+k+j),(maxr+l-j):(maxr+l+j)]*nbin/(1+2*j)^2)-0.0001))),
                                    sep="")
        }
      }
    }
  }
  dataframe <- data.frame(states=as.vector(microstates),count=rep(1,size[1]*size[2]))
  states <- aggregate(count ~ states, data=dataframe, FUN=sum)
  return(sum(states$count/(size[1]*size[2])*log(size[1]*size[2]/states$count))/logb)
}

entropypct <- function(physSpace,radii,nbin=3,logb=1){
  # physSpace : list of matrices at finest resolution representing data for each category (normalized to be between 0 and 1)
  # radii : vector of radii representing scales to be tested by increasing order
  # nbin : number of binning intervals for each category's characteristics
  # logb : base for the logarithm
  size <- c(nrow(physSpace[[1]]),ncol(physSpace[[1]]))
  microstates <- matrix("s",nrow=size[1],ncol=size[2])
  maxr <- max(radii)
  result <- matrix(0,nrow=size[1],ncol=size[2])
  for(i in 1:length(physSpace)){
    physSpace[[i]] <- rbind(physSpace[[i]][(maxr+1):2,],physSpace[[i]],physSpace[[i]][(size[1]-1):(size[1]-maxr),])
    physSpace[[i]] <- cbind(physSpace[[i]][,(maxr+1):2],physSpace[[i]],physSpace[[i]][,(size[2]-1):(size[2]-maxr)])
  }
  for(i in 1:length(physSpace)){
    for(j in radii){
      for(k in 1:size[1]){
        for(l in 1:size[2]){
          microstates[k,l] <- paste(microstates[k,l],
                                    as.character(max(0,floor(sum(physSpace[[i]][(maxr+k-j):(maxr+k+j),(maxr+l-j):(maxr+l+j)]*nbin/(1+2*j)^2)-0.0001))),
                                    sep="")
        }
      }
    }
  }
  dataframe <- data.frame(states=as.vector(microstates),count=rep(1,size[1]*size[2]))
  states <- aggregate(count ~ states, data=dataframe, FUN=sum)
  for(i in 1:length(states$count)){
    result[which(microstates == states$states[i])] <- 1/(size[1]*size[2])*log((size[1]*size[2])/states$count[i])
  }
  return(result)
}


#########################################
##### Classical inequality measures #####
#########################################


IncSegClemRear <- function(X,Y,Z,MX,MY,D){
  # D = degree of interpolating polynomial
  l <- length(Z)
  mx <- min(X)
  Mx <- max(X)
  my <- min(Y)
  My <- max(Y)
  intx <- seq(mx,Mx-(Mx-mx)/MX,by=(Mx-mx)/MX)
  inty <- seq(my,My-(My-my)/MY,by=(My-my)/MY)
  J <- NULL
  for(i in 1:MX){
    for(j in 1:MY){
      J <- c(J,paste(i,j,sep="."))
    }
  }
  j <- length(J)
  Zb <- data.frame(Z=Z,X=X,Y=Y,etiquette=rep(0,length(Z)),count=rep(1,l))
  Zb <- Zb[order(Z),]
  row.names(Zb) <- as.character(1:l)
  for(i in 1:length(Z)){
    X0 <- max(which(intx <= Zb$X[i]))
    Y0 <- max(which(inty <= Zb$Y[i]))
    Zb$etiquette[i] <- paste(X0,Y0,sep=".")
  }
  Zagg <- aggregate(Zb$count,by=list(Zb$Z),FUN="sum")
  k <- nrow(Zagg)
  Zagg$y <- rep(0,k)
  for(i in 1:k){
    Zagg$y[i] <- sum(Zagg$x[1:i])/l
  }
  ctable <- matrix(0,nrow=j+1,ncol=k)
  ctable2 <- matrix(0,nrow=j+1,ncol=k)
  ctable3 <- matrix(0,nrow=j+1,ncol=k)
  for(p in 1:k){
    for(q in 1:j){
      pp <- Zagg$Group.1[p]
      ctable[q,p] <- sum(Zb$count[which(Zb$Z == pp & Zb$etiquette == J[q])])
    }
    ctable[j+1,p] <- sum(ctable[1:j,p])
  }
  hn <- rep(0,j+1)
  for(q in 1:(j+1)){
    hn[q] <- sum(ctable[q,])
    if(hn[q] > 0){
      for(p in 1:k){
        ctable2[q,p] <- sum(ctable[q,1:p])/hn[q]
        if(ctable2[q,p] > 0 & ctable2[q,p] < 1){
          ctable3[q,p] <- -ctable2[q,p]*log(ctable2[q,p])/log(2)+(ctable2[q,p]-1)*log(1-ctable2[q,p])/log(2)
        }
        ctable2[q,p] <- 4*ctable2[q,p]*(1-ctable2[q,p])
      }
    }
  }
  inter2 <- rep(0,k-1)
  inter3 <- rep(0,k-1)
  for(p in 1:(k-1)){
    inter2[p] <- sum((ctable2[j+1,p]-ctable2[1:j,p])*hn[1:j])/length(Z)/ctable2[j+1,p]
    inter3[p] <- sum((ctable3[j+1,p]-ctable3[1:j,p])*hn[1:j])/length(Z)/ctable3[j+1,p]
  }
  polyfit2 <- lm(inter2 ~ poly(Zagg$y[1:(k-1)],D,raw=T))
  polyfit3 <- lm(inter3 ~ poly(Zagg$y[1:(k-1)],D,raw=T))
  Rcoef2 <- 6/(0:D+2)/(0:D+3)
  Rcoef3 <- 1
  for(i in 1:D){
    Rcoef3 <- c(Rcoef3,2/(i+2)^2+2*sum((-1)^(i:0)*choose(i,0:i)/(i:0+2)^2))
  }
  return(c(sum(polyfit3$coefficients*Rcoef3),sum(polyfit2$coefficients*Rcoef2)))
}

classIneqMeas <- function(X,Y,Z,MX,MY,D)
{
  avg <- mean(Z)
  l <- length(Z)
  #
  MM <- sum(abs(avg-Z))/(l*avg)
  GG <- ineq(Z,type="Gini")
  TT <- sum(Z[which(Z>0)]/l*log(Z[which(Z>0)]))
  #
  spDistVal <- newrandtoGrid(X,Y,Z,MX,MY,"mean")
  spDistRef <- newrandtoGrid(X,Y,rep(1,l),MX,MY,"sum")
  NSI <- sqrt(sum(spDistRef*(spDistVal-avg)^2)/sum((Z-avg)^2))
  REAR <- IncSegClemRear(X,Y,Z,MX,MY,D)
  return(c(MM,GG,TT,NSI,REAR[1],REAR[2]))
}

################
##### Misc #####
################

mult_centered_distrib <- function(dist,nrow,ncol,xcoords,ycoords,power=1,weights=rep(1,length(xcoords)))
  # nrow, ncol = total space size
  # xcoords, ycoords = coordinates of centers of attraction
  # /!\ dist size = nrow*ncol /!\
{
  M <- matrix(0,ncol=ncol,nrow=nrow)
  for(s in 1:length(xcoords)){
    for(i in 1:nrow){
      for(j in 1:ncol){
        M[i,j] <- M[i,j]+weights[s]/sqrt((i-xcoords[s])^2+(j-ycoords[s])^2)
      }
    }
  }
  Mrank <- matrix(rank(M,ties.method="first"),ncol=ncol,nrow=nrow)
  dist <- sort(dist,decreasing=T)
  dist <- sample(dist,length(dist),prob=rev(dist^power)/sum(dist^power))
  for(l in 1:length(dist)){
    M[which(Mrank == l)] <- dist[l]
  }
  return(M)
}




