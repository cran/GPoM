#' @title Periodic solution test
#'
#' @description Tests if a trajectory is periodic.
#'
#' @inheritParams  gPoMo
#' @inheritParams  autoGPoMoSearch
#'
#' @param wthresh Threshold used to detect the limit cycle.
#' @param fxPtThresh Threshold used to detect fixed points.
#'
#' @author Sylvain Mangiarotti, Flavie Le Jean
#'
#' @return \code{periodic}  An integer classifying the models:
#' diverging or unclassified trajectory (0),
#' period-1 trajectory (-1), period-2 trajectory (-2)
#' and fixed Point (2).
#'
#' @seealso \code{\link{autoGPoMoTest}}, \code{\link{gPoMo}}
#'
#' @examples
#' #Example
#' # Load data:
#' data('P1FxChP2')
#' # Test a period-1 trajectory
#' testP(P1FxChP2[,1:2], wthresh=0.1, fxPtThresh = 1e-6, show=0)
#' # Test a Fixed Point trajectory
#' testP(P1FxChP2[,3:4], wthresh=0.1, fxPtThresh = 1e-6, show=0)
#' # Test a chaotic trajectory
#' testP(P1FxChP2[,5:6], wthresh=0.1, fxPtThresh = 1e-6, show=0)
#' # Test a period-2 trajectory
#' testP(P1FxChP2[,7:8], wthresh=0.1, fxPtThresh = 1e-6, show=0)
#'
#' @export
testP <- function(data, wthresh=0.1, fxPtThresh=1E-4, show=0) {
  
  
  lg=length(data[,2])
  periodic = list()
  periodic$status <- 1
  if(is.finite(data[lg,2]) == FALSE  ) {
    periodic$status <- 0
    return(periodic)
  }
  
  if(lg < 500  ) {periodic$status <- 0}    #### if iter<500 don't test it yet
  
  if((lg >= 500)  && (lg < 1000)  ) {    #### only test p1
    
    vx=data[,1] 
    vy=data[,2]
    vx=vx[length(vx):1]           ##### inversion of the model to start from the end
    vy=vy[length(vy):1]
    
    px=mean(vx)        
    py=mean(vy)
    pc=c(px,py)                   ##### 'center point' of the model
    
    loop1=c()
    loop2=c()
    theta1=c()
    theta2=c()
    
    loops=0
    for(ii in 1:length(vx) ){
      pv=c(vx[ii],vy[ii])       ####### go throught the model 
      
      
      if (  (theta(pc, ii, vx, vy) >= 0 ) && (ii>1) && (theta(pc, (ii-1), vx, vy) < 0 ) )  {   ### detect when the last loop starts
        
        l1=1
        pp1=c(vx[ii],vy[ii])
        for(i1 in ii: length(vx)){
          pv=c(vx[i1],vy[i1])
          loop1[l1]=sqrt(sum((pv - pc) ^ 2))  ### distance from the point to the center
          theta1[l1]=theta(pc, i1, vx, vy)
          l1=l1+1
          
          if( (l1>2)  && (theta(pc, i1, vx, vy) >= 0 )  &&  (theta(pc, (i1-1), vx, vy) < 0 ) ) {loops=1; break}  ### ends loop1 starts loop2
        }
        
        l2=1
        pp2=c(vx[i1],vy[i1])
        for (i2 in i1: length(vx) ){  
          pv=c(vx[i2],vy[i2])
          loop2[l2]=sqrt(sum((pv - pc) ^ 2))
          theta2[l2]=theta(pc, i2, vx, vy)
          l2=l2+1
          if( (l2>2) &&  (theta(pc, i2, vx, vy) >= 0 )  &&  (theta(pc, (i2-1), vx, vy) < 0 ) ) {loops=2; break }##### ends loop2
        }
        
        if(  loops==2 ) {break}#### break loop '(ii in 1: length(vx))'
      }
      
    }
    
    lg1=length(loop1)
    lg2=length(loop2)
    
    if ( lg1==0 | lg2==0 ) {
      #### break if any 'initial point (theta=0)' of the loop has not found
      periodic$status <- 0
      return(periodic)
      if (show==1){
        plot(0,0)
        text(0, 0.7, "(TESTED FOR P1)")
        text(0, 0.9, "no cicles found")
        }
      }
    
    n=max(lg1,lg2)
    s1=spline(loop1, y=NULL, n=n, method='fmm')
    s2=spline(loop2, y=NULL, n=n, method='fmm') 
    
    errorp1=sqrt(mean((s2$y-s1$y)^2))
    periodic$errP1 <- errorp1
    lf=c(loop2, loop1)
    threshold=mean(lf)*wthresh
    
    distance=0
    for (ifp in 1: 100){          ### calculate the distance between the last 100 points
      plast=c(vx[ifp],vy[ifp])      # to check it if the model is a fixed point
      plast1=c(vx[ifp+1],vy[ifp+1])
      distance=distance+(sqrt(sum((plast1 - plast) ^ 2)))
      if (distance > fxPtThresh) {break}
    }
    
    if (distance < fxPtThresh) {   ### fixed point
      periodic$status <- 2 
    }
    
    else if (errorp1 < threshold)  {
      periodic$status <- -1
      if (show==1){
        dev.new()
        plot(loop1, type='l')
        lines(loop2, col='red')
        text("TEST FOR P1 \n")
        message("threshold = ", threshold, "\n")
        message("errorp1 = ", errorp1, "\n")
        message("DISTANCE < = ", fxPtThresh, "\n")
      }
    }
  }
  
  if(lg >=1000  ) {
    ###### test p1 or p2
    vx=data[,1] 
    vy=data[,2]
    vx=vx[length(vx):1]           ##### inversion of the model to start from the end
    vy=vy[length(vy):1]
    
    px=mean(vx)        
    py=mean(vy)
    pc=c(px,py)#### 'center point' of the model
    
    loop1=c()
    loop2=c()
    loop3=c()
    loop4=c()
    loops=0
    
    for(ii in 1:length(vx) ){
      pv=c(vx[ii],vy[ii]) ############# go throught the model 

      if (  (theta(pc, ii, vx, vy) >= 0 )
            && (ii>1)
            && (theta(pc, (ii-1), vx, vy) < 0 ) )  {
        ### detect when the last loop starts
        l1=1
        pp1=c(vx[ii],vy[ii])
        for(i1 in ii: length(vx)){
          pv=c(vx[i1],vy[i1])
          loop1[l1]=sqrt(sum((pv - pc) ^ 2))  ### distance from the point to the center
          l1=l1+1
          if( (l1>2)  && (theta(pc, i1, vx, vy) >= 0 )  &&  (theta(pc, (i1-1), vx, vy) < 0 ) ) {loops=1; break}  ### ends loop1 starts loop2
        }
        
        l2=1
        pp2=c(vx[i1],vy[i1])
        for (i2 in i1: length(vx) ){  
          pv=c(vx[i2],vy[i2])
          loop2[l2]=sqrt(sum((pv - pc) ^ 2))
          l2=l2+1
          if( (l2>2) &&  (theta(pc, i2, vx, vy) >= 0 )  &&  (theta(pc, (i2-1), vx, vy) < 0 ) ) {loops=2; break }##### ends loop2
        }
        
        l3=1
        for (i3 in i2: length(vx) ){  
          pv=c(vx[i3],vy[i3])
          loop3[l3]=sqrt(sum((pv - pc) ^ 2))
          l3=l3+1
          if(  (l3>2) &&  (theta(pc, i3, vx, vy) >= 0 )  &&  (theta(pc, (i3-1), vx, vy) < 0 ) ) {loops=3; break }##### ends loop3
        }
        
        l4=1
        for (i4 in i3: length(vx) ){  
          pv=c(vx[i4],vy[i4])
          loop4[l4]=sqrt(sum((pv - pc) ^ 2))
          l4=l4+1
          if(  (l4>2) &&  (theta(pc, i4, vx, vy) >= 0 )  &&  (theta(pc, (i4-1), vx, vy) < 0 ) ) {loops=4; break}##### ends loop4
        }
        
        if(  loops==4 ) {break}#### break loop '(ii in 1: length(vx))'
      }
      
    }
    
    lg1=length(loop1)
    lg2=length(loop2)
    lg3=length(loop3)
    lg4=length(loop4)
    
    if ( lg1==0 | lg2==0 | lg3==0 | lg4==0 )    {#break if any 'initial point (theta=0)' has not found
      periodic$status <- 0
      return(periodic)
      
      if (show==1){
        plot(0,0)
        text(0, 0.7, "(TESTED FOR P2)")
        text(0, 0.9, "no cicles found")
      }} 
    
    n=max(lg1,lg2)
    sp1.1=spline(loop1, y=NULL, n=n, method='fmm')
    sp1.2=spline(loop2, y=NULL, n=n, method='fmm') 
    
    errorp1=sqrt(mean((sp1.2$y-sp1.1$y)^2))
    periodic$errP1 <- errorp1
    
    n=max(lg1,lg3)
    s1=spline(loop1, y=NULL, n=n, method='fmm')
    s3=spline(loop3, y=NULL, n=n, method='fmm') 
    
    n=max(lg2,lg4)
    s2=spline(loop2, y=NULL, n=n, method='fmm')
    s4=spline(loop4, y=NULL, n=n, method='fmm')
    
    error1=sqrt(mean((s3$y-s1$y)^2))      # RMSE between LOOP1 AND LOOP3
    error2=sqrt(mean((s4$y-s2$y)^2))      # RMSE between LOOP2 AND LOOP4
    
    lf=c(loop1, loop2, loop3, loop4)
    threshold=mean(lf)*wthresh
    
    distance=0
    for (ifp in 1: 100){          ### calculate the distance between the last 100 points
      plast=c(vx[ifp],vy[ifp])      # to check if the model is a fixed point
      plast1=c(vx[ifp+1],vy[ifp+1])
      distance=distance+(sqrt(sum((plast1 - plast) ^ 2)))
      if (distance > fxPtThresh) {break}
    }
    
    
    if (distance < fxPtThresh) {   ### fixed point
      periodic$status <- 2
    }  
    
    else if (errorp1 < threshold) {
      periodic$status <- -1
      periodic$errP1 <- errorp1
      
      }
    else if (   (error1 < threshold) && (error2 < threshold)  ) {
      periodic$status <- -2
      periodic$errP2 <- error2
      }
    
    
    if (show==1){
      dev.new()
      plot(loop1, type='b')
      lines(loop2, type='b', col='red')
      lines(loop3, type='b', col='green')
      lines(loop4, type='b', col='blue')
      message("threshold = ", threshold, "\n")
      message("fxPtThresh = ",fxPtThresh, "\n" )
      message("errorp1 = ", errorp1, "\n")
      message("error1 = ", error1, "\n")
      message("error2 = ", error2, "\n")
      message("distance fp (last 100p) = ", distance, "\n")
      if (periodic$status == -1) message("-1: The trajectory is period-1", "\n")
      if (periodic$status == -2) message("-2: The trajectory is period-2", "\n")
      if (periodic$status == 0) message("0: The trajectory is diverging or undetermined", "\n")
      if (periodic$status == 2) message("2: The trajectory is a Fixed Point", "\n")
    }
    
  }  
  return (periodic) 

  
}  

theta <- function (pc=c(2), ii, vx, vy)
# calculates the angle of each point respect the reference line defined by the center point

{ pv=c(vx[ii],vy[ii])

if( (pv[2]-pc[2]) >= 0 ){
  cat=(pc[1]-pv[1])
  hip=sqrt( (pc[1]-pv[1])^2 + (pc[2]-pv[2])^2  )
  theta=(acos(cat/hip))
  if (cat==0 && hip==0) { theta=pi/2}
  
}

else {
  cat=(pv[1]-pc[1])
  hip=sqrt( (pc[1]-pv[1])^2 + (pc[2]-pv[2])^2  )
  theta=-(acos(cat/hip))
  if (cat==0 && hip==0) { theta=-(pi/2)}
}
return(theta)
}
