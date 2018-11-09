pbEqBelt<-function(n, k, interval, start_phi,
                   end_phi, start_theta){

  delta_phi <- interval
  phi <- start_phi-delta_phi
  delta_theta <- interval
  theta <- start_theta-delta_theta
  P <- numeric(0)
  index <- 0
  while (phi < end_phi){
    tmp <- pbNextEqBelt(n, k, phi, theta,
                        delta_phi, delta_theta)
    index <- index+1
    phi <- tmp[[1]]
    theta <- tmp[[2]]
    dTheta_dPhi<-tmp[[3]]
    a_temp <- cbind(theta, phi)
    P<- rbind(P,a_temp)
    dL_dPhi <- sqrt(dTheta_dPhi^2+sin(theta)^2)
    delta_phi <- interval/dL_dPhi
    delta_theta <- dTheta_dPhi*delta_phi
    # Add the end edge point
    if (phi+delta_phi > end_phi) delta_phi <- end_phi-phi
  }
  return(P)
}



pbNextEqBelt <- function(n, k, prev_phi, prev_theta,
                         delta_phi, delta_theta){

  phi <- prev_phi+delta_phi
  a <- 2/3-4*k/(3*n)
  b <- 8/(3*pi)
  t1 <- a+b*phi
  t2 <- a-b*phi

  #With consideration of the complex values
  if ( (a+b*phi) >= -1 && (a+b*phi) <= 1 )
  {
    theta0 <- acos(a+b*phi)
    if ( (a-b*phi) >= -1 && (a-b*phi) <= 1 )
    {
      theta1 <- acos(a-b*phi)
    }
    else
    {
      theta1 <- Im(acos(a-b*phi))
    }
  }
  else
  {
    theta0 <- Im(acos( a + b*phi ))
    if ( (a - b*phi) >= -1 && (a - b*phi) <= 1 ){
      theta1 <- acos(a-b*phi)
    }else{
      theta1 <- Im(acos( a - b*phi ))
    }
  }

  exp_theta <- prev_theta + delta_theta
  if (abs(theta0-exp_theta) < abs( theta1 - exp_theta )) {
    theta <- Re(theta0)
    dTheta_dPhi <- -b/sin(theta)
  } else {
    theta <- Re(theta1)
    dTheta_dPhi <- b/sin(theta)
  }
  return( list(phi, theta, dTheta_dPhi) )
}


############################
pbPolarCap <- function(n, k, interval,
                       start_phi, end_phi, flag){

  delta_phi <- interval
  phi <- start_phi - delta_phi
  P <- numeric(0)

  while (phi < end_phi){
    tmp <- pbNextPolarCap(n, k, phi, delta_phi, flag)
    phi <- tmp[[1]]
    theta <- tmp[[2]]
    dTheta_dPhi <- tmp[[3]]
    a_temp <- cbind(theta, phi)
    P <- rbind(P,a_temp)
    dL_dPhi <- sqrt( dTheta_dPhi^2 + sin(theta)^2 )
    delta_phi<- interval/dL_dPhi

    if (phi+delta_phi > end_phi) delta_phi <- end_phi-phi   #add the end edge point
  }
  return(P)
}


pbNextPolarCap <- function(n, k, prev_phi, delta_phi, flag){

  phi <- prev_phi + delta_phi
  b <- -k^2*pi^2/(12*n^2)

  if (flag==0){
    theta <- acos( 1 + b/phi^2 )           #Eq. 19, here cos(theta)=z;
    dTheta_dPhi <- 2*b/(phi^3*sin(theta))
  }else{
    theta <- acos( 1 + b/( phi-pi/2)^2 )   #Eq. 20, here cos(theta)=z;
    dTheta_dPhi <-  2*b/( (phi-pi/2)^3*sin(theta) )
  }

  return(list(phi,theta,dTheta_dPhi))
}




