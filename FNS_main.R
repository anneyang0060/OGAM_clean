# online gam
ogam <- function(K, X, y, n, m, delta_inner, delta_outer, Max_iter, band_select, K_band, C1=NULL, L=1, L_theta=0, L_sigma=0, pd1=1, pd2=3){
  
  eval_vec <<- seq(0.1, 0.9, length.out = m)
  d <- ncol(X)
  G <- rep(1,d)
  
  if(K==1){
    
    ## some initial values
    N <<- 0; h <<- rep(1,d)
    ## stored statistics for the main regression
    U_ <<- array(0, dim = c(pd1*d+1,m^d,L))
    V_ <<- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L)) 
    Q_ <<- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L))
    if(link!='identity'){
      R_ <<- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L))
    }else{
      R_ <<- 0
    }
    ## stored statistics for bandwidth selection
    if(band_select==TRUE){
      U_theta <<- array(0, dim = c(pd2*d+1,m^d,L_theta))
      V_theta <<- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta))
      Q_theta <<- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta))
      U_sigma <<- array(0, dim = c(pd1*d+1,m^d,L_sigma))
      V_sigma <<- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma))
      Q_sigma <<- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma))
      R_sigma <<- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L_sigma))
      R_theta <<- array(0, dim = c(pd2*d+1,pd2*d+1,pd2*d+1,m^d,L_theta))
      if(link!='identity'){
        R_sigma <<- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L_sigma))
        R_theta <<- array(0, dim = c(pd2*d+1,pd2*d+1,pd2*d+1,m^d,L_theta))
      }else{
        R_sigma <<- 0; R_theta <<- 0
      }
    }
    
  }  
  
  N <<- N + n
  if(K > K_band){band_select <- FALSE}
  if(K <= K_band & band_select==FALSE){Ch <<- C1}
  
  # bandwidth selection
  if(band_select==TRUE){
    
    h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
    eta_theta <- sapply(1:L_theta, function(l){((L_theta-l+1)/L_theta)^(1/7) * h_theta})#dim: d*L
    eta_sigma <- sapply(1:L_sigma, function(l){((L_sigma-l+1)/L_sigma)^(1/5) * h_sigma})
    
    if(K==1){
      
      # generate initial estimate and bandwidth
      initial_res <- initial(X,y,h_theta,eval_vec,pd2)
      beta0_theta_old <- initial_res$beta0; beta_theta_old <- initial_res$beta
      beta0_theta <- initial_res$beta0; beta_theta <- initial_res$beta
      idx_theta <- 1:L_theta
      centrds_theta <<- eta_theta
      rm(initial_res)
      
      initial_res <- initial(X,y,h_sigma,eval_vec,pd1)
      beta0_sigma_old <- initial_res$beta0; beta_sigma_old <- initial_res$beta
      beta0_sigma <- initial_res$beta0; beta_sigma <- initial_res$beta
      idx_sigma <- 1:L_sigma
      centrds_sigma <<- eta_sigma
      rm(initial_res)
      
    }else{
      
      # update index for pseudo-bandwidth and the centroids
      idx_theta <- sapply(1:L_theta,function(l){
        which.min(abs(eta_theta[1,l] - centrds_theta[1,]))
      })# idx of all d dimensions are the same, it is sufficient to compuute the first dim
      for(i in 1:d){
        centrds_theta[i,] <<- (centrds_theta[i,idx_theta] * (N-n) + eta_theta[i,] * n) / N
      }
      
      idx_sigma <- sapply(1:L_sigma,function(l){
        which.min(abs(eta_sigma[1,l] - centrds_sigma[1,]))
      })
      for(i in 1:d){
        centrds_sigma[i,] <<- (centrds_sigma[i,idx_sigma] * (N-n) + eta_sigma[i,] * n) / N
      }
    }
    
    # backfitting
    {
      # tt0 <- Sys.time()
      res_beta <- backfitting(beta0_theta, beta_theta,U_theta[,,idx_theta[1]],V_theta[,,,idx_theta[1]],
                              Q_theta[,,,idx_theta[1]],R_theta[,,,,idx_theta[1]], n, N, h_theta, pd2)
      if(K==1){rm(beta0_theta, beta_theta)}
      if(res_beta$delta < delta_outer){
        beta0_theta <<- res_beta$beta0
        beta_theta <<- res_beta$beta
      }
      rm(res_beta)
      # tt1 <- Sys.time()
      # print(paste0('theta backfitting:',round(difftime(tt1,tt0,units = 'mins'),3),'mins'))
      
      # tt0 <- Sys.time()
      res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma[,,idx_sigma[1]],V_sigma[,,,idx_sigma[1]],
                              Q_sigma[,,,idx_sigma[1]],R_sigma[,,,,idx_sigma[1]], n, N, h_sigma, pd1)
      if(K==1){rm(beta0_sigma, beta_sigma)}
      if(res_beta$delta < delta_outer){
        beta0_sigma <<- res_beta$beta0
        beta_sigma <<- res_beta$beta
      }
      rm(res_beta)
      # tt1 <- Sys.time()
      # print(paste0('sigma backfitting:',round(difftime(tt1,tt0,units = 'mins'),3),'mins'))
    }
    
    # compute the constants
    {
      beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
      for(i in 1:d){
        beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
      }
      if(link=='log'){
        for(i in 1:d){
          sigma2[i] <- mean(exp(beta_hat[,(1:d)[-i]]))*mean(exp(-beta_sigma[,i]))*exp(-beta0_sigma)
        }
      }
      if(link=='identity'){
        sigma2 <- mean((y-rowSums(beta_hat[,(1:d)])-beta0_sigma)^2)
      }
      if(link=='logit'){
        sigma2 <- 4*mean(1/(1-exp(-beta0_sigma-rowSums(beta_hat)))*(1-1/(1+exp(-beta0_sigma-rowSums(beta_hat)))))
      }
      Ch <<- (15 * sigma2/colMeans(beta_theta[,(2*d+1):(3*d)]^2))^(1/5)
    }
    
    # update statistics
    {
      res_update <- update_stats(U_theta,V_theta,Q_theta,R_theta, beta0_theta,beta_theta,eta_theta,
                                 idx_theta,n,N,pd2,L_theta)
      U_theta <<- res_update$U_
      V_theta <<- res_update$V_
      Q_theta <<- res_update$Q_
      R_theta <<- res_update$R_
      rm(res_update)
      
      res_update <- update_stats(U_sigma,V_sigma,Q_sigma,R_sigma, beta0_sigma,beta_sigma,eta_sigma,
                                 idx_sigma,n,N,pd1,L_sigma)
      U_sigma <<- res_update$U_
      V_sigma <<- res_update$V_
      Q_sigma <<- res_update$Q_
      R_sigma <<- res_update$R_
      rm(res_update)
    }
    
  }
  
  # compute bandwidth and candidates
  {
    h <<- sapply(1:d, function(i){min(Ch[i]*N^(-1/5),h[i])})
    eta <- sapply(1:L, function(l){((L-l+1)/L)^(1/5) * h}) #dim = d*L
  }
  
  
  # initial parametric estimate and combination rule
  if(K==1){
    
    initial_res <- initial(X,y,h,eval_vec,pd1)
    beta0_est <<- initial_res$beta0; beta_est <<- initial_res$beta
    
    idx <- 1:L
    centrds <<- eta
    
    par(mfrow=c(2,2))
    plot(beta_est[,1],type='l',main=paste0('K=',0)); plot(beta_est[,2],type='l')
    plot(beta_est[,3], type='l'); plot(beta_est[,4],type='l')
    
  }else{
    
    idx<-sapply(1:L,function(l){
      which.min(abs(eta[1,l] - centrds[1,]))
    })
    for(i in 1:d){
      centrds[i,] <- (centrds[i,idx] * (N-n) + eta[i,] * n) / N
    }
    centrds <<- centrds
    
  }
  
  # backfitting
  {
    res_beta <- backfitting(beta0_est, beta_est,U_[,,idx[1]],V_[,,,idx[1]],Q_[,,,idx[1]],R_[,,,,idx[1]], n, N, h, pd1)
    if(res_beta$delta < delta_outer){
      beta0_est <<- res_beta$beta0; beta_est <<- res_beta$beta
      print(paste0('backfitting modified : delta=', res_beta$delta))
    }else{
      print('backfitting modified failed')
    }
    rm(res_beta)
  }
  # print('backfitting OK')
  
  # update statistics
  {
    res_update <- update_stats(U_,V_,Q_,R_, beta0_est, beta_est, eta, idx, n, N, pd1, L)
    U_ <<- res_update$U_
    V_ <<- res_update$V_
    Q_ <<- res_update$Q_
    R_ <<- res_update$R_
    rm(res_update)
  }
  # print('update statistics OK')
  
}
