################################################################################
library(MASS)
library(nlme)
library(glmmLasso)
library(xtable)


################################################################################
#Global variable of spatial variable selection, contains optimization tuning parameters
control.default=list(maxIter=200,iwls=10^(-4),tol1=10^(-3),tol2=10^(-3))

################################################################################
#SpatialVS is the major function for this paper. It performs variable selection for the spatial Poisson regression
#model under adaptive elastic net penalty.
#dat.obj: List, input data. Must contains:
# X: numeric matrix, the covariates.
# y: integer vector, the response in counts.
# dist: numeric matrix, the distance matrix.
# offset: numeric vector, the offset item.
#alpha.vec: numeric vector, a vector of possible values of regularization parameter. The range is [0,1].
#lambda.vec: numeric vector, a vector of positive values of regularization parameter.
#method: string, the method to be used. Options are:
#"PQL": penalized quasi-likelihood method that considers spatial correlation.
#"PQL.nocor": penalized quasi-likelihood method that ignores spatial correlation.
#"APL": approximate penalized loglikelihood method that considers spatial correlation.
#"APL.nocor": approximate penalized loglikelihood method that ignores spatial correlation.
#plots: bool, if True, contour plot of AIC/BIC values is generated.
#intercept: bool, if True, an intercept item will be included in model.
#verbose: bool, if True, various updates are printed during each iteration of the algorithm.

SpatialVS=function(dat.obj=lyme.svs.eco0.dat, alpha.vec=seq(0.6, 1, by=0.05), lambda.vec=seq(0.15, 1, len=50), method="PQL", plots=F, intercept=T, verbose=T)
{

  if(method=="PQL")
  {
    SpatialVS_P_FUNCTION=SpatialVS_P_testPQL
    SpatialVS.out.FUNCTION=SpatialVS.out
  }

  if(method=="PQL.nocor")
  {
    SpatialVS_P_FUNCTION=SpatialVS_P_testPQL.nocor
    SpatialVS.out.FUNCTION=SpatialVS.out.nocor
  }

  if(method=="APL")
  {
    SpatialVS_P_FUNCTION=SpatialVS_P_testLP
    SpatialVS.out.FUNCTION=SpatialVS.out
  }

  if(method=="APL.nocor")
  {
    SpatialVS_P_FUNCTION=SpatialVS_P_testLP.nocor
    SpatialVS.out.FUNCTION=SpatialVS.out.nocor
  }

  if(verbose)
  {
    cat("Finding initial values \n")
  }

  start.val.ini<-try(spatialVS_simu_ini(simu.data=dat.obj, intercept=intercept, verbose=verbose))


  if(class(start.val.ini)=="try-error")
  {
    res=NULL
    return(res)
  }

  if(verbose)
  {
    cat("ini.beta:", start.val.ini$beta, "ini.theta:", start.val.ini$theta, "\n")
  }


  ###no adapative, i.e., EN

  if(verbose)
  {
     cat("Starting Major looping for alpha and lambda with EN. \n")
  }

  L.EN.obj<-try(SpatialVS_P_FUNCTION(data=dat.obj, start=start.val.ini, alpha.vec=alpha.vec, lambda.vec=lambda.vec, cc=0, adaptive=F, intercept=intercept, verbose=verbose), silent=T)

  if(class(L.EN.obj)=="try-error")
  {
    res=NULL
    return(res)
  }


  Lout.EN.obj<-try(SpatialVS.out.FUNCTION(obj=L.EN.obj, data=dat.obj, cc=0), silent=T)

  if(class(Lout.EN.obj)=="try-error")
  {
    res=NULL
    return(res)
  }

  ##########select best tuning parameters
  contour.out.EN.obj=SpatialVS.obj.contour.plot(SpatialVS.obj=L.EN.obj, SpatialVS.out.obj=Lout.EN.obj, type="BIC", plots=plots)

  alpha.EN.best.val=contour.out.EN.obj$alpha.best.val
  lambda.EN.best.val=contour.out.EN.obj$lambda.best.val

  if(verbose)
  {
    cat("Computing best alpha and lambda with EN. \n")
  }

  ###
  L.EN.best.obj=SpatialVS_P_FUNCTION(data=dat.obj, start=start.val.ini, alpha.vec=alpha.EN.best.val, lambda.vec=lambda.EN.best.val, cc=0, adaptive=F, intercept=intercept, verbose=verbose)

  lasso.weight=as.vector(L.EN.best.obj$beta)
  nn=nrow(dat.obj$X)
  lasso.weight[lasso.weight==0]=(1/nn)

  if(verbose)
  {
    cat("start.ini beta:", round(start.val.ini$beta, 4),"\n")
    cat("EN lasso weight:", round(lasso.weight, 4),"\n")
  }

  if(verbose)
  {
    cat("Starting Major looping for alpha and lambda with AEN. \n")
  }

  L.obj<-try(SpatialVS_P_FUNCTION(data=dat.obj, start=start.val.ini, alpha.vec=alpha.vec, lambda.vec=lambda.vec, cc=0, adaptive=T, lasso.weight=lasso.weight, intercept=intercept, verbose=verbose), silent=T)

  if(class(L.obj)=="try-error")
  {
    res=NULL
    return(res)
  }


  if(verbose)
  {
      cat("Computing best alpha and lambda AEN. \n")
  }

  Lout.obj<-try(SpatialVS.out.FUNCTION(obj=L.obj, data=dat.obj, cc=0), silent=T)

  if(class(Lout.obj)=="try-error")
  {
    res=NULL
    return(res)
  }

  ##########select best tuning parameters
  contour.out.obj=SpatialVS.obj.contour.plot(SpatialVS.obj=L.obj, SpatialVS.out.obj=Lout.obj, type="BIC", plots=plots)

  alpha.best.val=contour.out.obj$alpha.best.val
  lambda.best.val=contour.out.obj$lambda.best.val

  #########refit the best model
  L.best.obj=SpatialVS_P_FUNCTION(data=dat.obj, start=start.val.ini, alpha.vec=alpha.best.val, lambda.vec=lambda.best.val, cc=0, adaptive=T, lasso.weight=lasso.weight, intercept=intercept, verbose=verbose)

  Lout.best.obj=SpatialVS.out(obj=L.best.obj, data=dat.obj, cc=0)

  res=list(dat.obj=dat.obj, start=start.val.ini, L.obj=L.obj, Lout.obj=Lout.obj, contour.out.obj=contour.out.obj, L.best.obj=L.best.obj, Lout.best.obj=Lout.best.obj, L.EN.obj=L.EN.obj, Lout.EN.obj=Lout.EN.obj, contour.out.EN.obj=contour.out.EN.obj, L.EN.best.obj=L.EN.best.obj, lasso.weight=lasso.weight, method=method)

  return(res)
}

################################################################################
#return the summarized parameter estimates from the SpatialVS procedure.
#obj: the output of SpatialVS
SpatialVS.summary=function(obj)
{
  dat.obj=obj$data.obj
  L.best.obj=obj$L.best.obj
  Lout.best.obj=obj$Lout.best.obj

  res=c(L.best.obj$beta.res, L.best.obj$theta.res)

  return(res)
}

################################################################################
#Performs penalized quasi-likelihood method that considers spatial correlation.
#a sub-routine used by the major function SpatialVS.
#arguments are similar to the SpatialVS
SpatialVS_P_testPQL=function(data, start=NULL, control=control.default, alpha.vec, lambda.vec, cc, cov.fun.type="exact", adaptive=F, lasso.weight=NULL, intercept=T, verbose=T)
{
  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist

  n=length(y)

  beta.res=theta.res=b.res=alpha.res=lambda.res=NULL

  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }

  ## the covariance fun
  switch(cov.fun.type,
         "exact"={
           cov.fun=covobj_PQLpenb
         })

  for(i in 1:length(alpha.vec))
  {
    for(j in 1:length(lambda.vec))
    {

      if(verbose)
      {
        cat("i=", i, "alpha=", alpha.vec[i], "j=", j, "lambda=", lambda.vec[j], "\n")
      }

      nIter=0

      beta=start$beta
      theta=start$theta
      D=covStruct(theta, cc=cc, dist=dist)
      convCov=convPar=1



      p=length(beta)

      if(is.null(lasso.weight))
      {
        lasso.weight=rep(1, p)
      }


      b=start$b

      ###main event
      while ((control$maxIter > nIter) & (convCov > control$tol1| convPar > control$tol2 ) )
      {


        #cat("nIter=", nIter, "convCov=", convCov, "\n")

        nIter_beta=1
        convPar=1

        while((nIter_beta<2) & (convPar > control$tol2))   #nIter_beta<30
        {
          ###first step: get the current estimates of beta and b using iwls
          res1=iwls(X,y,D,offset,beta,b,tol=control$iwls,dist=dist)

          ####second step: the elastic-net step
          beta.new=res1$beta;
          mu.new=res1$mu;
          b.new=res1$b

          #print(mean(b.new))

          Z=X%*%beta.new-rep(1,n)+y/mu.new

          if(adaptive)
          {

            pf1.wts=1/(abs(lasso.weight))
            pf2.wts=rep(1,p)

            if(intercept)
            {
              pf1.wts[1]=0
              pf2.wts[1]=0
            }

            res2<-weighted.ada.enet(xmat=X, y=Z, obs.wts=as.vector(mu.new), alpha=alpha.vec[i], lambda=lambda.vec[j], pf1=pf1.wts, pf2=pf2.wts)
            beta.new=res2$coef

          }else{
                  pf1.wts = rep(1,p)
                  pf2.wts=rep(1,p)

                  if(intercept)
                  {
                     pf1.wts[1]=0
                     pf2.wts[1]=0
                  }

                  res2<-weighted.ada.enet(xmat=X, y=Z, obs.wts=as.vector(mu.new), alpha=alpha.vec[i], lambda=lambda.vec[j], pf1=pf1.wts, pf2=pf2.wts)

                  beta.new=res2$coef

                }


          convPar <- sum(sqrt(crossprod(beta.new-beta))/(1+sqrt(crossprod(beta.new))))

          #cat("        nIter_beta=", nIter_beta, "convPar=", convPar, "\n")

          nIter_beta=nIter_beta+1

          #print(convPar)
          #print(beta.new)

          beta=beta.new
          b=b.new


        }#end of beta iteration

        ###third step: update the theta
        
        mu=as.vector(exp(X%*%beta + b + offset))
        ystar=as.vector(X%*%beta+b+(y-mu)/mu)


        nonzerobeta=beta[beta!=0]
        X_L=as.matrix(X[, beta!=0])

        if(adaptive)
        {

          nonzerolasso.weight=abs(lasso.weight)[beta!=0]

          Sigma_diag=alpha.vec[i]*lambda.vec[j]/(abs(nonzerobeta)*abs(nonzerolasso.weight))+2*lambda.vec[j]*(1-alpha.vec[i])


          if(intercept & beta[1]!=0)
          {
             Sigma_diag[1]=0
          }



        }else{
               Sigma_diag=alpha.vec[i]*lambda.vec[j]/(abs(nonzerobeta))+2*lambda.vec[j]*(1-alpha.vec[i])

               if(intercept & beta[1]!=0)
               {
                 Sigma_diag[1]=0
               }


             }

        Sigma_pen=diag(Sigma_diag, nrow=length(Sigma_diag), ncol=length(Sigma_diag))


        system.time(res3<-optim(theta, cov.fun, mu=mu, nonzerobeta=nonzerobeta, n=n, X_L=X_L, Sigma_pen=Sigma_pen,  dist=dist, lower=-5, upper=10, method="L-BFGS-B", cc=cc, offset=offset, X=X, y=y, b.start=b, tol=control.default$iwls, beta=beta))

        theta.new=res3$par

        ###covergence check
        convCov <- sum(sqrt(crossprod(theta.new-theta))/(1+sqrt(crossprod(theta.new))))

        nIter=nIter+1

        #print(paste("Number of iterations:", nIter))

        ###update
        theta=theta.new;
        D=covStruct(theta,cc=cc,dist=dist)

      }#end of main iteration


      res_b=iwls_b(X=X,y=y,D=D,offset=offset,beta=beta,b=b,tol=control$iwls,dist=dist)
      b=res_b$b

      beta.res=rbind(beta.res, as.vector(beta))
      b.res=rbind(b.res, t(as.matrix(b))) # change b to t(b)
      theta.res=rbind(theta.res, as.vector(theta))
      alpha.res=cbind(alpha.res, alpha.vec[i])
      lambda.res=cbind(lambda.res, lambda.vec[j])

      if(verbose)
      {
       print(beta)
       print(theta)
       print(nIter)
      }

    }
  }

  ori.theta.res=cbind(exp(theta.res[,1]), exp(theta.res[,2]))

  return(list(beta.res=beta.res, b.res=b.res, theta.res=theta.res, ori.theta.res=ori.theta.res, alpha.res=alpha.res, lambda.res=lambda.res, alpha.vec=alpha.vec, lambda.vec=lambda.vec))

}


################################################################################
#Performs penalized quasi-likelihood method that ignores spatial correlation.
#a sub-routine used by the major function SpatialVS.
#arguments are similar to the SpatialVS
SpatialVS_P_testPQL.nocor=function(data, start=NULL, control=control.default, alpha.vec, lambda.vec, cc, cov.fun.type="exact", adaptive=F, lasso.weight=NULL, intercept=T, verbose=T)
{
  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist
  
  n=length(y)
  
  beta.res=theta.res=b.res=alpha.res=lambda.res=NULL
  
  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }
  
  ## the covariance fun
  switch(cov.fun.type,
         "exact"={
           cov.fun=covobj_PQLpen.nocor
         })
  
  for(i in 1:length(alpha.vec))
  {
    for(j in 1:length(lambda.vec))
    {
      if(verbose)
      {
        cat("i=", i, "alpha=", alpha.vec[i], "j=", j, "lambda=", lambda.vec[j], "\n")
      }
      
      
      nIter=0
      
      beta=start$beta    
      theta=start$theta[1]     
      D=covStruct.nocor(theta, n=n)      
      convCov=convPar=1
      p=length(beta)
      
      if(is.null(lasso.weight))
      {
        lasso.weight=rep(1, p)
      }
      
      b=start$b
      
      ###main event
      while ( (control$maxIter > nIter) & (convCov > control$tol1| convPar > control$tol2 ) )
      {
        

        convPar=1
        while(convPar > control$tol2)
        {
          ###first step: get the current estimates of beta and b using iwls
          res1=iwls(X,y,D,offset,beta,b,tol=control$iwls,dist=dist)
          
          ####second step: the elastic-net step
          beta.new=res1$beta 
          mu.new=res1$mu 
          b.new=res1$b
          
          Z=X%*%beta.new-rep(1,n)+y/mu.new

          if(adaptive)
          {
            # adaptive elastic-net
            
            pf1.wts=1/(abs(lasso.weight))
            pf2.wts=rep(1,p)
            
            if(intercept)
            {
              pf1.wts[1]=0
              pf2.wts[1]=0 
            }
            
            res2<-weighted.ada.enet(xmat=X, y=Z, obs.wts=as.vector(mu.new), alpha=alpha.vec[i], lambda=lambda.vec[j], pf1=pf1.wts, pf2=pf2.wts)
            
            beta.new=res2$coef
          
          }else{
          
                 pf1.wts = rep(1,p)
                 pf2.wts=rep(1,p)
            
                 if(intercept)
                 {
                   pf1.wts[1]=0
                   pf2.wts[1]=0 
                 }
            
                 res2<-weighted.ada.enet(xmat=X, y=Z, obs.wts=as.vector(mu.new), alpha=alpha.vec[i], lambda=lambda.vec[j], pf1=pf1.wts, pf2=pf2.wts)
                  
                  beta.new=res2$coef
          
          }  
          
          
          convPar <- sum(sqrt(crossprod(beta.new-beta))/(1+sqrt(crossprod(beta.new))))
          
          beta=beta.new 
          b=b.new
          ###third step: update the theta
        }
        
        
        mu=as.vector(exp(X%*%beta + b + offset))        
        ystar=as.vector(X%*%beta+b+(y-mu)/mu)
        
        nonzerobeta=beta[beta!=0]    
        X_L=X[, beta!=0]
        
        if(adaptive)
        {
          nonzerolasso.weight=abs(lasso.weight)[beta!=0]
          Sigma_diag=alpha.vec[i]*lambda.vec[j]/(abs(nonzerobeta)*abs(nonzerolasso.weight))+2*lambda.vec[j]*(1-alpha.vec[i])
                    
          if(intercept & beta[1]!=0)
          {
             Sigma_diag[1]=0
          }
          
        }else{
               Sigma_diag=alpha.vec[i]*lambda.vec[j]/(abs(nonzerobeta))+2*lambda.vec[j]*(1-alpha.vec[i])   
               
               if(intercept & beta[1]!=0)
               {
                 Sigma_diag[1]=0
               }           
               
             }
                
        Sigma_pen=diag(Sigma_diag, nrow=length(Sigma_diag), ncol=length(Sigma_diag))          

        res3=optim(theta, cov.fun, mu=mu, nonzerobeta=nonzerobeta, n=n, X_L=X_L, Sigma_pen=Sigma_pen, ystar=ystar, dist=dist, lower=-5, upper=10, method="L-BFGS-B", cc=cc) 
        
        theta.new=res3$par
        
        ###covergence check
        convCov <- sum(sqrt(crossprod(theta.new-theta))/(1+sqrt(crossprod(theta.new))))
        
        nIter=nIter+1
        
        #print(paste("Number of iterations:", nIter))
        
        ###update
        theta=theta.new; 
        D=covStruct.nocor(theta,n=n)
      }
      
      
      res_b=iwls_b(X=X,y=y,D=D,offset=offset,beta=beta,b=b,tol=control$iwls,dist=dist)
      b=res_b$b
      
      beta.res=rbind(beta.res, as.vector(beta))
      b.res=rbind(b.res, t(as.matrix(b))) # change b to t(b)
      theta.res=rbind(theta.res, as.vector(theta))
      alpha.res=cbind(alpha.res, alpha.vec[i])
      lambda.res=cbind(lambda.res, lambda.vec[j])
      
      if(verbose)
      {
        print(beta)
        print(theta)
      }
      
    }
  }
  
  ori.theta.res=exp(theta.res)
  return(list(beta.res=beta.res, b.res=b.res, theta.res=theta.res, ori.theta.res=ori.theta.res,
              alpha.res=alpha.res, lambda.res=lambda.res, alpha.vec=alpha.vec, lambda.vec=lambda.vec)) 
}

################################################################################
#Performs approximate penalized loglikelihood method that considers spatial correlation.
#a sub-routine used by the major function SpatialVS.
#arguments are similar to the SpatialVS
SpatialVS_P_testLP=function(data, start=NULL, control=control.default, alpha.vec, lambda.vec, cc, cov.fun.type="exact", adaptive=F, lasso.weight=NULL, intercept=T, verbose=T)
{

  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist
  

  
  n=length(y)
  
  beta.res=theta.res=b.res=alpha.res=lambda.res=NULL
  
  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }
  
  ## the covariance fun
  switch(cov.fun.type,
         "exact"={
           cov.fun=covobj_LP
         })
  
  for(i in 1:length(alpha.vec))
  {
    for(j in 1:length(lambda.vec))
    {
 
      if(verbose)
      {
        cat("i=", i, "alpha=", alpha.vec[i], "j=", j, "lambda=", lambda.vec[j], "\n")
      }      
      
      nIter=0
      
      beta=start$beta    
      theta=start$theta     
      D=covStruct(theta,cc=cc, dist=dist)      
      convCov=convPar=1
      p=length(beta)
      
      if(is.null(lasso.weight))
      {
        lasso.weight=rep(1, p)
      }
      
      pf1.wts=1/abs(lasso.weight)
      
      if(intercept)
      {
        pf1.wts[1]=0
      }

      b=start$b
      
      ###main event
      while ( (control$maxIter > nIter) & (convCov > control$tol1| convPar > control$tol2))
      {
        

        
        nIter_beta=1
        convPar=1
        
        # estimate beta given theta
        while( (convPar > control$tol2)&(20 > nIter_beta) )
        {
          nIter_beta=nIter_beta+1
          
          beta.ini=beta
          
          for(s in 1:p)
          {
            e_s=rep(0, p)
            e_s[s]=1
            
            res1=iwls_b_new(X=X,y=y,D=D,offset=offset,beta=beta,b=b,tol=control$iwls,dist=dist)
            
            beta=res1$beta 
            mu=res1$mu 
            b=res1$b 
            eta=log(mu)
            
            x_s=as.numeric(X[,s])
            DW=D%*%diag(mu)
            
            if(intercept & s==1)
            {              
               beta_sm=sum(x_s*(mu-y))+0.5*sum(diag(solve(DW+diag(rep(1, n)))%*%DW%*%diag(x_s)))           
               h_sm=sum(x_s*x_s*mu)    
            
            }else{
                   beta_sm=sum(x_s*(mu-y))+2*lambda.vec[j]*(1-alpha.vec[i])*beta[s]+0.5*sum(diag(solve(DW+diag(rep(1, n)))%*%DW%*%diag(x_s)))           
                   h_sm=sum(x_s*x_s*mu)+2*lambda.vec[j]*(1-alpha.vec[i])            
                 }

            if(intercept & s==1)
            {              
               d_sm=-beta_sm/h_sm
            
               
            }else{
                    w_s=pf1.wts[s]
                    d_sm=median(c((w_s*lambda.vec[j]*alpha.vec[i]-beta_sm)/h_sm, -beta[s], (-w_s*lambda.vec[j]*alpha.vec[i]-beta_sm)/h_sm))
                 }
            
            beta=beta+d_sm*e_s
            
            
            
          }
          
          beta.new=beta          
          convPar <- sum(sqrt(crossprod(beta.new-beta.ini))/(1+sqrt(crossprod(beta.new))))
          #print(convPar)
 
        }

        # estimate covariance parameters
        res3=optim(theta, cov.fun, beta=beta.new, X=X, y=y, cc=cc, b.start=b, offset=offset,  dist=dist, lower=-5, upper=10, method="L-BFGS-B")
        
        theta.new=res3$par
        
        ###covergence check
        convCov <- sum(sqrt(crossprod(theta.new-theta))/(1+sqrt(crossprod(theta.new))))

        #print(convCov)
      
        nIter=nIter+1
        
        #print(paste("Number of iterations:", nIter))
        
        ###update
        beta=beta.new 
        theta=theta.new 
        D=covStruct(theta,cc=cc,dist=dist)
        
      }
      
      beta.res=rbind(beta.res, as.vector(beta))
      b.res=rbind(b.res, t(as.matrix(b))) # change b to t(b)
      theta.res=rbind(theta.res, as.vector(theta))
      alpha.res=cbind(alpha.res, alpha.vec[i])
      lambda.res=cbind(lambda.res, lambda.vec[j])
      
      if(verbose)
      {
        print(beta)
        print(theta)
        print(nIter)
      }
      
    }
    
  }
  
  ori.theta.res=cbind(exp(theta.res[,1]), exp(theta.res[,2]))
  
  return(list(beta.res=beta.res, b.res=b.res, theta.res=theta.res, ori.theta.res=ori.theta.res, alpha.res=alpha.res, lambda.res=lambda.res, alpha.vec=alpha.vec, lambda.vec=lambda.vec)) 

}

################################################################################
#Performs approximate penalized loglikelihood method that ignores spatial correlation.
#a sub-routine used by the major function SpatialVS.
#arguments are similar to the SpatialVS
SpatialVS_P_testLP.nocor=function(data, start=NULL, control=control.default, alpha.vec, lambda.vec, cc, cov.fun.type="exact", adaptive=F, lasso.weight=NULL, intercept=T, verbose=T)
{

  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist
  
  n=length(y)
  
  beta.res=theta.res=b.res=alpha.res=lambda.res=NULL
  
  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }
  
  ## the covariance fun
  switch(cov.fun.type,
         "exact"={
           cov.fun=covobj_LP.nocor
         })
  
  for(i in 1:length(alpha.vec))
  {
    for(j in 1:length(lambda.vec))
    {
      if(verbose)
      {
        cat("i=", i, "alpha=", alpha.vec[i], "j=", j, "lambda=", lambda.vec[j], "\n")
      } 
      
      nIter=0
      
      beta=start$beta    
      theta=start$theta[1]     
      D=covStruct.nocor(theta,n=n)      
      convCov=convPar=1
      p=length(beta)
      
      b=start$b

      if(is.null(lasso.weight))
      {
        lasso.weight=rep(1, p)
      }
      
      pf1.wts=1/abs(lasso.weight)
      
      if(intercept)
      {
        pf1.wts[1]=0
      }
      
      ###main event
      while((control$maxIter > nIter) & (convCov > control$tol1| convPar > control$tol2 ) )
      {
        
        # iterate beta given theta
        # estimate beta given theta
        nIter_beta=1
        convPar=1
        
        # estimate beta given theta
        while( (convPar > control$tol2)&(20 > nIter_beta) )
        {
          nIter_beta=nIter_beta+1
          beta.ini=beta
          
          for(s in 1:p)
          {
            e_s=rep(0, p)
            e_s[s]=1
            
            res1=iwls_b_new(X=X,y=y,D=D,offset=offset,beta=beta,b=b,tol=control$iwls,dist=dist)
            
            beta=res1$beta 
            mu=res1$mu 
            b=res1$b 
            eta=log(mu)
            
            x_s=as.numeric(X[,s])
            DW=D%*%diag(mu)
            
            if(intercept & s==1)
            {              
               beta_sm=sum(x_s*(mu-y))+0.5*sum(diag(solve(DW+diag(rep(1, n)))%*%DW%*%diag(x_s)))           
               h_sm=sum(x_s*x_s*mu)    
            
            }else{
                   beta_sm=sum(x_s*(mu-y))+2*lambda.vec[j]*(1-alpha.vec[i])*beta[s]+0.5*sum(diag(solve(DW+diag(rep(1, n)))%*%DW%*%diag(x_s)))           
                   h_sm=sum(x_s*x_s*mu)+2*lambda.vec[j]*(1-alpha.vec[i])            
                 }

            if(intercept & s==1)
            {              
               d_sm=-beta_sm/h_sm
            
            }else{
                    w_s=pf1.wts[s]
                    d_sm=median(c((w_s*lambda.vec[j]*alpha.vec[i]-beta_sm)/h_sm, -beta[s], (-w_s*lambda.vec[j]*alpha.vec[i]-beta_sm)/h_sm))
                 }

            
            beta=beta+d_sm*e_s
            
          }
          beta.new=beta          
          convPar <- sum(sqrt(crossprod(beta.new-beta.ini))/(1+sqrt(crossprod(beta.new))))
          #print(convPar)
          
        }
        
        

        # estimate covariance parameters
        res3=optim(theta, cov.fun, beta=as.vector(beta.new), X=X, y=y, cc=cc, b.start=as.vector(b), offset=offset,  dist=dist, lower=-5, upper=10, method="L-BFGS-B")
        
        theta.new=res3$par
        
        ###covergence check
        convCov <- sum(sqrt(crossprod(theta.new-theta))/(1+sqrt(crossprod(theta.new))))


        
        nIter=nIter+1
        
        #print(paste("Number of iterations:", nIter))
        
        ###update
        beta=beta.new 
        theta=theta.new 
        D=covStruct.nocor(theta,n=n)
        
      }
      
      #print(paste("Number of iterations:", nIter))
      
      beta.res=rbind(beta.res, as.vector(beta))
      b.res=rbind(b.res, t(as.matrix(b))) # change b to t(b)
      theta.res=rbind(theta.res, as.vector(theta))
      alpha.res=cbind(alpha.res, alpha.vec[i])
      lambda.res=cbind(lambda.res, lambda.vec[j])
      
      if(verbose)
      {
        print(beta)
        print(theta)
        print(nIter)
      }
          
    }
  }
  
  ori.theta.res=cbind(exp(theta.res))
  
  return(list(beta.res=beta.res, b.res=b.res, theta.res=theta.res, ori.theta.res=ori.theta.res, alpha.res=alpha.res, lambda.res=lambda.res, alpha.vec=alpha.vec, lambda.vec=lambda.vec))
   
}

################################################################################
#Used by SpatialVS in the estimation of covariance parameters under spatial correlation
covobj_LP=function(theta, cc, dist, X, y, offset=offset, beta, b.start, tol=control.default$iwls)
{
  D=covStruct(theta,cc=cc,dist=dist)
  n=length(y)
  
  b=b.start  
  mu=exp(X%*%as.vector(beta)+b+offset)
  mu=as.vector(mu)
  
  
  detextraterm=determinant(D%*%diag(as.vector(mu))+diag(rep(1,n)))
  
  res=sum(mu)-sum(y*(X%*%as.vector(beta)+b))+0.5*t(b)%*%solve(D)%*%b+0.5*(detextraterm$modulus)*(detextraterm$sign)
  
  return(as.numeric(res))
  
}

################################################################################
#Used by SpatialVS in the estimation of covariance parameters under no spatial correlation
covobj_LP.nocor=function(theta, cc, dist, X, y, offset=offset, beta, b.start, tol=control.default$iwls)
{
  n=length(y)
  
  D=covStruct.nocor(theta,n=n)
  
  b=b.start
  eta=X%*%as.vector(beta)+b+offset
  eta=as.vector(eta)
  mu=exp(eta)
  
  detextraterm=determinant(D%*%diag(as.vector(mu))+diag(rep(1,n)))
  
  res=0.5*t(b)%*%solve(D)%*%b +0.5*(detextraterm$modulus)*(detextraterm$sign)
  res=as.numeric(res)
  
  return(res)
  
}

################################################################################
#Function used to compute the AIC/BIC of spatial fitting from functions such SpatialVS_P_testPQL
#obj is the output of SpatialVS_P_testPQL
#data is the data matrix
#cc is a variable used in computing covariance matrix
SpatialVS.out=function(obj, data, cc)
{
  
  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist
  
  n=length(y)
  
  lld=bic=aic=dev=ERIC=rep(0, dim(obj$beta.res)[1])
  
  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }
  
  for(i in 1:length(dev))
  {

    tt=as.numeric(X%*%obj$beta.res[i,]) + as.numeric(obj$b.res[i,]) + offset
    D=covStruct(obj$theta.res[i,], cc=cc, dist=dist)
    
    detextraterm=determinant(D%*%diag(as.vector(exp(tt)))+diag(rep(1,n)))
    
    if(sum(log(factorial(y)))==Inf)
    {
      lld[i]=-sum(exp(tt)) + as.numeric( t(y) %*% (tt) ) -  
        as.numeric(0.5*t(obj$b.res[i,])%*%solve(D)%*%(obj$b.res[i,]))
        -0.5*(detextraterm$modulus)*(detextraterm$sign)
    }else{
      lld[i]=-sum(exp(tt)) + as.numeric( t(y) %*% (tt) ) - sum(log(factorial(y))) - 
        as.numeric(0.5*t(obj$b.res[i,])%*%solve(D)%*%(obj$b.res[i,]))
        -0.5*(detextraterm$modulus)*(detextraterm$sign)
    }
    
    k=sum(obj$beta.res[i,]!=0)+2
    
    bic[i]=-2*lld[i]+k*log(n)
    
    aic[i]=-2*lld[i]+k*2
    
    dev[i]=sum(2*(y[y!=0]*log(y[y!=0]/exp(tt[y!=0]))))+sum(2*(-(y-exp(tt))))
    

    
    ERIC[i]=-2*lld[i]+2*k*log(n/((obj$alpha.res[,i])*(obj$lambda.res[,i])))
  }
  
  return(list(lld=lld, bic=bic, dev=dev, aic=aic, ERIC=ERIC))
  
}

################################################################################
#Function used to compute the AIC/BIC of spatial fitting from functions such SpatialVS_P_testPQL.nocor
#obj is the output of SpatialVS_P_testPQL
#data is the data matrix
#cc is a variable used in computing covariance matrix
SpatialVS.out.nocor=function(obj, data, cc)
{
  
  y=data$y
  X=data$X
  offset=data$offset
  dist=data$dist
  
  n=length(y)
  
  lld=bic=aic=dev=rep(0, dim(obj$beta.res)[1])
  
  ### if there is no offset term in the model
  if(is.null(offset))
  {
    offset=rep(0, n)
  }
  
  for(i in 1:length(dev))
  {
    tt=as.numeric(X%*%obj$beta.res[i,]) + as.numeric(obj$b.res[i,]) + offset
    D=covStruct.nocor(obj$theta.res[i,],n=n)
    
    detextraterm=determinant(D%*%diag(as.vector(exp(tt)))+diag(rep(1,n)))
    
    lld[i]=-sum(exp(tt)) + as.numeric( t(y) %*% (tt-offset) ) - 
      as.numeric(0.5*t(obj$b.res[i,])%*%solve(D)%*%(obj$b.res[i,]))-
      0.5*(detextraterm$modulus)*(detextraterm$sign)
    
    k=sum(obj$beta.res[i,]!=0)+1
    
    bic[i]=-2*lld[i]+k*log(n)
    
    aic[i]=-2*lld[i]+k*2
    
    dev[i]=sum(2*(y[y!=0]*log(y[y!=0]/exp(tt[y!=0]))))+sum(2*(-(y-exp(tt))))
  }
  
  return(list(lld=lld, bic=bic, dev=dev, aic=aic))
  
}

################################################################################
#intermediate function, performs iterated weighted least squares
iwls=function(X, y, D, offset, beta, b, tol, dist)
{
  
  eta = X%*%beta + b + offset
  mu = exp(eta)
  

  
  convEta=10
  
  while(convEta>tol)
  {
    ##define the working response and weight
    
    mu=as.vector(mu)
    Y = eta + (y-mu)/mu - offset
    V = diag(1/mu) + D
    
    solveV=solve(V)
    
    
    ###update beta and b 
    XVX = crossprod(X, solveV)%*%X
    XVY = crossprod(X, solveV)%*%Y
    DVX = crossprod(D, solveV)%*%X
    DVY = crossprod(D, solveV)%*%Y
    
    
    beta.new=solve(XVX, XVY)
    b.new=DVY-DVX%*%beta.new
    
    eta.new=X%*%beta.new + b.new + offset
    
    ###covergence 
    convEta <- sum((eta.new-eta)^2)/sum(eta.new^2)
    
    ###updata the eta and mu
    eta = eta.new
    mu = exp(eta)
    
    #print(convEta)
  } 
  
  return(list(beta=beta.new, b=b.new, eta=eta, mu=mu))
  
}

################################################################################
#Used by SpatialVS in the estimation of covariance parameters under penalty
covobj_PQLpenb=function(theta, cc, mu, nonzerobeta, beta, b.start, tol=control.default$iwls, n, X_L, X, y, Sigma_pen,  dist, offset)
{ 
  D=covStruct(theta,cc=cc,dist=dist)
  
  res1=iwls_b(X=X, y=y, D=D, offset=offset, beta=beta, b=b.start, tol=tol,dist=dist)
  b=res1$b
  mu=res1$mu


  V=diag(1/mu)+D
  solveV=solve(V)
  detV=determinant(V)
  detXVX=determinant(t(X_L)%*%solveV%*%X_L+Sigma_pen)

  ystar=X%*%beta+b+(y-mu)/mu

  res=as.numeric(0.5*detV$modulus*detV$sign+0.5*detXVX$modulus*detXVX$sign+0.5*t(ystar-X_L%*%nonzerobeta)%*%solveV%*%(ystar-X_L%*%nonzerobeta))

  return(res)
}

################################################################################
#Used by SpatialVS in the estimation of covariance parameters under penalty when ignores spatial correlation
covobj_PQLpen.nocor=function(theta, cc, mu, nonzerobeta, n, X_L, Sigma_pen, ystar, dist)
{

  V=diag((1/mu)+exp(theta))

  solveV=diag(1/((1/mu)+exp(theta)))

  detV=sum(log((1/mu)+exp(theta)))
  detXVX=determinant(t(X_L)%*%solveV%*%X_L+Sigma_pen)

  return(as.numeric(0.5*detV+0.5*detXVX$modulus*detXVX$sign+
                      0.5*t(ystar-X_L%*%nonzerobeta)%*%solveV%*%(ystar-X_L%*%nonzerobeta)))
}

################################################################################
###get the current estimates of  b using iwls
iwls_b=function(X, y, D, offset, beta, b, tol, dist)
{  

  
  eta = X%*%beta + b + offset
  mu = exp(eta) 
  convEta=10  
  solveD=solve(D)
    

    
  while(convEta>tol)
  {
    ##define the working response and weight
  
    mu=as.vector(mu)
    fir_devb = (mu-y) + solveD%*%b
    sec_devb = diag(mu) + solveD 
    
    ###b 

    b.new=b-solve(sec_devb)%*%fir_devb
    eta.new=X%*%beta + b.new + offset
    
    ###covergence 
    convEta <- sum((eta.new-eta)^2) /sum(eta^2)
    
    ###updata the eta and mu
    eta = eta.new

    
    mu = exp(eta)
    b = b.new
    
    #print(convEta)
  } 
  
  return(list(beta=as.numeric(beta), b=as.numeric(b), eta=as.numeric(eta), mu=as.numeric(mu)))
  
}

################################################################################
####construct the covariance matrix
covStruct=function(theta,cc,dist)
{
  sigma_sq=exp(theta[1])
  alpha=exp(theta[2])
  
  cov=exp(-dist/alpha)
  
  if(cc==0)
  {
    res=sigma_sq*cov

    
    return(res)
  
  } else{
    tt=Wendland(dist, theta=cc, k=1, dimension=2)
    cov_new=cov*tt
    cov_new=as(cov_new, "sparseMatrix")
    return((sigma_sq)*cov_new)
  } 

}

################################################################################
####construct the covariance matrix
covStruct.nocor=function(theta, n)
{
  sigma_sq=exp(theta)
  
  return(sigma_sq*diag(rep(1, n)))
}

################################################################################
#make a contour plot for the AIC or BIC for penalty tuning parameters.
SpatialVS.obj.contour.plot=function(SpatialVS.obj, SpatialVS.out.obj, type="BIC", nlevels=5, sep=5, plots=T)
{
  
  x=SpatialVS.obj$alpha.vec
  y=SpatialVS.obj$lambda.vec
  
  switch(type,
         "BIC"={
           z=matrix(SpatialVS.out.obj$bic, nrow=length(x), ncol=length(y), byrow=T)
         },
         "AIC"={
           z=matrix(SpatialVS.out.obj$aic, nrow=length(x), ncol=length(y), byrow=T)
         },
         "Deviance"={
           z=matrix(SpatialVS.out.obj$dev, nrow=length(x), ncol=length(y), byrow=T)
         },
         "LogLikelihood"={
           z=matrix(SpatialVS.out.obj$lld, nrow=length(x), ncol=length(y), byrow=T)
         })
  
  id.x=order(x)
  id.y=order(y)
 

  tmp.idx.mat=(z==min(z))
  alpha.best.val=x[rowSums(tmp.idx.mat)>0]
  lambda.best.val=y[colSums(tmp.idx.mat)>0]
  alpha.best.val=alpha.best.val[1]
  lambda.best.val=lambda.best.val[1]
    
  if(plots)
  {
    contour(x=x[id.x], y=y[id.y], z=z[id.x, id.y],  xlab=expression(alpha), ylab=expression(lambda),levels=seq(min(z)+0.01, min(z)+3,,5), main=type)
  
    points(alpha.best.val, lambda.best.val, pch=16, col=2, lwd=2)
    

  }
  
  print(c(alpha.best.val, lambda.best.val))
  
  res=list(x=x, y=y, z=z, alpha.best.val=alpha.best.val, lambda.best.val=lambda.best.val)
  
  return(invisible(res))
}


####################################################################
# simulate model matrix in simulation study
spatialVS_simuX=function(ll, numofcor, rho, noeff_corr,  p.total)
{
  # number of sample
  n=ll^2 
  
  # model matrix
  X=NULL
  
  # if we consider correlation structure among noeff predictors
  tmpX=NULL
  
  Xsq=0.2
  #Xsq=0.5
  
  # independent predictors
  if(all(numofcor==0))
  {
    tmp=rnorm(n*p.total, mean=0, sd=sqrt(Xsq))
    tmp[abs(tmp)>2*sqrt(Xsq)]=2*sqrt(Xsq)
    X=cbind(X, matrix(tmp, n, p.total))    
  }else{
    # correlated predictors
    for(i in 1:length(numofcor))
    {
      cor.mat=matrix(0, numofcor[i], numofcor[i])
      for(j in 1:numofcor[i])
      {
        cor.mat[j,]=rho[i]^(abs(j-(1:numofcor[i])))
      }
      
      X=cbind(X, mvrnorm(n = n, mu=rep(0, numofcor[i]), Sigma=Xsq*cor.mat))
      
      if(noeff_corr){
        # same correlation structure for noeff predictors
        tmpX=cbind(tmpX, mvrnorm(n = n, mu=rep(0, numofcor[i]), Sigma=Xsq*cor.mat))       
      }
    }
    
    if(noeff_corr){
      pp=p.total-ncol(X)-ncol(tmpX)
      X=cbind(X, tmpX, matrix(rnorm(n*pp, mean=0, sd=sqrt(Xsq)), n, pp))
    }else{
      pp=p.total-ncol(X)
      X=cbind(X, matrix(rnorm(n*pp, mean=0, sd=sqrt(Xsq)), n, pp))
    }
  }
  
  
  return(scale(X, scale=TRUE))
}

####################################################################
# simulate data in simulation study
spatialVS_simu_data=function(ll, numofcor, rho, cov.pars, betaeff, noeff_corr=F, p.total)
{
  n=ll^2
  
  pp=length(betaeff)
  
  contX=TRUE
  while(contX)
  {
     X=spatialVS_simuX(ll=ll, numofcor=numofcor,rho=rho,noeff_corr=noeff_corr,  p.total=p.total)
    
    if(max(exp(X[,1:pp]%*%betaeff))<100)
    {
      contX=FALSE
    }
  
  }

  
  beta=c(betaeff, rep(0, p.total-pp))
  
  # distance matrix
  dist=dist1(as.matrix(expand.grid(x = seq(1, 10, l = ll), y = seq(1, 10, l = ll))))
  
  D=covStruct(theta=log(cov.pars), cc=0, dist=dist)
  b=mvrnorm(1, mu=rep(0, n), Sigma=D)
  
  y=rpois(n, lambda=exp(X%*%beta+b))
  
  return(simu.data=list(y=y, X=X, dist=dist, offset=rep(0, n)))
  
}

####################################################################
#function compute the distance matrix
dist1=function(mat)
{
  nn=nrow(mat)
  
  xmat1=kronecker(mat[,1], t(rep(1, nn)), "*")
  xmat2=kronecker(rep(1, nn), t(mat[,1]), "*")
  
  a1=(xmat1-xmat2)^2

  ymat1=kronecker(mat[,2], t(rep(1, nn)), "*")
  ymat2=kronecker(rep(1, nn), t(mat[,2]), "*")
  
  a2=(ymat1-ymat2)^2
  
  res=sqrt(a1+a2)
  
  return(res)
  
}
################################################################################
# calculate initial values for the simulation study
spatialVS_simu_ini=function(simu.data, intercept=T, verbose=T)
{
  y=simu.data$y
  X=simu.data$X
  dist=simu.data$dist
  offset=simu.data$offset
  
  ID=1:length(y)
  glmmPQL.ini=glmmPQL(y~X-1+offset(offset), random=~1|ID, family=poisson, verbose = F)

  sigma.ini=try(exp(attr(glmmPQL.ini$apVar, "Pars")["reStruct.ID"]), silent=T)
  
  b.ini=as.numeric(ranef(glmmPQL.ini)[,1])
  beta.ini=fixef(glmmPQL.ini)
  
  if(class(sigma.ini)=="try-error"|sigma.ini<0.001)
  {
    sigma.ini=0.3
    b.ini=rnorm(length(y), mean=0, sd=sigma.ini)
    
    D=covStruct(theta=log(c(sigma.ini^2, 4)), cc=0, dist=dist)  
    b.ini=as.numeric(iwls_b(X=X, y=y, D=D, offset=offset, beta=beta.ini, b=b.ini, tol=0.001, dist=dist)$b)
  }
  
  dd.ini=optim(par=log(1.5), fn=cov.fun.fixsigma, method="Brent", lower=-5, upper=10, sigma.ini=sigma.ini, dist=dist, glmmPQL.ini=glmmPQL.ini,beta.ini=beta.ini, b.ini=b.ini, y=y, X=X)$par
  
  start.ini=list(beta=beta.ini, theta=c(log(sigma.ini^2), dd.ini), b=b.ini)
  
  tmp=SpatialVS_P_testPQL(data=simu.data, start=start.ini, alpha.vec=0, lambda.vec=0, cc=0, adaptive=F, intercept=intercept, verbose=verbose)

  start=list(beta=as.numeric(start.ini$beta), theta=as.vector(start.ini$theta), b=as.numeric(start.ini$b), lasso.weight=as.numeric(start.ini$beta))
  
  
  return(start=start)
}

####################################################################
#the main function used in the simulation study
#mm is the number of repeats
#ll is the size of the spatial region
#numofcor: number of correlated covariates
#betaeff: beta values
#cov.pars: covariance parameters
#p.total: total number of covariates
#other variables are similar to spatialVS
spatialVS_simu=function(mm, ll=15, numofcor=5, rho=0.8, betaeff=c(-0.5, 0.75, 1, -0.75, -1), cov.pars=c(0.1, 5), p.total=15, alpha.vec=seq(0.6, 1, by=0.1), lambda.vec=seq(0.01, 0.5, by=0.02), cc=0, adaptive=F, noeff_corr=F, nocor=F, lambda.vec.ext=NULL, realX=F, realX.file=NULL)
{  

  #simulate data
  if(!realX)
  {
  
    simu.data=spatialVS_simu_data(ll=ll, numofcor=numofcor, rho=rho, betaeff=betaeff,cov.pars=cov.pars, noeff_corr=noeff_corr, p.total=p.total)
  
  }
  
  if(realX)
  {
    load(file=realX.file)
    simu.data=spatialVS_simu_data_LymeX(dat.obj=lyme.eco0.boot.dat$dat.obj, fit.obj=lyme.eco0.boot.dat$fit.obj, subset.no=ll^2, intercept=F)  
  }
  
  
  # initial values
  start=try(spatialVS_simu_ini(simu.data=simu.data, intercept=F, verbose=F), silent=T)
  
  if(class(start)=="try-error")
  {
    return(simu.data=simu.data)
  }
  

  ##########################
  # select the range of lambda parameter
  # if the initial range of lambda parameter is not big enough, increase
  PQLmax.obj=try(SpatialVS_P_testPQL(data=simu.data,start=start, alpha.vec=1, lambda.vec=c(max(lambda.vec), 1, 10), cc=cc, adaptive=adaptive, lasso.weight=start$beta, intercept=F, verbose=F), silent=T)
  
  PQLmaxout.obj=try(SpatialVS.out(obj=PQLmax.obj, data=simu.data, cc=0), silent=T)
  
  if((class(PQLmax.obj)!="try-error")&(class(PQLmaxout.obj)!="try-error"))
  {

    bic.max=PQLmaxout.obj$bic
    
    if((bic.max[1]>bic.max[3])|(bic.max[1]>bic.max[2]))
    {
      lambda.vec=lambda.vec.ext
    }
    
    if(bic.max[2]>bic.max[3])
    {
      lambda.vec=c(lambda.vec, seq(max(lambda.vec.ext), 10, by=1))
    }
    
  }

  
  tmp.PQL=try(SpatialVS(dat.obj=simu.data, alpha.vec=alpha.vec, lambda.vec=lambda.vec, method="PQL", plots=F, intercept=F, verbose=F), silent=T)
   
  PQL.obj=tmp.PQL$L.best.obj 
  PQLout.obj=tmp.PQL$Lout.best.obj
  start=tmp.PQL$start
  lasso.weight=tmp.PQL$lasso.weight
  
  if(nocor)
  {
    tmp.PQL.nocor=try(SpatialVS(dat.obj=simu.data, alpha.vec=alpha.vec, lambda.vec=lambda.vec, method="PQL.nocor", plots=F, intercept=F, verbose=F), silent=T)
   
    PQL_nocor.obj=tmp.PQL.nocor$L.best.obj 
    PQL_nocorout.obj=tmp.PQL.nocor$Lout.best.obj
    start.nocor=tmp.PQL.nocor$start
    lasso.weight.nocor=tmp.PQL.nocor$lasso.weight
    
    return(list(simu.data=simu.data, start=start, start.nocor=start.nocor, PQL.obj=PQL.obj, PQLout.obj=PQLout.obj, PQL_nocor.obj=PQL_nocor.obj, PQL_nocorout.obj=PQL_nocorout.obj, lasso.weight=lasso.weight, lasso.weight.nocor=lasso.weight.nocor))
  
  }

  return(list(simu.data=simu.data, start=start, PQL.obj=PQL.obj, PQLout.obj=PQLout.obj, lasso.weight=lasso.weight))

}


####################################################################
#the main simulation function for the APL method
#arguments are similar to spatialVS
spatialVS_simu_LP=function(mm, alpha.vec=seq(0.6, 1, by=0.1), lambda.vec=seq(0.1, 1.5, by=0.1), cc=0, adaptive=F, name, nocor=F)
{  
  load(name)
  simu.data=result[[mm]]$simu.data
  start=result[[mm]]$start
  lasso.weight=result[[mm]]$lasso.weight
  
  if(is.null(start))
  {
    return(simu.data=simu.data)
  }

  LP.obj=try(SpatialVS_P_testLP(data=simu.data, start=start, alpha.vec=alpha.vec, lambda.vec=lambda.vec, cc=cc, adaptive=adaptive, lasso.weight=lasso.weight, intercept=F, verbose=F), silent=T)
  
  if(class(LP.obj)=="try-error")
  {
    return(list(simu.data=simu.data))
  }
  
  LPout.obj=try(SpatialVS.out(obj=LP.obj, data=simu.data, cc=0), silent=T)
  
  if(class(LPout.obj)=="try-error")
  {
    return(list(simu.data=simu.data))
  }
  
  
  if(nocor)
  {    
    
    start.nocor=result[[mm]]$start.nocor
    lasso.weight.nocor=result[[mm]]$lasso.weight.nocor
    
    if(is.null(start.nocor))
    {
      return(simu.data=simu.data)
    }
    
    LP_nocor.obj=try(SpatialVS_P_testLP.nocor(data=simu.data,start=start.nocor, alpha.vec=alpha.vec, lambda.vec=lambda.vec, cc=cc, adaptive=adaptive, lasso.weight=lasso.weight.nocor, intercept=F, verbose=F), silent=T)
    
    
    if(class(LP_nocor.obj)=="try-error")
    {
      return(list(simu.data=simu.data))
    }
    
    LP_nocorout.obj=try(SpatialVS.out.nocor(obj=LP_nocor.obj, data=simu.data, cc=0), silent=T)
    
    if(class(LP_nocorout.obj)=="try-error")
    {
      return(list(simu.data=simu.data))
    }
    
    return(list(simu.data=simu.data, start=start, LP.obj=LP.obj, LPout.obj=LPout.obj, LP_nocor.obj=LP_nocor.obj, LP_nocorout.obj=LP_nocorout.obj))  
  }
   
  return(list(simu.data=simu.data, start=start, LP.obj=LP.obj, LPout.obj=LPout.obj))
}

####################################################################
#simulation study for existing methods: LASSO, p-value based methods.
spatialVS_simu_existing_methods=function(mm, name, lasso.formula)
{  

  load(name)
  simu.data=result[[mm]]$simu.data
 

  pvalue.res<-try(spatialVS.glmmPQL.fit(dat.obj=simu.data, intercept=F, verbose=F))
  
  if(class(pvalue.res)=="try-error")
  {
    pvalue.res=NULL
  }

  backward.res<-try(spatialVS.glmmPQL.fit.backward(dat.obj=simu.data, intercept=F, verbose=F))

  if(class(backward.res)=="try-error")
  {
    backward.res=NULL
  }
  
  glmmLasso.res<-try(spatialVS.glmmLasso.fit(dat.obj=simu.data, formula=lasso.formula, intercept=F, verbose=F))

  if(class(glmmLasso.res)=="try-error")
  {
    glmmLasso.res=NULL
  }


  res=list(pvalue.res=pvalue.res, backward.res=backward.res, glmmLasso.res=glmmLasso.res)
  
  return(res)
   
}

####################################################################
# collect results in the simulation study
spatialVS_sim_coll=function(betaeff=c(-0.5, 0.75, 1, -0.75, -1),cov.pars=c(0.1, 5), name, pdfname, nocor=T, method="PQL", p.zero=10)
{  
  load(name)
  beta.res=theta.res=thres.res=NULL
  beta_nocor.res=theta_nocor.res=thres_nocor.res=NULL
  
  if (method=="LP")
  {
    for(i in 1:length(result))
    {
      #print(i)
      tmp=result[[i]]
      

      index=which.min(tmp$LPout.obj$bic)
      
      thres.res=rbind(thres.res, c(tmp$LP.obj$alpha.res[,index], tmp$LP.obj$lambda.res[,index]))
      beta.res=rbind(beta.res, tmp$LP.obj$beta.res[index,])
      theta.res=rbind(theta.res, tmp$LP.obj$ori.theta.res[index,])
      
      
      if(nocor)
      {
        index=which.min(tmp$LP_nocorout.obj$bic)
        
        thres_nocor.res=rbind(thres_nocor.res, c(tmp$LP_nocor.obj$alpha.res[,index], tmp$LP_nocor.obj$lambda.res[,index]))
        beta_nocor.res=rbind(beta_nocor.res, tmp$LP_nocor.obj$beta.res[index,])
        theta_nocor.res=rbind(theta_nocor.res, tmp$LP_nocor.obj$ori.theta.res[index,])  
      }
      
    }
    
  }else{
      for(i in 1:length(result))
      {
      #print(i)
      tmp=result[[i]]
      
      index=which.min(tmp$PQLout.obj$bic)
      
      if(length(tmp$PQLout.obj$bic)==1)
      {
        index=1
      }
      
      thres.res=rbind(thres.res, c(tmp$PQL.obj$alpha.res[,index], tmp$PQL.obj$lambda.res[,index]))
      beta.res=rbind(beta.res, tmp$PQL.obj$beta.res[index,])
      theta.res=rbind(theta.res, tmp$PQL.obj$ori.theta.res[index,])
      
      
      if(nocor)
      {
        index=which.min(tmp$PQL_nocorout.obj$bic)
        
        if(length(tmp$PQL_nocorout.obj$bic)==1)
        {
          index=1
        }
        
        thres_nocor.res=rbind(thres_nocor.res, c(tmp$PQL_nocor.obj$alpha.res[,index], tmp$PQL_nocor.obj$lambda.res[,index]))
        beta_nocor.res=rbind(beta_nocor.res, tmp$PQL_nocor.obj$beta.res[index,])
        theta_nocor.res=rbind(theta_nocor.res, tmp$PQL_nocor.obj$ori.theta.res[index,])
        
      }
      
    }
  }
  

  

  
  # remove the not coverge cases
  id=(theta.res[,1]>exp(-5))&(theta.res[,1]<100)&
     (theta.res[,2]>exp(-5))&(theta.res[,2]<100)
  
  
  beta.res=beta.res[id, ]
  theta.res=theta.res[id, ]
  thres.res=thres.res[id, ]
  

  
  if(nocor){
    id_nocor=(theta_nocor.res>exp(-5))&(theta_nocor.res<100)
    
    beta_nocor.res=beta_nocor.res[id_nocor, ]
    theta_nocor.res=theta_nocor.res[id_nocor]
    thres_nocor.res=thres_nocor.res[id_nocor, ]
  }

  beta=c(betaeff, rep(0, p.zero))
  
  p=length(beta)
  p1=length(betaeff)
  

  
  # histograms of parameters
  pdf(file=pdfname, height=6, width=6)
  
  for(i in 1:p)
  {
    hist(beta.res[,i], xlab="", main=eval(expression(substitute(beta[i], list(i=i)))))
    abline(v=beta[i], col=2)
  }
  
  for(i in 1:dim(theta.res)[2])
  {
    hist(theta.res[,i], xlab="", main=eval(expression(substitute(theta[i], list(i=i)))))
    abline(v=cov.pars[i], col=2)
  }
  
  for(i in 1:dim(thres.res)[2])
  {
    hist(thres.res[,i], xlab="", main=eval(expression(substitute(lambda[i], list(i=i)))))
  }
  
  if(nocor)
  {
    for(i in 1:p)
    {
      hist(beta_nocor.res[,i], xlab="", main=eval(expression(substitute(beta[i], list(i=i)))))
      abline(v=beta[i], col=2)
    }
    
    hist(theta_nocor.res, xlab="", main=expression(theta[1]))
    abline(v=cov.pars[1], col=2)
    
    for(i in 1:dim(thres_nocor.res)[2])
    {
      hist(thres_nocor.res[,i], xlab="", main=eval(expression(substitute(lambda[i], list(i=i)))))
    }
  }
  dev.off()
  
  function1=function(x)
  {
    sum(x!=0)
  }
  
  function2=function(x)
  {
    sum(x[(p1+1):p]==0)
  }
  
  function3=function(x)
  {
    sum(x[1:p1]==0)
  }
  
  modelsize=apply(beta.res, 1, function1)
  corr0=apply(beta.res, 1, function2)
  inc0=apply(beta.res, 1, function3)
  
  perc.freq=apply(beta.res, 2, function1)/(dim(beta.res)[1])
  
  if(nocor)
  {
    modelsize_nocor=apply(beta_nocor.res, 1, function1)
    corr0_nocor=apply(beta_nocor.res, 1, function2)
    inc0_nocor=apply(beta_nocor.res, 1, function3)
    
    perc.freq_nocor=apply(beta_nocor.res, 2, function1)/(dim(beta_nocor.res)[1])
    
    return(rbind(c(mean(modelsize), mean(corr0), mean(inc0), perc.freq),
                 c(mean(modelsize_nocor), mean(corr0_nocor), mean(inc0_nocor), perc.freq_nocor)))
  }
  
  return(c(mean(modelsize), mean(corr0), mean(inc0), perc.freq) )
}

#######################################################################
#collect the results in the simulation study for the existing methods.
spatialVS_sim_coll.ExM=function(betaeff=c(-0.5, 0.75, 1, -0.75, -1),cov.pars=c(0.1, 5), name, p.zero=10)
{  
  load(name)
  
  nn=length(result)
  
  pvalue.res.mat=NULL
  backward.res.mat=NULL
  glmmLasso.res.mat=NULL
  
  
  for(i in 1:nn)
  {
    tmp=result[[i]]
    pvalue.res.mat=rbind(pvalue.res.mat, tmp$pvalue.res)
    backward.res.mat=rbind(backward.res.mat, tmp$backward.res)
    glmmLasso.res.mat=rbind(glmmLasso.res.mat, tmp$glmmLasso.res)
    
  }
  

  
  beta=c(betaeff, rep(0, p.zero))
  
  p=length(beta)
  p1=length(betaeff)
  


  function1=function(x)
  {
    sum(x!=0)
  }
  
  function2=function(x)
  {
    sum(x[(p1+1):p]==0)
  }
  
  function3=function(x)
  {
    sum(x[1:p1]==0)
  }
  
  sum.fun=function(mat)
  {
    if(!is.null(mat))
    {
     modelsize=apply(mat, 1, function1)
     corr0=apply(mat, 1, function2)
     inc0=apply(mat, 1, function3)
  
     tres=c(mean(modelsize), mean(corr0), mean(inc0))
    }else{
           tres=NULL
         }
         
    return(tres)
    
  }
  
  pvalue.sum=sum.fun(mat=pvalue.res.mat)
  back.sum=sum.fun(mat=backward.res.mat)
  glmmLasso.sum=sum.fun(mat=glmmLasso.res.mat)

  res=rbind(pvalue.sum, back.sum, glmmLasso.sum)
  
  return(res)
  
}


################################################################################
#function for adaptive elastic net, used in SpatialVS
weighted.ada.enet=function(xmat, y, obs.wts, alpha, lambda, pf1, pf2)
{
  max.iter=100
  eps=1e-4
   
  wts.sqrt=sqrt(obs.wts)
  xmat.wtd=sweep(xmat, 1, wts.sqrt, "*")
  y.wtd=y*wts.sqrt
  
  xx.mat=t(xmat.wtd)%*%xmat.wtd
  xx.mat.inv=solve(xx.mat)
  beta0=xx.mat.inv%*%t(xmat.wtd)%*%y.wtd
  
  beta0=as.vector(beta0)
  beta1=beta0

  pp=length(beta0)
  
  aa=colSums(xmat.wtd^2)+2*lambda*(1-alpha)*pf2
  aa=as.vector(aa)
  cc=pf1*lambda*alpha

  max.dd=1
  iter=1
  
  
  while(iter<max.iter & max.dd>eps)
  {  
    for(j in 1:pp)
    {
      rr=y.wtd-xmat.wtd%*%beta1
      dyy=rr+xmat.wtd[,j]*beta1[j]      
      bb=sum(dyy*xmat.wtd[,j])
      btmp=soft.thresholding(z=bb, r=cc[j])     
      beta1[j]=btmp/aa[j]  
    }
      
    iter=iter+1
    max.dd=max(abs(beta1-beta0))
    beta0=beta1
    
    #cat(max.dd, iter, "\n")
    
  }
  

 
  res=list(coef=beta1)
  
  return(res)
}
################################################################################
#function used in adaptive enet
soft.thresholding=function(z, r)
{
   res=0
   
   if(z>0 &r<abs(z))
   {
     res=z-r
   }

   if(z<0 &r<abs(z))
   {
     res=z+r
   }

  return(res)
}
################################################################################
#simulation setting that uses the Lyme disease covariates.
spatialVS_simu_data_LymeX=function(dat.obj, fit.obj, subset.no=NULL, intercept=T)
{
  X=dat.obj$X
  offset=dat.obj$offset
  dist=dat.obj$dist   

  if(is.null(subset.no))
  {
    nn=nrow(X)
    idx=1:nn
  }else{
         nn=subset.no
         idx=sample(1:(nrow(X)), nn)
         idx=sort(idx)
       }
  
  beta.res=as.vector(fit.obj$beta.res)
  theta.res=fit.obj$theta.res
  
  dist.sim=dist[idx, idx]
  
  D=covStruct(theta=theta.res, cc=0, dist=dist.sim)

  b=mvrnorm(1, mu=rep(0, nn), Sigma=D)
  
  X.sim=X[idx,]
  offset.sim=offset[idx]
  
  if(!intercept)
  {
    X.sim=X.sim[,-1]
    offset.sim=offset.sim+beta.res[1]
    beta.res=beta.res[-1]
  }
  
  lambda=exp(X.sim%*%beta.res+b+offset.sim)
  y=rpois(nn, lambda=lambda)
  

  
  simu.data=list(y=y, X=X.sim, dist=dist.sim, offset=offset.sim, b.true=b)
  
  return(simu.data)

}
################################################################################
#bootstrap for the SpatialVS output for the Lyme disease data
spatialVS_boot_data_LymeX=function(dat.obj, fit.obj, subset.no=NULL)
{
  X=dat.obj$X
  offset=dat.obj$offset
  dist=dat.obj$dist   

  if(is.null(subset.no))
  {
    nn=nrow(X)
    idx=1:nn
  }else{
         nn=subset.no
         idx=sample(1:(nrow(X)), nn)
         idx=sort(idx)
       }
  
  beta.res=as.vector(fit.obj$beta.res)
  b.res=as.vector(fit.obj$b.res)
  theta.res=fit.obj$theta.res
  
  dist.sim=dist[idx, idx]
  
  D=covStruct(theta=theta.res, cc=0, dist=dist.sim)
  
  Omega=D/exp(theta.res[1])
  xx.mat=diag(nn)-1/nn
  sigma2.hat.adj=var(b.res)*(nn-1)/sum(xx.mat*Omega)
  
  sigma2.hat.adj=max(sigma2.hat.adj, exp(theta.res[1]))
  
  D.adj=sigma2.hat.adj*Omega
 
  b=mvrnorm(1, mu=rep(0, nn), Sigma=D.adj)
  
  X.sim=X[idx,]
  offset.sim=offset[idx]
  
  lambda=exp(X.sim%*%beta.res+b+offset.sim)
  y=rpois(nn, lambda=lambda)
  

  
  simu.data=list(y=y, X=X.sim, dist=dist.sim, offset=offset.sim, b.true=b)
  
  return(simu.data)

}

################################################################################
#bootstrap for the SpatialVS output for the Lyme disease data
lyme.fit.boot=function(sf.id, boot.dat, alpha.vec=seq(0.05, 1,, 10), lambda.vec=seq(0.05, 10, len=10), method="PQL")
{
   dat.obj=boot.dat$dat.obj
   fit.obj=boot.dat$fit.obj
   
   tmp.sim.dat=spatialVS_boot_data_LymeX(dat.obj=dat.obj, fit.obj=fit.obj)
   tmp.sim.dat.obj<-SpatialVS(dat.obj=tmp.sim.dat, alpha.vec=alpha.vec, lambda.vec=lambda.vec, method=method, intercept=T)
   res=SpatialVS.summary(obj=tmp.sim.dat.obj)
   return(res)
   
}

################################################################################
#collect bootstrap results
lyme.boot.res.collection=function(filename, seqs)
{

  
  aa=paste(filename,".",seqs[1],sep="")
  bb=load(file=aa)
  tmp=get(bb)  
  
  par.est=NULL
  
  mm=length(tmp)


  
  for(j in 1:mm)
  {
    xtmp=tmp[[j]]
    if(!is.null(xtmp))
    {
      par.est=rbind(par.est, xtmp)
    }
  }
  nn=length(seqs)
  if(nn>1)
  {
   for(i in 2:nn)
   {
    aa=paste(filename,".",seqs[i],sep="")
    bb=load(file=aa)
    tmp=get(bb)

    mm=length(tmp)
    for(j in 1:mm)
    {
      xtmp=tmp[[j]]
      if(!is.null(xtmp))
      {
       par.est=rbind(par.est, xtmp)              
      }     
    }
   }
  }
   
  res=par.est 

  print(dim(res))
    
  return(res)
}

################################################################################
#function for CI
tmp.ci.fun=function(x)
{
  alpha=0.05
  res=quantile(x,c(alpha/2, 1-alpha/2),na.rm=T)
  res=as.vector(res)
  return(res)
}
################################################################################
#function for CI
tmp.BC.ci.fun=function(x)
{
  alpha=0.05
  zz=x[1]
  x=x[-1]
  bb=mean(x<=zz)
  alpha.adj=pnorm(2*qnorm(bb)+qnorm(c(alpha/2, 1-alpha/2)))
  res=quantile(x,alpha.adj,na.rm=T)
  res=as.vector(res)
  return(res)
}
################################################################################
#function for CI
lyme.boot.est.ci=function(boot.dat, est.mat) 
{
  fit.obj=boot.dat$fit.obj
  
  beta.est=as.vector(fit.obj$beta.res)
  theta.est=as.vector(fit.obj$theta.res)
  
  
  est.mat.idx=rowSums(est.mat==10)+rowSums(est.mat==(-5))
  est.mat=est.mat[est.mat.idx==0, ]
  
  #print(dim(est.mat))
  #bias correction
  est.mat[,1]=est.mat[,1]-mean(est.mat[,1])+beta.est[1]
  est.mat[,length(beta.est)+1]=est.mat[,length(beta.est)+1]-mean(est.mat[,length(beta.est)+1])+theta.est[1]
  est.mat[,length(beta.est)+2]=est.mat[,length(beta.est)+2]-mean(est.mat[,length(beta.est)+2])+theta.est[2]
  
  var.name=c("Intercept","X.Dvlpd_NLCD06", "X.Forest_NLCD06", "X.Scrub_NLCD06", "X.Tract_Frag06", "X.FragPerim06", "CWED_DF06", "TECI_DF06", "CWED_FH06", "TECI_FH06", "CWED_HD06", "TECI_HD06", "pop_den_census10", "median_age", "mean_income")
  
  ci.mat=apply(est.mat, 2, tmp.ci.fun)
  ci.mat=t(ci.mat)

  xtmp=rbind(c(beta.est, theta.est), est.mat)
  ci.bc.mat=apply(xtmp, 2, tmp.BC.ci.fun)
  ci.bc.mat=t(ci.bc.mat)
  
  res=cbind(c(beta.est, theta.est), ci.mat, ci.bc.mat)
  
  pp=length(beta.est)
  pp.theta=length(theta.est)
  
  res[(pp+1):(pp+pp.theta),]=exp(res[(pp+1):(pp+pp.theta),])
  
  theta.var.name=c("sigma2", "dd")
  var.name=c(var.name, theta.var.name[1:pp.theta])
  
  res=data.frame(res)
  row.names(res)=as.vector(var.name)
  colnames(res)=c("Est", "CI.lower", "CI.upper", "CI.bc.lower", "CI.bc.upper")
  
  return(res)  
  
}
################################################################################
#function for glmmPQL, used in comparison
spatialVS.glmmPQL.fit=function(dat.obj, intercept=T, verbose=T)
{
  y=dat.obj$y
  X=dat.obj$X
  dist=dat.obj$dist
  offset=dat.obj$offset
    
  ID=1:length(y)
  
  tfit=glmmPQL(y~X-1+offset(offset), random=~1|ID, family=poisson, verbose = F)
  
  tab=summary(tfit)
  
  if(verbose)
  {
   print(tab)
  }
  
  res=as.numeric(tab$tTable[,"p-value"]<=.05)

  if(intercept)
  {
    res=res[-1]
  }
  
  return(res)

}  
################################################################################
#functions for backward selection, used in comparison
spatialVS.glmmPQL.fit.backward=function(dat.obj, intercept=T, verbose=T)
{
  y=dat.obj$y
  X=dat.obj$X
  dist=dat.obj$dist
  offset=dat.obj$offset
    
  ID=1:length(y)
    
  pp=ncol(X)
  var.id=1:pp
  
  
  keep.mat=matrix(0, pp+1, pp)
  keep.mat[1,]=1
  
  i=1
  
  nIter=pp
  if(intercept)
  {
    nIter=pp-1
  }
  
  vs.flag=T
  
  while(i<=nIter & vs.flag)
  {
    
    xmat=X[,keep.mat[i,]==1]
    cur.ids=var.id[keep.mat[i,]==1]
    
    i=i+1
    keep.mat[i,]=keep.mat[i-1,]
    
    tfit=glmmPQL(y~xmat-1+offset(offset), random=~1|ID, family=poisson, verbose = F)
    tab=summary(tfit)
    xres=tab$tTable[,"p-value"]
    
    if(intercept)
    {
      xres=xres[-1]
      cur.ids=cur.ids[-1]
    }
    
    if(any(xres>=0.05))
    {
       drop.id=cur.ids[xres==max(xres)]
       drop.id=drop.id[1]
       keep.mat[i, drop.id]=0
    
    }else{
           vs.flag=F
         }

  }

  res=keep.mat[i,]
  
  if(intercept)
  {
    res=res[-1]
  }
  
  return(res)
}  
#################################################################################
#function for glmmLasso, used in comparison
spatialVS.glmmLasso.fit=function(dat.obj, formula=y~x1+x2+x3+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+offset(offset), intercept=F, verbose=T)
{
  y=dat.obj$y
  X=dat.obj$X
  
  dist=dat.obj$dist
  offset=dat.obj$offset
  
  ID=1:length(y)
  
  tmp.dat=cbind(y, X, offset, ID)
  tmp.dat=as.data.frame(tmp.dat)
  tmp.dat[,"ID"]=as.factor(ID)
    
  pp=ncol(X)
  nn=length(ID)
  
  
  if(intercept)
  {
    clt=list(print.iter=verbose, start=c(1, rep(0, pp-1+nn)), q.start=0.7)

  }else{
         clt=list(print.iter=verbose, start=c(1, rep(0, nn+pp)), q.start=0.7, center=F)
       }
  
  PQL.fit<-try(glmmPQL(formula, random = ~1|ID, family = poisson(link = log),  data=tmp.dat, niter = 10, verbose = F), silent=T)
  
  if((class(PQL.fit)!="try-error")[1])
  { 
    Delta.start<-c(as.numeric(PQL.fit$coef$fixed),as.numeric(t(PQL.fit$coef$random$ID)))
    Q.start<-as.numeric(VarCorr(PQL.fit)[1,1])
    
    if(is.na(Q.start))
    {
      Q.start=0.7
    }
    
    clt=list(print.iter=verbose, start=Delta.start, q.start=Q.start, center=F)
  }
  
  #lambda.vec=seq(0.01, 100,, 100)
  lambda.vec=seq(500, 200,, 100)

  bic.vec=double(length(lambda.vec))
  
  for(i in 1:length(lambda.vec))
  {
    if(verbose)
    {
     cat("i=", i, "lambda=", lambda.vec[i], "\n")
    }
    

    tfit <- glmmLasso(formula, rnd = list(ID=~1), family = poisson(link = log), lambda=lambda.vec[i], data=tmp.dat, control=clt)
    bic.vec[i]=tfit$bic
    
    tmp.pp=length(as.vector(coef(tfit)))
    clt$start[1:tmp.pp]=as.vector(coef(tfit))

  }
  
  lamba.best=lambda.vec[bic.vec==min(bic.vec)]
  lamba.best=lamba.best[1]
  
  tfit <- glmmLasso(formula, rnd = list(ID=~1), family = poisson(link = log), lambda=lamba.best, data=tmp.dat, control=clt)
  
  res=as.numeric(round(coef(tfit),3)!=0)

  res=res[-1]
  
  return(res)
    
}
#################################################################################
#collect simulation results for Lyme disease covariate
spatialVS_sim_coll.realX=function(beta, name, method="PQL")
{  
  load(name)
    
  beta.res=theta.res=thres.res=NULL
  
  if (method=="LP")
  {
    for(i in 1:length(result))
    {
      #print(i)
      tmp=result[[i]]
      
      index=which.min(tmp$LPout.obj$bic)
      
      thres.res=rbind(thres.res, c(tmp$LP.obj$alpha.res[,index], tmp$LP.obj$lambda.res[,index]))
      beta.res=rbind(beta.res, tmp$LP.obj$beta.res[index,])
      theta.res=rbind(theta.res, tmp$LP.obj$ori.theta.res[index,])
      
     
    }
    
  }else{
      
      for(i in 1:length(result))
      {
      #print(i)
      tmp=result[[i]]
      
      index=which.min(tmp$PQLout.obj$bic)
      
      if(length(tmp$PQLout.obj$bic)==1)
      {
        index=1
      }
      
      thres.res=rbind(thres.res, c(tmp$PQL.obj$alpha.res[,index], tmp$PQL.obj$lambda.res[,index]))
      beta.res=rbind(beta.res, tmp$PQL.obj$beta.res[index,])
      theta.res=rbind(theta.res, tmp$PQL.obj$ori.theta.res[index,])
      
      
    }
  }
  
  # remove the not coverge cases
  id=(theta.res[,1]>exp(-5))&(theta.res[,1]<100)&
     (theta.res[,2]>exp(-5))&(theta.res[,2]<100)

  
  beta.res=beta.res[id, ]
  theta.res=theta.res[id, ]
  thres.res=thres.res[id, ]
  
  aa=abs(beta)
  idx=order(aa, decreasing=T)
  beta=aa[idx]
  beta.res=beta.res[ ,idx]
   
  p=length(beta)
  p1=sum(beta!=0)
  
  function1=function(x)
  {
    sum(x!=0)
  }
  
  function2=function(x)
  {
    sum(x[(p1+1):p]==0)
  }
  
  function3=function(x)
  {
    sum(x[1:p1]==0)
  }
  
  modelsize=apply(beta.res, 1, function1)
  corr0=apply(beta.res, 1, function2)
  inc0=apply(beta.res, 1, function3)
  
  return(c(mean(modelsize), mean(corr0), mean(inc0)) )
}
################################################################################ 
#collect simulation results for existing methods using Lyme covariate
spatialVS_sim_coll.ExM.realX=function(beta, name)
{  
  load(name)
  
  nn=length(result)
  
  pvalue.res.mat=NULL
  backward.res.mat=NULL
  glmmLasso.res.mat=NULL
  
  
  for(i in 1:nn)
  {
    tmp=result[[i]]
    pvalue.res.mat=rbind(pvalue.res.mat, tmp$pvalue.res)
    backward.res.mat=rbind(backward.res.mat, tmp$backward.res)
    glmmLasso.res.mat=rbind(glmmLasso.res.mat, tmp$glmmLasso.res)
    
  }
  
  beta=beta[-1]
  pvalue.res.mat=pvalue.res.mat[,-1]
  backward.res.mat=backward.res.mat[,-1]
  
  aa=abs(beta)
  idx=order(aa, decreasing=T)
  beta=aa[idx]
  
  pvalue.res.mat=pvalue.res.mat[ ,idx]
  backward.res.mat=backward.res.mat[ ,idx]
  glmmLasso.res.mat=glmmLasso.res.mat[ ,idx] 
   
  p=length(beta)
  p1=sum(beta!=0)
  


  function1=function(x)
  {
    sum(x!=0)
  }
  
  function2=function(x)
  {
    sum(x[(p1+1):p]==0)
  }
  
  function3=function(x)
  {
    sum(x[1:p1]==0)
  }
  
  sum.fun=function(mat)
  {
    if(!is.null(mat))
    {
     modelsize=apply(mat, 1, function1)
     corr0=apply(mat, 1, function2)
     inc0=apply(mat, 1, function3)
  
     tres=c(mean(modelsize), mean(corr0), mean(inc0))
    }else{
           tres=NULL
         }
         
    return(tres)
    
  }
  
  pvalue.sum=sum.fun(mat=pvalue.res.mat)
  back.sum=sum.fun(mat=backward.res.mat)
  glmmLasso.sum=sum.fun(mat=glmmLasso.res.mat)

  res=rbind(pvalue.sum, back.sum, glmmLasso.sum)
  
  return(res)
  
}
################################################################################
#print the simulation results into table
spatialVS_sim_coll.ExM.print.tab=function(file.name.ext, betaeff=c(-0.5, 0.75, 1, -0.75, -1), p.zero=10)
{
  case1.adp.res=spatialVS_sim_coll.ExM(name=paste("resultExM_case1_adaptive", file.name.ext, sep=""), betaeff=betaeff, p.zero=p.zero)
  case2.adp.res=spatialVS_sim_coll.ExM(name=paste("resultExM_case2_adaptive", file.name.ext, sep=""), betaeff=betaeff, p.zero=p.zero)
  case3.adp.res=spatialVS_sim_coll.ExM(name=paste("resultExM_case3_adaptive", file.name.ext, sep=""), betaeff=betaeff, p.zero=p.zero)
  case4.adp.res=spatialVS_sim_coll.ExM(name=paste("resultExM_case4_adaptive", file.name.ext, sep=""), betaeff=betaeff, p.zero=p.zero)
  case5.adp.res=spatialVS_sim_coll.ExM(name=paste("resultExM_case5_adaptive", file.name.ext, sep=""), betaeff=betaeff, p.zero=p.zero)

  xx=cbind(rbind(case1.adp.res[1, 1:3], case2.adp.res[1, 1:3],case3.adp.res[1, 1:3],case4.adp.res[1, 1:3],case5.adp.res[1, 1:3]),
         rbind(case1.adp.res[2, 1:3], case2.adp.res[2, 1:3],case3.adp.res[2, 1:3],case4.adp.res[2, 1:3],case5.adp.res[2, 1:3]),
         rbind(case1.adp.res[3, 1:3], case2.adp.res[3, 1:3],case3.adp.res[3, 1:3],case4.adp.res[3, 1:3],case5.adp.res[3, 1:3]))
  library(xtable)
  colnames(xx)=c("Modelsize","Correctnum", "Wrongnum", "Modelsize","Correctnum", "Wrongnum", "Modelsize","Correctnum", "Wrongnum")
  rownames(xx)=paste("Case", 1:5)
  print(xtable(round(xx, 2)),sanitize.text.function = function(x) { x })
}
################################################################################
#perform simulation study using real covariate from the Lyme data
spatialVS_simu_realX=function(mm, ll=15, numofcor=5, rho=0.8, betaeff=c(-0.5, 0.75, 1, -0.75, -1), cov.pars=c(0.1, 5), p.total=15, alpha.vec=seq(0.6, 1, by=0.1), lambda.vec=seq(0.01, 0.5, by=0.02), cc=0, adaptive=F, noeff_corr=F, nocor=F, lambda.vec.ext=NULL, realX=F, realX.file=NULL)
{  

  #simulate data
  if(!realX)
  {
  
    simu.data=spatialVS_simu_data(ll=ll, numofcor=numofcor, rho=rho, betaeff=betaeff,cov.pars=cov.pars, noeff_corr=noeff_corr, p.total=p.total)
  
  }
  
  if(realX)
  {
    load(file=realX.file)
    simu.data=spatialVS_simu_data_LymeX(dat.obj=lyme.eco0.boot.dat$dat.obj, fit.obj=lyme.eco0.boot.dat$fit.obj, subset.no=ll^2, intercept=F)  
  }
  
  ##################################################
  #print(lambda.vec)

  alpha.vec=seq(0.05, 1,,10)
  lambda.vec=seq(0.05, 7, len=20)
  
  tmp.PQL=try(SpatialVS(dat.obj=simu.data, alpha.vec=alpha.vec, lambda.vec=lambda.vec, method="PQL", plots=F, intercept=F, verbose=F), silent=T)
   
  PQL.obj=tmp.PQL$L.best.obj 
  PQLout.obj=tmp.PQL$Lout.best.obj
  start=tmp.PQL$start
  lasso.weight=tmp.PQL$lasso.weight
  
  
  if(nocor)
  {
    tmp.PQL.nocor=try(SpatialVS(dat.obj=simu.data, alpha.vec=alpha.vec, lambda.vec=lambda.vec, method="PQL.nocor", plots=F, intercept=F, verbose=F), silent=T)
   
    PQL_nocor.obj=tmp.PQL.nocor$L.best.obj 
    PQL_nocorout.obj=tmp.PQL.nocor$Lout.best.obj
    start.nocor=tmp.PQL.nocor$start
    lasso.weight.nocor=tmp.PQL.nocor$lasso.weight
     
    return(list(simu.data=simu.data, start=start, start.nocor=start.nocor, PQL.obj=PQL.obj, PQLout.obj=PQLout.obj, PQL_nocor.obj=PQL_nocor.obj, PQL_nocorout.obj=PQL_nocorout.obj, lasso.weight=lasso.weight, lasso.weight.nocor=lasso.weight.nocor))
  
  }

  return(list(simu.data=simu.data, start=start, PQL.obj=PQL.obj, PQLout.obj=PQLout.obj, lasso.weight=lasso.weight))

}
################################################################################
#perform simulation study for APL using covariate from real Lyme data
spatialVS_simu_LP_realX=function(mm, alpha.vec=seq(0.6, 1, by=0.1), lambda.vec=seq(0.1, 1.5, by=0.1), cc=0, adaptive=F, name, nocor=F)
{  
  alpha.vec=seq(0.05, 1,,10)
  lambda.vec=seq(0.05, 7, len=20)

  load(name)
  simu.data=result[[mm]]$simu.data
  start=result[[mm]]$start
  lasso.weight=result[[mm]]$lasso.weight
  
  if(is.null(start))
  {
    return(simu.data=simu.data)
  }


  LP.obj=try(SpatialVS_P_testLP(data=simu.data, start=start, alpha.vec=alpha.vec, lambda.vec=lambda.vec, cc=cc, adaptive=adaptive, lasso.weight=lasso.weight, intercept=F, verbose=F), silent=T)
  
  
  
  if(class(LP.obj)=="try-error")
  {
    return(list(simu.data=simu.data))
  }
  
  LPout.obj=try(SpatialVS.out(obj=LP.obj, data=simu.data, cc=0), silent=T)
  
  if(class(LPout.obj)=="try-error")
  {
    return(list(simu.data=simu.data))
  }
  
  
  if(nocor)
  {    
    
    start.nocor=result[[mm]]$start.nocor
    lasso.weight.nocor=result[[mm]]$lasso.weight.nocor
    
    if(is.null(start.nocor))
    {
      return(simu.data=simu.data)
    }
    
    LP_nocor.obj=try(SpatialVS_P_testLP.nocor(data=simu.data,start=start.nocor, alpha.vec=alpha.vec, lambda.vec=lambda.vec, cc=cc, adaptive=adaptive, lasso.weight=lasso.weight.nocor, intercept=F, verbose=F), silent=T)
    
    
    if(class(LP_nocor.obj)=="try-error")
    {
      return(list(simu.data=simu.data))
    }
    
    LP_nocorout.obj=try(SpatialVS.out.nocor(obj=LP_nocor.obj, data=simu.data, cc=0), silent=T)
    
    if(class(LP_nocorout.obj)=="try-error")
    {
      return(list(simu.data=simu.data))
    }
    
    return(list(simu.data=simu.data, start=start, LP.obj=LP.obj, LPout.obj=LPout.obj, LP_nocor.obj=LP_nocor.obj, LP_nocorout.obj=LP_nocorout.obj))  
  }
   
  return(list(simu.data=simu.data, start=start, LP.obj=LP.obj, LPout.obj=LPout.obj))
}
################################################################################
#print simulation results into table
spatialVS_sim_coll.APL.PQL.print.tab=function(file.name.ext, betaeff=c(-0.5, 0.75, 1, -0.75, -1), cov.pars=c(0.1, 5), p.zero=10)
{
  #######################################################################
  # adaptive elastic net, LP
  print("APL")
  case1.adpLP.res=spatialVS_sim_coll(name=paste("resultLP_case1_adaptive", file.name.ext, sep=""), pdfname="case1LP_adp_paras.pdf", method="LP",betaeff=betaeff, cov.pars=cov.pars, p.zero=p.zero)
  case2.adpLP.res=spatialVS_sim_coll(name=paste("resultLP_case2_adaptive", file.name.ext, sep=""), pdfname="case2LP_adp_paras.pdf", method="LP",betaeff=betaeff, cov.pars=cov.pars, p.zero=p.zero)
  case3.adpLP.res=spatialVS_sim_coll(name=paste("resultLP_case3_adaptive", file.name.ext, sep=""), pdfname="case3LP_adp_paras.pdf", method="LP",betaeff=betaeff, cov.pars=cov.pars, p.zero=p.zero)
  case4.adpLP.res=spatialVS_sim_coll(name=paste("resultLP_case4_adaptive", file.name.ext, sep=""), pdfname="case4LP_adp_paras.pdf", method="LP",betaeff=betaeff, cov.pars=cov.pars, p.zero=p.zero)
  case5.adpLP.res=spatialVS_sim_coll(name=paste("resultLP_case5_adaptive", file.name.ext, sep=""), pdfname="case5LP_adp_paras.pdf", method="LP",betaeff=betaeff, cov.pars=cov.pars, p.zero=p.zero)

  # consider/not consider spatial correlation
   xx=cbind(rbind(case1.adpLP.res[1, 1:3], case2.adpLP.res[1, 1:3],case3.adpLP.res[1, 1:3],case4.adpLP.res[1, 1:3],case5.adpLP.res[1, 1:3]), rbind(case1.adpLP.res[2, 1:3], case2.adpLP.res[2, 1:3],case3.adpLP.res[2, 1:3],case4.adpLP.res[2, 1:3],case5.adpLP.res[2, 1:3]))
   library(xtable)
   colnames(xx)=c("Modelsize","Correctnum", "Wrongnum", "Modelsize","Correctnum", "Wrongnum")
   rownames(xx)=paste("&Case", 1:5)
   print(xtable(round(xx, 2)),sanitize.text.function = function(x) { x })

   print("PQL")
   
   case1.adp.res=spatialVS_sim_coll(name=paste("result_case1_adaptive", file.name.ext, sep=""), pdfname="case1_adp_paras.pdf",betaeff=betaeff, cov.pars=cov.pars, p.zero=p.zero)
   case2.adp.res=spatialVS_sim_coll(name=paste("result_case2_adaptive", file.name.ext, sep=""), pdfname="case2_adp_paras.pdf",betaeff=betaeff, cov.pars=cov.pars, p.zero=p.zero)
   case3.adp.res=spatialVS_sim_coll(name=paste("result_case3_adaptive", file.name.ext, sep=""), pdfname="case3_adp_paras.pdf",betaeff=betaeff, cov.pars=cov.pars, p.zero=p.zero)
   case4.adp.res=spatialVS_sim_coll(name=paste("result_case4_adaptive", file.name.ext, sep=""), pdfname="case4_adp_paras.pdf",betaeff=betaeff, cov.pars=cov.pars, p.zero=p.zero)
   case5.adp.res=spatialVS_sim_coll(name=paste("result_case5_adaptive", file.name.ext, sep=""), pdfname="case5_adp_paras.pdf",betaeff=betaeff, cov.pars=cov.pars, p.zero=p.zero)

   # consider/not consider spatial correlation
   xx=cbind(rbind(case1.adp.res[1, 1:3], case2.adp.res[1, 1:3],case3.adp.res[1, 1:3],case4.adp.res[1, 1:3],case5.adp.res[1, 1:3]), rbind(case1.adp.res[2, 1:3], case2.adp.res[2, 1:3],case3.adp.res[2, 1:3],case4.adp.res[2, 1:3],case5.adp.res[2, 1:3]))
   library(xtable)
   colnames(xx)=c("Modelsize","Correctnum", "Wrongnum", "Modelsize","Correctnum", "Wrongnum")
   rownames(xx)=paste("&Case", 1:5)
   print(xtable(round(xx, 2)),sanitize.text.function = function(x) { x })

}
################################################################################
#print simulation results into table
spatialVS_sim_coll.realX.print.tab=function()
{
  case1.adp.res=spatialVS_sim_coll.realX(beta=c(0.000, -0.161, 0.009, 0.000, 0.000, 0.000, 0.064, 0.000, 0.503, 0.000, 0.000,-0.048, 0.000, 0.185), name="result_realX_adaptive", method="PQL")

  case2.adp.res=spatialVS_sim_coll.realX(beta=c(0.000, -0.161, 0.009, 0.000, 0.000, 0.000, 0.064, 0.000, 0.503, 0.000, 0.000,-0.048, 0.000, 0.185), name="resultLP_realX_adaptive", method="LP")

  case3.adp.res=spatialVS_sim_coll.ExM.realX(beta=c(0.000, -0.161, 0.009, 0.000, 0.000, 0.000, 0.064, 0.000, 0.503, 0.000, 0.000,-0.048, 0.000, 0.185), name="resultExM_realX_adaptive")

  xx=rbind(case2.adp.res, case1.adp.res, case3.adp.res)
  library(xtable)
  colnames(xx)=c("Modelsize","Correctnum", "Wrongnum")
  rownames(xx)=c("APL", "PQL", "Pvalue", "Back", "glmmLasso")
  print(xtable(round(xx, 2)),sanitize.text.function = function(x) { x })

}
################################################################################
#Used by SpatialVS in the estimation of covariance parameters under spatial correlation
cov.fun.fixsigma=function(par, sigma.ini, dist, glmmPQL.ini, beta.ini,  b.ini, y, X)
{

  D=covStruct(theta=c(log(sigma.ini^2), par),cc=0, dist=dist)

  mu=exp(glmmPQL.ini$fitted[,"fixed"]+b.ini)

  V=diag(1/mu)+D

  Vinv=solve(V)

  ytrans=b.ini+(y-mu)/mu

  return(as.numeric(determinant(V)$modulus+determinant(t(X)%*%Vinv%*%X)$modulus+t(ytrans-X%*%beta.ini)%*%Vinv%*%(ytrans-X%*%beta.ini)))
}
################################################################################
#intermediate function, performs iterated weighted least squares
iwls_b_new=function(X, y, D, offset, beta, b, tol, dist)    #penalty for zero sum
{
  #browser()

  eta = X%*%beta + b + offset
  mu = exp(eta)
  convEta=10
  solveD=solve(D)

  nn=length(b)
  zmat=matrix(1, ncol=nn, nrow=nn)
  #b=rep(0, nn)

  nIter=0

  while(convEta>tol)
  {
    ##define the working response and weight

    nIter=nIter+1

    mu=as.vector(mu)
    fir_devb = (mu-y) + solveD%*%b+zmat%*%b
    sec_devb = diag(mu) + solveD+ zmat

    ###b
    #browser()
    b.new=b-solve(sec_devb)%*%fir_devb
    eta.new=X%*%beta + b.new + offset

    ###covergence
    convEta <- sum((eta.new-eta)^2) /sum(eta^2)

    #cat("convEta=", convEta, "nIter=", nIter, "mean.b=", mean(b.new), "\n")

    ###updata the eta and mu
    eta = eta.new
    #browser()

    mu = exp(eta)
    b = b.new

    #print(convEta)
  }

  #cat("\n")

  return(list(beta=as.numeric(beta), b=as.numeric(b), eta=as.numeric(eta), mu=as.numeric(mu)))

}
################################################################################
#plot correlations in real data
cov.corr.plot=function(dat0=lyme.svs.eco0.dat, dat1=lyme.svs.eco1.dat)
{
  X0=dat0$X
  X1=dat1$X
  xmat=rbind(X0, X1)
  xmat=xmat[,-1]

  xmat0=X0[,-1]
  xmat1=X1[,-1]

  cmat0=cor(xmat0)
  cmat1=cor(xmat1)

  cmat=cor(xmat)


  print(round(cmat,2))

  rr=cmat[row(cmat)>col(cmat)]

  rr1=abs(rr)

  hist(rr1, xlab="Absolute Values of Pairwise Correlations", main="", col="grey")

  print(round(range(rr1), 3))

  #browser()

}
################################################################################
#plot observed counts vs expected counts from the model
PQL.exp.vs.obs.plot=function(obj)
{
  dat=obj$dat.obj
  obj=obj$fit.obj
    
  y=dat$y
  X=dat$X
  
  offset=dat$offset
  dist=dat$dist
  
  beta=obj$beta.res
  b=obj$b.res
  
  beta=as.vector(beta)
  b=as.vector(b)
 
  tmp1=X%*%beta+b+offset
  mu=exp(tmp1)

  res=cbind(y, mu)
  return(invisible(res)) 
}

################################################################################
#plot observed counts vs expected counts from the model
overall.exp.vs.obs=function(obj0, obj1)
{
  a0=PQL.exp.vs.obs.plot(obj=obj0)
  a1=PQL.exp.vs.obs.plot(obj=obj1)
  aa=rbind(a0, a1)

  y=aa[,1]
  mu=aa[,2]
  
  plot(y, mu, pch=16, xlim=c(0, 150), ylim=c(0, 150), xlab="Observed Counts", ylab="Estimated Counts")
  abline(0,1, lwd=2, col=2)

  print(1-sum((y-mu)^2)/sum((y-mean(y))^2))
  
}

################################################################################







 