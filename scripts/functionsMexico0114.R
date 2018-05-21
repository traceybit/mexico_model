################################################
## Pella-Tomlinson Recovery Projections Model for Mexico
## Functions file
## Gavin McDonald (gmcdonald@bren.ucsb.edu)
## Sustainable Fisheries Group / EcoAnalytics
## May 6, 2015
################################################

## Load packages
library(truncnorm)

################################################
### Implementation vector
delayVec <- c(2:20)

################################################
## Define functions

################################################
## g_legend function
## Grab legend from ggplot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

################################################
## bioModel function
## Biological model in terms of reference points b and f
## b defined as B/BMSY and f defined as F/FMSY
## Uses Pella-Tomlinson model

bioModel = function(b,phi,g,f_intervention,f_nonintervention)
{
  b_next = b + ((phi + 1) / phi ) * g * b * (1 -  (b ^ phi)  / (phi + 1)) - g * f_intervention * b - g * f_nonintervention * b
  bmax = (phi+1)^(1/phi) - 0.1
  return(max(0,min(bmax,b_next)))
}

################################################
## MSY function
## Calculates MSY using life history information
MSY = function(g,K,phi)
{
  msy = g * K / (phi + 1) ^ (1 / phi)
  return(msy)
}

################################################
## BMSY function
## Calculates BMSY using life history information
BMSY = function(K,phi)
{
  bmsy = K / (phi + 1) ^ (1 / phi)
  return(bmsy)
}

################################################
## econModel function
## Economic model in terms of reference points b and f
## b defined as B/BMSY and f defined as F/FMSY

econModel = function(g,K,phi,p,f,b,c,beta)
{
  revenue = p * f * b * MSY(g,K,phi)
  cost = c * (g * f) ^ beta
  pi =  revenue - cost
  return(list(pi=pi,
              revenue=revenue,
              cost=cost))
}

################################################
## profitMSY function
## Determine profit at MSY and BMSY

profitMSY = function(g,K,phi,p,b,c,beta)
{
  revenue = p * 1 * 1 * MSY(g,K,phi)
  cost = c * (g * 1) ^ beta
  pi =  revenue - cost
  return(pi)
}

################################################
## parameterPull function
## Pulls random value from truncated normal distribution

parameterPull = function(upper,lower,expected)
{
  if (lower == expected & expected == upper) parameter = expected else parameter = rtruncnorm(1, a=lower, b=upper, mean=expected, sd = (upper - lower)/4)
  return(parameter)
}

################################################
## recoveryTime function
## Calculates time to recovery (first time period when vector value of bProjectionVec > recoveryCutoff)
## bProjectionVec defined in terms of B/BMSY and recoveryCutoff defined in terms of B/BMSY

recoveryTime = function(bProjectionVec,recoveryCutoff,t)
{
  if (max(bProjectionVec) < recoveryCutoff) time = NA else time = min(which(bProjectionVec >= recoveryCutoff))
  return(time)
}

################################################
## profitOptim function
## Internal Optimization Function 
## To be used with dynamicPolicy function
## Gives (negative) value function value for value function iteration code

profitOptim = function(f_intervention,f_nonintervention,b,p,K,c,g,phi,beta,V,bvec,delta)
{  
  
  f_intervention = max(0,f_intervention)
  
  f_nonintervention = max(0,f_nonintervention)
  
  profit = econModel(g,K,phi,p,f_intervention,b,c,beta)$pi
  
  bnext = bioModel(b,phi,g,f_intervention,f_nonintervention)
  
  if (bnext < bvec[1]) Vnext = 0 else Vnext = approx(bvec,V,bnext)$y #spline(bvec,V,xout=bnext)
  
  negout= -(profit + delta*Vnext)

  if(!is.numeric(negout)) browser()
  if(is.na(negout)) browser()

  return(list(negout=negout))
}

################################################
## dynamicPolicy function
## Run Dynamic Optimization
## Solves for optimal policy function f (as function of bvec) given model parameters

dynamicPolicy = function(K,g,phi,p,c,beta,disc,bvec,f_nonintervention)
{
  
  tol=.01 ## Convergence tolerance
  
  delta= 1/(1+disc) ## Discount parameter
  t=0
  
  f1= matrix(1,length(bvec),1)
  Vnew= matrix(0,length(bvec),1)
  diff= 10*tol
  
  while (t<4 | diff>tol)
  {
    t= t+1
    V= Vnew
    oldf1= f1
    
    for (i in 1:length(bvec))
    {
      b= bvec[i]
      if(i==1)
      {guess= 1}
      else
      {guess= f1[i-1]}

      FishOut= optim(par=guess,fn=profitOptim,gr=NULL,lower=0,upper=((phi+1)/phi-0.01-f_nonintervention),f_nonintervention=f_nonintervention,b=b,p=p,K=K,c=c,g=g,phi=phi,beta=beta,V=V,bvec=bvec,delta=delta,method="L-BFGS-B")
          
      Vnew[i]= -FishOut$value
      f1[i]= FishOut$par
      
      
    } ## Close bvec loop
    
    diff= sum(abs(f1-oldf1))
    #print(diff)
    
  }## Close while loop
  
  
  return(list(f1=f1))
  
} ## Close function

################################################
## policy function
## Generates four policy functions, one for each scenario

policy = function(s,g,K,phi,p,f0Int,f0NonInt,c,beta,disc,bvec)
{

 ## Policy vector for status quo; maintain f0_Intervention
 if (s==1) f1 = rep(f0Int,length(bvec))
 
 ## Policy vector that maximimzes food production; set at f=1 for all conditions, including f_NonIntervention
 if (s==2) f1 = rep(max(0,1-f0NonInt),length(bvec))
 
 ## Policy vector that minimizes recovery time; set at 0 until b=1 then set at f=1, including f_NonIntervention
 zeros = which(bvec<1)
 ones = which(bvec>=1)
 if (s==3) f1 = c(rep(0,length(zeros)),rep(max(0,1-f0NonInt),length(ones)))
 
 ## Policy vector that dynamically maximizes NPV
 if (s==4) f1 = dynamicPolicy(K,g,phi,p,c,beta,disc,bvec,f0NonInt)$f1
 
 ## Policy vector for closing the fishery; f0_Intervention = 0
 if (s==5) f1 = 0
 
 ## Policy vector for open access; maintain f0_Intervention. 
 ## Mortality will be modified in projection loop depending on lambda
 if (s==6) f1 = rep(f0Int,length(bvec))
 
 return(list(f1=f1))
}

################################################
## NPV function
## Calculate NPV given a discount rate and vector of profits

NPV = function(discount,P)
{
  pv=vector()
  for (i in 1:length(P))
  {
    pv[i] = P[i] / (1+discount)^i
  }
  return(sum(pv))
}

################################################
## annuity function
## Calculate annuity given a discount rate and vector of profits

annuityFunc = function(discount,P)
{
  n = length(P)
  if (discount == 0) annuity = NPV(discount,P)/n else{
    annuity = NPV(discount,P) / ((1 - 1/(1+discount)^n)/discount)
  }
  return(annuity)
}

################################################
## breakEvenNPV function
## Calculates NPV break-even point [year]

breakEvenNPV = function(discount,P)
{
  for (i in 1:length(P))
  {
    npvCurrent = NPV(discount,P[1:i])
    if (npvCurrent >= 0) pointer = i-1
    if (npvCurrent >= 0) break
    if (i==length(P) & npvCurrent <= 0) pointer=NA
  }
  return(pointer)
}

################################################
## projectionModel

projectionModel = function(params,S,Thetas,CatchShareLoop,LegalLoop)
{
  ##Set up loops
  N = params$N
  T = params$T
  if (CatchShareLoop == "yes") L = 2 else L = 1
  if (LegalLoop == "yes") M = 2 else M = 1
  
  ## Set recovery definition cutoff in terms of B/BMSY
  cutoff = 0.8
  
  maxb = (params$phi_expected+1)^(1/params$phi_expected)
  
  ## initialize vector that contains range of b values from 0 to maxbN
  bVEC = seq(0,maxb+.1,.1) 
  
  ## initialize policies array
  policies = array(NA,dim=c(length(S),N,Thetas,L,M,length(bVEC)))
  
  ## initialize biomass reference point projection array
  bProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize absolute biomass projection array
  BProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize intervention fishing mortality reference point projection array
  fIntProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize non-intervention fishing mortality reference point projection array
  fNonIntProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize total fishing mortality reference point projection array
  fTotalProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize intervention absolute harvest projection array
  HIntProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize non-intervention absolute harvest projection array
  HNonIntProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize total absolute harvest projection array
  HTotalProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize revenue projection array
  revenueProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize cost projection array
  costProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize profit projection array
  profitProjections = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec),T))
  
  ## initialize time to biological recovery matrix
  timeToRecovery = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec)))
  
  ## initialize NPV matrix
  npv = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec)))
  
  ## initialize annuity matrix
  annuity = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec)))
  
  ## initialize NPV break-even point matrix
  breakEven = array(NA,dim=c(length(S),N,Thetas,L,M,length(delayVec)))
  
  ## Loop over all delay scenarios
  for (h in 1:length(delayVec)){  ## added this loop
  
    ## Loop over all policy scenarios
    for (i in 1:length(S)){
    
      ## Loop over all Monte-Carlo iterations
      for (j in 1:N){
      
        ## Loop over legal theta values
        for (k in 1:Thetas){
        
          ## Loop over catch share price and cost scalars
          for (l in 1:L){
          
            ## Loop over illegal harvest elimination
            for (m in 1:M){
            
              ## For the first monte-carlo iteration, use the expected values of each parameter
              if (j == 1){
              discN = params$disc_expected
              gN = params$g_expected
              phiN = params$phi_expected
              maxFN = (phiN+1)/phiN
              KN = params$K_expected
              b0N = min(maxb,params$b0_expected)
              f0_totalN = min(maxFN,params$f0_total_expected)
              thetaN_legal = params$theta_legal_expected
              thetaN_mexico = params$theta_mexico
              pN = params$p_expected  
              betaN = params$beta_expected	
              cN = params$c_expected
              gamma_pN = params$gamma_p
              gamma_cN = params$gamma_c
              lambdaN = params$lambda
              
              ## For all other iterations, pull a parameter from a truncated normal distribution
              ## defined by the upper bound, lower bound, and expected mean
            } else {
              discN = parameterPull(params$disc_upper,params$disc_lower,params$disc_expected)
              gN = parameterPull(params$g_upper,params$g_lower,params$g_expected)
              phiN = parameterPull(params$phi_upper,params$phi_lower,params$phi_expected)
              maxFN = (phiN+1)/phiN
              KN = parameterPull(params$K_upper,params$K_lower,params$K_expected)
              b0N = min(maxb,parameterPull(params$b0_upper,params$b0_lower,params$b0_expected))
              f0_totalN = min(maxFN,parameterPull(params$f0_total_upper,params$f0_total_lower,params$f0_total_expected))
              thetaN_legal = parameterPull(params$theta_legal_upper,params$theta_legal_lower,params$theta_legal_expected)
              thetaN_mexico = params$theta_mexico
              pN = parameterPull(params$p_upper,params$p_lower,params$p_expected)
              betaN = parameterPull(params$beta_upper,params$beta_lower,params$beta_expected)
              cN = parameterPull(params$c_upper,params$c_lower,params$c_expected)
              gamma_pN = params$gamma_p
              gamma_cN = params$gamma_c
              lambdaN = params$lambda
            }
            
            ## Set policies for current Monte-Carlo iteration parameters - assume perfect information
            f0IntN = f0_totalN * thetaN_legal * thetaN_mexico
            
            f0NonIntN = f0_totalN * (1 - thetaN_legal * thetaN_mexico)
            
            policies[i,j,k,l,m,] = policy(S[i],gN,KN,phiN,pN,f0IntN,f0NonIntN,cN,betaN,discN,bVEC)$f1 
            
            ## Set non-intervention f0 depending on illegal harvest elimination loop
            if (m == 1) f0NonIntN = f0_totalN * (1 - thetaN_legal * thetaN_mexico) 
            if (m == 2) f0NonIntN = f0_totalN * (1 - thetaN_mexico)
            
            ## Set price and cost depending on catch share loop
            if (l == 2) { 
                pN = pN * gamma_pN
                cN = cN * gamma_cN
            }
            
            # policies[i,j,k,l,m,] = policy(S[i],gN,KN,phiN,pN,f0IntN,f0NonIntN,cN,betaN,discN,bVEC)$f1
            
            ## Loop over all time steps
            for (n in 1:T){
              
              #print(paste("Scenario ",i,"of",length(S),"; Monte Carlo Iteration",j,"of",N,"; Theta", k,"of",Thetas,"; Catch Share Loop",l,"of",L,"; Illegal Elimination Loop",m,"of",M,"; Time Step",n,"of",T,sep=" "))
              
                if (n == 1) {
                    bProjections[i,j,k,l,m,h,n] = b0N
                    fIntProjections[i,j,k,l,m,h,n] = f0IntN
                    fNonIntProjections[i,j,k,l,m,h,n] = f0_totalN * (1 - thetaN_legal * thetaN_mexico)
                    pN = params$p_expected  
                    cN = params$c_expected
                } else {
                    if (l == 2 & n >= delayVec[h]) {  ## added code here
                        pN = params$p_expected * gamma_pN
                        cN = params$c_expected * gamma_cN
                    }
                    if (m == 1) fNonIntProjections[i,j,k,l,m,h,n] = f0_totalN * (1 - thetaN_legal * thetaN_mexico) 
                    if (m == 2 & n < delayVec[h]) fNonIntProjections[i,j,k,l,m,h,n] = f0_totalN * (1 - thetaN_legal * thetaN_mexico)
                    if (m == 2 & n >= delayVec[h]) fNonIntProjections[i,j,k,l,m,h,n] = f0_totalN * (1 - thetaN_mexico)
                    ## Open access policy scenario
                    if (S[i] == 6){
                        if (n >= delayVec[h]) { # Only implement improved policy if implementation year is reached
                            fIntProjections[i,j,k,l,m,h,n] = max(0,min(maxFN,fIntProjections[i,j,k,l,m,h,n-1] + 
                                                                           lambdaN * 
                                                                           profitProjections[i,j,k,l,m,h,n-1] / profitMSY(gN,KN,phiN,pN,bProjections[i,j,k,l,m,h,n-1],cN,betaN)))
                        }
                    } else fIntProjections[i,j,k,l,m,h,n] = f0IntN
                    ## All other policy scenarios
                    if (S[i] < 6) {
                        if (n >= delayVec[h]) { # Only implement improved policy if implementation year is reached
                            if (bProjections[i,j,k,l,m,h,n] < bVEC[1]) fIntProjections[i,j,k,l,m,h,n] = approx(bVEC,policies[i,j,k,l,m,],bVEC[1])$y else{
                                fIntProjections[i,j,k,l,m,h,n] = approx(bVEC,policies[i,j,k,l,m,],bProjections[i,j,k,l,m,h,n])$y
                            }
                        } else fIntProjections[i,j,k,l,m,h,n] = f0IntN
                    }
                }
            
              
              fTotalProjections[i,j,k,l,m,h,n] = fIntProjections[i,j,k,l,m,h,n] + fNonIntProjections[i,j,k,l,m,h,n]
              
              BProjections[i,j,k,l,m,h,n] = bProjections[i,j,k,l,m,h,n] * BMSY(KN,phiN)
              
              HIntProjections[i,j,k,l,m,h,n] = fIntProjections[i,j,k,l,m,h,n] * gN * bProjections[i,j,k,l,m,h,n] * BMSY(KN,phiN)
              
              HNonIntProjections[i,j,k,l,m,h,n] = fNonIntProjections[i,j,k,l,m,h,n] * gN * bProjections[i,j,k,l,m,h,n] * BMSY(KN,phiN)
              
              HTotalProjections[i,j,k,l,m,h,n] = HIntProjections[i,j,k,l,m,h,n] + HNonIntProjections[i,j,k,l,m,h,n]
              
              revenueProjections[i,j,k,l,m,h,n] = econModel(gN,KN,phiN,pN,fIntProjections[i,j,k,l,m,h,n],bProjections[i,j,k,l,m,h,n],cN,betaN)$revenue
              
              costProjections[i,j,k,l,m,h,n] = econModel(gN,KN,phiN,pN,fIntProjections[i,j,k,l,m,h,n],bProjections[i,j,k,l,m,h,n],cN,betaN)$cost
              
              profitProjections[i,j,k,l,m,h,n] = econModel(gN,KN,phiN,pN,fIntProjections[i,j,k,l,m,h,n],bProjections[i,j,k,l,m,h,n],cN,betaN)$pi
              
              if (n < T) bProjections[i,j,k,l,m,h,n+1] = bioModel(bProjections[i,j,k,l,m,h,n],phiN,gN,fIntProjections[i,j,k,l,m,h,n],fNonIntProjections[i,j,k,l,m,h,n])
              
            } ## End loop over time steps
            
            timeToRecovery[i,j,k,l,m,h] = recoveryTime(bProjections[i,j,k,l,m,h,],cutoff,T)
            
            npv[i,j,k,l,m,h] = NPV(discN,profitProjections[i,j,k,l,m,h,])
            
            annuity[i,j,k,l,m,h] = annuityFunc(discN,profitProjections[i,j,k,l,m,h,])
            
            breakEven[i,j,k,l,m,h] = breakEvenNPV(discN,profitProjections[i,j,k,l,m,h,])
      
          } ## End loop over illegal harvest
            
        } ## End loop over catch schares
            
      } ## End loop over thetas
      
    } ## End loop over Monte-Carlo iterations
    
  } ## End loop over policy scenarios
  
  } ## End loop over delay vec
  
  return(list(bVEC=bVEC,
              policies=policies,
              bProjections=bProjections,
              BProjections=BProjections,
              fIntProjections=fIntProjections,
              fNonIntProjections=fNonIntProjections,
              fTotalProjections=fTotalProjections,
              HIntProjections=HIntProjections,
              HNonIntProjections=HNonIntProjections,
              HTotalProjections=HTotalProjections,
              revenueProjections=revenueProjections,
              costProjections=costProjections,
              profitProjections=profitProjections,
              timeToRecovery=timeToRecovery,
              npv=npv,
              breakEven=breakEven,
              cutoff=cutoff,
              annuity=annuity))
  
}