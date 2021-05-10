##########################################################################
##########################################################################
# Project:
# Script purpose: 
# Usage example: 
# original code from the paper Economou et al., 2020
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon May 10 13:42:50 2021
##########################################################################
##########################################################################
##########################################################################
##########################################################################
# Project: organoid patterning 
# Script purpose: functions for RD of 2 components
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon May 10 15:59:40 2021
##########################################################################
##########################################################################

# Functions and scripts to:
# 1). Generate parameterisations of Jacobian matrix which satisfy the criteria for diffusion driven instability 
# (as outlined in supplementary note) for a 2-component RD system
# 2). Convert parameterisations of Jacobian matrix to the form of RD system used in analyses by including degradation rates, 
# constant production rates and upper and lower reaction rate limits, and scale rates to form a spatial waves in the simulation grid
# 3). Run simulations of inhibiting each component in RD system (this is the perturbation response)

# Functions are equivalent to those for 3-component simulations (see 3-compenent simulations for details)

#For simulating inhibitions there are two functions: multi_perturb_2c which inhibits both species a and b at the level of response, and multi_perturb_2c_type_2 which inhibits both species at the level of production

#####################

#For running simulations

library(deSolve) #Version: 1.13	Depends: R (>= 2.15.0)
library(plyr) #Version: 1.8.4	Depends: R (>= 3.1.0)
library(doMC) #Version: 1.3.4	Depends: R (>= 2.14.0), foreach(>= 1.2.0), iterators(>= 1.0.0), parallel
registerDoMC(cores=6)

#######################
#Functions

#Unlike for 3-component system (or any system with 2> components) only statable (non-oscillating) waves can for for a 2-component system, therefore all parameter sets are type 1 (see 3-component simulations for details)

generate_parameters_2c<-function(x,top){
  
  par<-runif(4,0,1) 
  
  diff<-runif(2,1,10000)
  
  Da<-500/diff[1]
  Aa<-par[1]*top[1,1]
  Ab<-par[2]*top[1,2]
  Db<-500/diff[2]
  Ba<-par[3]*top[2,1]
  Bb<-par[4]*top[2,2]
  
  Y<-rep(NA,12)
  
  
  Tr<-Aa+Bb
  
  Det<-(Aa*Bb)-(Ab*Ba)
  
  
  if((Tr<0)&(Det>0)){
    
    D1<-(Da*Bb)+(Db*Aa)
    
    D2<-2*sqrt(Da*Db*(Det))
    
    if((D1>D2)&(D2>0)){
      
      k2_min<-(D1)/(2*Da*Db)
      
      jacmat<-matrix(c(Aa-Da*k2_min,Ab,Ba,Bb-Db*k2_min),byrow=TRUE,ncol=2)
      
      eig<-eigen(jacmat)
      
      maxeig<-(1:2)[Re(eig$values)==max(Re(eig$values))]
      
      eigvec<-Re(eig$vectors[,maxeig])
      
      phase<-sign(c(eigvec[1]/eigvec[2]))
      
      Y<-c(1,Da,Aa,Ab,Db,Ba,Bb,phase,k2_min,max(Re(eig$values)),abs(eigvec))
      
    }
    
  }
  
  return(Y)	
  
}



rdmod2c<-function(time,state,parms){
  with(as.list(parms),{
    
    n<-length(state)/2
    
    a<-state[1:n]
    b<-state[(n+1):(2*n)]	
    
    z<-seq(0,1,length=n)
    dz<-1/(n-1)	
    L<-L0+S*time
    dL<-S
    
    Th1<-G0*GZ
    Th2<-G0
    
    Ph1<-Th1/L
    Ph2<-Th2/L
    
    #diffusion
    Diffa<-(Da/(L^2))*diff(diff(c(a[1],a,a[n])))/(dz^2)
    Diffb<-(Db/(L^2))*diff(diff(c(b[1],b,b[n])))/(dz^2)
    
    #reaction
    F<-Caa*a+Cab*b+Ba
    G<-Cba*a+Cbb*b+Bb
    
    if(bounds==TRUE){
      F[F>maxA]<-maxA
      F[F<0]<-0
      G[G>maxB]<-maxB
      G[G<0]<-0		
    }
    
    
    #flow
    A<-numeric(length=n)
    A[z<Ph1]<-0
    A[(Ph1<=z)&(z<Ph2)]<-S*(z[(Ph1<=z)&(z<Ph2)]*L-Th1)/(Th2-Th1)
    A[Ph2<=z]<-S
    
    #advection
    Adva<-((z*dL-A)/L)*diff(c(a[1],a,a[n]),2)/(2*dz)
    Advb<-((z*dL-A)/L)*diff(c(b[1],b,b[n]),2)/(2*dz)
    
    #dilution due to growth
    Sdil<-numeric(length=n)
    Sdil[z<Ph1]<-0
    Sdil[(Ph1<=z)&(z<Ph2)]<-S/(Th2-Th1)
    Sdil[Ph2<=z]<-0
    
    
    #rate of change	
    da<-F+Diffa-Ra*a-Sdil*a+Adva
    db<-G+Diffb-Rb*b-Sdil*b+Advb
    
    return (list(c(da,db)))	
  })
  
}



parameters2c<-function(par,S,L0,G0,GZ,bounds){
  parset<-with(as.list(par),{
    
    #first add degradation rates to self interactions
    if(Aa>0){
      Ra<-sample(seq(0.01,1,by=0.01),1)
      Caa<-Aa+Ra
    }else{
      Ra<--Aa
      Caa<-0
    }	
    if(Bb>0){
      Rb<-sample(seq(0.01,1,by=0.01),1)
      Cbb<-Bb+Rb
    }else{
      Rb<--Bb
      Cbb<-0
    }	
    
    #set equilibrium levels to (1,1) and calculate background rates
    X<-matrix(c(Aa,Ab,Ba,Bb),byrow=TRUE,ncol=2)
    Xinv<-solve(X)
    eqm<-c(1,1)
    bgX<-solve(Xinv,eqm)
    bg<-c(bgX,eqm)
    
    #set max rates at random between 1.5 and 3 fold equilibrium reaction rate
    lim<-runif(1,-1,1)
    maxA<-(1+2^lim)*Ra*bg[3]
    maxB<-(1+2^lim)*Rb*bg[4]		
    
    
    #add parameters to vector, with zeros for b/g production and max
    parms<-c(Da=Da,
             Caa=Caa,
             Cab=Ab,
             Ba=-bg[1],
             Ra=Ra,
             maxA=maxA,
             Db=Db,
             Cba=Ba,
             Cbb=Cbb,
             Bb=-bg[2],
             Rb=Rb,
             maxB=maxB,
             S=S,
             L0=L0,
             G0=G0,
             GZ=GZ,
             bounds=bounds)
    
    return(parms)
  })
  parset
}




pattern_search_2c<-function(param,target_k2,target_l){
  
  test_k2<-as.numeric(param["k2"])
  
  test_l<-as.numeric(param["l"])
  
  scale_k2<-(sqrt(target_k2)/sqrt(test_k2))^2 #remember, these are squares!
  
  scale_l<-target_l/test_l
  
  parms<-parameters2c(param,0,50,50,0.5,FALSE)
  
  #prelimimarily scale according to tp k2 so si is in right area
  scale_parms<-parms*scale_l
  scale_parms["L0"]<-scale_parms["L0"]/scale_l #this is not a rate so don't scale with above
  
  scale_parms["Da"]<-scale_parms["Da"]/scale_k2
  scale_parms["Db"]<-scale_parms["Db"]/scale_k2
  
  
  #run a preliminary simulation from which to calculate the period
  
  parms<-scale_parms
  
  #equilibrium initial levels
  X<-matrix(c((as.list(parms)$Caa-as.list(parms)$Ra),as.list(parms)$Cab,as.list(parms)$Cba,(as.list(parms)$Cbb-as.list(parms)$Rb)),byrow=TRUE,ncol=2)
  Y<-c(-as.list(parms)$Ba,-as.list(parms)$Bb)
  
  Z<-solve(X,Y)
  
  n<-250
  
  a<-rep(Z[1],n)
  b<-rep(Z[2],n)
  
  #as doing unbounded to measure period initiate with perturbation at one end
  a[1]<-a[1]+.1
  
  ## Initial conditions
  state    <- c(a,b)
  
  ## RUNNING the model:
  
  times  <- seq(0, 2500, by = 10)   # output wanted at these time intervals
  
  out <- ode.1D(y = state, times = times, func = rdmod2c, parms = parms,nspec = 2)
  
  
  #measure period
  
  L0<-as.list(parms)$L0
  
  #first check enough waves (ie >=3) to measure, if not scale diffusion by 2/3
  
  wavesp<-length((1:250)[diff(sign(diff(out[251,2:251])))==2])
  wavesm<-length((1:250)[diff(sign(diff(out[251,2:251])))==-2])
  
  if((wavesp<3)|(wavesm<3)){#print("hello")
    
    parms["Da"]<-parms["Da"]*(2/3)^2
    parms["Db"]<-parms["Db"]*(2/3)^2
    
    out <- ode.1D(y = state, times = times, func = rdmod2c, parms = parms,nspec = 2)
    
  }
  
  stop<-FALSE
  
  i<-5
  
  while(stop==FALSE){
    
    start<-((1:n)[sign(diff(out[i,2:(n+1)]))==0])
    
    if(length(start)==0){
      
      start<-n
      
    }else{
      
      start<-start[1]
      
    }
    
    L<-((1:n)[diff(sign(diff(out[i,2:(n+1)])))!=0])	
    
    L<-L[L<start]
    
    if((length(L)>4)&(2*L[2]==L[4])){ 
      
      stop<-TRUE
      
    }else{
      
      if(i==251){
        
        stop<-TRUE
        
      }
      
    }
    
    i<-i+1
    
  }
  
  if(length(L)<2){
    
    outputs<-matrix(ncol=12,dimnames=list(NULL,c("Type","Da","Aa","Ab","Db","Ba","Bb","ph","K2","l","va","vb")))
    outputs [1,]<-rep(NA,12)
    
  }else{
    
    
    k<-2*pi/(L0*L[2]/n)
    
    test_k2<-k^2 #this is the period of the system measured empirically
    
    
    Da<-parms["Da"]
    Aa<-parms["Caa"]-parms["Ra"]
    Ab<-parms["Cab"]
    Db<-parms["Db"]
    Ba<-parms["Cba"]
    Bb<-parms["Cbb"]-parms["Rb"]
    
    jacmat<-matrix(c(Aa-Da* test_k2,Ab,Ba,Bb-Db* test_k2),byrow=TRUE,ncol=2)
    
    eig<-eigen(jacmat)
    
    test_l<-max(Re(eig$values)) #this is the rate based on the measured period
    
    
    #scale based on empirical measures
    
    scale_k2<-(sqrt(target_k2)/sqrt(test_k2))^2 #remember, these are squares!
    
    scale_l<-target_l/test_l
    
    scale_parms<-parms*scale_l
    
    scale_parms["Da"]<-scale_parms["Da"]/scale_k2
    scale_parms["Db"]<-scale_parms["Db"]/scale_k2
    
    
    #now to run a growing simulation need to change the domain size, add a growth rate and bounds
    
    scale_parms["L0"]<-12.5
    scale_parms["S"]<-0.0075
    scale_parms["bounds"]<-TRUE
    
    parms<-scale_parms
    
    
    #known k2
    
    k2<-target_k2
    
    outputs<-with(as.list(parms),{
      
      
      Y<-matrix(ncol=12,dimnames=list(NULL,c("Type","Da","Aa","Ab","Db","Ba","Bb","ph","K2","l","va","vb")))
      
      Da<-Da
      Aa<-Caa-Ra
      Ab<-Cab
      Db<-Db
      Ba<-Cba
      Bb<-Cbb-Rb
      
      
      jacmat<-matrix(c(Aa-Da*k2,Ab,Ba,Bb-Db*k2),byrow=TRUE,ncol=2)
      
      eig<-eigen(jacmat)
      
      maxeig<-(1:2)[Re(eig$values)==max(Re(eig$values))]
      
      eigvec<-Re(eig$vectors[,maxeig])
      
      phase<-sign(c(eigvec[1]/eigvec[2]))
      
      if(length(maxeig)==1){
        
        Y[1,]<-c(1,Da,Aa,Ab,Db,Ba,Bb,phase,k2,max(Re(eig$values)),abs(eigvec))
        
      }else{
        
        Y[1,]<-rep(NA,12)
        
      }
      
      return(Y)
      
    })
    
  }
  
  paroutput<-c(as.numeric(param["sets"]),parms[1:12],outputs[,c(1,8:12)])
  
  out_table<-numeric()
  
  out_table
  
}



multi_perturb_2c<-function(out_table){
  
  parms<-c(out_table[2:13],S=00,L0=375,G0=250,GZ=0.2,bounds=TRUE)
  
  #equilibrium initial levels
  X<-matrix(c((as.list(parms)$Caa-as.list(parms)$Ra),as.list(parms)$Cab,as.list(parms)$Cba,(as.list(parms)$Cbb-as.list(parms)$Rb)),byrow=TRUE,ncol=2)
  Y<-c(-as.list(parms)$Ba,-as.list(parms)$Bb)
  
  Z<-solve(X,Y)
  
  n<-100
  a<-rep(Z[1],n)
  b<-rep(Z[2],n)
  
  #noise
  a<-a+runif(n,-0.01,0.01)
  b<-b+runif(n,-0.01,0.01)
  
  ## Initial conditions
  state    <- c(a,b)
  
  ## RUNNING the model:
  
  gen<-50000
  step<-50
  
  times  <- seq(0, gen, by = step)   # output wanted at these time intervals
  
  max_time<-(gen/step)+1
  
  check_time<-max_time-50
  
  maxtimes<-max(times)
  
  out <- ode.1D(y = state, times = times, func = rdmod2c, parms = parms, nspec = 2)
  
  final<-out[max_time,-1]
  
  minmaxA<-find_waves(final[1:n])
  minmaxB<-find_waves(final[(n+1):(2*n)])
  
  if((minmaxA==minmaxB)){
    minmax<-minmaxA
    
    if(minmax>0){
      
      if((minmax>4)&(minmax<9)){
        
        peaks<-TRUE
        
        #steady<-out[151,-1]
        
        finalmeanA<-mean(final[1:n])
        finalmeanB<-mean(final[(n+1):(2*n)])
        
        finalheight<-diff(range(final))
        
        
        state    <- final
        
        out <- ode.1D(y = state, times = times, func = rdmod2c, parms = parms, nspec = 2)
        
        unperturbed<-out[max_time,-1]
        
        steadymeanA<-mean(unperturbed[1:n])
        steadymeanB<-mean(unperturbed[(n+1):(2*n)])
        
        steadyminmaxA<-find_waves(unperturbed[1:n])
        steadyminmaxB<-find_waves(unperturbed[(n+1):(2*n)])
        
        steadyheight<-diff(range(unperturbed))
        
        threshA<-abs((steadymeanA-finalmeanA)/steadymeanA)
        threshB<-abs((steadymeanB-finalmeanB)/steadymeanB)
        
        
        if(((threshA<0.01)&(threshB<0.01))&((steadyminmaxA==minmax)&(steadyminmaxB==minmax))&(abs((steadyheight-finalheight)/steadyheight)<0.01)){ #ie does not drift over second half of sim
          
          parms2<-parms		
          
          
          ##INHIBIT 'A'
          
          inh<-.25
          stop<-FALSE
          j<-0 #power raised to
          k<-0 #number of 'good' responses recovered
          
          
          outputsA<-numeric()
          
          while(stop==FALSE){
            
            inh1<-1-inh^j
            inh2<-1
            
            parms["Caa"]<-as.numeric(parms2["Caa"])*inh1
            parms["Cab"]<-as.numeric(parms2["Cab"])*inh2
            parms["Cba"]<-as.numeric(parms2["Cba"])*inh1
            parms["Cbb"]<-as.numeric(parms2["Cbb"])*inh2
            
            
            state    <- final			
            
            out <- ode.1D(y = state, times = times, func = rdmod2c, parms = parms, nspec = 2)
            
            finalrun<-out[max_time,-1]
            
            finalminmaxA<-find_waves(finalrun[1:n])
            finalminmaxB<-find_waves(finalrun[(n+1):(2*n)])
            
            find_peaks<-FALSE
            
            if((finalminmaxA==minmax)&(finalminmaxB==minmax)){	#does it have the same number peaks?
              
              find_peaks<-TRUE
              
              steadyrun<-out[check_time,-1]
              
              finalheight<-diff(range(finalrun[1:n]))
              steadyheight<-diff(range(steadyrun[1:n]))
              
              if(abs((finalheight-steadyheight)/finalheight)<0.01){	#is the pattern relatively stable
                
                outputsA<-rbind(outputsA,finalrun)
                
                k<-k+1
                
              }
              
            }  	
            
            if((k==3)|(j==10)){
              
              stop<-TRUE
              
            }
            
            j<-j+1
            
          }
          
          outputsA<-rbind(unperturbed,outputsA)
          
          if(length(outputsA[,1])<4){
            
            outputsA<-rbind(outputsA,rep(NA,(2*n)))
            
          }
          
          if(NA%in%outputsA){
            
            if(find_peaks==FALSE){
              
              if((finalminmaxA>0)|(finalminmaxB>0)){
                
                summaryA<-c(2,j,rep(NA,18))	# peak number changes
                
              }else{
                
                summaryA<-c(3,j,rep(NA,18))	# waves died
                
              }
            }else{
              
              summaryA<-c(4,j,rep(NA,18))	# waves still growing/shrinking
              
            }
            
          }else{
            summaryA<-c(1,j,perturb_summary_2c(outputsA))	# works
          }
          
          
          
          ##INHIBIT 'B'
          
          inh<-.25
          stop<-FALSE
          j<-0 #power raised to
          k<-0 #number of 'good' responses recovered
          
          
          outputsB<-numeric()
          
          while(stop==FALSE){
            
            inh1<-1
            inh2<-1-inh^j
            
            parms["Caa"]<-as.numeric(parms2["Caa"])*inh1
            parms["Cab"]<-as.numeric(parms2["Cab"])*inh2
            parms["Cba"]<-as.numeric(parms2["Cba"])*inh1
            parms["Cbb"]<-as.numeric(parms2["Cbb"])*inh2
            
            
            state    <- final			
            
            out <- ode.1D(y = state, times = times, func = rdmod2c, parms = parms, nspec = 2)
            
            finalrun<-out[max_time,-1]
            
            finalminmaxA<-find_waves(finalrun[1:n])
            finalminmaxB<-find_waves(finalrun[(n+1):(2*n)])
            
            find_peaks<-FALSE
            
            if((finalminmaxA==minmax)&(finalminmaxB==minmax)){	#does it have the same number peaks?
              
              find_peaks<-TRUE
              
              steadyrun<-out[check_time,-1]
              
              finalheight<-diff(range(finalrun[1:n]))
              steadyheight<-diff(range(steadyrun[1:n]))
              
              if(abs((finalheight-steadyheight)/finalheight)<0.01){	#is the pattern relatively stable
                
                outputsB<-rbind(outputsB,finalrun)
                
                k<-k+1
                
              }
              
            }  	
            
            if((k==3)|(j==10)){
              
              stop<-TRUE
              
            }
            
            j<-j+1
            
          }
          
          outputsB<-rbind(unperturbed,outputsB)
          
          if(length(outputsB[,1])<4){######DIAGNOSE TYPE DEPENDING ON WHETHER STABILISED
            outputsB<-rbind(outputsB,rep(NA,(2*n)))
          }
          
          if(NA%in%outputsB){
            
            if(find_peaks==FALSE){
              
              if((finalminmaxA>0)|(finalminmaxB>0)){
                
                summaryB<-c(2,j,rep(NA,18))	# peak number changes
                
              }else{
                
                summaryB<-c(3,j,rep(NA,18))	# waves died
                
              }
              
            }else{
              
              summaryB<-c(4,j,rep(NA,18))	# waves still growing/shrinking
              
            }
            
          }else{
            
            summaryB<-c(1,j,perturb_summary_2c(outputsB))	# works
            
          }
          
          
          
          summary<-as.numeric(c(threshA,threshB,summaryA,summaryB))
          
          
        }else{
          
          summary<-c(rep(NA,2),rep(c(6,rep(NA,19)),2))	# waves growing/shrinking
          
          
        }	
        
      }else{
        
        summary<-c(rep(NA,2),rep(c(7,rep(NA,19)),2))	# wrong number of waves
        
      }
      
    }else{
      
      summary<-c(rep(NA,2),rep(c(5,rep(NA,19)),2))	# no waves
      
    }
    
  }else{
    
    summary<-c(rep(NA,2),rep(c(8,rep(NA,19)),2))	# weird - diff numbers of waves in diff components
    
  }
  
  as.numeric(summary)	
  
}


multi_perturb_2c_type_2<-function(out_table){
  
  parms<-c(out_table[2:13],S=00,L0=375,G0=250,GZ=0.2,bounds=TRUE)
  
  #equilibrium initial levels
  X<-matrix(c((as.list(parms)$Caa-as.list(parms)$Ra),as.list(parms)$Cab,as.list(parms)$Cba,(as.list(parms)$Cbb-as.list(parms)$Rb)),byrow=TRUE,ncol=2)
  Y<-c(-as.list(parms)$Ba,-as.list(parms)$Bb)
  
  Z<-solve(X,Y)
  
  n<-100
  a<-rep(Z[1],n)
  b<-rep(Z[2],n)
  
  #noise
  a<-a+runif(n,-0.01,0.01)
  b<-b+runif(n,-0.01,0.01)
  
  ## Initial conditions
  state    <- c(a,b)
  
  ## RUNNING the model:
  
  gen<-50000
  step<-50
  
  times  <- seq(0, gen, by = step)   # output wanted at these time intervals
  
  max_time<-(gen/step)+1
  
  check_time<-max_time-50
  
  maxtimes<-max(times)
  
  out <- ode.1D(y = state, times = times, func = rdmod2c, parms = parms, nspec = 2)
  
  final<-out[max_time,-1]
  
  minmaxA<-find_waves(final[1:n])
  minmaxB<-find_waves(final[(n+1):(2*n)])
  
  if((minmaxA==minmaxB)){
    minmax<-minmaxA
    
    if(minmax>0){
      
      if((minmax>4)&(minmax<9)){
        
        peaks<-TRUE
        
        #steady<-out[151,-1]
        
        finalmeanA<-mean(final[1:n])
        finalmeanB<-mean(final[(n+1):(2*n)])
        
        finalheight<-diff(range(final))
        
        
        state    <- final
        
        out <- ode.1D(y = state, times = times, func = rdmod2c, parms = parms, nspec = 2)
        
        unperturbed<-out[max_time,-1]
        
        steadymeanA<-mean(unperturbed[1:n])
        steadymeanB<-mean(unperturbed[(n+1):(2*n)])
        
        steadyminmaxA<-find_waves(unperturbed[1:n])
        steadyminmaxB<-find_waves(unperturbed[(n+1):(2*n)])
        
        steadyheight<-diff(range(unperturbed))
        
        threshA<-abs((steadymeanA-finalmeanA)/steadymeanA)
        threshB<-abs((steadymeanB-finalmeanB)/steadymeanB)
        
        
        if(((threshA<0.01)&(threshB<0.01))&((steadyminmaxA==minmax)&(steadyminmaxB==minmax))&(abs((steadyheight-finalheight)/steadyheight)<0.01)){ #ie does not drift over second half of sim
          
          parms2<-parms		
          
          
          ##INHIBIT 'A'
          
          inh<-.25
          stop<-FALSE
          j<-0 #power raised to
          k<-0 #number of 'good' responses recovered
          
          
          outputsA<-numeric()
          
          while(stop==FALSE){
            
            inh1<-1-inh^j
            inh2<-1
            
            parms["Caa"]<-as.numeric(parms2["Caa"])*inh1
            parms["Cab"]<-as.numeric(parms2["Cab"])*inh1
            parms["Ba"]<-as.numeric(parms2["Ba"])*inh1
            
            
            state    <- final			
            
            out <- ode.1D(y = state, times = times, func = rdmod2c, parms = parms, nspec = 2)
            
            finalrun<-out[max_time,-1]
            
            finalminmaxA<-find_waves(finalrun[1:n])
            finalminmaxB<-find_waves(finalrun[(n+1):(2*n)])
            
            find_peaks<-FALSE
            
            if((finalminmaxA==minmax)&(finalminmaxB==minmax)){	#does it have the same number peaks?
              
              find_peaks<-TRUE
              
              steadyrun<-out[check_time,-1]
              
              finalheight<-diff(range(finalrun[1:n]))
              steadyheight<-diff(range(steadyrun[1:n]))
              
              if(abs((finalheight-steadyheight)/finalheight)<0.01){	#is the pattern relatively stable
                
                outputsA<-rbind(outputsA,finalrun)
                
                k<-k+1
                
              }
              
            }  	
            
            if((k==3)|(j==10)){
              stop<-TRUE
            }
            
            j<-j+1
            
          }
          
          outputsA<-rbind(unperturbed,outputsA)
          
          if(length(outputsA[,1])<4){
            
            outputsA<-rbind(outputsA,rep(NA,(2*n)))
            
          }
          
          if(NA%in%outputsA){
            
            if(find_peaks==FALSE){
              
              if((finalminmaxA>0)|(finalminmaxB>0)){
                
                summaryA<-c(2,j,rep(NA,18))	# peak number changes
                
              }else{
                
                summaryA<-c(3,j,rep(NA,18))	# waves died
                
              }
              
            }else{
              
              summaryA<-c(4,j,rep(NA,18))	# waves still growing/shrinking
              
            }
            
          }else{
            
            summaryA<-c(1,j,perturb_summary_2c(outputsA))	# works
            
          }
          
          
          
          ##INHIBIT 'B'
          
          inh<-.25
          stop<-FALSE
          j<-0 #power raised to
          k<-0 #number of 'good' responses recovered
          
          
          outputsB<-numeric()
          
          while(stop==FALSE){
            
            inh1<-1
            inh2<-1-inh^j
            
            parms["Cba"]<-as.numeric(parms2["Cba"])*inh2
            parms["Cbb"]<-as.numeric(parms2["Cbb"])*inh2
            parms["Bb"]<-as.numeric(parms2["Bb"])*inh2
            
            
            state    <- final			
            
            out <- ode.1D(y = state, times = times, func = rdmod2c, parms = parms, nspec = 2)
            
            finalrun<-out[max_time,-1]
            
            finalminmaxA<-find_waves(finalrun[1:n])
            finalminmaxB<-find_waves(finalrun[(n+1):(2*n)])
            
            find_peaks<-FALSE
            
            if((finalminmaxA==minmax)&(finalminmaxB==minmax)){	#does it have the same number peaks?
              
              find_peaks<-TRUE
              
              steadyrun<-out[check_time,-1]
              
              finalheight<-diff(range(finalrun[1:n]))
              steadyheight<-diff(range(steadyrun[1:n]))
              
              if(abs((finalheight-steadyheight)/finalheight)<0.01){	#is the pattern relatively stable
                
                outputsB<-rbind(outputsB,finalrun)
                
                k<-k+1
                
              }
              
            }  
            
            if((k==3)|(j==10)){
              
              stop<-TRUE
              
            }
            
            j<-j+1
            
          }
          
          outputsB<-rbind(unperturbed,outputsB)
          
          if(length(outputsB[,1])<4){######DIAGNOSE TYPE DEPENDING ON WHETHER STABILISED
            
            outputsB<-rbind(outputsB,rep(NA,(2*n)))
            
          }
          
          if(NA%in%outputsB){
            
            if(find_peaks==FALSE){
              
              if((finalminmaxA>0)|(finalminmaxB>0)){
                
                summaryB<-c(2,j,rep(NA,18))	# peak number changes
                
              }else{
                
                summaryB<-c(3,j,rep(NA,18))	# waves died
                
              }
              
            }else{
              
              summaryB<-c(4,j,rep(NA,18))	# waves still growing/shrinking
              
            }
            
          }else{
            
            summaryB<-c(1,j,perturb_summary_2c(outputsB))	# works
            
          }
          
          summary<-as.numeric(c(threshA,threshB,summaryA,summaryB))
          
          
        }else{
          
          summary<-c(rep(NA,2),rep(c(6,rep(NA,19)),2))	# waves growing/shrinking
          
        }	
        
      }else{
        
        summary<-c(rep(NA,2),rep(c(7,rep(NA,19)),2))	# wrong number of waves
        
      }
      
    }else{
      
      summary<-c(rep(NA,2),rep(c(5,rep(NA,19)),2))	# no waves
      
    }
    
  }else{
    
    summary<-c(rep(NA,2),rep(c(8,rep(NA,19)),2))	# weird - diff numbers of waves in diff components
    
  }
  
  as.numeric(summary)	
  
}


#Same as equivalent function in '3-component simulations' file

find_waves<-function(x){
  
  w<-x-mean(x)
  
  length(w[diff(sign(w))!=0])
  
}


#Same as equivalent function in '3-component simulations' file

peak_value<-function(x){
  
  x_reindexed<-x[2:99] #99 because in these sims n=100, should make it 2:(n-1)
  
  peaktrough<-x_reindexed-mean(x_reindexed)
  
  minmax<-diff(sign(diff(x)))
  
  maxvalues<-x[(minmax<0)&(peaktrough>0)] #identify peaks and ensure actual peak (not local maximum in long flat trough)
  
  mean(maxvalues)
  
}


#Same as equivalent function in '3-component simulations' file

trough_value<-function(x){
  
  x_reindexed<-x[2:99]
  
  peaktrough<-x_reindexed-mean(x_reindexed)
  
  minmax<-diff(sign(diff(x)))
  
  maxvalues<-x[(minmax>0)&(peaktrough<0)]
  
  mean(maxvalues)
  
}


#Same as equivalent function in '3-component simulations' file

perturb_stats<-function(response){
  
  l<-length(response[1,])
  
  row_trends<-rowMeans(response)
  
  N<-length(row_trends)
  
  global<-(row_trends[2:N]-row_trends[1])/(row_trends[1])
  
  peak_trends<-t(apply(response,1,function(x) peak_value(x)))
  peak<-(peak_trends[2:N]-peak_trends[1])/(peak_trends[1])
  
  trough_trends<-t(apply(response,1,function(x) trough_value(x)))
  trough<-(trough_trends[2:N]-trough_trends[1])/(trough_trends[1])	
  
  data.frame(global=global,peak=peak,trough=trough)
  
}


perturb_summary_2c<-function(outputs){
  
  n<-length(outputs[1,])/2
  
  responseA<-outputs[,1:n]
  responseB<-outputs[,(n+1):(2*n)]
  
  statsA<-perturb_stats(responseA)
  statsB<-perturb_stats(responseB)
  
  
  globalA<-as.numeric(statsA$global)
  globalB<-as.numeric(statsB$global)
  
  peakA<-as.numeric(statsA$peak)
  peakB<-as.numeric(statsB$peak)
  
  troughA<-as.numeric(statsA$trough)
  troughB<-as.numeric(statsB$trough)
  
  l<-length(globalA)
  
  m<-round(l/2)
  
  output<-c(globalAmin=globalA[1],peakAmin=peakA[1],troughAmin=troughA[1],
            globalBmin=globalB[1],peakBmin=peakB[1],troughBmin=troughB[1],
            globalAmid=globalA[m],peakAmid=peakA[m],troughAmid=troughA[m],
            globalBmid=globalB[m],peakBmid=peakB[m],troughBmid=troughB[m],
            globalAmax=globalA[l],peakAmax=peakA[l],troughAmax=troughA[l],
            globalBmax=globalB[l],peakBmax=peakB[l],troughBmax=troughB[l])
  
  output
  
}
