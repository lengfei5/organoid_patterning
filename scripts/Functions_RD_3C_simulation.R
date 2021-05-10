#Functions and scripts to:
#1. Generate parameterisations of Jacobian matrix which satisfy the criteria for diffusion driven instability (as outlined in supplementary note) for a 3-component RD system
#2. Convert parameterisations of Jacobian matrix to the form of RD system used in analyses by including degradation rates, constant production rates and upper and lower reaction rate limits, and scale rates to form a spatial waves in the simulation grid
#3. Run simulations of inhibiting each component in RD system


#####################
#For running simulations
library(deSolve) #Version: 1.13	Depends: R (>= 2.15.0)
library(plyr) #Version: 1.8.4	Depends: R (>= 3.1.0)
library(doMC) #Version: 1.3.4	Depends: R (>= 2.14.0), foreach(>= 1.2.0), iterators(>= 1.0.0), parallel
registerDoMC(cores=6)

#######################

#Functions


#Function to generates a random parameter set and return it if it satisfies the criteria for a DDI (otherwise NAs are returned)
#In output vector, for parameter sets that will give a stable pattern of waves, the first term is 1, for parameter sets that will give an oscillating pattern the first terms is -1, and for parameter sets where based on the criteria used the two cannot be distinguished the first term is 0
#The next 12 terms in the output vector are the reaction rate and diffusion parameters
#The next two terms give the phases of components B and C relative to A respectively
#The next term give the wavelength parameter k^2 at the turning point on which conditions for DDI are evaluated (see supplementary notes)
#The final four terms give the largest eigen value and associated eigen vectors

#Note that essentially this same function was used to identify parameter sets for topology search (see emthods) except that rather than generating random parameter sets, parameter values were systematically varied as described

generate_parameters<-function(x,top){
	
	par<-runif(9,0,1) 

	diff<-runif(3,1,10000)

	Da<-500/diff[1]
	Aa<-par[1]*top[1,1]
	Ab<-par[2]*top[1,2]
	Ac<-par[3]*top[1,3]
	Db<-500/diff[2]
	Ba<-par[4]*top[2,1]
	Bb<-par[5]*top[2,2]
	Bc<-par[6]*top[2,3]
	Dc<-500/diff[3]
	Ca<-par[7]*top[3,1]
	Cb<-par[8]*top[3,2]
	Cc<-par[9]*top[3,3]

	Y<-rep(NA,20)

	#is the system stable in the absence of diffusion?

	jacmd0<-matrix(c(Aa,Ab,Ac,Ba,Bb,Bc,Ca,Cb,Cc),byrow=TRUE,ncol=3)

	a1_0<--Aa-Bb-Cc

	a2_0<-Aa*Cc+Aa*Bb+Bb*Cc-Cb*Bc-Ab*Ba-Ac*Ca

	a3_0<--det(jacmd0)
	
	

	b_a3<-Da*Db*Dc

	c_a3<--(Db*Da*Cc+Da*Dc*Bb+Db*Dc*Aa)

	d_a3<-Da*(Bb*Cc-Cb*Bc)+Db*(Aa*Cc-Ca*Ac)+Dc*(Aa*Bb-Ab*Ba)

	h_a3<-a3_0



	b_a123<-(Db+Dc)*((Da^2)+Db*Dc+Da*Db+Da*Dc)

	c_a123<--Aa*(Db+Dc)*(2*Da+Db+Dc)-Bb*(Da+Dc)*(Da+2*Db+Dc)-Cc*(Da+Db)*(Da+Db+2*Dc)

	d_a123<-Da*(2*Aa*Cc+2*Aa*Bb+2*Bb*Cc+(Cc^2)+(Bb^2)-Ab*Ba-Ac*Ca)+Db*(2*Aa*Cc+2*Aa*Bb+2*Bb*Cc+(Aa^2)+(Cc^2)-Ba*Ab-Bc*Cb)+Dc*(2*Aa*Cc+2*Aa*Bb+2*Bb*Cc+(Aa^2)+(Bb^2)-Ac*Ca-Bc*Cb)

	h_a123<-a1_0*a2_0-a3_0


	if((a1_0>0)&(a3_0>0)&(h_a123>0)){
		#print("reaction terms stable")
	
		DDIa3<-F
		DDIa123<-F
		
		#for a3		
		if((d_a3<0)|((c_a3<0)&(c_a3^2>3*b_a3*d_a3))){
			
			a3_k2_TP<-(-c_a3+sqrt((c_a3^2)-3*b_a3*d_a3))/(3*b_a3)
			
			y_TP_a3<-b_a3*a3_k2_TP^3+c_a3*a3_k2_TP^2+d_a3*a3_k2_TP+h_a3	
			
			if(y_TP_a3<0){
			
				DDIa3<-T
			
			}	
		}
			
		#for a123	
		if((d_a123<0)|((c_a123<0)&(c_a123^2>3*b_a123*d_a123))){
			
			a123_k2_TP<-(-c_a123+sqrt((c_a123^2)-3*b_a123*d_a123))/(3*b_a123)
			
			y_TP_a123<-b_a123*a123_k2_TP^3+c_a123*a123_k2_TP^2+d_a123*a123_k2_TP+h_a123
			
			if(y_TP_a123<0){
			
				DDIa123<-T
		
			}		
		
		}
		
		if((DDIa3==T)|(DDIa123==T)){
			#print("diffusion driven instability")
		
			c_a2<-Db*Dc+Db*Da+Da*Dc
	
			d_a2<-Dc*Bb+Dc*Aa+Db*Cc+Db*Aa+Da*Bb+Da*Cc

			h_a2<-Aa*Cc+Aa*Bb+Bb*Cc-Cb*Bc-Ab*Ba-Ac*Ca
		
		
			if((DDIa3==T)&(DDIa123==F)){
				
				#criteria for stable waves
				
				d_a1<-Da+Db+Dc
				
				h_a1<--Aa-Bb-Cc
				
				
				c_a12<-c_a2+d_a1^2
				
				d_a12<-2*d_a1*h_a1+d_a2
				
				h_a12<-h_a1^2+h_a2
				
				a12_k2_TP<--d_a12/(2*c_a12)
					
				#estimate phase from jacobian evaluated at k^2 of turning point
				jacmat<-matrix(c(Aa-Da*a3_k2_TP,Ab,Ac,Ba,Bb-Db*a3_k2_TP,Bc,Ca,Cb,Cc-Dc*a3_k2_TP),byrow=TRUE,ncol=3)
					
				eig<-eigen(jacmat)
					
				maxeig<-(1:3)[Re(eig$values)==max(Re(eig$values))]
					
				eigvec<-Re(eig$vectors[,maxeig])
					
				phase<-sign(c(eigvec[1]/eigvec[2],eigvec[1]/eigvec[3]))

				if((c_a12*a12_k2_TP^2+d_a12*a12_k2_TP+h_a12)>0){

					Y<-c(1,Da,Aa,Ab,Ac,Db,Ba,Bb,Bc,Dc,Ca,Cb,Cc,phase,a3_k2_TP,max(Re(eig$values)),abs(eigvec))
					
				}
				
				else{	
				
					Y<-c(0,Da,Aa,Ab,Ac,Db,Ba,Bb,Bc,Dc,Ca,Cb,Cc,phase,a3_k2_TP,max(Re(eig$values)),abs(eigvec))	
						
				}
					
			}
						
			if((DDIa123==T)&(DDIa3==F)){
				
				#criteria for oscillating waves
					
				a2_k2_TP<--d_a2/(2*c_a2)

				#estimate phase from jacobian evaluated at k^2 of turning point
				jacmat<-matrix(c(Aa-Da*a123_k2_TP,Ab,Ac,Ba,Bb-Db*a123_k2_TP,Bc,Ca,Cb,Cc-Dc*a123_k2_TP),byrow=TRUE,ncol=3)
					
				eig<-eigen(jacmat)
					
				maxeig<-(1:3)[Re(eig$values)==max(Re(eig$values))]
					
				eigvec<-Re(eig$vectors[,maxeig])
					
				phase<-sign(c(eigvec[1]/eigvec[2],eigvec[1]/eigvec[3]))

				if((c_a2*a2_k2_TP^2+d_a2*a2_k2_TP+h_a2)>0){
				
					Y<-c(-1,Da,Aa,Ab,Ac,Db,Ba,Bb,Bc,Dc,Ca,Cb,Cc,0,0,0,0,0,0,0)	
							
				}
				
				else{
		
					Y<-c(0,Da,Aa,Ab,Ac,Db,Ba,Bb,Bc,Dc,Ca,Cb,Cc,0,0,0,0,0,0,0)	
					
				}
				
		
			}
			
			if((DDIa123==T)&(DDIa3==T)){
	
				Y<-c(0,Da,Aa,Ab,Ac,Db,Ba,Bb,Bc,Dc,Ca,Cb,Cc,0,0,0,0,0,0,0)
					
			}
	
		}
	}	
	
	Y	
	
}



#Model for simulating 3-component RD system on a growing domian for three species a, b and c
#Implementation of domain growth based on apical growth in Crampin et al (2002) Bull Math Biol 64, 747-769, except the growth zone of size G0, is split into an apical non growing part (of size GZ*G0, where GZ<1) with growth occurring in the remainder of the growt zone
#Note that for simulations used in inhibition analyses, there is no domain growth (S=0)

rdmod<-function(time,state,parms){
	with(as.list(parms),{
	
	n<-length(state)/3
	
	a<-state[1:n]
	b<-state[(n+1):(2*n)]	
	c<-state[(2*n+1):(3*n)]
	
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
	Diffc<-(Dc/(L^2))*diff(diff(c(c[1],c,c[n])))/(dz^2)
	
	#reaction
	F<-Caa*a+Cab*b+Cac*c+Ba
	G<-Cba*a+Cbb*b+Cbc*c+Bb
	H<-Cca*a+Ccb*b+Ccc*c+Bc

	if(bounds==TRUE){

		F[F>maxA]<-maxA
		F[F<0]<-0
		G[G>maxB]<-maxB
		G[G<0]<-0
		H[H>maxC]<-maxC
		H[H<0]<-0		

	}


	#flow
	A<-numeric(length=n)
	A[z<Ph1]<-0
	A[(Ph1<=z)&(z<Ph2)]<-S*(z[(Ph1<=z)&(z<Ph2)]*L-Th1)/(Th2-Th1)
	A[Ph2<=z]<-S

	#advection
	Adva<-((z*dL-A)/L)*diff(c(a[1],a,a[n]),2)/(2*dz)
	Advb<-((z*dL-A)/L)*diff(c(b[1],b,b[n]),2)/(2*dz)
	Advc<-((z*dL-A)/L)*diff(c(c[1],c,c[n]),2)/(2*dz)

	#dilution due to growth
	Sdil<-numeric(length=n)
	Sdil[z<Ph1]<-0
	Sdil[(Ph1<=z)&(z<Ph2)]<-S/(Th2-Th1)
	Sdil[Ph2<=z]<-0


	#rate of change	
	da<-F+Diffa-Ra*a-Sdil*a+Adva
	db<-G+Diffb-Rb*b-Sdil*b+Advb
	dc<-H+Diffc-Rc*c-Sdil*c+Advc
		
	return (list(c(da,db,dc)))	
	
	})
	
}


#Function 'parameters' generates parameter set from output of 'generate_parameters' for simulation using 'rdmod'
#Function adds degrdation rates, background production rates and upper limits to reaction rates
#For simulations of growing domains (as described above for rdmod), linear growth rate of the domain is added 'S' alongside and initial length for the domain 'L0'
#Length of 'growth zone' is specified in 'G0', and a non-growing 'apical' proportion specified by 'GZ'
#For simulation of systems without growth, set S=0
#If 'bound=TRUE' upper and lower limits on reaction rate will be used, resulting in stable waves forming, otherwise if 'bounds=FALSE' waves will grow to infinity

parameters<-function(par,S,L0,G0,GZ,bounds){
	
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
		if(Cc>0){
			Rc<-sample(seq(0.01,1,by=0.01),1)
			Ccc<-Cc+Rc
		}else{
			Rc<--Cc
			Ccc<-0
		}	

		#set equilibrium levels to (1,1,1) and calculate background rates
		X<-matrix(c(Aa,Ab,Ac,Ba,Bb,Bc,Ca,Cb,Cc),byrow=TRUE,ncol=3)
		Xinv<-solve(X)
		eqm<-c(1,1,1)
		bgX<-solve(Xinv,eqm)
		bg<-c(bgX,eqm)

		#set max rates at random between 1.5 and 3 fold equilibrium reaction rate
		lim<-runif(1,-1,1)
		maxA<-(1+2^lim)*Ra*bg[4]
		maxB<-(1+2^lim)*Rb*bg[5]
		maxC<-(1+2^lim)*Rc*bg[6]		
		
		#add parameters to vector, with zeros for b/g production and max
		parms<-c(Da=Da,
				Caa=Caa,
				Cab=Ab,
				Cac=Ac,
				Ba=-bg[1],
				Ra=Ra,
				maxA=maxA,
				Db=Db,
				Cba=Ba,
				Cbb=Cbb,
				Cbc=Bc,
				Bb=-bg[2],
				Rb=Rb,
				maxB=maxB,
				Dc=Dc,
				Cca=Ca,
				Ccb=Cb,
				Ccc=Ccc,
				Bc=-bg[3],
				Rc=Rc,
				maxC=maxC,
				S=S,
				L0=L0,
				G0=G0,
				GZ=GZ,
				bounds=bounds)
		
		return(parms)

	})
	
	parset

}



#Function 'pattern_search' takes a parameterisation of the Jacobian matrix 'param' which gives a series of stable waves (ie not oscillating - type 1 in 'generate_parameters'), as output by function 'generate_parameters'
#First, background production rates, degradation rates and limits to reaction rates are added as according the linear RD model used for simulations
#Second, the rate and diffusion parameters are secaled as according to an empirically derived target wavelength 'target_k2' and growth rate 'target_l' so that when the parameter set is simulated a series of non-discretised waves will form on the simulation grid, and the waves will grow and stablise within the duration of a simulation 
#'target_k2' and 'target_l' are derived from simulating an RD system on an apically growing domain (as described above), and scaling rate parameters and diffusion to give a pattern of stripe appearance in the 'growth zone' resembling that seen for the rugae - 'target_k2' is derived from the observed wavelength and 'target_l' from the largest eigen value of the system
#Scaling is a two-step process: 
#First, the system is scaled according to the wavelength parameter and eigen value output from 'generate_parameters', which should allow approximately the correct number of waves to grow at the desired rate
#Second, the scaled parameter set is simulated and the wavelength determined (and from this the associated eigen value derived)
#Based on these empirically derived measures the parameter set is rescaled to give the final parameter set 


pattern_search<-function(param,target_k2,target_l){
	
	test_k2<-as.numeric(param["k2"])

	test_l<-as.numeric(param["l"])

	scale_k2<-(sqrt(target_k2)/sqrt(test_k2))^2 #remember, these are squares!

	scale_l<-target_l/test_l

	#use 'parameters' function for additional rates as in RD model
	#set S=0 so no domian growth and bounds=FALSE as saturating the reaction rates to prevent waves growing infinitely appeared to slightly affect wavelength
	parms<-parameters(param,0,50,50,0.5,FALSE) 

	#prelimimarily scale according to tp k2 so si is in right area
	scale_parms<-parms*scale_l
	scale_parms["L0"]<-scale_parms["L0"]/scale_l #this is not a rate so don't scale with above

	scale_parms["Da"]<-scale_parms["Da"]/scale_k2
	scale_parms["Db"]<-scale_parms["Db"]/scale_k2
	scale_parms["Dc"]<-scale_parms["Dc"]/scale_k2


	#run a preliminary simulation from which to calculate the period

	parms<-scale_parms 

	#for initial conditions set level of each species to that of the equilibrium for the reaction terms and initiate waves in first point of simultaion grid so waves grow sequentially from boundary

	#equilibrium initial levels
	X<-matrix(c((as.list(parms)$Caa-as.list(parms)$Ra),as.list(parms)$Cab,as.list(parms)$Cac,as.list(parms)$Cba,(as.list(parms)$Cbb-as.list(parms)$Rb),as.list(parms)$Cbc,as.list(parms)$Cca,as.list(parms)$Ccb,(as.list(parms)$Ccc-as.list(parms)$Rc)),byrow=TRUE,ncol=3)
	Y<-c(-as.list(parms)$Ba,-as.list(parms)$Bb,-as.list(parms)$Bc)

	Z<-solve(X,Y)

	n<-250

	a<-rep(Z[1],n)
	b<-rep(Z[2],n)
	c<-rep(Z[3],n)

	#as doing unbounded to measure period initiate with perturbation at one end
	a[1]<-a[1]+.1

	## Initial conditions
	state    <- c(a,b,c)
	
	## RUNNING the model:

	times  <- seq(0, 2500, by = 10)   # output wanted at these time intervals

	out <- ode.1D(y = state, times = times, func = rdmod, parms = parms,nspec = 3)
                

	#measure period

	L0<-as.list(parms)$L0

	#first check enough waves (ie >=3) to measure, if not scale diffusion by 2/3

	wavesp<-length((1:250)[diff(sign(diff(out[251,2:251])))==2]) #peaks
	wavesm<-length((1:250)[diff(sign(diff(out[251,2:251])))==-2]) #troughs

	if((wavesp<3)|(wavesm<3)){
		
		parms["Da"]<-parms["Da"]*(2/3)^2
		parms["Db"]<-parms["Db"]*(2/3)^2
		parms["Dc"]<-parms["Dc"]*(2/3)^2
	
		#rerun simulation
	
		out <- ode.1D(y = state, times = times, func = rdmod, parms = parms,nspec = 3)

	}


	#searching where waves fill whole simultion field and first two waves are same length

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

	if(length(L)<2){ #ie not find two waves of same length in field over course of simulation
		
		outputs<-matrix(ncol=20,dimnames=list(NULL,c("Type","Da","Aa","Ab","Ac","Db","Ba","Bb","Bc","Dc","Ca","Cb","Cc","phAB","phAC","K2","l","va","vb","vc")))
		
		outputs [1,]<-rep(NA,20)
		
	}else{
	
		k<-2*pi/(L0*L[2]/n)

		test_k2<-k^2 #this is the period of the system measured empirically


		Da<-parms["Da"]
		Aa<-parms["Caa"]-parms["Ra"]
		Ab<-parms["Cab"]
		Ac<-parms["Cac"]
		Db<-parms["Db"]
		Ba<-parms["Cba"]
		Bb<-parms["Cbb"]-parms["Rb"]
		Bc<-parms["Cbc"]
		Dc<-parms["Dc"]
		Ca<-parms["Cca"]
		Cb<-parms["Ccb"]
		Cc<-parms["Ccc"]-parms["Rc"]

		jacmat<-matrix(c(Aa-Da* test_k2,Ab,Ac,Ba,Bb-Db* test_k2,Bc,Ca,Cb,Cc-Dc* test_k2),byrow=TRUE,ncol=3)
			
		eig<-eigen(jacmat)

		test_l<-max(Re(eig$values)) #this is the rate based on the measured period


		#scale based on empirical measures

		scale_k2<-(sqrt(target_k2)/sqrt(test_k2))^2 #remember, these are squares!

		scale_l<-target_l/test_l

		scale_parms<-parms*scale_l

		scale_parms["Da"]<-scale_parms["Da"]/scale_k2
		scale_parms["Db"]<-scale_parms["Db"]/scale_k2
		scale_parms["Dc"]<-scale_parms["Dc"]/scale_k2


		#now to run a growing simulation need to change the domain size, add a growth rate and bounds

		scale_parms["L0"]<-12.5
		scale_parms["S"]<-0.0075
		scale_parms["bounds"]<-TRUE

		parms<-scale_parms

		
		#parameters are now scaled
		
		#for phase of system, evaluate eigen vectors of system Jacobian matrix at target k2

		k2<-target_k2 #k2 is known from above

		outputs<-with(as.list(parms),{

			Y<-matrix(ncol=20,dimnames=list(NULL,c("Type","Da","Aa","Ab","Ac","Db","Ba","Bb","Bc","Dc","Ca","Cb","Cc","phAB","phAC","K2","l","va","vb","vc")))


			Da<-Da
			Aa<-Caa-Ra
			Ab<-Cab
			Ac<-Cac
			Db<-Db
			Ba<-Cba
			Bb<-Cbb-Rb
			Bc<-Cbc
			Dc<-Dc
			Ca<-Cca
			Cb<-Ccb
			Cc<-Ccc-Rc
			
			
			jacmat<-matrix(c(Aa-Da*k2,Ab,Ac,Ba,Bb-Db*k2,Bc,Ca,Cb,Cc-Dc*k2),byrow=TRUE,ncol=3)
			
			eig<-eigen(jacmat)
					
			maxeig<-(1:3)[Re(eig$values)==max(Re(eig$values))]
					
			eigvec<-Re(eig$vectors[,maxeig])
								
			phase<-sign(c(eigvec[1]/eigvec[2],eigvec[1]/eigvec[3]))
					
			if(length(maxeig)==1){
		
				Y[1,]<-c(1,Da,Aa,Ab,Ac,Db,Ba,Bb,Bc,Dc,Ca,Cb,Cc,phase,k2,max(Re(eig$values)),abs(eigvec))
		
			}else{
		
				Y[1,]<-rep(NA,20)
		
			}
			
			return(Y)
	
		})

	}

	paroutput<-c(as.numeric(param["sets"]),parms[1:21],outputs[,c(1,14:20)]) #return scaled parameters and associated rate and phase measures
	
	out_table<-numeric()
	
	out_table<-paroutput

	out_table

}



#Function 'multi_perturb' takes a parameter set as output by 'pattern_search' and simulates a series of waves
#Starting from this pattern of waves, each component (A, B and C) is then inhibited and the effect on waves of each of A, B and C is recorded
#Output is given as a change in the mean level of each component after inhibition of each pathway relative to an unperturbed state
#As system is modelling the biological network of Wnt, BMP and Hh (respectively A, B and C) production of A is inhibited, while response to B and C is inhibited
#Inhibition is initially a complete block of the pathway but if the waves are lost then the inhibition strength is successively reduced until waves are maintained after inhibtion (as the simulations are modelling a biological process where small molecule inhibitors alter the width and intensity of a series of stripes of expression, waves need to be maintained after inhibition)
#As well as recording the effect of the strongest inhibition that maintains waves, the effect of two successively smaller inhibitions are also recorded
#As well as recording the effect of inhibition on the mean level of each component, the change in the level of the wave peaks and troughs relative the the unperturbed states are also recorded (to ensure waves are not flattening out) 

multi_perturb<-function(out_table){
	
	parms<-c(out_table[2:22],S=00,L0=375,G0=250,GZ=0.2,bounds=TRUE)

	#initial values for species are small perturbations around the equilibrium value of the reaction terms 
	
	#equilibrium initial levels
	X<-matrix(c((as.list(parms)$Caa-as.list(parms)$Ra),as.list(parms)$Cab,as.list(parms)$Cac,as.list(parms)$Cba,(as.list(parms)$Cbb-as.list(parms)$Rb),as.list(parms)$Cbc,as.list(parms)$Cca,as.list(parms)$Ccb,(as.list(parms)$Ccc-as.list(parms)$Rc)),byrow=TRUE,ncol=3)
	Y<-c(-as.list(parms)$Ba,-as.list(parms)$Bb,-as.list(parms)$Bc)

	Z<-solve(X,Y)

	n<-100
	a<-rep(Z[1],n)
	b<-rep(Z[2],n)
	c<-rep(Z[3],n)
	
	#noise
	a<-a+runif(n,-0.01,0.01)
	b<-b+runif(n,-0.01,0.01)
	c<-c+runif(n,-0.01,0.01)

	## Initial conditions
	state    <- c(a,b,c)

	## RUNNING the model:

	gen<-50000 #longer simulation length used so higher proportion of parameterisations produce a stable pattern of waves
	step<-50

	times  <- seq(0, gen, by = step)   # output wanted at these time intervals

	max_time<-(gen/step)+1
	
	#timepoint to use to ensure that end of simulation is stable
	check_time<-max_time-50

	maxtimes<-max(times)	
	
	#run simulation
	out <- ode.1D(y = state, times = times, func = rdmod, parms = parms, nspec = 3)

	#extract final pattern of waves
	final<-out[max_time,-1]

	#find number of successive peaks and troughs
	minmaxA<-find_waves(final[1:n])
	minmaxB<-find_waves(final[(n+1):(2*n)])
	minmaxC<-find_waves(final[(2*n+1):(3*n)])

	if((minmaxA==minmaxB)&(minmaxA==minmaxC)){ #ensure that all species have same number of waves (ie nothing strange happening at edges of simulation field)
		
		minmax<-minmaxA
	
		if(minmax>0){ #ie series of waves do indeed form 
	
			if((minmax>4)&(minmax<9)){ #ie 2 to 4.5 complete waves
			
				peaks<-TRUE

				#record mean levels for each species
				finalmeanA<-mean(final[1:n])
				finalmeanB<-mean(final[(n+1):(2*n)])
				finalmeanC<-mean(final[(2*n+1):(3*n)])

				finalheight<-diff(range(final)) #amplitude of largest wave
		
				
				#first run a simulation of the unperturbed system
		
				state    <- final
	
				out <- ode.1D(y = state, times = times, func = rdmod, parms = parms, nspec = 3)
				
				#final waves of unpertured system
				unperturbed<-out[max_time,-1]

				steadymeanA<-mean(unperturbed[1:n])
				steadymeanB<-mean(unperturbed[(n+1):(2*n)])
				steadymeanC<-mean(unperturbed[(2*n+1):(3*n)])
	
				steadyminmaxA<-find_waves(unperturbed[1:n])
				steadyminmaxB<-find_waves(unperturbed[(n+1):(2*n)])
				steadyminmaxC<-find_waves(unperturbed[(2*n+1):(3*n)])
	
				steadyheight<-diff(range(unperturbed))
			
				threshA<-abs((steadymeanA-finalmeanA)/steadymeanA)
				threshB<-abs((steadymeanB-finalmeanB)/steadymeanB)
				threshC<-abs((steadymeanC-finalmeanC)/steadymeanC)
				
				#check that the unperturbed system is stable (ie does not drift over second half of sim) by confirming that there is no major change in the mean values of each wave and amplitude of biggest wave
				if(((threshA<0.01)&(threshB<0.01)&(threshC<0.01))&((steadyminmaxA==minmax)&(steadyminmaxB==minmax)&(steadyminmaxC==minmax))&(abs((steadyheight-finalheight)/steadyheight)<0.01)){ 
					
					parms2<-parms		

					
					#Successively inhibit each component
					#first completely block interaction, then successively reduce inhibition strength by 0.25

					##INHIBIT 'A'

					inh<-.25
					stop<-FALSE
					j<-0 #power raised to
					k<-0 #number of 'good' responses recovered


					outputsA<-numeric()
	
					while(stop==FALSE){
					
						inh1<-1-inh^j
						inh2<-1
						inh3<-1
	
						parms["Caa"]<-as.numeric(parms2["Caa"])*inh1
						parms["Cab"]<-as.numeric(parms2["Cab"])*inh1
						parms["Cac"]<-as.numeric(parms2["Cac"])*inh1
						parms["Ba"]<-as.numeric(parms2["Ba"])*inh1

		
						state    <- final	#series of waves from initial simulation		

						out <- ode.1D(y = state, times = times, func = rdmod, parms = parms, nspec = 3)

						finalrun<-out[max_time,-1]

						finalminmaxA<-find_waves(finalrun[1:n])
						finalminmaxB<-find_waves(finalrun[(n+1):(2*n)])
						finalminmaxC<-find_waves(finalrun[(2*n+1):(3*n)])

						find_peaks<-FALSE

						if((finalminmaxA==minmax)&(finalminmaxB==minmax)&(finalminmaxC==minmax)){	#does it have the same number peaks?
		
							find_peaks<-TRUE
		
							steadyrun<-out[check_time,-1]
		
							finalheight<-diff(range(finalrun[1:n]))
							steadyheight<-diff(range(steadyrun[1:n]))
	
							if(abs((finalheight-steadyheight)/finalheight)<0.01){	#is the pattern relatively stable
			
								outputsA<-rbind(outputsA,finalrun)
			
								k<-k+1
			
							}

						}  	
	
						if((k==3)|(j==10)){ #stop simulations when have recorded three successive inhibitions, or if before then inhibition strength drops to 1-0.25^10
							stop<-TRUE
						}
	
						j<-j+1
	
					}
	
					outputsA<-rbind(unperturbed,outputsA)

					if(length(outputsA[,1])<4){ #could not perforem three successive inhibitions
					
						outputsA<-rbind(outputsA,rep(NA,(3*n)))
					
					}

					#based on what is recorded in outputsA record output statistics only if simulations 'worked' 
					#record different types of 'failure' by first term in summary vector, where '1' is worked

					if(NA%in%outputsA){
						
						if(find_peaks==FALSE){
						
							if((finalminmaxA>0)|(finalminmaxB>0)|(finalminmaxC>0)){
		
								summaryA<-c(2,j,rep(NA,27))	# peak number changes
		
							}else{
			
								summaryA<-c(3,j,rep(NA,27))	# waves died
			
							}
						}else{
		
							summaryA<-c(4,j,rep(NA,27))	# waves still growing/shrinking
		
						}
		
					}else{
						summaryA<-c(1,j,perturb_summary(outputsA))	# works! record summary statistics (see function 'perturb_summary')
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
						inh3<-1
		
						parms["Caa"]<-as.numeric(parms2["Caa"])*inh1
						parms["Cab"]<-as.numeric(parms2["Cab"])*inh2
						parms["Cac"]<-as.numeric(parms2["Cac"])*inh3
						parms["Cba"]<-as.numeric(parms2["Cba"])*inh1
						parms["Cbb"]<-as.numeric(parms2["Cbb"])*inh2
						parms["Cbc"]<-as.numeric(parms2["Cbc"])*inh3
						parms["Cca"]<-as.numeric(parms2["Cca"])*inh1
						parms["Ccb"]<-as.numeric(parms2["Ccb"])*inh2
						parms["Ccc"]<-as.numeric(parms2["Ccc"])*inh3

						state    <- final	

						out <- ode.1D(y = state, times = times, func = rdmod, parms = parms, nspec = 3)
		
						finalrun<-out[max_time,-1]

						finalminmaxA<-find_waves(finalrun[1:n])
						finalminmaxB<-find_waves(finalrun[(n+1):(2*n)])
						finalminmaxC<-find_waves(finalrun[(2*n+1):(3*n)])

						find_peaks<-FALSE

						if((finalminmaxA==minmax)&(finalminmaxB==minmax)&(finalminmaxC==minmax)){	#does it have the same number peaks?
		
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

					if(length(outputsB[,1])<4){
						
						outputsB<-rbind(outputsB,rep(NA,(3*n)))
					
					}
	
					if(NA%in%outputsB){
					
						if(find_peaks==FALSE){
					
							if((finalminmaxA>0)|(finalminmaxB>0)|(finalminmaxC>0)){
		
								summaryB<-c(2,j,rep(NA,27))
		
							}else{
			
								summaryB<-c(3,j,rep(NA,27))
			
							}
						}else{
		
							summaryB<-c(4,j,rep(NA,27))
		
						}
	
					}else{
				
						summaryB<-c(1,j,perturb_summary(outputsB))
				
					}



					##INHIBIT 'C'

					inh<-.25
					stop<-FALSE
					j<-0 #power raised to
					k<-0 #number of 'good' responses recovered


					outputsC<-numeric()
	
					while(stop==FALSE){
				
						inh1<-1
						inh2<-1
						inh3<-1-inh^j
		
						parms["Caa"]<-as.numeric(parms2["Caa"])*inh1
						parms["Cab"]<-as.numeric(parms2["Cab"])*inh2
						parms["Cac"]<-as.numeric(parms2["Cac"])*inh3
						parms["Cba"]<-as.numeric(parms2["Cba"])*inh1
						parms["Cbb"]<-as.numeric(parms2["Cbb"])*inh2
						parms["Cbc"]<-as.numeric(parms2["Cbc"])*inh3
						parms["Cca"]<-as.numeric(parms2["Cca"])*inh1
						parms["Ccb"]<-as.numeric(parms2["Ccb"])*inh2
						parms["Ccc"]<-as.numeric(parms2["Ccc"])*inh3
	
						state    <- final

						out <- ode.1D(y = state, times = times, func = rdmod, parms = parms, nspec = 3)

						finalrun<-out[max_time,-1]

						finalminmaxA<-find_waves(finalrun[1:n])
						finalminmaxB<-find_waves(finalrun[(n+1):(2*n)])
						finalminmaxC<-find_waves(finalrun[(2*n+1):(3*n)])

						find_peaks<-FALSE

						if((finalminmaxA==minmax)&(finalminmaxB==minmax)&(finalminmaxC==minmax)){	#does it have the same number peaks?
		
							find_peaks<-TRUE
		
							steadyrun<-out[check_time,-1]
		
							finalheight<-diff(range(finalrun[1:n]))
							steadyheight<-diff(range(steadyrun[1:n]))

							if(abs((finalheight-steadyheight)/finalheight)<0.01){	#is the pattern relatively stable
			
								outputsC<-rbind(outputsC,finalrun)
			
								k<-k+1
			
							}

						}  
	
						if((k==3)|(j==10)){
							stop<-TRUE
						}
		
						j<-j+1
	
					}

					outputsC<-rbind(unperturbed,outputsC)

					if(length(outputsC[,1])<4){
				
						outputsC<-rbind(outputsC,rep(NA,(3*n)))
				
					}

					if(NA%in%outputsC){
					
						if(find_peaks==FALSE){
					
							if((finalminmaxA>0)|(finalminmaxB>0)|(finalminmaxC>0)){
			
								summaryC<-c(2,j,rep(NA,27))
			
							}else{
			
								summaryC<-c(3,j,rep(NA,27))
			
							}
				
						}else{
		
							summaryC<-c(4,j,rep(NA,27))
		
						}
	
					}else{
				
						summaryC<-c(1,j,perturb_summary(outputsC))
				
					}


					summary<-as.numeric(c(threshA,threshB,threshC,summaryA,summaryB,summaryC)) #as well as perturbation statistics, also record 'threshX' as a measure of drift in unperturbed system
		
	
				}else{
			
					summary<-c(rep(NA,3),rep(c(6,rep(NA,28)),3))	# waves growing/shrinking in unperturbed state
			
			
				}	
		
			}else{

			summary<-c(rep(NA,3),rep(c(7,rep(NA,28)),3))	# wrong number of waves (too many or few for simulation field)
	
			}
	
		}else{

			summary<-c(rep(NA,3),rep(c(5,rep(NA,28)),3))	# no waves after initial simulation
	
		}
	
	}else{

		summary<-c(rep(NA,3),rep(c(8,rep(NA,28)),3))	# weird - diff numbers of waves in diff components

	}
		
	as.numeric(summary)	
	
}


#Function 'find_waves' takes a vector 'x' describing a wave and measures the number of successive peaks and troughs

find_waves<-function(x){
	
	w<-x-mean(x)
	
	length(w[diff(sign(w))!=0])
	
}


#Function 'peak_value' takes a vector 'x' describing a wave, identifies peaks, and evaluates the mean value of x at the peaks

peak_value<-function(x){
	
	x_reindexed<-x[2:99] #99 because in these sims n=100, should make it 2:(n-1)
	
	peaktrough<-x_reindexed-mean(x_reindexed)
	
	minmax<-diff(sign(diff(x)))
	
	maxvalues<-x[(minmax<0)&(peaktrough>0)] #identify peaks and ensure actual peak (not local maximum in long flat trough)

	mean(maxvalues)
	
}


#Function 'trough_value' takes a vector 'x' describing a wave, identifies troughs, and evaluates the mean value of x at the troughs

trough_value<-function(x){
	
	x_reindexed<-x[2:99]
	
	peaktrough<-x_reindexed-mean(x_reindexed)
	
	minmax<-diff(sign(diff(x)))
	
	maxvalues<-x[(minmax>0)&(peaktrough<0)]

	mean(maxvalues)
	
}


#Function 'perturb_stats' takes a matrix where the each row is the wave for a particular species, first the unperturbed state, then the three inhibitions, smallest to largest, and returns the response of the mean, peak and trough level for each component, relative to the unperturbed state, for the three inhibition strengths

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


#Function 'perturb_summary' takes a matrix where the first row is the concatenated unperturbed waves for species a, b and c, and the second to fourth rows are the waves for the three inhibitions in order of increasing strength (as generated in 'multi_perturb') and returns the response of the mean, peak and trough level for each component, relative to the unperturbed state, for the three inhibition strengths
#'perturb_summary' uses the function 'perturb_stats' to return the statistics for each species

perturb_summary<-function(outputs){

	n<-length(outputs[1,])/3

	responseA<-outputs[,1:n]
	responseB<-outputs[,(n+1):(2*n)]
	responseC<-outputs[,(2*n+1):(3*n)]

	statsA<-perturb_stats(responseA)
	statsB<-perturb_stats(responseB)
	statsC<-perturb_stats(responseC)

	globalA<-as.numeric(statsA$global)
	globalB<-as.numeric(statsB$global)
	globalC<-as.numeric(statsC$global)

	peakA<-as.numeric(statsA$peak)
	peakB<-as.numeric(statsB$peak)
	peakC<-as.numeric(statsC$peak)

	troughA<-as.numeric(statsA$trough)
	troughB<-as.numeric(statsB$trough)
	troughC<-as.numeric(statsC$trough)
	
	l<-length(globalA)

	m<-round(l/2)

	output<-c(globalAmin=globalA[1],peakAmin=peakA[1],troughAmin=troughA[1],
			globalBmin=globalB[1],peakBmin=peakB[1],troughBmin=troughB[1],
			globalCmin=globalC[1],peakCmin=peakC[1],troughCmin=troughC[1],
			globalAmid=globalA[m],peakAmid=peakA[m],troughAmid=troughA[m],
			globalBmid=globalB[m],peakBmid=peakB[m],troughBmid=troughB[m],
			globalCmid=globalC[m],peakCmid=peakC[m],troughCmid=troughC[m],
			globalAmax=globalA[l],peakAmax=peakA[l],troughAmax=troughA[l],
			globalBmax=globalB[l],peakBmax=peakB[l],troughBmax=troughB[l],
			globalCmax=globalC[l],peakCmax=peakC[l],troughCmax=troughC[l])

	output
	
}


###################

#Running scripts

#define topologies

topo<-list(topo_3_3_f=matrix(c(1,0,-1,-1,-1,0,0,-1,-1),byrow=T,ncol=3),
	topo_4_3_f=matrix(c(1,1,0,0,-1,-1,1,0,-1),byrow=T,ncol=3),
	topo_5_3_f=matrix(c(-1,-1,-1,-1,-1,0,1,0,-1),byrow=T,ncol=3),
	topo_6_3_f=matrix(c(-1,-1,-1,-1,-1,0,0,-1,-1),byrow=T,ncol=3),

	topo_3_3_r=matrix(c(-1,-1,0,0,1,1,1,0,-1),byrow=T,ncol=3),
	topo_4_3_r=matrix(c(-1,0,1,1,1,0,0,-1,-1),byrow=T,ncol=3),
	topo_5_3_r=matrix(c(-1,-1,0,-1,-1,1,0,-1,-1),byrow=T,ncol=3),
	topo_6_3_r=matrix(c(-1,-1,0,-1,-1,1,1,0,-1),byrow=T,ncol=3),

	topo_AI_3_1=matrix(c(1,0,-1,-1,-1,0,1,1,-1),byrow=T,ncol=3),
	topo_AI_5_2=matrix(c(1,1,-1,-1,-1,0,1,0,-1),byrow=T,ncol=3),
	topo_AI_4_1=matrix(c(1,-1,-1,0,-1,-1,1,0,-1),byrow=T,ncol=3),
	topo_AI_7_1=matrix(c(1,0,-1,0,-1,-1,1,1,-1),byrow=T,ncol=3),
	topo_AI_7_2=matrix(c(1,0,-1,0,-1,-1,1,-1,-1),byrow=T,ncol=3),

	topo_SD_3_1=matrix(c(-1,-1,0,0,1,1,-1,-1,-1),byrow=T,ncol=3),
	topo_SD_5_2=matrix(c(-1,-1,0,1,1,1,0,-1,-1),byrow=T,ncol=3),
	topo_SD_4_1=matrix(c(-1,0,1,-1,1,1,0,-1,-1),byrow=T,ncol=3),
	topo_SD_7_1=matrix(c(-1,0,1,0,1,1,-1,-1,-1),byrow=T,ncol=3),
	topo_SD_7_2=matrix(c(1,0,1,0,1,1,1,-1,-1),byrow=T,ncol=3))


#X number of parameterisations to try for DDI
#Z number of parameterisations to run for inhibition simulation

X<-1:500000

Z<-1000

for(i in 1:18){
	
	fileparm<-paste("~/Desktop/Simulation output/parameters ",names(topo[i]),".txt",sep="")

	fileoutput<-paste("~/Desktop/Simulation output/output ",names(topo[i]),".txt",sep="")

	print(i)
	print("parameters generate")

	##########PART 1 
	##########generate parameters
	
	Y<-numeric()
	
	print(system.time(Y<-aaply(X,1,generate_parameters,top=topo[[i]])))

	colnames(Y)<-c("Type","Da","Aa","Ab","Ac","Db","Ba","Bb","Bc","Dc","Ca","Cb","Cc","phAB","phAC","k2","l","va","vb","vc")

	Y<-Y[!is.na(Y[,1]),]

	parsets<-Y[Y[,1]==1,]


	##########PART 2 SCALE DATASETS
	##########scale datasets
	
	print("parameters scale")
	
	parsetl<-length(parsets[,1])

	number<-Z+50 #increase over Z incase method for scaling fails on some parameter sets 

	sets<-sample(1:parsetl,number)

	par<-cbind(parsets[sets,],sets)

	par<-as.data.frame(par)

	print(system.time(out_table<-daply(par,.(par$sets),pattern_search,target_k2= 0.314719528,target_l= 0.00454855,.parallel=F)))

	#check that there are the required number of parameter sets
	l<-(1:number)[is.na(out_table[,"phAB"])]

	if(length(l)>0){
		out_table<-out_table[-l,]
		print("need some more please...")
		print(names(topo[[i]]))
		print(length(l))
	}
	
	#needed to scale pattern between scaling simulations and inhibition simulations
	out_table[,"Da"]<-out_table[,"Da"]*10*10
	out_table[,"Db"]<-out_table[,"Db"]*10*10
	out_table[,"Dc"]<-out_table[,"Dc"]*10*10
	
	
	OT<-as.data.frame(out_table[1:Z,])
	
	write.table(OT,fileparm,row.names=F,sep="\t")


	##########PART 3 
	##########run sim
	
	print("simulations")

	
	range<-1:100
	print(system.time(pert<-daply(OT[range,],.(OT[range,]$V1),multi_perturb,.parallel=T)))
	
	colnames(pert)<-c("threshA","threshB","threshC",
	
	"A_type","A_strength","A_globalAmax","A_peakAmax","A_troughAmax","A_globalBmax","A_peakBmax","A_troughBmax","A_globalCmax","A_peakCmax","A_troughCmax","A_globalAmid","A_peakAmid","A_troughAmid","A_globalBmid","A_peakBmid","A_troughBmid","A_globalCmidx","A_peakCmid","A_troughCmid","A_globalAmin","A_peakAmin","A_troughAmin","A_globalBmin","A_peakBmin","A_troughBmin","A_globalCmin","A_peakCmin","A_troughCmin",
	
	"B_type","B_strength","B_globalAmax","B_peakAmax","B_troughAmax","B_globalBmax","B_peakBmax","B_troughBmax","B_globalCmax","B_peakCmax","B_troughCmax","B_globalAmid","B_peakAmid","B_troughAmid","B_globalBmid","B_peakBmid","B_troughBmid","B_globalCmid","B_peakCmid","B_troughCmid","B_globalAmin","B_peakAmin","B_troughAmin","B_globalBmin","B_peakBmin","B_troughBmin","B_globalCmin","B_peakCmin","B_troughCmin",
	
	"C_type","C_strength","C_globalAmax","C_peakAmax","C_troughAmax","C_globalBmax","C_peakBmax","C_troughBmax","C_globalCmax","C_peakCmax","C_troughCmax","C_globalAmid","C_peakAmid","C_troughAmid","C_globalBmid","C_peakBmid","C_troughBmid","C_globalCmid","C_peakCmid","C_troughCmid","C_globalAmin","C_peakAmin","C_troughAmin","C_globalBmin","C_peakBmin","C_troughBmin","C_globalCmin","C_peakCmin","C_troughCmin")
	
	write.table(pert,fileoutput,row.names=T,sep="\t")
	
	for(j in 1:((Z/100)-1)){ #print number simulated every 100
		
		range_next<-range+(j*100)
		print(system.time(pert<-daply(OT[range_next,],.(OT[range_next,]$V1),multi_perturb,.parallel=F)))
	
		write.table(pert,fileoutput,row.names=T,sep="\t",append=T,col.names=F)
			
	}
	
}

