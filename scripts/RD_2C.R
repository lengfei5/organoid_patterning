

###################
rm(list=ls())
#Running scripts

#define topologies

topo<-list(topo_AI=matrix(c(1,-1,1,-1),byrow=T,ncol=2),
	topo_SD=matrix(c(1,1,-1,-1),byrow=T,ncol=2))


resDir = paste0('../results/RD_topology_test/')
if(!dir.exists(resDir)) dir.create(resDir)

#X number of parameterisations to try for DDI
#Z number of parameterisations to run for inhibition simulation

X<-1:50000
Z<-50

source('Functions_RD_2components_simulation.R')

for(i in 1:2){
	
  # i = 1
  
	fileparm<-paste(resDir,names(topo[i]),".txt",sep="")
	fileoutput<-paste(resDir,names(topo[i]),".txt",sep="")

	print(i)
	print("parameters generate")
  
	##########PART 1 
	##########generate parameters
	
	Y<-numeric()
	
	print(system.time(Y<- aaply(X,1, generate_parameters_2c, top=topo[[i]])))

	colnames(Y)<-c("Type","Da","Aa","Ab","Db","Ba","Bb","ph","k2","l","va","vb")

	Y<-Y[!is.na(Y[,1]),]

	parsets<-Y[Y[,1]==1,]
  
  
	##########PART 2 SCALE DATASETS
	##########scale datasets
	print("parameters scale")
	
	parsetl<-length(parsets[,1])

	number<-Z+50

	sets<-sample(1:parsetl, size = number, replace = FALSE)

	par<-cbind(parsets[sets,],sets)

	par<-as.data.frame(par)

	print(system.time(out_table <- daply(par,.(par$sets), pattern_search_2c, target_k2 = 0.314719528, 
	                                     target_l= 0.00454855, .parallel=T)))
	
	
	l<-(1:number)[is.na(out_table[,"ph"])]
  
	if(length(l)>0){
		out_table<-out_table[-l,]
		print("need some more please...")
		print(names(topo[[i]]))
		print(length(l))
	}
  
	out_table[,"Da"]<-out_table[,"Da"]*10*10
	out_table[,"Db"]<-out_table[,"Db"]*10*10
	
	
	OT<-as.data.frame(out_table[1:Z,])
	
	write.table(OT,fileparm,row.names=F,sep="\t")


	##########PART 3 
	##########run sim
	
	print("simulations")

	
	range<-1:100
	#for inhibition of response run this
	print(system.time(pert<-daply(OT[range,],.(OT[range,]$V1),multi_perturb_2c,.parallel=T)))

	#for inhibition of production run this
	#print(system.time(pert<-daply(OT[range,],.(OT[range,]$V1),multi_perturb_2c_type_2,.parallel=T)))


	colnames(pert)<-c("threshA","threshB",
	
	"A_type","A_strength","A_globalAmax","A_peakAmax","A_troughAmax","A_globalBmax","A_peakBmax","A_troughBmax","A_globalAmid","A_peakAmid","A_troughAmid","A_globalBmid","A_peakBmid","A_troughBmid","A_globalAmin","A_peakAmin","A_troughAmin","A_globalBmin","A_peakBmin","A_troughBmin",
	
	"B_type","B_strength","B_globalAmax","B_peakAmax","B_troughAmax","B_globalBmax","B_peakBmax","B_troughBmax","B_globalAmid","B_peakAmid","B_troughAmid","B_globalBmid","B_peakBmid","B_troughBmid","B_globalAmin","B_peakAmin","B_troughAmin","B_globalBmin","B_peakBmin","B_troughBmin")
	
	write.table(pert,fileoutput,row.names=T,sep="\t")
	
	for(j in 1:((Z/100)-1)){
		
		range_next<-range+(j*100)
		print(system.time(pert<-daply(OT[range_next,],.(OT[range_next,]$V1),multi_perturb_2c,.parallel=T)))
	
		write.table(pert,fileoutput,row.names=T,sep="\t",append=T,col.names=F)
	
		
	}
	
}