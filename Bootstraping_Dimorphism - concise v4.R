####################################################
####Bootstrap Analysis for Sexual Dimorphism
####Developed by Caio Kenup and Marcelo Segall
####20120-09-10
####################################################

###########
###Setup###
###########
{
  
##Loading packages necessary to run the analysis
require(reshape2)
require(plyr)
require(lubridate)
require(dplyr)
require(data.table)
require(magrittr)
require(combinat)
require(stringr)
require(devtools)
require(beepr)
require(readxl)
require(parallel)
options(stringsAsFactors=FALSE)


###Import custom functions from GitHub
source_url("https://raw.githubusercontent.com/KenupCF/Kenup-Repo/master/Quick%20Functions.R")

}


###################################
###Loading and manipulating data###
###################################
{
  
#Define the working directory (where the dataset is)
setwd("C:\\Users\\caiok\\Dropbox\\03-Work\\01-Science\\01-Papers\\00-Sendo Escritos\\Bootstrap_Owl\\Analysis - Lean Version")

#Importing dataset 
dat<- readRDS("dados_para_artigo.rds")

#The dataset has the following variables
#Individual: 
#Pair: The identifier of what pair the individual belongs to 
#(this is filled even if the individual did not belong to a pair)
#Species: The species of the individual
#One column for each song feature measured

######################
####Data wrangling####
######################

#Create a copy of original data set
dat0 <- dat
#Correct double spacing in couple ID
dat$Pair<-gsub(x=dat$Pair,"  "," ")
#Save original individual ID
dat$OrigIndividual<-dat$Individual
#Create new Individual ID (base on Species and Pair)
dat%<>%group_by(Species,Pair)%>%
  mutate(Individual=paste(unique(Species),unique(Pair),1:n(),sep="_"))

###Defining song variables column names
##All the columns, except:
##"Individual","Pair" e "Species"
	songVars<-colnames(dat)[
		!colnames(dat)%in%
		  c("Individual","Pair","Species")]

##Convert all song variables to numeric (this assumes all measured song features are continuous)
dat[,songVars]<-apply(dat[,songVars],2,as.numeric)

##Defining individuals without pairs
dat%<>%						
	group_by(Pair)%>%		
	mutate(lonely=n()<2)	

##Filtering out individuals without pairs
dat<-dat%>%
	filter(!lonely)
}


#############################
###Boostraping simulations###	
#############################

## Parameters for the simulations

nsimu<- 10 #Number of iterations

counter<-30  #Counter variable (after how many iterations will a result be printed)

#Spliting data by species (therefore the analysis can be run for different species at the same time)

dat.split<-split(dat,dat$Species) #Splitting dataset by species

#--------------------------------------------#
#----------Preparing parallelization---------#
#--------------------------------------------#
  noCores<-min(detectCores(),length(dat.split)) #define number of cluster to run paralellization
  cl<-makeCluster(noCores) #make cluster of processing
  split.sims<-split(1:length(dat.split),sort(1:length(dat.split)%%noCores)) #split simulation into 'roughly' equal chunks
  clusterEvalQ(cl, library(plyr))#load package on each cluster
  
  #Loading packages on clusters
  clusterEvalQ(cl, library(dplyr))#load package on each cluster
  clusterEvalQ(cl, library(magrittr))#load package on each cluster
  clusterEvalQ(cl, library(MASS))#load package on each cluster
  clusterEvalQ(cl, library(reshape2))#load package on each cluster
  clusterEvalQ(cl, library(stringr))#load package on each cluster
  clusterEvalQ(cl, library(gtools))#load package on each cluster
  
  #Import all objects so far into the new clusters
  clusterExport(cl,ls())  #Export all objects from main enviroment to each cluster enviroment

start.<-now() #Start time of simulations

#For each core, run a set of species (as defined in split.sims)
simResu<-parLapply(cl, split.sims, function(s){

    
###Create an empty list to hold results (by species)  
all.resu<-list() 

#For each species:
for (sp in s){
		d=dat.split[[sp]] #define data set (the data regarding that particular species)
		spName<-names(dat.split)[sp] #extract species names
		
	  ###Define which variables are available for each species ]
		##### (Keep variables without at least one value that is not an NA)
		songVars_Sp<-songVars[apply(d[,songVars],2,function(x){any(!is.na(x))})] 
		
		#Create empty object to hold results (by variable)
  	resu<-list()
  
  #For each variable (in each species)
  for(v in 1:length(songVars_Sp)){
		
		#Define current variable
		var<-songVars_Sp[v]
		
    #Split data set by couples
		d_pairs<-split(d,d$Pair)

		###Apply for each couple the following criteria:
		d2<-lapply(d_pairs,function(x){		
		  
		  #If the couple has a missing value on the current variable
			if(any(is.na(x[,var]))){ 
				x<-NULL #remove the couple
			#If not
			}else{ 
			  class(x)<-'data.frame'
			  #Arrange individuals from decreasing value of the song variable
			  x<-x[order(x[,var],decreasing = T),] 
			  #Define 'assumed gender' 
			  #(2 for the sex with the greatest value, 
			  #1 for the sex with smallest value) 
        x$Assumed_sex<-2:1
			}
			return(x)
			})
		d2<-rbind.fill(d2)
			
		#Calculate the number of couples remaining
		no.couples<-nrow(d2)/2 
	
		#Empty list for results assuming dimorphism
		dimorfLoop<-list()
		
		#Empty list for results assuming monomorphism
		monomorfLoop<-list() 
	
	##For each iteration (of each variable, in each species)
	for(i in 1:nsimu){ 
	
		###Assembling simulated couples
		
			###Assuming dimorphism:
			
	    #Sample all individuals assumed as sex '1' (with replacement)
			sex1<-sample((d2%>%filter(Assumed_sex==1))$Individual,replace=TRUE) 
	    #Sample all individuals assumed as sex '2' (with replacement)
			sex2<-sample((d2%>%filter(Assumed_sex==2))$Individual,replace=TRUE)
			
			#pairs of individuals will be assigned based on their order
			#(first element of sex1 will be paired with first element of sex2, etc)
			
			###Assuming monomorphism:
			
			#Empty vector to keep individuals
			unk_sex<-rep(NA,nrow(d2))
			
			#first individual in the pair
			#generate a vector with any n/2 individuals from the population, 
			#sampled with replacement
			fst_coup<-sample(d2$Individual,length(unk_sex)/2,replace=TRUE)
				
			#second individual in the pair 
			#any n/2 individuals from the population, sampled with replacement			
			#AS LONG AS they are not already the first individual in the pair
			scd_coup<-sapply(
				fst_coup,function(x){
					sample(d2$Individual[!d2$Individual%in%x],1)})
					
			###Bind individuals by couple
			#(Pair1, Pair1, Pair2, Pair2, PairN, PairN...)
			unk_sex[seq(from=1,to=(length(unk_sex)-1),by=2)]<-fst_coup
			unk_sex[seq(from=2,to=length(unk_sex),by=2)]<-scd_coup			

			
		##Creating the dataset for iterated dimorphism
			
			dimorfLoop[[i]]<-
			   # for each individual of sex '1'
				lapply(1:length(sex1),function(j){
					e<-d2%>%
					#filter the dataframe for his pair
					dplyr::filter(Individual%in%c(sex1[j],sex2[j]))%>% 
          #create iteration and pair identifiers
					dplyr::mutate(Iteration=i,Pair=j)})%>% #depois
				rbind.fill() #bind all individual couples pairs
				
			#reduce number of columns
			dimorfLoop[[i]]%<>%
				select_(.dots=list("Individual","Pair","Iteration",var))
			
		##Creating the dataset for iterated monomorphism
			
			#new dataframe
			monomorfLoop[[i]]<-data.frame( 
			    #With all individuals sampled sorted in pairs
					Individual=unk_sex, 
					Pair=rep(1:(length(unk_sex)/2),each=2)) 
					 
     #Adding song variables of each individual through a merge
			monomorfLoop[[i]]<-merge(
				monomorfLoop[[i]],
				d2%>%select_(.dots=list("Individual",var)),by='Individual')
			
      #create iteration identifiers
			monomorfLoop[[i]]%<>%
				mutate(Iteration=i)%>%
				arrange(Pair)
				
	
	###Print counter of iterations
	if(i%%counter==0){print(paste(spName," - ",var,": iteration no. ",i,sep=""))}
	
	} #end of iteration ('i') loop
	

		##Bind iterations
		dimorfLoop.indiv<-rbind.fill(dimorfLoop)
		monomorfLoop.indiv<-rbind.fill(monomorfLoop)

		###Create pair-iteration combination identifier
		dimorfLoop.indiv%<>%
			dplyr::mutate(Iter_Pair=paste(Iteration,Pair,sep="_"))
		monomorfLoop.indiv%<>%
			dplyr::mutate(Iter_Pair=paste(Iteration,Pair,sep="_"))
		
	###Temporary object splitting by Pair-iteration combination
	temp1<-split(dimorfLoop.indiv,dimorfLoop.indiv$Iter_Pair)
	temp2<-split(monomorfLoop.indiv,monomorfLoop.indiv$Iter_Pair)
	
	
	###CALCULATING RATIO AND DIFFERENCES FOR EACH ITERATION
	
	###Apply for each data.frame on split the following operation
	temp1<-lapply(temp1,function(x){
				x<-x[order(x[,var]),] #sort data.frame by song variable value
				Ratio<-x[,var][2]/x[,var][1] #calculate ratio (greater/smaller)
				Diff<-x[,var][2]-x[,var][1]  #calculate difference (greater - smaller)
				Pair<-unique(x$Pair)		     #define pair
				Iteration<-unique(x$Iteration) #define iteration
				#aggregate all calculated values
				resu<-data.frame(Pair=Pair,Iteration=Iteration,
				                 Ratio=Ratio,Diff=Diff,Var=var)
				return(resu)
			})
	temp2<-lapply(temp2,function(x){
				x<-x[order(x[,var]),] #sort data.frame by song variable value
				Ratio<-x[,var][2]/x[,var][1] #calculate ratio (greater/smaller)
				Diff<-x[,var][2]-x[,var][1]  #calculate difference (greater - smaller)
				Pair<-unique(x$Pair)		     #define pair
				Iteration<-unique(x$Iteration) #define iteration
				#aggregate all calculated values
				resu<-data.frame(Pair=Pair,Iteration=Iteration,
				                 Ratio=Ratio,Diff=Diff,Var=var)
				return(resu)
			})
		
	#rebind the splitted dataframes
	dimorfLoop.couple<-rbind.fill(temp1) 
	monomorfLoop.couple<-rbind.fill(temp2)	
	
	###Summarising per iteration
	dimorfLoop.var<-dimorfLoop.couple%>%
		group_by(Iteration)%>%
		summarise(
			Diff=mean(Diff), 		
			Ratio=mean(Ratio))%>%		
		mutate(Var=var)				
	monomorfLoop.var<-monomorfLoop.couple%>%
		group_by(Iteration)%>%
		summarise(
			Diff=mean(Diff),		
			Ratio=mean(Ratio))%>%	
		mutate(Var=var)					
	
	##save all results for the current variable 
	##(at the couple, individual and variable level)
	resu[[v]]<-
		list(
		  #Couple level
			couple.df=rbind.fill(
				dimorfLoop.couple%>%mutate(Hypothesis="alt"), 
				monomorfLoop.couple%>%mutate(Hypothesis="null"))%>%
				mutate(Var=var,SP=spName),
			#Individual level
			indiv.df=rbind.fill(
				dimorfLoop.indiv%>%mutate(Hypothesis="alt"),
				monomorfLoop.indiv%>%mutate(Hypothesis="null"))%>%
				mutate(Var=var,SP=spName),
			#Variable level
			var.df=rbind.fill(
				dimorfLoop.var%>%mutate(Hypothesis="alt"),
				monomorfLoop.var%>%mutate(Hypothesis="null"))%>%
				mutate(Var=var,SP=spName,no.couples=no.couples)
				)
	class(resu[[v]]$indiv.df)<-"data.frame"
	resu[[v]]$indiv.df$Value<-resu[[v]]$indiv.df[,var]
	resu[[v]]$indiv.df[,var]<-NULL
	
} #end of variable ('v') loop
	
	###naming the list by variable name
	names(resu)<-songVars_Sp
	
	###saving results by species
	all.resu[[sp]]<-
		list(
			indiv.df=rbind.fill(lapply(resu,function(r){r$indiv.df})),
			couple.df=rbind.fill(lapply(resu,function(r){r$couple.df})),
			var.df=rbind.fill(lapply(resu,function(r){r$var.df})))
			
} #end of species ('sp' loop)

###naming list by name of the species
all.resu<-all.resu[!sapply(all.resu,is.null)]
names(all.resu)<-names(dat.split)[s]
return(all.resu)
})
stopCluster(cl) #interrupt parallelization

##Calculate how many hours were ellapsed since start of analysis
hoursElapsed<-abs(as.numeric(difftime(start.,now(),units="hours")))

###Reassembling split object from parallelization
if(paralleling){  
all.resu<-list()
for(i in 1:length(simResu)){
  slots<-seq(from=length(all.resu)+1,by=1,length.out=length(simResu[[i]]))
  all.resu[slots]<-simResu[[i]]
  print(i)
  }
}

#Save image with results
save.image("run.RData")


###################################################################
###Create summary dataframes including all species and variables###
###################################################################
{
bootstrap.indiv<-lapply(all.resu,function(r){
	r$indiv.df
	})%>%rbind.fill
	
bootstrap.couple<-lapply(all.resu,function(r){
	r$couple.df
	})%>%rbind.fill
	
bootstrap.var<-lapply(all.resu,function(r){
	r$var.df
	})%>%rbind.fill	
}	

################################
###Generating summary metrics###
################################

summ1<-bootstrap.couple%>%
	dplyr::group_by(SP,Var,Hypothesis,Iteration)%>%
	dplyr::summarise(Diff=mean(Diff),
	                 Ratio=mean(Ratio),
	                  no.couples=n())%>%
	dplyr::group_by(SP,Var,Hypothesis)%>%
	dplyr::mutate(
		     meanDiff=mean(Diff),
	       meanRatio=mean(Ratio))%>%
  ungroup()
  
finalSummary<-summ1%>%
  #For each Variable in each studied Species
	dplyr::group_by(SP,Var)%>%
	dplyr::summarise(
	 # Calculate a p-value for the observed ratios and differences between sexes
	 # (proportion of iterations where 
	 # the mean observed value for dimorphism is greater than 
	 # the expected under lack of dimorphism)
	 p_ratio_mean=1-((sum(meanRatio[Hypothesis=="alt"]
	                      >=Ratio[Hypothesis=="null"])+1) /
	               (sum(Hypothesis=="null")+1)),
	  p_diff_mean=1-((sum(meanDiff[Hypothesis=="alt"]
	                      >=Diff[Hypothesis=="null"])+1) /
	              (sum(Hypothesis=="null")+1)),
	  no.couples=unique(no.couples))%>%
  arrange(p_diff_mean)

head(finalSummary)

# Export final outputs in a .csv file
write.csv(finalSummary,"results.csv")