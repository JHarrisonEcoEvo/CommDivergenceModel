# install.packages("MCMCpack")
# install.packages("ecodist")
# install.packages("VGAM")
library(MCMCpack) 
library(ecodist)
library(VGAM)

#setting seed to ensure consistency
set.seed(14567)
options(scipen = 99)

#Stuff to revisit when all said and done and make sure we like
#1. the divergence function computes divergences on summed abundances for communities, not an average for all pairwise distances between individuals
#this may be ideal, but it is worth thinking about more.
#2. We should think more about the zero inflated Poisson used for H in the Dirichlet process. 



#--------------------------------#
#-----BIG BLOCK OF FUNCTIONS-----#
#--------------------------------#


#functions to generate communities using a Dirichlet process. Parameters are:
# 	indiv =  number of individuals in each community, 
# 	m = number of microbe taxa in each individual
#	numcom = the number of communities in total
#	abund_microbe =the maximum abundance of a microbe
#	parameter = clustering parameter for Dirichlet process


#generate with identical values for all individuals in all communities

#for debugging
#m = microbes

generateSame = function(indiv, m, numcom, abund_microbe,iterations,parameter){
  x = list()
  y = list()
  
  #outside loop uses a zero-inflated Poisson to build a probability distribution sampled by the Dirichlet process
  #inside loops call the Dirichlet process, and then assign it to all individuals in a community y, which is a list element in x the meta-community
  
  for(j in 1:numcom){
  	
  	#pstr0 is the probability of a structural zero, so higher numbers result in sparser data
   # hostzero = round(rdirichlet(1, rzipois(m, lambda = 1, pstr0 = .5))*(abund_microbe))
    
    for(k in iterations){
    	Hout = 	dirichletprocess(m, abund_microbe, parameter)
    }
	for(i in 1:indiv){
		y[[i]] = Hout
	}
	
	#assign new community to a list element in the meta-community
	x[[j]] = y
		 	 
  }
	return(x)
}


generateDiff = function(indiv, m, numcom, abund_microbe,iterations,parameter){
	  x = list()
  	  y = list()
  for(j in 1:numcom){
  	for(i in 1:indiv){
	  	hostzero = round(rdirichlet(1, rzipois(m, lambda = 1, pstr0 = .5))*(abund_microbe))
	
	    for(k in iterations){
	    	Hout = 	dirichletprocess(hostzero,parameter)
	    	hostzero=Hout
	    }	
		y[[i]] = Hout
	}
	x[[j]] = y	 
  }
	return(x)
}


#Function to build microbial communities using a Dirichlet process
#	H is the probability base distribution ( we shall have it be the zero inflated poisson)
#	parameter is the scaling parameter, when it is larger there is a higher probability of sitting at a new table...adding a draw to a different microbe. When lower, more likely to add draw to most abundant microbe
#inspired by this article (https://en.wikipedia.org/wiki/Dirichlet_process)

#parameters are m = the number of microbial taxa, abund_microbe = summed abundance of microbes, paramter = conc. parameter.
#these are inherited from the generateSame function


dirichletprocess = function( m, abund_microbe, parameter){
	
	#define concentration parameter
  	alpha = parameter
  	
  	#Make a uniform probability distribution of the same length as the number of microbes
  	H = rep(1/m, m)

	#Make a vector for the new community we are building. 
	  #Its length is determined by the length of the input distribution
	  D = vector(mode="numeric", length= length(H))
	  
	  #Make a vector of probabilities from the input distribution (H). 
	  #take the maximum value of this vector and assign a one to the corresponding position in D
	  
	  newH  = as.vector(rdirichlet(1,H))
	  D[which(newH == max(newH))] = 1
	  
	  for(i in 1:length(H)){
		  newMicrobe = alpha/(alpha + i-1)
		  oldMicrobe = D/(alpha + i-1)
		  
		  #As a reminder, the above sum to 1, so it works for a probability distribution. 
		  #sum(newMicrobe, oldMicrobe)
		  
		  #bind those vectors together and sample that Dirichlet distribution
		  dir_draws = rdirichlet(1,append(oldMicrobe, newMicrobe))
		  
		  #find max of the draws again. If the max element is a previously present microbe add 1 to its abundance, if not, then randomly add 1 to a previously unobserved microbe. 
		  MicrobeIncreasing = which(dir_draws == max(as.vector(dir_draws)))
		  
		  #This conditional works on what element of MicrobeIncreasing is greatest. This is the microbe getting an additional
		  #read. If the element chosen is the last element in dir_draws, then it means we need to assign a 1 to a new
		  #element in D that was previously 0. A new microbe is now being observed. If the element chosen is NOT
		  #the last element, then we add one to a previously observed microbe
		  
		  if(MicrobeIncreasing == (1+length(H))){
		  	D[sample(which(D == 0),1)] = 1
		  }else{
		  	D[MicrobeIncreasing] = D[MicrobeIncreasing] + 1	
		  }
	  }
	  return(D)
}
 	
# # dirichletprocess = function(H, parameter){
  # D=NA
  # for (i in 1:length(H)){
    # D[[i]]=0
    # #sit at new table; slash a new microbe is born....new element in vector
    
    # #there is a problem with the first part of this conditional. It takes the value from H and dumps it into D
    # #it should instead assign the value of 1 to that element of D, or add a 1 to whatever value was in D at that position
    # #it should also only overwrite elements of D that equaled zero, so we need another conditional
    # #to enforce this
    
    # #the runif calculates a random probability. If this value is less than 
    # #parameter/parameter + i-1 then add a one to a random value that previously equaled zero
    # if(runif(1, 0,1)< (parameter/(parameter+i-1))){
      
      # #pick a random integer from 1 through the length of the input probability distribution
      # pick = round(runif(1, 1, length(H)))
      
      # #replace the 0 in D with the value from H corresponding to the randomly selected pick
      # D[[i]] = H[[pick]]  
    # }
    # #add to an existing vector...assign the next read to a preexisting microbe
    
    # #the runif command here is redundant...because this should always run whenever we get a no from the first conditional
    # #we also need to update H here
    # else if(runif(1, 0,1)< (H[[i]]/(parameter+i-1))){
      # D[[i]] = H[[i]] +1 
      
    # }
  # }
  # return(D)
# }

#function to do dirichlet/other replacement of individuals with new microbes
#Parameters are: 
#community = a list representing a community that we are replacing individuals inside

#community = sandbox
replacement = function(community){  	
  	
	#choose a community to replace
  	replaced_comm=round(runif(1, 1, length(community)))
  	
  	#choose an individual at random to replace
	replaced=round(runif(1, 1,individuals))
	  
	 #Replace individual in community with zeros to facilitate summing at next step
	 community[[replaced_comm]][[replaced]] = rep(0,length(community[[replaced_comm]][[replaced]]))			
	  	  
	 #summing function (calculate abundance of each microbe in a community and then divide by number of hosts)
	 #consider how this average calculation may affect things when lots of individuals are present
	 dirichletVector = Reduce("+",community[[replaced_comm]])/individuals
	  
	 community[[replaced_comm]][[replaced]] = round(rdirichlet(1, dirichletVector)*(abund_microbe))

	 return(community)
}

#Function to compute pairwise divergence between two communities.
#IMPORTANT: this sums abundances across individuals from each community, it does not compute all possible pairwise differences
#and extract an average

#inputs are lists of lists

#for debugging
# comm1 = sandbox[[1]]
# comm2 = sandbox[[2]]
# method2 = "bray"
divergence = function(comm1, comm2, method2){
  
	  out = distance(rbind(Reduce("+",comm1), Reduce("+",comm2)), method=method2)
  
  
  #compute distance metric between two communities
  # p=1
  # for(m in 1:length(comm1)){
  #   for(n in 1:length(comm2)){
  #     out[[p]] = distance(rbind(as.data.frame(comm1[[m]]), as.data.frame(comm2[[n]])), method=method2)
  #     p = p+ 1
  #   }
  # }
  # 
  return(out)
}


model = function(sandbox){
  require("MCMCpack")
  #clearing output list(holds mean divergence)
  divOut = NULL
  z = max(plotpoints)
  k=0
 	while(k<z){
   		sandbox = replacement(sandbox) 
    	k=k+1
    	print(k)
	if (k %in% plotpoints){  	
	    div=NA
	    m=1
	    for(i in 1:length(sandbox)){
	        for(j in 1:length(sandbox)){
	          if(i != j){
	            div[m] = divergence(sandbox[[i]],sandbox[[j]], "bray")
	            m=m+1
	          }else{next}
	        }
	     }
	      #outputs average divergence in ascending order
	      divOut[length(divOut)+1] = mean(div)
	    }
    }
  return(list(divOut, sandbox))
}


#----------------------------------#
#-----VERY IMPORTANT VARIABLES-----#
#----------------------------------#

communities = 10
individuals = 2
microbes = 1000
abund_microbe = 10000	
parameter = 0.5
iterations = 100

#points at which we calculate the divergence
plotpoints = seq(from = 0, to = 4000, by=50)

#---------------------------------------#
#-----GENERATING INITIAL CONDITIONS-----#
#---------------------------------------#

sandbox = generateSame(individuals, microbes, communities, abund_microbe, iterations, parameter)
#sandbox = generateDiff(individuals, microbes, communities, abund_microbe, iterations, parameter)

#save original communities
sandbox_original = sandbox

#save communities post running model
out = model(sandbox)

#plot divergence versus time
plot(plotpoints[2:length(plotpoints)], out[[1]], ylab="Divergence", xlab = "Time step")


length(out[[1]])