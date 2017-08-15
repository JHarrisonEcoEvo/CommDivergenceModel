#------------------------------------------#
#-----BOILERPLATE FOR LOADING PACKAGES-----#
#------------------------------------------#

# install.packages("MCMCpack")
# install.packages("ecodist")
# install.packages("VGAM")
library(MCMCpack) 
library(ecodist)
#library(VGAM)

#setting seed to ensure consistency
set.seed(14567)
options(scipen = 99)

#Stuff to revisit when all said and done and make sure we like
#1. the divergence function computes divergences on summed abundances for communities, not an average for all pairwise distances between individuals
#this may be ideal, but it is worth thinking about more.

#--------------------------------#
#-----BIG BLOCK OF FUNCTIONS-----#
#--------------------------------#


#functions to generate communities using a Dirichlet process. Parameters are:
# 	indiv =  number of individuals in each community, 
# 	m = number of microbe taxa in each individual
#	  numcom = the number of communities in total
# 	abund_microbe =the maximum abundance of a microbe
#	  parameter = clustering parameter for Dirichlet process


#GenerateSame makes a series of communities (k locations/sites/communities)
#each with n individuals. Individuals within a community start with the same microbial assemblage.
#But, communities differ.

#for debugging
# m = 5
# abund_microbe

generateSame = function(indiv, microbes, abund_microbe,parameter){
  x = list()
  y = list()
  
  Hout = 	dirichletprocess(microbes, abund_microbe, parameter)
  
  for(i in 1:indiv){
    y[[i]] = Hout
  }
  #assign new community to a list element in the meta-community
  x[[1]] = y
  return(x)
}

############################################################################
#Function to build microbial communities using a Dirichlet process
#inspired by this article (https://en.wikipedia.org/wiki/Dirichlet_process)
############################################################################

#parameters inherited from the generateSame function are:
#				m = the number of microbial taxa, 
#				abund_microbe = summed abundance of microbes, paramter = conc. parameter.
#       m=100
#       abund_microbe = 1000
#       parameter=10

dirichletprocess = function( m, abund_microbe, parameter){
  
  #define concentration parameter
  alpha = parameter
  
  #Make a uniform probability distribution of the same length as the number of microbes
  H = rep(1/m, m)
  
  #Make a vector for the new community we are building. 
  #Its length is determined by the length of the input distribution
  D = vector(mode="numeric", length = length(H))
  
  #Make a vector of probabilities from the input distribution (H). 
  #take the maximum value of this vector and assign a one to the corresponding position in D
  #This adds a one to a random microbe to start the community
  newH  = as.vector(rdirichlet(1,H))
  D[which(newH == max(newH))] = 1
  
  #Run the process for time specified by "abund_microbe"
  #Abund_microbe is the maximum abundance of any given microbe.
  #this implements the Dirichlet process on D, building the comm. for the number of steps specified by "abund_microbe"
  
  for(i in 1:abund_microbe){
    newMicrobe = alpha/(alpha + i-1)
    oldMicrobe = D/(alpha + i-1)
    
    #As a reminder, the above sum to 1, so it works for a probability distribution. 
    #sum(newMicrobe, oldMicrobe)
    
    #bind those vectors together and sample that Dirichlet distribution
    dir_draws = rdirichlet(1,append(oldMicrobe, newMicrobe))
    
    #find max of the draws again. If the max element is a previously present microbe add 1 to its abundance, if not, then randomly add 1 to a previously unobserved microbe. 
    MicrobeIncreasing = which(dir_draws == max(as.vector(dir_draws)))
    
    #This conditional works on what element of MicrobeIncreasing is greatest. This is the microbe getting an additional observation. 
    #If the element chosen is the last element in dir_draws, then it means we need to assign a 1 to a new
    #element in D that was previously 0. A new microbe is now being observed. If the element chosen is NOT
    #the last element, then we add one to a previously observed microbe
    
    if(MicrobeIncreasing == (1+length(H))){
      #this adds a one to a random microbe that was previously unobserved. But, only operates if there is a
      #prev. unobserved microbe
      if(length(which(D == 0)) > 1 ){
        D[sample(which(D == 0),1)] = 1
      }else{next}}
    else{
      #add a microbe to whichever got choosen by the rdirichlet function above
      D[MicrobeIncreasing] = D[MicrobeIncreasing] + 1	
    }
  }
  
  return(D)
}

#function to do dirichlet/other replacement of individuals with new microbes
#Parameters are: 
#community = a list representing a community that we are replacing individuals inside

#community = sandbox
replacement = function(community){  	
  
  #choose a community to replace
  replaced_comm=round(runif(1, 1, length(community)))
  
  #choose an individual at random to replace
  replaced=round(runif(1, 1,length(community[[replaced_comm]])))
  
  #Replace individual in community with zeros to facilitate summing at next step
  community[[replaced_comm]][[replaced]] = rep(0,length(community[[replaced_comm]][[replaced]]))			
  
  #summing function (calculate abundance of each microbe in a community and then divide by number of hosts)
  #consider how this average calculation may affect things when lots of individuals are present
  dirichletVector = Reduce("+",community[[replaced_comm]])/individuals
  
  community[[replaced_comm]][[replaced]] = round(rdirichlet(1, dirichletVector)*(abund_microbe))
  
  return(community)
}

#Function to compute pairwise divergence between two communities.
#IMPORTANT: this sums abundances across individuals from each community, 
#it does not compute all possible pairwise differences and extract an average

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

#######################################################################
#RUN MODEL. This wrapper function incorporates all the above functions
#######################################################################
model = function(sandbox,distancemetric){
  require("MCMCpack")
  #clearing output list(holds mean divergence)
  divOut = NULL
  z = max(plotpoints)
  k=0
  m=0
  while(k<=z){
    sandbox = replacement(sandbox) 
    if (k %in% plotpoints){
      m=m+1
      print(paste("Sampling Divergence at Step ",k,"/",max(plotpoints)))
      divOut[[m]] = divergence(sandbox[[1]],sandbox[[2]], as.character(distancemetric))
      #outputs average straight divergence
    }
    k=k+1
  }
  return(list(divOut, sandbox))
  #corresponds to out
}


#----------------------------------#
#-----VERY IMPORTANT VARIABLES-----#
#----------------------------------#

#num of communities
communities = 1

#number of individuals in each community
individuals = 100

#number of microbe slots per individual
microbes = 1000

#max amt of microbes per microbe slot 
abund_microbe = 1000000	

#parameter for the Dirichlet function. Higher numbers create less 
parameter = 10

#points at which we calculate the divergence
plotpoints = seq(from = 0, to = 2000, by=10)


sandbox = generateSame(10,1,10000,1)
sandbox[[2]] = NA
sandbox[[2]] = sandbox[[1]]
out = model(sandbox,"bray-curtis")
plot.new()
title(main="Drift Only",xlab="Time Steps", ylab="Drift ")
plot.window(xlim = c(0,max(plotpoints)), ylim = c(0,1), xaxt="n", yaxt="n")
box()
axis(2, at = c(0,0.5,1), labels=c(0,0.5,1))
axis(1, at = c(0,max(plotpoints)/2,max(plotpoints)), labels = c(0,max(plotpoints)/2,max(plotpoints)))
lines(plotpoints,out[[1]], type = "s", col = "red")


#use a riemann sum with column 1 of out