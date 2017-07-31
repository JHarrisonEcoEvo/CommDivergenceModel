
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

generateSame = function(indiv, microbes, numcom, abund_microbe,parameter){
  x = list()
  y = list()
  
  for(j in 1:numcom){
    Hout = 	dirichletprocess(microbes, abund_microbe, parameter)
    
    for(i in 1:indiv){
      y[[i]] = Hout
    }
    #assign new community to a list element in the meta-community
    x[[j]] = y
  }
  return(x)
}

#GenerateDiff makes a series of communities (k locations/sites/communities)
#each with n individuals. Individuals within a community start with DIFFERENT microbial assemblage.
#And, communities differ.
generateDiff = function(indiv, microbes, numcom, abund_microbe,parameter){
  x = list()
  y = list()
  for(j in 1:numcom){
    for(i in 1:indiv){
      y[[i]] = dirichletprocess(microbes, abund_microbe, parameter)
    }
    #assign new community to a list element in the meta-community
    x[[j]] = y
  }
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
model = function(sandbox,mode1,sensitivity,stoptype,distancemetric){
  running = TRUE
  require("MCMCpack")
  #clearing output list(holds mean divergence)
  divOut = NULL
  z = max(plotpoints)
  k=0
  counter = 0
  endstep = NA
  #if we aren't doing it the smart way, set the endstep to the max plotpoints
  if(mode1 == "normal"){
    endstep = max(plotpoints)
  }
  while(k<=z && running == TRUE){
    sandbox = replacement(sandbox) 
    if (k %in% plotpoints){  	
      print(paste("Sampling Divergence at Step ",k,"/",max(plotpoints)))
      div=NA
      m=1
      for(i in 1:length(sandbox)){
        for(j in 1:length(sandbox)){
          if(i != j){
            div[m] = divergence(sandbox[[i]],sandbox[[j]], as.character(distancemetric))
            m=m+1
            #First if statement resets the counter to zero if we totally diverged, then converged
            if(div[m-1] != 1 && counter > 1){
              counter = 0
            }
            else if(as.character(mode1) == "smart" && div[m-1] == 1){
              #iterate the counter if we have reached total divergence. 
              counter = 1 + counter
              #run until we hit sensitivity (num of steps post total divergence)
              #save when total divergence happened
              if(counter < sensitivity){
                break
              }
              else{
                running = FALSE
                print(paste("autostopping at ",endstep))
                endstep = length(div) - sensitivity
                if(mode == "assume"){
                  divOut[div[m]:max(plotpoints)] = 1
                  #ADD A BUNCH OF 1S BECAUSE ASSUME
                  running = FALSE
                }
                if(mode == "break"){
                  #do nothing
                  running = FALSE
                }
              }
            }
          }
        }
      }
      #outputs average divergence in ascending order
      divOut[length(divOut)+1] = mean(div)
    }
    k=k+1
  }
  return(list(divOut, sandbox, endstep))
  #corresponds to out
}


#----------------------------------#
#-----VERY IMPORTANT VARIABLES-----#
#----------------------------------#

#num of communities
communities = 10

#number of individuals in each community
individuals = 100

#number of microbe slots per individual
microbes = 1000

#max amt of microbes per microbe slot 
abund_microbe = 10000	

#parameter for the Dirichlet function. Higher numbers create less 
parameter = 10

#points at which we calculate the divergence
plotpoints = seq(from = 0, to = 200000, by=1000)
parameterSet=c(1,10,100,300, 500, 1000, 2000)
indiv = c(10,20,50)
microbes = c(50,200,500)
colors = list("red","orange","green","aquamarine","blue","purple","black","burlywood","cadetblue","chartreuse","chocolate","coral","cornflowerblue","cornsilk","cyan","darkblue","darkgoldenrod","darkolivegreen")
#---------------------------------------#
#-----GENERATING INITIAL CONDITIONS-----#
#---------------------------------------#

colorcount=0
#doing iterating parameters on ONE graph
p=0
#-------------------------------------#
#----GRAPHING/SIMULATION OF MODEL-----#
#-------------------------------------#
k=1
j=1
pdf(file="Output.pdf", width=8.5, height=11)
cols = c("black", "blue", "green", "orange", "red", "cadetblue")
par(mfrow=c(length(indiv),length(microbes)))
for(j in 1:length(microbes)){
  for(k in 1:length(indiv)){
    #making new plot for the parameter changing
    plot.new()
    plot.window(xlim = c(0,max(plotpoints)), ylim = c(0,1), xaxt="n", yaxt="n")
    title(main=paste("Individuals ", indiv[[k]], "Microbes ", microbes[[j]]),xlab="Time Steps", ylab="Divergence")
    box()
    axis(2, at = c(0,0.5,1), labels=c(0,0.5,1))
    axis(1, at = c(0,max(plotpoints)/2,max(plotpoints)), labels = c(0,max(plotpoints)/2,max(plotpoints)))
    for (p in 1:length(parameterSet)){
      colorcount = colorcount + 1
      sandbox = generateSame(indiv[[k]], microbes[[j]], communities, abund_microbe, parameterSet[p])
      out = model(sandbox, "smart","assume" ,5, "bray")
      print(paste("Finished model with individuals: ", indiv[[k]],", microbes: ", microbes[[j]], ", parameter: ", parameterSet[p]))
      lines(plotpoints[1:length(plotpoints)], out[[1]],col = paste(colors[[p]]),type = "s", lwd=2)
    }
  }
}
dev.off()
