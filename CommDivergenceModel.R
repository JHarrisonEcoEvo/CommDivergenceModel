
require(MCMCpack) 
require(ecodist)

#setting seed to ensure consistency
set.seed(14567)
options(scipen = 99)

#read life table data

lifehistorytable = read.csv("./data/hypotheticalShortlived.csv")
lifehistorytable = read.csv("./data/humanfemale.csv")


#Stuff to revisit when all said and done and make sure we like
#1. the divergence function computes divergences on summed abundances for communities, not an average for all pairwise distances between individuals
#this may be ideal, but it is worth thinking about more.
#Does it make sense to try and use a D. process to replace the community in a replaced individual

#---------------- -------------#
#-----High level functions-----#
#------------------------------#
#These call lower level functions that are specified below
#"model" runs the model, the highest level function
#generateSame and generateDiff generate initial simulated communities
#replacement - replaces individuals at appropriate time steps

#debugging stuff, will delete when done
# commSame = TRUE
# distancemetric="bray-curtis"
# numComm=2
# numIndiv=5
# minInd = 2
# numMicrobes=10
# microbeAbund=10000
# conc.par=10
# timesteps = 1000

model <- function(lifehistorytable, 
                 commSame = TRUE,
                 dispersal = TRUE,
                 minInd = 8,
                 distancemetric="bray-curtis", 
                 numComm=10, 
                 numIndiv=10, 
                 numMicrobes=100, 
                 microbeAbund=10000, 
                 conc.par=10, 
                 timesteps = 1000){
  #Parameters are:
  #commSame - Boolean to specify if communities should be made up of identical individuals or different individuals
  #distancemetric - any distance metric accepted by the distance function of ecodist
  #numComm - number of communities
  #numIndiv - number of individuals in each community
  #numMicrobes - number of microbe slots per individual
  #microbeAbund, abundance of microbes, summed across taxa
  #conc.par - Concentration parameter
  #timesteps - how long to run model
  
  
  # calculate timestep to calculate divergence for plotting
  plotpoints=seq(0,timesteps, by=timesteps/10) 
  
  #obtain life history data
  agesThetas <- formatLifeTable(lifehistorytable)
  
  #vector to hold output
  divOut  <- NA
  divOutVar <- NA
  #counter
  k <- 0
  
  #generate initial simulated communities, of specified parameters
  if(commSame == TRUE | missing(commSame)){
    community <- generateSame(numIndiv, numMicrobes, numComm, microbeAbund,conc.par, agesThetas)
    print("Generating meta-community using commSame == TRUE")
    print("Individuals within a community will start with identical microbiomes, but communities will differ")
  }
  if(commSame == FALSE){
    community <- generateDiff(numIndiv, numMicrobes, numComm, microbeAbund,conc.par, agesThetas)
    print("Generating meta-community using commSame == FALSE")
    print("Individuals within a community will start with different microbiomes")
  }
  
  #run model for appropriate number of time steps
  repeat{
    k <- k+1
    #replace members of communities at each iteration
    community <- replacement(community, numComm, numMicrobes, numIndiv, microbeAbund, agesThetas); community
    
    #let member(s) of communities disperse at each iteration
    if(dispersal == TRUE){
      community <- disperse(community,minInd)
    }
    
    #calculate divergence at time steps specified by plotpoints
      if (k %in% plotpoints){ 
        #vector to hold distance indices
        div <- NULL
        m <- 1
        print(paste("Sampling Divergence at Step ",k,"/",max(plotpoints)))
        for(i in 1:length(community)){
          for(j in 1:length(community)){
            if(i != j){
              #compute pairwise distance between all communities
              div[m] <- divergence(community[[i]][[1]],community[[j]][[1]], as.character(distancemetric))
              m  <- m+1
            }
          }
        }
        #Calculate average pairwise distance among all communities
        #also calculate variance in pairwise distance among communities
        divOut[length(divOut)+1] <- mean(div)
        divOutVar[length(divOutVar)+1] <- var(div)
      }
    if (k >= max(plotpoints)){
      return(list(divOut,divOutVar,community))
      break
    }
  }
}

######################################################################################################
#GenerateSame makes a series of communities (k locations/sites/communities) each with n individuals. #
#Individuals within a community start with the same microbial assemblage, but communities differ.    #
######################################################################################################
#for debugging
# m = 5
# abund_microbe

generateSame <- function(numIndiv, numMicrobes, numComm, microbeAbund,conc.par, agesThetas){
  # 	indiv =  number of individuals in each community, 
  # 	m = number of microbe taxa in each individual
  #	  numcom = the number of communities in total
  # 	abund_microbe =the maximum abundance of a microbe
  #	  parameter = clustering parameter for Dirichlet process
  
  x = list()
  y = list()
  z = list()
  for(j in 1:numComm){
    Hout <- dirichletprocess(numMicrobes, microbeAbund, conc.par)
    ages <- round(runif(numIndiv, min = 0,max = max(agesThetas[[1]])))
    for(i in 1:numIndiv){
      y[[i]] <- Hout
    }
    #assign new community to a list element in the meta-community
    x[[j]] <- y
    z[[j]] <- list(x[[j]],ages)
  }
  return(z)
}

##################################################################################################
#GenerateDiff makes a series of communities (k locations/sites/communities)each with n individuals. 
#Individuals within a community start with DIFFERENT microbial assemblages.
##################################################################################################

generateDiff = function(numIndiv, numMicrobes, numComm, microbeAbund,conc.par, agesThetas){
  x = list()
  y = list()
  z = list()

  for(j in 1:numComm){
    for(i in 1:numIndiv){
      y[[i]] <- dirichletprocess(numMicrobes, microbeAbund, conc.par)
    }
    #assign new community to a list element in the meta-community
    x[[j]] <- y
    ages <- round(runif(numIndiv, min = 0,max = max(agesThetas[[1]])))
    z[[j]] <- list(x[[j]],ages)
  }
  return(z)
}

##########################################################################################
#Replacement function - replaces individuals in each community as a function of their age#
##########################################################################################
#for debugging
  # sandbox = generateSame(10, 100, 10, 10000, 3)
  # community = sandbox

replacement = function(community, numComm, numMicrobes, numIndiv, microbeAbund, agesThetas){  	
  
  #increase age for all individuals in all communities by 1
  for(l in 1:numComm){
    community[[l]][[2]] <-  community[[l]][[2]]+1
  }
  
  #replace individuals in each community as a function of their age, do for all communities
  for(l in 1:numComm){
    #First, compute probability of death for each individual in the community 
    #sample from a binomial distribution parameterized via the input life history table
    dpcount <- 1
    deathProb <- 1
    for(h in community[[l]][[2]]){
      #this conditional sets any individuals that age past the max known age to have a 100% probability of death
      if(h >= max(agesThetas[[1]])){
        deathProb[dpcount] <- 1
        dpcount <- dpcount+1
      }else{
        deathProb[dpcount] <- agesThetas[[2]][h == agesThetas[[1]]]
        dpcount <- dpcount+1
      }
    }
    
    #make a vector of individuals that should be replaced and pass to the remainder of this function
    deadAlive = NA
    for(h in 1:length(deathProb)){
      deadAlive[h] <- rbinom(1, 1, prob = deathProb[h])
    }
    replaced <- which(deadAlive == 1)
    
    #Replace individual in community with zeros to facilitate summing at next step
    for(rp in replaced){
      community[[l]][[1]][[rp]] <- rep(0,numMicrobes)
      community[[l]][[2]][[rp]] <- 0 #Reset the age
    }
    #summing function (calculate abundance of each microbe in a community and then divide by number of hosts)
    #consider how this average calculation may affect things when lots of individuals are present
    if(length(replaced) < numIndiv){
      dirichletVector = Reduce("+",community[[l]][[1]])/numIndiv
      for(rp in replaced){
        #note that rounding can make this vector sum to < numMicrobe, but not by much
        community[[l]][[1]][[rp]] <- as.vector(round(rdirichlet(1, dirichletVector)*(microbeAbund)))
      }
    }else{
      print("Error: all individuals in a community died. Run with more individuals and this might not happen.")
    }
  }
  return(community)
}

#function that implements dispersal. 
  #Assumptions: 1. individuals are more likely to disperse OUT of larger populations
  #             2. There is a uniform probability regarding WHERE these individuals go
  #             3. Communities don't have to be the same size obviously, but do have a minimum size
  #             4. There is a unform prob. distribution for which indidivual in the choosen community disperses. Might be worth playing with this more.

disperse <- function(community, minInd){
  #remember that sapply returns a vector, while lapply returns a list
  numindiv_eachComm <- sapply(community, FUN= function(community){length(community[[1]])})
  
  #take a deviate from a Dirichlet distribution parameterized using the number of individuals in each community.
  #this generates higher probabilities for communities with more individuals.
  disperseProbs <- rdirichlet(1, numindiv_eachComm)
  #pick the largest probability. This will be the community that has the dispersing individual
  CommWithDisperser <- which(disperseProbs == max(disperseProbs))
  
  #choose the destination community and individual to replace
  destinationComm <- round(runif(1,min=1, max=length(community)))
  
  #make sure we actually disperse (instead of replacing an indiv. in the dispersal community)
  if(destinationComm == CommWithDisperser){
    repeat{
      destinationComm <- round(runif(1,min=1, max=length(community)))
      #test that there are enough individuals in this community. 
      #we don't want dispersal to cause a community to go extinct eventually, so we specifiy a min. number of individuals allowable
      if(minInd >= length(community[[CommWithDisperser]][[1]])){
        repeat{
          disperseProbs <- rdirichlet(1, numindiv_eachComm)
          CommWithDisperser <- which(disperseProbs == max(disperseProbs))
          if(minInd >= length(community[[CommWithDisperser]][[1]])){
            break 
          }
        }
      }
      if(destinationComm != CommWithDisperser){
        break
      }
    }
  }
  
  #pick the actual individual within the community to disperse. This individual is choosen using a uniform prob. distribution.
  #extending from 1 to how ever many individuals are in the community
  disperser <- round(runif(1,min=1, max=length(community[[CommWithDisperser]][[1]])))
    #The disperser: community[[CommWithDisperser]][[1]][[disperser]]
  
  #calculate 1 plus the number of individuals in the destination community, so we can append the disperser to that community
  destIndiv <- length(community[[destinationComm]][[1]])+1
  community[[destinationComm]][[1]][[destIndiv]] <- community[[CommWithDisperser]][[1]][[disperser]]
  #dont forget to add its age
  community[[destinationComm]][[2]][[destIndiv]] <- community[[CommWithDisperser]][[2]][[disperser]]
  print(paste("Individual ", disperser, " from community ", CommWithDisperser, " is dispersing. Yay, go forth!"))
  
  #remove the dispersering individual from the original community
  community[[CommWithDisperser]][[1]][[disperser]] <- NULL
  community[[CommWithDisperser]][[2]] <- community[[CommWithDisperser]][[2]][-disperser]
  
  return(community)
  }  

#---------------- -------------#
#-----Low level functions------#
#------------------------------#

############################################################################
#Function to build microbial communities using a Dirichlet process
#inspired by this article (https://en.wikipedia.org/wiki/Dirichlet_process)
############################################################################

dirichletprocess = function(numMicrobes, microbeAbund, conc.par){
  #numMicrobes - is the number of microbial taxa
  #microbeAbund - is the total microbial abundance
  #conc.par - is the concentration parameter that controls the evenness of the simulated community
  
  #Make a vector for the new community we are building. 
  D = vector(mode="numeric", length = numMicrobes)
  
  #Add a one to a random microbe to start the community
  D[round(runif(1, min=1, max = numMicrobes))] = 1
  
  #Run the Dirichlet process for time specified by "microbeAbund"
  for(i in 1:microbeAbund){
    newMicrobe = conc.par/(conc.par + i-1)
    oldMicrobe = D/(conc.par + i-1)
    
    #bind those vectors together and sample that Dirichlet distribution
    dir_draws = rdirichlet(1,append(oldMicrobe, newMicrobe))
    
    #find max of the draws again. If the max element is a previously present microbe add 1 to its abundance, if not, then randomly add 1 to a previously unobserved microbe. 
    MicrobeIncreasing = which(dir_draws == max(as.vector(dir_draws)))
    
    #This conditional works on what element of MicrobeIncreasing is greatest. This is the microbe getting an additional observation. 
    #If the element chosen is the last element in dir_draws, then it means we need to assign a 1 to a new
    #element in D that was previously 0. A new microbe is now being observed. If the element chosen is NOT
    #the last element, then we add one to a previously observed microbe
    
    if(MicrobeIncreasing == (1+numMicrobes)){
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

#Function to compute pairwise divergence between two communities.
#IMPORTANT: this sums abundances across individuals from each community, 
#it does not compute all possible pairwise differences and extract an average

divergence <- function(comm1, comm2, distancemetric){
  out <- distance(rbind(Reduce("+",comm1), Reduce("+",comm2)), method=distancemetric)
  return(out)
}

#format life history table
#Initial format assumes two columns, the first with ranges, the second with percent dying within that range
formatLifeTable <- function(filename){
  #expand the range and assign the probability of death for that range to all integers
  ages <- NULL
  thetas <- NULL #vector of survival probabilities to go with each age
  for(i in 1:length(filename[,1])){
    start <- as.numeric(gsub("(\\d+)-\\d+", "\\1",filename[i,1]))
    end <- as.numeric(gsub("\\d+-(\\d+)", "\\1",filename[i,1]))
    Range <- seq(start, end, by=1)
    ages <- c(ages,Range)
    #Note that here we divide the probability of death by the age range, so that the summed probability
    #of an individual dying over a age range matches the life table
    thetas <- c(thetas, rep(filename[i,2]/length(Range), length(Range)))
  }
  return(list(ages, thetas))
}


#########
#Testing#
#########

# #simulate a community using defaults, except for conc. parameter
# out <- model(dat, conc.par=4)
# 
# #this is what community one looks like
# #ten host individuals are present, each with 100 different possible taxa
# out[3][[1]][[2]][[1]]


#-----------------#
#-----Plotting-----#
#-----------------#

#points at which we calculate the divergence
plotpointsv = seq(0,10000, by=10000/10)
minsize = 8
parameterSet=c(3,30,100)
indiv = c(50,100)
microbes = c(50,100)
colors = c("cadetblue1","cadetblue4","darkolivegreen1","darkolivegreen4","indianred1","indianred4","slateblue1","slateblue2","plum1","plum4")
colorcount=0
p=0
k=1
j=1


#pdf(file="Output.pdf", width=8.5, height=11)
par(mfrow=c(length(indiv),length(microbes)))
for(j in 1:length(microbes)){
  for(k in 1:length(indiv)){
    #making new plot for the parameter changing
    plot.new()
    plot.window(xlim = c(0,max(plotpointsv)), ylim = c(0,1), xaxt="n", yaxt="n")
    title(main=paste("Individuals ", indiv[[k]], "Microbes ", microbes[[j]]),xlab="Time Steps", ylab="Divergence")
    box()
    axis(2, at = c(0,0.5,1), labels=c(0,0.5,1))
    axis(1, at = c(0,max(plotpointsv)/2,max(plotpointsv)), labels = c(0,max(plotpointsv)/2,max(plotpointsv)))
    for (p in 1:length(parameterSet)){
      colorcount = colorcount + 1
       out = model(dat, dispersal=FALSE, minInd = 8, numIndiv = indiv[[k]], numMicrobes = microbes[[j]], microbeAbund = 1000, conc.par = parameterSet[p], timesteps = 10000)
      # print(paste("Finished model with individuals (no dispersal): ", indiv[[k]],", microbes: ", microbes[[j]], ", parameter: ", parameterSet[p]))
       lines(plotpointsv[1:length(plotpointsv)], out[[1]],col = paste(colors[[p]]),lwd=2, lty=1)
      # 
      out = model(dat, dispersal=TRUE, numComm = 2, minInd = 8, numIndiv = indiv[[k]], numMicrobes = microbes[[j]], microbeAbund = 1000, conc.par = parameterSet[p], timesteps = 10000)
      lines(plotpointsv[1:length(plotpointsv)], out[[1]],col = paste(colors[[p]]),type = "s", lwd=2, lty=3)
      print(paste("Finished model with individuals (dispersal): ", indiv[[k]],", microbes: ", microbes[[j]], ", parameter: ", parameterSet[p]))

    }
  }
}
#dev.off()
