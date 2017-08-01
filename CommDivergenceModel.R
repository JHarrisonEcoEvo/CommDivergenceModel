
#PROBLEM: the model calculates a prob of death for each integer, since the prob is across a RANGE, this means things die way more often then they should
#I am changing a 10% chance of dying over a decade, into a 10% chance of dying every year for a decade!!!!
# 
# #simulate some data. 
# #Will need to use this to show users what input data should like
# #Remember to put this, or something like it, in the help file example. clean up first 
# age_cat  <-  seq(0,101, by=10)
# for(i in 1:length(age_cat)){
#   if(i<=99){
#   age_cat[i] <- paste(age_cat[i], "-", as.numeric(age_cat[i+1])-1, sep="")
#   }
# }
# probDeath_cat  <-  runif(11, min = 0,max = 0.01)
# lifehistorytable  <-  data.frame(as.character(age_cat), probDeath_cat)
# lifehistorytable[,1] <- as.character(lifehistorytable[,1])
# lifehistorytable[11,1] <- "100-102"
# dat <-lifehistorytable

#example data can be found in ./data

#Annual phlox, time in days: phloxdrummondii.csv
#human from the 2014 Social Security area population https://www.ssa.gov/oact/STATS/table4c6.html
#humanmale.csv
#humanfemale.csv
#--------------------------#
#-----LOADING PACKAGES-----#
#--------------------------#
require(MCMCpack) 
require(ecodist)

#setting seed to ensure consistency
set.seed(14567)
options(scipen = 99)

#Stuff to revisit when all said and done and make sure we like
#1. the divergence function computes divergences on summed abundances for communities, not an average for all pairwise distances between individuals
#this may be ideal, but it is worth thinking about more.

#---------------- -------------#
#-----High level functions-----#
#------------------------------#
#These call lower level functions that are specified below
#"model" runs the model, the highest level function
#generateSame and generateDiff generate initial simulated communities
#replacement - replaces individuals at appropriate time steps
dat= read.csv("./data/humanfemale.csv")
model = function(lifehistorytable, 
                 commSame = TRUE,
                 distancemetric="bray-curtis", 
                 numComm=10, 
                 numIndiv=10, 
                 numMicrobes=100, 
                 microbeAbund=10000, 
                 conc.par=10, 
                 timesteps = 1000){
  #Parameters are:
  #commSame - Boolean to specify if communities should be made up of identical individuals or different individuals
  #distancemetric - is any distance metric accepted by the distance function of ecodist
  #numComm - number of communities
  #numIndiv - number of individuals in each community
  #numMicrobes - number of microbe slots per individual
  #microbeAbund, abundance of microbes, summed across taxa
  #conc.par - Concentration parameter
  #timesteps - how long to run model
  
  plotpoints=seq(0,timesteps, by=timesteps/10) # timesteps to calculate divergence for plotting
  
  #obtain life history data
  agesThetas <- formatLifeTable(lifehistorytable)
  
  #vector to hold output
  divOut  <- NA
  divOutVar <- NA
  #counter
  k <- 0
  
  #generate initial simulated communities, of specified parameters
  if(commSame == TRUE | missing(commSame)){
    community <- generateSame(numIndiv, numMicrobes, numComm, microbeAbund,conc.par)
    print("Generating meta-community using commSame == TRUE")
    print("Individuals within a community will start with identical microbiomes, but communities will differ")
  }
  if(commSame == FALSE){
    community <- generateDiff(numIndiv, numMicrobes, numComm, microbeAbund,conc.par)
    print("Generating meta-community using commSame == FALSE")
    print("Individuals within a community will start with different microbiomes")
  }
  
  #run model for appropriate number of time steps
  repeat{
    k <- k+1
    #replace members of communities at each iteration
    community <- replacement(community, numComm, numMicrobes, numIndiv, microbeAbund) 
    
    #calculate divergence at time steps specified by plotpoints
      if (k %in% plotpoints){ 
        #vector to hold distance indices
        div <- NA
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
##########################################################
#GenerateSame makes a series of communities (k locations/sites/communities) each with n individuals. 
#Individuals within a community start with the same microbial assemblage, but communities differ.

#for debugging
# m = 5
# abund_microbe

generateSame <- function(numIndiv, numMicrobes, numComm, microbeAbund,conc.par){
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

##########################################################
#GenerateDiff makes a series of communities (k locations/sites/communities)each with n individuals. 
#Individuals within a community start with DIFFERENT microbial assemblages.

generateDiff = function(numIndiv, numMicrobes, numComm, microbeAbund,conc.par){
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

##########################################################

#for debugging
  # sandbox = generateSame(10, 100, 10, 10000, 3)
  # community = sandbox

replacement = function(community, numComm, numMicrobes, numIndiv, microbeAbund){  	
  
  #age all individuals in all communities by 1
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
#---------------- -------------#
#-----Low level functions------#
#------------------------------#

############################################################################
#Function to build microbial communities using a Dirichlet process
#inspired by this article (https://en.wikipedia.org/wiki/Dirichlet_process)
############################################################################

dirichletprocess = function(numMicrobes, microbeAbund, conc.par){
  #Make a vector for the new community we are building. 
  D = vector(mode="numeric", length = numMicrobes)
  
  #Add a one to a random microbe to start the community
  D[round(runif(1, min=1, max = numMicrobes))] = 1
  
  #Run the Dirichlet process for time specified by "microbeAbund"
  for(i in 1:microbeAbund){
    newMicrobe = conc.par/(conc.par + i-1)
    oldMicrobe = D/(conc.par + i-1)
    
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

#-----------------#
#-----Testing-----#
#-----------------#

#points at which we calculate the divergence
plotpoints = seq(from = 0, to = 200000, by=1000)
parameterSet=c(1,10,100,300, 500, 1000, 2000)
indiv = c(10,20,50)
microbes = c(50,200,500)
colors = list("red","orange","green","aquamarine","blue","purple","black","burlywood","cadetblue","chartreuse","chocolate","coral","cornflowerblue","cornsilk","cyan","darkblue","darkgoldenrod","darkolivegreen")
colorcount=0
p=0
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
