
rm(list = ls())


library(subplex)            # for model optimization, faster than optim
library(tidyverse)
library(deSolve)
library(magrittr)
library(bbmle)
library(pomp)
require(pomp)
library(coda)
library(doParallel)         # for parallel computing
library(pbmcapply)          # for parallel computing
library(foreach)            # faster than a for loop


#### Used for: Max likelihood analysis to fit mechanistic transmission model to data

# The first part of code combines data from various spreadsheets (and we will strive to make a single master spreadsheet for all future experiments!)


##Feeding Rate data
data <- read_csv("data/AJV2_6Nov2019_JLH.csv")
names(data)

#Life table data
life <- read_csv("data/LifeTable_AJV2_3Dec2019.csv")
names(life)

#check data
head(data)


head(life)

#merge data together
left_join(data, life, by=c("Genotype", "Spore_level", "Animal", "Infected_Uninfected")) %>%
  mutate(Individual_ID = paste(Genotype, 
                               Spore_level, 
                               Animal, 
                               sep = "_")) %>%
  as.data.frame() -> all.data.joined

names(all.data.joined)

##removing unnecessary columns
drops <- c("Sex.y", "Length_Day_7", "Length_Day_20", "DOB", "Total.Offspring")
all.data.joined <-all.data.joined[ , !(names(all.data.joined) %in% drops)]

# colnames(all.data.joined)<-c("Genotype", "Spore_level", "Plate", "Animal", "Sex", "Technical_Replicate", "Fluorometry_Reading", "Length_mm", "Size_mm2", "Final_length_mm", "Final_Size_mm2", "Growth_mm2", "Infected_Uninfected", "Susceptible", "Infected", "Spores_ul", "Total_Fecundity", "Control_Number", "Control_Technical_Replicate", "Control_Fluor_Reading", "USE_Control_Fluor_Reading", "Plate_Control_Numeric", "Time.till.death.days", "Censured", "Growth_mm", "Age_1st_clutch", "Indivdual_ID")

##Remove males from dataframe
all.female <- subset(all.data.joined, all.data.joined$Sex == "F")








#### LOAD PO DATA ####
##adding AJ2 active phenoloxidase data
all_po_data <- read_csv("data/AJ2_PO_for_R.csv")
names(all_po_data)
head(all_po_data)


##cleaning PO data 

##Removing time point 0
subset(all_po_data, all_po_data$Time !="0") -> po_data
##Removing NAs
subset(po_data, po_data$Abs != "na") -> po_data
##Making all the data numeric
as.character(po_data$Abs) -> po_data$Abs
as.numeric(po_data$Abs) -> po_data$Abs
as.character(po_data$Prop_inf) -> po_data$Prop_inf
as.numeric(po_data$Prop_inf) -> po_data$Prop_inf

names(all_po_data)
head(all_po_data)

length(all_po_data)
all_po_data
#### po_data includes:Genotypes (control=PBS), Spores (exposed to), ####
## Sample (sample number), Replicate (2 reps from each sample),
## time (in hours), Abs (absorbance value), and 
## Prop_inf (proportion of infected individuals in the sample)


# GOALS:
# (1) take the average for each technical replicate
# (2) calculate delta PO: take the mean_po at time point 4.5 and subtract it from mean_po at time point 1.0
# (3) add that data frame to the dataframe above (all.data)
clone.names <- as.character(unique(all_po_data$Genotype))
clone.names

# here the mean we want to calculate doesn't work because the na in the original dataframe
# are not real R's NA, but they are just character strings na. Two possible solutions:
# you delete the rows where those na are, or you change those na into proper NA. 
# I chose the second solution, at line 154/155 (transformed into numeric first, because it was a factor)

# avg_po1 <- all_po_data %>% 
#   mutate(Abs = ifelse(Abs == "na", NA, Abs)) %>% 
#   mutate(Abs = as.numeric(Abs)) %>% 
#   arrange(Genotype, Spores, Time, Replicate, Sample) %>% 
#   group_by(Sample) %>%
#   mutate(mean_po = mean(Abs, na.rm = T)) %>% 
#   select(Genotype, Spores, Time, Sample, mean_po) %>% 
#   distinct() %>% as.data.frame()
# avg_po1 

#### CALCULATE GENOTYPE PO ####


geno_po <- all_po_data %>% 
  mutate(Abs = ifelse(Abs == "na", NA, Abs)) %>% 
  mutate(Abs = as.numeric(Abs)) %>% 
  arrange(Genotype, Spores) %>% 
  group_by(Genotype, Spores) %>%
  mutate(mean_po = mean(Abs, na.rm = T)) %>% 
  select(Genotype, Spores, mean_po) %>% 
   distinct() %>% as.data.frame()
geno_po



#### Means of the controls  - these means for each treatment (i.e., plate_control) are then used in the calculation below ####
 

summary_df <- all.female %>%
  group_by(Plate_Control) %>%
  mutate(Control_mean = mean(USE_Control_Flour_Reading, na.rm = T)) %>%
  select(Plate_Control, Control_mean) %>%
  distinct() %>% as.data.frame()



#match control mean with plate treatment
map = setNames(summary_df$Control_mean,summary_df$Plate_Control)
all.female$Control_mean <- map[as.character(all.female$Plate)]
all.female$Control_mean

v <- 10 ## volume, mL
t <- 7  ## time, hours

k <- all.female$Control_mean/all.female$Flourometry_Reading  ##difference in feeding compared to control

##calculate feeding rate (Sarnelle and Wilson)
all.female$fr_sw <- log(k) * (v/t)

#calculate mean fluor reading per animal (keep rate > 0)
data <- all.female %>%
  #keep fr_sw values > 0
  filter(fr_sw > 0) %>%
  #group based on these variables
  group_by(Genotype, Spore_level, Animal, Plate, Control_mean) %>%
  dplyr::mutate(ID = Genotype) %>% 
  #means based on groups above
  mutate(mm = mean(Length_mm, na.rm = TRUE)) %>% 
  mutate(A1 = mean(Flourometry_Reading, na.rm = TRUE)) %>% 
  mutate(A0 = mean(Control_mean, na.rm = TRUE)) %>% 
  mutate(Susceptible = mean(Susceptible, na.rm = TRUE)) %>% 
  mutate(Infected = mean(Infected, na.rm = TRUE)) %>% 
  mutate(P0 = mean(Spores_ul, na.rm = TRUE)) %>% 
  select(ID, Genotype, Spore_level, Animal, Plate, Control_mean, mm, A1, A0, Susceptible, Infected, P0) %>%
  distinct() %>% as.data.frame()
data


## 'data' is a dataframe with columns:
## 'P0' parasite/pathogen spore dose
## 'Animal' biological replicate number
## 'ID' a numeric identifier for genotype treatment (added below)
## 'A0' initial algae concentration (control values - these are plate-specific)
## 'mm' gives the size of the animal (convert to surface area: L^2)
## 'A1' final algae concentration (final algal concentration - fluorometry reading)
## 'Susceptible' number of susceptible hosts (0  = infected, 1 = susceptible: for this case, all treatments had 1 animal/vial)
## 'Infected' infection outcome (0 = uninfected, 1 = infected)

sum(is.na(data))
sum(is.na(data$Infected))
sum(is.na(data$mm))
colSums(is.na(data))
sum(is.na(data$A1))
sum(is.na(data$A0))


# remove missing data from the S/I columns...this was causing a problem below in the NLL function.
data <- data[!is.na(data$Infected), ]
data

# NEXT, remove missing data from mm
data <- data[!is.na(data$mm), ]
data


#### All of this code needs to loop through multiple genotypes, and needs to match each individual's feeding rate with that animal's length, 
#### infection status (at the end of the experiment), spore exposure and plate-specific control (A0) value
#### Bolker's book has this cool trick which will help with some of this.

## Which model to use? (currently, 1-4)
model <- 'mod1'


## 'params' is a named vector of parameters f, sd, and u
## For the purposes of the fitting, these parameters should be
## able to take any values from -Inf to Inf, so they must be
## transformed using the log (for f and sd) or logit (for u)


## It is useful to sweep across a large number of parameter
## combinations before attempting any ML fitting. The code below
## sets up 'nseq' initial parameter estimates spread across
## a hyperdimensional volume bounded by lower and upper estimates
## of each parameter.

## With the current model formulation, the
## genotype-specific feeding parameters (f1-f1) need to be transformed
## in order to vary between -Inf and Inf
sobol_design(lower = c(u1 = log(0.0001),         ## genotype-specific per-spore susceptibility
                      u2 = log(0.0001),
                      u3 = log(0.0001),
                      u4 = log(0.0001),
                      f1 = log(0.001),          ## genotype feeding rate
                      f2 = log(0.001),
                      f3 = log(0.001),
                      f4 = log(0.001),
                      a  = -0.01,               ## per-spore affect on feeding rate/strength of anorexia
                      t  = 1.0,
                      sd = log(0.1)),           ## observation error in feeding
            upper = c(u1 = log(0.0009),
                      u2 = log(0.0009),
                      u3 = log(0.0009),
                      u4 = log(0.0009),
                      f1 = log(30000),
                      f2 = log(30000),
                      f3 = log(30000),
                      f4 = log(30000),
                      a = 0.01,
                      t  = 1.0,
                      sd = log(0.1)),
            nseq = 500) -> pars


## compute the negative log-likelihood
## 'params' is a named vector of parameter values
## 'data' is the full dataset
## 'mod' is a character string giving the name of the feeding model to fit

nll2 <- function(params, data, mod) {
  T <- 7                                                         ## length of feeding/exposure in hrs.
  v <- 10                                                        ## volume 
  
  nll <- 0                                                       ## for storage of the final NLL
  for (id in unique(as.numeric(data$ID))) {                                  ## For each GENOTYPE value, compute the NLL of its set of parameter values and add this to the existing total
    
    fhat <- c(exp(params["f1"]),                                 ## transform age/day-specific feeding rate back to the                                                                        natural scale and choose the fhat value (per genotype)
              exp(params["f2"]),
              exp(params["f3"]),
              exp(params["f4"]))[id] %>% unname
    a  <- params["a"] %>% unname                                 ## per-spore effect on feeding can already vary between -Inf (a<0 implies a decrease in feeding rate with spores) and Inf (a>0 implies an increase in feeding rate with spores)
    sd <- exp(params["sd"]) %>% unname
    
    u  <- c(exp(params["u1"]),                                   ## transform genotype-specific susceptibility back to the natural scale and choose the u value for this genotype
            exp(params["u2"]),
            exp(params["u3"]),
            exp(params["u4"]))[id] %>% unname
    
    data2 <- subset(data, ID == id)                              ## extract the data for this genotype
    
    ## compute the size- and spore-dependent feeding rate
    ## You can specify as many different feeding models as you want
    ## I'll just use the two hypothetical models above (one with an
    ## exponential change in feeding with spores, and one with a linear
    ## change) as an example:
    f_pred <- switch(mod,
                     mod1 = fhat*(data2$mm),
                     mod2 = fhat*(data2$mm)*exp(-a*data2$P0*data2$mm),
                     mod3 = fhat*(data2$mm)*exp(-1* a*data2$P0),                         
                     mod4 = ifelse(data2$Time_treatment == 1,fhat*(data2$mm)*exp(-1.5*a*data2$P0), 
                                   fhat*(data2$mm)*exp(-1.5*a*data2$P0))
    )
    
## if any of the feeding rates are negative, this can generate a negative prob of infection,  breaking dbinom and generating errors
## these are bad parameter sets, so simply assign them very high -logLik values. Note: you shouldn't assign a log-likelihood of 0, as although 
## the likelihood must be a positive number, the log-likelihood can be negative if the likelihood is small, which it often will be when the fitting 
## algorithm gets started. Thus, assigning a value of 0 might actually be saying that this parameter set isn't terrible whereas assigning it a value of 
## Inf says that this parameter set is the worst.
    
    if (any (f_pred < 0))
      nll <- nll + Inf
    else {                                                                          ## Otherwise, estimate the amount of algae that is remaining
      exp.A1 <- log(data2$A0) - f_pred*(T/v)                                        ## Compute the expected amount of food remaining, assuming that food in the vial changes according to the ODE dA/dt = -f_pred*A, which has the solution A(t) = A(0)*exp(-f_pred*t). Then, after T time units, A(T) = A(0)*exp(-f_pred*T). Taking logarithms, you have log(A(T)) = log(A(0)) - f_pred*T
      lambda <- 1 - exp(-u*f_pred*data2$P0)                                         ## compute the expected probability of infection, given feeding and per-spore susceptibility
      
      
## Compute the negative log-liklihood of observing log(data2$A1) algae remaining, under the expectation exp.A1, and the negative log-likelihood of observing the binary infection outcome data2$Infected, 
## under the expected probability of infection lambda; sum these two NLLs, and add that sum to the total NLL across all genotypes
      nll <- -sum(dnorm(log(data2$A1), 
                        mean = exp.A1,
                        sd = sd,
                        log = TRUE) %>% sum,
                  dbinom(data2$Infected,
                         size = 1,
                         prob = lambda,
                         log = TRUE) %>% sum) + nll
    }
  }
  return(nll)
}



## To do this we will do the following
## the function apply(pars,1,function(p)...) says:
## 1. Take the information in the dataframe 'pars' and
## 2. For each row of that dataframe (specify performing an option on the row by setting the second value passed to apply equal to 1)
## 3. Pass the values in that row to the function 'nll2'. Thus 'p' in function(p) is the current row in pars. You have to do the unlist(p) because pars is a dataframe, whereas nll2 is expected 'pars' to be a named numeric vector. 
## (You can see this by running class(pars[1,]) versus class(unlist(pars[1,]))

## Conceptually, what apply is doing is doing the following for loop:
## ll <- vector(length = nrow(pars))
## for (i in 1:nrow(pars)) {
##      p <- unlist(pars[i,])
##      ll[i] <- nll2(p, data, mod = 'mod1')
## }


apply(pars,
      1,
      function(p) nll2(unlist(p), data, mod = 'mod4')) -> ll 



####################
## Just as a note: when  taking guesses at parameters for fitting a model, 
## I actually want to see a lot of guesses that result in Inf for the likelihood. 
## This is because, if I have no idea what the parameter values actually are, I want to 
## be sure that I bracket the full range of the plausible. If most of the estimates 
## "work" and return a value for the loglik (and in particular if those logliks don't vary a lot), 
## then you should broaden the parameter ranges for the initial guesses. This, of course,
## also means that you will need to make more guesses!
####################



## remove all parameter sets whose -logLik was Inf
pars2 <- pars[which(ll!=Inf), ]



## for each of the parameter guesses that produced a finite NLL,
## start a Nelder-Mead optimizer to find the parameter values that minimize
## the NLL. Again, you need to use unlist to turn the parameter dataframe into
## a numeric vector.


fits  <- vector(mode = 'list', length = nrow(pars2))                                          ## create a list to store the output from each Nelder-Mead optimization




# using subplex instead of optim, it is much much faster  
# older code - without parallel it was running insanely slowly!
# for(i in 1:nrow(pars2)) {
#     print(i)
#     subplex(par = unlist(pars2[i,]), 
#           fn = nll2,
#           data = data,
#           mod ='mod1',
#           hessian = FALSE,
#           control = list(maxit = 10) # set to 10 just for testing code
#           ) -> fits[[i]]
# }


no_cores <- detectCores() - 1  # Detect cores number and use detectCores-1 because we don't need to use all cores. 
cl       <- makeCluster(no_cores)
registerDoParallel(cl)

getDoParWorkers()             # sanity check to make sure we are actually running in parallel

# use foreach and dopar as loop using cores = Total core -1
fits <- foreach(i = 1:nrow(pars2)) %dopar%{
  library(subplex) # has to be included here, not sure why!
  library(doParallel)         # for parallel computing
  library(foreach)            # faster than a for loop
  print(i)
  subplex(par = unlist(pars2[i,]), 
          fn = nll2,
          data = data,
          mod ='mod4',
          hessian = FALSE,
          control = list(maxit = 10) # set to 10 just for testing code
  ) -> fits[[i]]
}

stopCluster(cl)



## before proceeding, there's something else we need to do here, which is check for convergence. For each output from optim stored in fits, we need to look at the convergence value. 
## If convergence is anything other than 0, it means that the optimizer had not yet converged. For example, if fits[[1]]$convergence = 1, it means that the optimizer hit the max. 
## number of function iterations and gave up without converging. When maxit was set to 200 previously, I was getting this convergence value for every parameter set.
## When I increased maxit to 5000 (which obviously really slows everything down), I started to see convergence. Of course, this greatly, greatly increased how long the 
## for loop takes to run!


fits[[1]]$convergence


## The following line tells which of the fits actually converged. Those are the only ones to pay attention to.
good.ests <- which(unlist(lapply(fits, function(f) f$convergence)) == 0)
good.ests
fits2     <- fits[good.ests]
fits2

## Extract the parameter estimates from fits2
## We'll do this using the function lapply, since fits2 is a list
## What lapply(fits2, function(f) f$pars) is doing is essentially
## the following for loop:
## f <- vector(mode = 'list', length = length(fits))
## for (i in 1:length(fits2))
##      f[[i]] <- fits2[[i]]$pars

## What I want to do with those parameter estimates is put them into a data.frame where each row is the set of optimized parameter values
## I am doing this by taking the list f and unlisting it, which turns it from a list into a vector with a length equal to (the number of estimated parameters) * (the length of fits). 
## I want to just put this into a matrix with a number of columns that is equal to the number of parameters that I am estimating. The 'matrix' function will essentially take this 
## really long vector and "wrap" it into the matrix by taking the first (number of parameters) values and putting them into the first row, then taking values (number of parameters)+1 to 2*(number of parameters) 
## and putting them into the second row, and so on. I then turn the matrix into a dataframe.## I could do this with the extra line
## parsets <- as.data.frame(matrix(unlist(f), ncol = ncol(pars), byrow = T))
## I've just put all of that into the following more concise form

lapply(fits2, function(f) f$par) %>% 
  unlist %>%
  matrix(., ncol = ncol(pars), byrow = T) %>% 
  as.data.frame -> parsets
colnames(parsets) <- names(fits[[1]]$par)


library(plyr)
## transform the parameter estimates back to the natural scale
mutate(parsets,
       u1 = exp(u1)/(1 + exp(u1)),
       u2 = exp(u2)/(1 + exp(u2)),
       u3 = exp(u3)/(1 + exp(u3)),
       u4 = exp(u4)/(1 + exp(u4)),
       f1 = exp(f1),
       f2 = exp(f2),
       f3 = exp(f3),
       f4 = exp(f4),
       sd = exp(sd)) -> parsets



## One last thing to do: pull out the NLL for each parameter set. Each initial guess converges to a different final spot, and these may differ a lot in NLL. We're only interested in parameter sets with high likelihood.
parsets$nll <- lapply(fits2, function(f) f$value) %>% unlist
range(parsets$nll) ## you can see that there is a *huge* range here

bestpars <- unlist(parsets[which.min(parsets$nll), ])                                            ## pull out just the best parameter set




####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####### ********************************** PLOTTING ****************************************########
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


## Using this best-fitting parameter set, let's compare the predicted and observed algae levels.


pred.A1 <- vector(length = nrow(data))                                                        ## storage for the predicted algae measurements

for (i in 1:nrow(data)) {
  g     <- data$ID[i]                                                                         ## select the current genotype
  fhat  <- unname(bestpars[paste0("f1", g)])                                                  ## pull the genotype-specific parameter estimate for fhat
  a     <- unname(bestpars["a"])                                                              ## pull the estimate of a
  l     <- data$mm[i]                                                                         ## pull the observed length
  A0    <- data$A0[i]                                                                         ## pull the observed initial algae measurement
  A1    <- data$A1[i]                                                                         ## pull the observed final algae measurement
  P0    <- data$P0[i]                                                                         ## pull the observed spore load
  ## compute the predicted algae
  pred.A1[i] <- log(A0) - (fhat*(l^2)*exp(a*P0))*(7)
}




#### NEW PLOT ####
data %>%
  group_by(ID, P0, Time_treatment) %>%
  ggplot()+
  geom_point(aes(x = data$mm^2, y = log(data$A1), color = as.factor(P0)))+
  geom_line(aes(x = data$mm^2,  y = pred.A1, color = as.factor(P0)))+
  labs (x = "Host surface area", y = "Intake rate") +
  theme(axis.text = element_text(size  = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank()) +
  theme_classic()+
  facet_wrap(~ ID + Time_treatment)







## Plot observed against predicted
plot(c(prev.g1.p150,
       prev.g1.p300,
       prev.g2.p150,
       prev.g2.p300,
       prev.g3.p150,
       prev.g3.p300,
       prev.g4.p150,
       prev.g4.p300),
     c(lambda.g1.p150,
       lambda.g1.p300,
       lambda.g2.p150,
       lambda.g2.p300,
       lambda.g3.p150,
       lambda.g3.p300,
       lambda.g4.p150,
       lambda.g4.p300),
     pch = 21,
     cex = 2,
     bg = c(1,1,2,2,3,3,4,4),
     xlab = "Obseved prevalence",
     ylab = "predicted prevalence")
abline(0, 1)
#### PLOT PROBABILITY OF INFECTION ACROSS SIZES: #### 
# plot a line showing how the probability of infection
# changes with size, assuming the best-fitting parameters
pred.prev <- expand.grid(SA = seq(0.5,2.5,0.1), ID = 1:6, spores = c(0, 300))
mutate(pred.prev, 
       ## predicted feeding rate: fhat*SA*exp(a*spores)
       f = bestpars[paste0("f1", ID)]*SA*exp(bestpars["a"]*spores),
       ## predicted probability of infection: 1 - exp(-u*f*spores)
       infProb = 1-exp(-bestpars[paste0("u",ID)]*f*spores)
) -> pred.prev




data %>%
  group_by(ID, P0) %>%
  ggplot()+
  geom_point(aes(x = data$mm, y = data$Infected, color = as.factor(P0)))+
  geom_line(data = pred.prev, mapping = aes(x = SA, y = infProb, color = as.factor(spores))) + 
  #ggPredict(pred.prev, data = data, se = TRUE, family = binomial)+
  #geom_line_interactive(seq(0.5, 3, by = 0.5), eval(pred.prev), lwd = 2, col = grey(0.85))
  labs (x = "Surface area", y = "Frequency of Infection") +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        legend.title = element_blank()) +
  theme_classic()+
  facet_wrap(~ ID)





