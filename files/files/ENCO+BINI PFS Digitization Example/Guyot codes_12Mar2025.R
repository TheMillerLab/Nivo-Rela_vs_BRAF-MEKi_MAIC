

#============================ Introduction ==================================#

# This program allows users to input a digitalized version of an extracted Kaplan-Meier Curve and associated 
# number of risk information to reconstruct individual patient level data using a method developed by Guyot et al., 2012 

# This program requires 2 file inputs, both in .csv format: 1) digitalized curve file and 2) number of risk file 

# 1) Digitalized curve file:
  # Should be formatted into 2 columns: 
  # "x": the x-axis coordinates of the digitized curve
  # "Curve1": the corresponding y-axis coordinates of the digitized curve


# 2) The number of risk file 
  # Should be formatted into 6 columns:
  # "Time" and "Nrisk": pairs of time and number of patients at risk at the corresponding time
  # "Median Survival", "CI", and "Number of Events": corresponding information regarding the survival curve directly extracted from the published literature
  # "Unit": the unit used for time


#============================ Load packages =================================#

rm(list = ls())
library(ggplot2)
library(grid)
library(dplyr)
library(gridExtra)
library(data.table)
library(MASS)
library(splines)
library(survival)
library(zoo)
library(scales)
library(survminer)



#============================== File inputs  ================================#

# please adjust set working directory path accordingly:
setwd()

# Read in digitalized curve file 
km_file <- here::here("files/MAIC Input Files for Dr.Miller/ENCO+BINI PFS Digitization Example","PFS_Dummer2018_02.28.2023.csv") # Adjust to name of digitized curve file
km_file
file.exists(km_file)
digC <- fread(km_file)%>% mutate(Curve1 = y) %>% dplyr::select(x, Curve1)

nrisk_file <- here::here("files/MAIC Input Files for Dr.Miller/ENCO+BINI PFS Digitization Example","PFS_Dummer2018_NumAtRisk_02.26.2023.csv")  # Adjust to name of number of risk file
nRisk <- fread(nrisk_file)%>%filter(Nrisk != 0)
#========================== Guyot Method functions =============================#

###---- 'FUN_KM_RECON' is a function that reconstructs individual patient-level data (IPD) from digitized survival curves and number-at-risk tables.

FUN_KM_RECON <- function(rawkmfile, rawnriskfile, totev, totp = 0){
  
  # Rename column names
  colnames(rawkmfile) <- c("t","p")  
  
  # Ensure time = 0 exists in the dataset 
  if(sum(rawkmfile$t==0)==0){
    new_row <- data.table(t = 0, p = 1) # create a new row with time = 0 and survival probability = 1
    rawkmfile <- rbind(new_row,rawkmfile)  # Append to the dataset
  } else{
    rawkmfile$p[which(rawkmfile$t==0)] <- 1 # If time = 0 exists, set survival probability to 1
  }
  
  # Reconstruction function adapted from the Guyot paper
  if (!is.null(rawnriskfile)) {
    # Compute lower and upper index bounds for risk table alignment
    nrisk_lower_upper <- Nrisk_low_up(t.risk0 = rawnriskfile$Time,
                                         t.S1    = rawkmfile$t,
                                         n.risk0 = rawnriskfile$Nrisk)
    # Reconstruct IPD using the Guyot method
    ipd_recon <-   reconKMGuyot(tot.events = totev, 
                                t.S        = rawkmfile$t, 
                                S          = rawkmfile$p, 
                                t.risk     = nrisk_lower_upper$t.risk1, 
                                lower      = nrisk_lower_upper$lower1, 
                                upper      = nrisk_lower_upper$upper1, 
                                n.risk     = nrisk_lower_upper$n.risk1, 
                                tol        = .01)
    return(list(surv = rawkmfile, nrisk = rawnriskfile, IPD = ipd_recon$ipd))
  } else {
    
    # If number-at-risk file is NOT provided, use total patients (totp) instead
    ipd_recon <-   reconKMGuyot(tot.events = totev, 
                                t.S        = rawkmfile$t, 
                                S          = rawkmfile$p, 
                                t.risk     = 0, 
                                lower      = 1, 
                                upper      = length(rawkmfile$t), 
                                n.risk     = totp, 
                                tol        = .01)
    return(list(surv = rawkmfile, nrisk = rawnriskfile, IPD = ipd_recon$ipd))
  }
}

# Algorithm to create a raw dataset from digitized Kaplan-Meier curve
# tot.events = total number of events reported. If not reported, then use "NA";
# t.S        = times from the digitized curve;  
# S          = survivor probabilities values extracted from digitized curve;

# If a KM curve reports a time horizon and it's n at risk, then
# Read in published numbers at risk, n.risk, at time, t.risk, lower and upper
# lower, upper, n.risk are all VECTORS

# lower      = index of times in t.S corresponding to each at-risk time
# upper      = index of times in t.S corresponding to each at-risk time; lagged version of lower
# n.risk     = at risk patients
# tol        = tolerance

reconKMGuyot <- function(tot.events, t.S, S, t.risk, lower, upper, n.risk, tol) {
  # Define key variables
  n.int <-length(n.risk) # Number of risk intervals
  n.t <- upper[n.int] # Number of time points in the KM curve
  
  # Initialize vectors for event counts, at-risk population, and censored data
  n.censor <- rep(0, (n.int-1))      # Estimated number of censored observations in each supplied interval
  n.hat    <- rep(n.risk[1]+1, n.t)  # Estimated at-risk counts at all x-axis points at which KM is digitized
  cen      <- rep(0, n.t)            # Estimated censored counts at all x-axis points at which KM is digitized
  d        <- rep(0, n.t)       # Estimated events counts at all x-axis points derived from censorship assumptions
  KM.hat   <- rep(1, n.t)       # Vector that tracks the time of last event (we need this because while the KM curve
  last.i   <- rep(1, n.int)     # guarantees at least one event occurs at the start of each interval, this may not be
  sumdL    <- 0                 # the case for our digitized curve.)
  
  # Iterate over each risk interval to estimate censored and event counts
  if (n.int > 1){
    for (i in 1:(n.int-1)){
      
      # Estimate number of censored observations in the interval
      n.censor[i]<- round(n.risk[i]*S[lower[i+1]]/S[lower[i]] - n.risk[i+1])

      # Adjust censorship count until estimated at-risk population matches the known number-at-risk
      while((n.hat[lower[i+1]]>n.risk[i+1])||((n.hat[lower[i+1]]<n.risk[i+1])&&(n.censor[i]>0))){
        if (n.censor[i]<=0){
          cen[lower[i]:upper[i]]<-0 # If no censorship, set to 0
          n.censor[i]<-0
        }
        
        if (n.censor[i]>0){
          # Distribute censored observations evenly over the interval
          cen.t<-rep(0,n.censor[i])
          for (j in 1:n.censor[i]){
            cen.t[j]<- t.S[lower[i]] + j*(t.S[lower[(i+1)]]-t.S[lower[i]])/(n.censor[i]+1)
          }
          cen[lower[i]:upper[i]]<-hist(cen.t,breaks=t.S[lower[i]:lower[(i+1)]], plot=F)$counts
        }

        # Compute number of events and at-risk population
        n.hat[lower[i]]<-n.risk[i]
        last<-last.i[i]
        
        for (k in lower[i]:upper[i]){
          if (i==1 & k==lower[i]){
            d[k]<-0
            KM.hat[k]<-1
          }
          else {
            d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
            KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          }
          n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
          if (d[k] != 0) last<-k
        }
        n.censor[i]<- n.censor[i]+(n.hat[lower[i+1]]-n.risk[i+1]) 
      }

      if (n.hat[lower[i+1]]<n.risk[i+1]) n.risk[i+1]<-n.hat[lower[i+1]]
      last.i[(i+1)]<-last
    }
  }
  
  # Time interval n.int.
  if (n.int>1){
    #Assume same censor rate as average over previous time intervals.
    n.censor[n.int]<- min(round(sum(n.censor[1:(n.int-1)])*(t.S[upper[n.int]]-
                                                              t.S[lower[n.int]])/(t.S[upper[(n.int-1)]]-t.S[lower[1]])), n.risk[n.int])
  } 
  if (n.int==1){n.censor[n.int]<-0} #make the censorship zero for all times if condition is true.
  if (n.censor[n.int] <= 0){ #estimated number of censored might be negative in the n.int>0 case so we correct. 
    cen[lower[n.int]:(upper[n.int]-1)]<-0
    n.censor[n.int]<-0
  }
  if (n.censor[n.int]>0){
    cen.t<-rep(0,n.censor[n.int])
    for (j in 1:n.censor[n.int]){
      cen.t[j]<- t.S[lower[n.int]] +
        j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
    }
    cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]],
                                             plot=F)$counts
  }
  
  # Find no. events and no. at risk on each interval to agree with K-M estimates read from curves
  n.hat[lower[n.int]]<-n.risk[n.int]
  last<-last.i[n.int]
  for (k in lower[n.int]:upper[n.int]){
    if(KM.hat[last] !=0){
      d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))} else {d[k]<-0}
    KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
    n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
    # No. at risk cannot be negative
    if (n.hat[k+1] < 0) {
      n.hat[k+1]<-0
      cen[k]<-n.hat[k] - d[k]
    }
    if (d[k] != 0) last<-k
  }
  
  # If total no. of events reported, adjust no. censored so that total no. of events agrees.
  if (tot.events != "NA"){
    if (n.int>1){
      sumdL<-sum(d[1:upper[(n.int-1)]])
      # If total no. events already too big, then set events and censoring = 0 on all further time intervals
      if (sumdL >= tot.events){
        d[lower[n.int]:upper[n.int]]<- rep(0,(upper[n.int]-lower[n.int]+1))
        if(upper[n.int]!=lower[n.int]) cen[lower[n.int]:(upper[n.int]-1)]<- rep(0,(upper[n.int]-lower[n.int]))
        else cen[lower[n.int]:(upper[n.int]-1)]<- rep(0,1)
        n.hat[(lower[n.int]+1):(upper[n.int]+1)]<- rep(n.risk[n.int],(upper[n.int]+1-lower[n.int]))
      }
    }
    # Otherwise adjust no. censored to give correct total no. events
    if ((sumdL < tot.events)|| (n.int==1)){
      sumd<-sum(d[1:upper[n.int]])
      while ((sumd > tot.events)||((sumd< tot.events)&&(n.censor[n.int]>0))){
        n.censor[n.int]<- n.censor[n.int] + (sumd - tot.events)
        if (n.censor[n.int]<=0){
          cen[lower[n.int]:(upper[n.int]-1)]<-0
          n.censor[n.int]<-0
        }
        if (n.censor[n.int]>0){
          cen.t<-rep(0,n.censor[n.int])
          for (j in 1:n.censor[n.int]){
            cen.t[j]<- t.S[lower[n.int]] + j*(t.S[upper[n.int]]-t.S[lower[n.int]])/(n.censor[n.int]+1)
          }
          cen[lower[n.int]:(upper[n.int]-1)]<-hist(cen.t,breaks=t.S[lower[n.int]:upper[n.int]], plot=F)$counts
        }
        n.hat[lower[n.int]]<-n.risk[n.int]
        last<-last.i[n.int]
        for (k in lower[n.int]:upper[n.int]){
          d[k]<-round(n.hat[k]*(1-(S[k]/KM.hat[last])))
          KM.hat[k]<-KM.hat[last]*(1-(d[k]/n.hat[k]))
          if (k != upper[n.int]){
            n.hat[k+1]<-n.hat[k]-d[k]-cen[k]
            # No. at risk cannot be negative
            if (n.hat[k+1] < 0) {
              n.hat[k+1]<-0
              cen[k]<-n.hat[k] - d[k]
            }
          }
          if (d[k] != 0) last<-k
        }
        sumd<- sum(d[1:upper[n.int]])
      }
    }
  }
  
  ret1 <- cbind(t.S,S,n.hat[1:n.t],d,cen)

 # Create IPD
  t.IPD     <-rep(t.S[n.t],n.risk[1]) # Time points
  event.IPD <-rep(0,n.risk[1]) # Event indicators (0 = censored, 1 = event)

  # Assign event time and event indicator for each event, as separate row in t.IPD and event.IPD
  k=1
  for (j in 1:n.t){
    if(d[j]!=0){
      t.IPD[k:(k+d[j]-1)]<- rep(t.S[j],d[j])
      event.IPD[k:(k+d[j]-1)]<- rep(1,d[j])
      k<-k+d[j]
    }
  }
  t.risk <- t.risk
  t.S <- t.S
  d <- d
  n.hat<<-n.hat
  
  # Assign censor time and event indicator for each censor, as separate row in t.IPD and event.IPD
  for (j in 1:(n.t-1)){
    if(cen[j]!=0){
      t.IPD[k:(k+cen[j]-1)]<- rep(((t.S[j]+t.S[j+1])/2),cen[j])
      event.IPD[k:(k+cen[j]-1)]<- rep(0,cen[j])
      k<-k+cen[j]
    }
  }

  # Create the final IPD dataset
  IPD <- data.frame(t = t.IPD, ev = event.IPD)
  
  # Adjust for follow-up length difference between digitized curve and IPD
  if(IPD[nrow(IPD),"t"] <  t.S[length(t.S)]){
    IPD[nrow(IPD),"t"] <- t.S[length(t.S)]
  }
  
  return(list(dat = ret1, ipd = IPD))    
  
}

#-------------------------------------------------------------------------------#

###---- 'Nrisk_low_up' is a function to get upper and lower index for the number at risk file
Nrisk_low_up=function(t.risk0,t.S1,n.risk0) {
  
  # Iterate over each time point in the number-at-risk table (except the first and last).
  for (i in 2:(length(t.risk0)-1)) {
    
    # Identify time points in the survival curve (t.S1) that fall within the given at-risk interval.
    temp=which((t.risk0[i]<=t.S1)*(t.S1<t.risk0[i+1])==1)
    # If no valid time points exist within this interval, mark it as missing (-1)
    if (length(temp)==0) {
      t.risk0[i]=-1
      n.risk0[i]=-1
      print(paste0("There is no event/point in between interval ", t.risk0[i], " and ", t.risk0[i+1]))
    }
  }
  
  # Remove any marked (-1) intervals, keeping only valid time points
  t.risk1=t.risk0[which(t.risk0!=-1)]
  n.risk1=n.risk0[which(n.risk0!=-1)]
  
  # Initialize lower and upper index vectors
  lower1=t.risk1
  upper1=t.risk1
  
  # Identify the lower and upper index for each at-risk time
  for (i in 1:(length(t.risk1)-1)) {
    lower1[i]=min(which((t.risk1[i]<=t.S1)*(t.S1<t.risk1[i+1])==1))
    upper1[i]=max(which((t.risk1[i]<=t.S1)*(t.S1<t.risk1[i+1])==1))
  }
  
  lower1[length(lower1)]=min(which(t.risk1[length(t.risk1)]<=t.S1))
  upper1[length(upper1)]=max(which(t.risk1[length(t.risk1)]<=t.S1))
  
  result=list("lower1"=lower1,"upper1"=upper1,"t.risk1"=t.risk1,"n.risk1"=n.risk1)
  return(result)
}


#========================== Reconstruct IPD =============================#

###--- The 'process_KM_data' function uses the above Guyot functions to reconstruct IPD, 

process_KM_data <- function(digC, nRisk) {
  
  # Reconstruct IPD
  IPD <- FUN_KM_RECON(rawkmfile = digC, rawnriskfile = nRisk, totev = nRisk$'Number of Events'[1])

  # Return the reconstructed IPD 
  list(
    IPD = IPD$IPD
  )
}


#============================== Generate outputs =============================#

###--- The 'process_KM_data' function inputs (digC = digitized curve file, nRisk = number at risk file) 
# were already defined at the top of this program. Run this program to generate the IPD.

result <- process_KM_data(digC, nRisk)

#============================== Save Outputs ==============================#

# Save the reconstructed IPD as a CSV file
write.csv(result$IPD, "Reconstructed_IPD.csv", row.names=FALSE)

#Check curve
library(survival)
library(survminer)

dat_curve <- result$IPD
# Convert data to a survival object
surv_obj <- Surv(time = dat_curve$t, event = dat_curve$ev)

# Fit Kaplan-Meier survival curve
km_fit <- survfit(surv_obj ~ 1, data = dat_curve)
# Extract the median survival time
median_time <- summary(km_fit)$table["median"]
print(median_time)

ggsurvplot(
  km_fit,
  data = dat_curve,
  risk.table = TRUE,         # Show number at risk
  risk.table.y.text = TRUE,  # Show risk table labels
  risk.table.height = 0.2,   # Adjust risk table height
  break.time.by = 4,         # Set x-axis ticks (every 2 months)
  xlab = "Time (Months)", 
  ylab = "Survival Probability",
  surv.median.line = "hv",   # Show median survival line
  ggtheme = theme_minimal()
)



###--- Program ends here ---###


