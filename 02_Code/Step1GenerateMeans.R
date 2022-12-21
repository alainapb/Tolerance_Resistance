#### Project: Feeding rate methods paper
#### Calls data from AJ1, AJ2, and Juv/Adult from the Day/Night experiment


# plotting options http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html#Marginal%20Histogram%20/%20Boxplot

#### This code takes the raw flurometry data and: 
#### Removes whatever data is necessary for question at hand. In this case, spores levels > 0 and males
#### Step 1: Calculates the mean of the controls
#### Step 2: Matches the control mean with the paired per-plate treatment reading
#### Step 3: Calculates the Sarnell and Wilson feeding rate per animal
#### Step 4: Takes the per-individual (technical replicates) average for feeding rates
#### Step 5: Removes any negative feeding rates (as this represents technical errors or non-feeding inividuals)
#### Step 6: Calculates the mean feeding rate per treatment (in this case, genotype)
#### Last updated:  October 2019 JLH

#### I removed A45_0 # 13 from the AJ2 csv file....couldn't recall how to exclude based on multiple variables in R
#### I removed STD_0 # 11 from the AJ csv file....couldn't recall how to exclude based on multiple variables in R

# remove any global variables to start off with a clean slate
rm(list = ls())
### Make sure to clear all: Session > ClearWorkspace


library(tidyverse)
library(MESS)
library(here)

#### Load data


##########################
# 11 Nov 2019: JLH - This code all works with this spreadsheet
#df <- read_csv("AJV2_6Nov2019_JLH.csv",) #### Load the treatment data
##########################

##########################
# 11 Nov 2019: JLH  - Now try with the newly-created combo spreadsheet
# remove animals where age at first clutch  = NA; means these animals died on day 6
#df <- read_csv("all_data_combined_AJV2.csv") 
df <- read_csv(here("01_Data", "AJV2_6Nov2019_JLH.csv"))

names(df)

## rename columnes just to prevent any potential errors
#colnames(df) <- c("genotype",  "spore_level", "plate_treatment", "animal", "sex", "technical_rep", "fluor_reading",  "size_microscope", "size_mm", "surface_area", "final_length", "final_length_mm", "final_surace_area", "total_growth", "inf_status", "susceptible", "infected",  "spores_per_uL", "fecundity", "control_no.", "control_tech_rep" ,"control_fluor", "use_fluor_control", "Plate_Control")

#### Create this summary to get the means of the controls  - these means for each treatment (i.e., Plate_Control) are then used in the calculation below 
summary_df <- df %>%
  group_by(Plate_Control) %>%
  mutate(control_mean = mean(USE_Control_Flour_Reading, na.rm = T)) %>%
  select(Plate_Control, control_mean) %>%
  distinct() %>% as.data.frame()
summary_df

summary_df  <- summary_df[!(summary_df$Plate_Control == ""),] ### for some reasont here was some weird spore level NA - this bi to fcode removes that. 
summary_df


df  <- df[!(df$Genotype == ""),] ### for some reason there was some weird spore level NA - this bi to fcode removes that. 
df

# match control mean with plate treatment
map = setNames(summary_df$control_mean,summary_df$Plate_Control)
df$control_mean <- map[as.character(df$Plate)]

#calculate rate
#also, make sure that these numbers match those in Excel - just as a double check that all is ordered correctly. DONE: 21 Oct 2019 JLH.

# v and t change depending on the round of the experiment (V1 or V2)
v = 10
t = 7
df$fr_sw = log(df$control_mean/df$Flourometry_Reading) * (v/t)



animal_mean1 <- df %>%
  filter(fr_sw > 0) %>%
  filter(Sex=="F") %>%
  group_by(Genotype, Spore_level, Animal, Plate, Plate_Control, control_mean) %>%
  summarise(mean_feeding_rate = mean(fr_sw, na.rm = TRUE))
animal_mean1

a1 <- as.data.frame(animal_mean1, na.rm = T)
head(a1)
a1.data <- write.csv(a1, "MeansForFeedingRateAJV2.csv")  ## export file to excel spreadsheet

#### EXPORT THESE DATA FOR TREE REGRESSION 5 Nov 2019 JLH #### 




