rm(list = ls())


#### Phenoloxidase ####

library(tidyverse)
library(cowplot)  
library(here)
library(flextable)
library(ggplot2)
library(viridis)


#### LOAD PO DATA ####
##adding AJ2 active phenoloxidase data
## JLH My changes here are b/c it looked like these calculations weren't averaged over the technical replicates? Were they and I missed something? The code below should calculate the mean for each biological replicate, which can then be used to calculate delta PO or whatever.
## I have not looked yet to see whether or not this changed any of the results. 

all_po_data <- read.csv(here("01_Data", "AJ2_PO_for_R.csv"))
names(all_po_data)
head(all_po_data)

#all this cleaning can be done in a dplyr pipeline. Much nicer to see and to #read!
##cleaning PO data 

##Removing time point 0

all_po_data <- all_po_data %>% 
  filter(Time != 0) %>% 
  filter(!is.na(Abs)) %>%  #filter out NAs
  mutate(Abs = as.numeric(Abs), 
         Prop_inf = as.numeric(Prop_inf))


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



#### CALCULATE DELTA PO ####

# # subset the data to only include time "4.5"
final <- subset(all_po_data, all_po_data$Time == "4.5")
final$Abs
# # 
# # ### Now look at the change in abs which equivalent to the amount of active PO produced ####
# # #ADD CHANGE IN ABS (ACTIVE PHENOLOXIDASE)
initial <- subset(all_po_data, all_po_data$Time == 1)
initial$Abs
# # 
# 
# # if you run this, remember to convert factors to numeric (both final$abs and
# # intiial$abs are factors, you have to convert them to numeric), like
final$Abs <- as.numeric(final$Abs)
initial$Abs <- as.numeric(initial$Abs)

# #calculating the change in absorbance
delta_po <- final$Abs - initial$Abs
# 
# # multiplying by 1000 to get "activity units" (Mucklow & Ebert, 2003)
delta_po <- delta_po*1000
# # adding the delta_po values back to dataframe
initial$delta_po <- delta_po
# 

#average across replicates and divide by the number of individuals in the sample
po <- initial %>%
  dplyr::group_by(Genotype, Spores, Sample, Prop_inf, No_in_sample) %>%
  summarise(delta_po = mean(delta_po))

po$cor_delta_po <- po$delta_po/po$No_in_sample

po <- subset(po, Genotype != "Control")

names(po)[2]<-"Spore_level"



names(po)


avg_po1 <- po %>% 
  filter(Genotype != "Control") %>% 
  arrange(Genotype, Spore_level, Sample) %>% 
  group_by(Genotype, Spore_level) %>%
  mutate(mean_po = mean(cor_delta_po, na.rm = T)) %>% 
  mutate(se_deltapo = sd(cor_delta_po, na.rm = TRUE)/sqrt(n())) %>%
  dplyr::select(Genotype, Spore_level, mean_po, se_deltapo) %>% 
  distinct() %>% as.data.frame()
avg_po1 

#### PLOT PO and SPORE DOSE (can't plot it with pathogen load b/c it's averaged across) ####




singlefigtheme <-theme_classic()+theme(axis.title.x = element_text(size = 18),
                                       axis.title.y = element_text(size = 18),
                                       axis.text.x = element_text(size = 14),
                                       axis.text.y = element_text(size = 14))




genolabel <- c("Genotype 1", "Genotype 2", "Genotype 3", "Genotype 4") 
names(genolabel) <- c("A45", "Midland 252", "Midland 281", "Standard")



#### Fix = mean and SE are incorrect ####
plota <- avg_po1  %>% 
  ggplot(aes(x = as.factor(Spore_level), 
             y = as.numeric(mean_po), 
             #size = Prop_inf,
             #shape = as.factor(Spore_level),
             color = as.factor(Genotype)))+
  geom_point(alpha = 0.3, position = position_dodge(0.2), size = 5)+
  geom_errorbar(aes(ymin = mean_po - se_deltapo,
                    ymax = mean_po + se_deltapo),
                width = 0, 
                position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(0.2), size = 5) +
 singlefigtheme +
  theme(legend.position = c(0.10, 0.90),
        panel.grid.major = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.75), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))+
  labs(y = expression("Phenoloxidase activity"))+
  labs(x = expression("Pathogen dose (no. mL"^-1*")"))+
  labs(
    colour = "Host genotype",
    shape = ""
  )
plota 






#### STATS #### 
avg_po2 <- po %>% 
  filter(Genotype != "Control") %>% 
  arrange(Genotype, Spore_level, Sample) %>% 
  group_by(Genotype, Spore_level) %>%
  mutate(mean_po = mean(cor_delta_po, na.rm = T)) %>% 
  dplyr::select(Genotype, Spore_level,Prop_inf, mean_po) %>% 
  distinct() %>% as.data.frame()
avg_po2 
 

lm_po_1 <- lm(mean_po ~ Genotype + as.factor(Spore_level) + Prop_inf, data = avg_po2)
lm_po_2 <- lm(mean_po ~ Genotype + as.factor(Spore_level), data = avg_po2)
lm_po_3 <- glm(mean_po ~ Genotype 
               + as.factor(Spore_level), 
               family = gaussian(link = "log"),
               data = avg_po2)

bbmle::AICctab(lm_po_1, lm_po_2, lm_po_3)


plot(lm_po_3)
stats1 <- car::Anova(lm_po_3, type = 'III')
Term <- rownames(stats1)
stats1 <- cbind(Term, stats1)
autofit(colformat_num(x = flextable(stats1), j = "Pr(>Chisq)", digits = 4))







#### Fold Diferences #####

avg_po3 <- po %>% 
  filter(Genotype != "Control") %>% 
  arrange(Genotype, Spore_level) %>% 
  group_by(Genotype, Spore_level) %>%
  mutate(mean_po = mean(cor_delta_po, na.rm = T)) %>% 
  dplyr::select(Genotype, Spore_level, mean_po) %>% 
  distinct() %>% as.data.frame()
avg_po3
