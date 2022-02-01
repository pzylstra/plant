library(plant)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(frame)
library(assertthat)
library(extraDistr)

# Import species traits and settings
interval <- 5
tr <- read.csv("Sp_Inputs.csv")  %>%
  mutate(Species = as.character(Species))

# STEP 1: Grow the forest
system.time(res <- grow_forest(dat = tr, B_lf1 = 1)) 

# STEP 2: Collect fire behaviour parameters
system.time(fireDat <- frameDynTab(dat=res, tr, upper = 100, interval = interval, 
                                   sample = 0.5, transects = 10, propDead = 0, leafForm = "Flat", 
                                   lwRat = 3, leafA = 0.002547, ram = 5, ignitionTemp = 260, moist = 1, 
                                   G.C_rat = 3, C.C_rat = 0.1, deltaL = 0.46, lat = -35, map = 1000, mat = 20))

# STEP 3: Model flammability dynamics
system.time(dyn <- firePlant(dat=fireDat, db.path = "out.plant.db", reps = 10,
                             slope = 10, slopeSD = 5, slopeRange = 15, 
                             temp = 30, tempSD = 5, tempRange = 3,
                             DFMC = 0.1, DFMCSD = 0.01, DFMCRange = 2, 
                             wind = 20, windSD = 5, windRange = 20,
                             moistureMultiplier = 1, moistureSD = 0.01, moistureRange = 1.5,
                             fLine = 100, leafVar = 0.1, updateProgress = TRUE) %>%
              mutate(yr = step*interval,
                     dCan = sc4/100))


######################################################
# Plot forest
result <-  res%>% 
  tidy_patch() %>% 
  plant:::expand_state() 
tab <- result$species%>%
  drop_na() %>%
  left_join(tr, by = c("species" = "Species")) %>%
  mutate(species = name)

# Plot size distribution to pdf
pdf("forest.pdf", height = 8, width = 13)
tab %>% plot_size_distribution_patch()
dev.off()



# Summarise mortality
mort <- dyn %>%
  group_by(yr) %>%
  summarise_if(is.numeric,mean) %>%
  select(yr, dCan)


# Plot behaviour
ros <- ggplot(data = dyn, 
              aes(x = yr, y = ros_kph, group = yr)) +
  ggtitle("Rate of spread") +
  geom_boxplot(outlier.size = 1)+
  labs(y = "ROS (km/h)" ) +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))

fh <- ggplot(data = dyn, 
             aes(x = yr, y = fh, group = yr)) +
  ggtitle("Flame height") +
  geom_boxplot(outlier.size = 1)+
  labs(y = "Height (m)" ) +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))

fl <- ggplot(data = dyn, 
             aes(x = yr, y = fl, group = yr)) +
  ggtitle("Flame length") +
  geom_boxplot(outlier.size = 1)+
  labs(y = "Length (m)" ) +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))

sc <- ggplot(data = dyn, 
             aes(x = yr, y = Height, group = yr)) +
  ggtitle("Scorch height") +
  geom_boxplot(outlier.size = 1)+
  labs(y = "Height (m)" ) +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))

sev <- ggplot(data = dyn, 
              aes(x = yr, y = severity, group = yr)) +
  ggtitle("Fire severity") +
  geom_boxplot(outlier.size = 1)+
  ylim(0,5) +
  labs(y = "Severity" ) +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))

can <- ggplot(data = dyn, 
              aes(x = yr, y = dCan, group = yr)) +
  ggtitle("Canopy death") +
  geom_boxplot(outlier.size = 1)+
  labs(y = "Likelihood of death" ) +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))

mcan <- ggplot(data = mort, 
               aes(x = yr, y = dCan)) +
  ggtitle("Canopy death") +
  geom_line()+
  ylim(0,1) +
  labs(y = "Likelihood of death" ) +
  theme_bw() +
  theme(axis.text.x  = element_text(vjust=1.5, size=14, colour = "black"),
        axis.text.y  = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        plot.title = element_text(vjust=1.5, face="bold", size=14, colour = "black"))


windows(15,10)
ros+fh+fl+sc+sev+mcan+
  plot_layout(ncol = 3)

