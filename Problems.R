library(plant)
library(dplyr)
library(tidyr)
library(frame)
library(assertthat)
library(purrr)

dat <- grow_forest(dat = tr, B_lf1 = 1)
tr <- read.csv("Sp_Inputs.csv")  %>%
  mutate(Species = as.character(Species))
lat <- -35
map <- 1000
mat <- 20
sample <- 0.5
transects <- 10
age <- 30
rec <- 1

f <- function(age) {
  # 1. Modify plant outputs by adding in crown parameters and stratifying them
  comm <- stratify_community(dat, tr, age, lat, map, mat, sample, transects)
  # 2. Build the 'structure' table needed for modelling in frame
  Structure <- frame::buildStructureP(comm, age, rec)
}

# The function gives correct outputs when run for the age 30,as for other ages
age<-c(30)
purrr::map(age, f) %>%
  dplyr::bind_rows()

# If you run it after age 25, however, the midstorey and elevated strata get confused and there is an NA output
age<-c(25, 30)
purrr::map(age, f) %>%
  dplyr::bind_rows()

# Somehow it's carrying the previous age over, because if you repeat the age 30, it works fine the 2nd time
age<-c(25, 30, 30)
purrr::map(age, f) %>%
  dplyr::bind_rows()
