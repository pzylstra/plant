#' Calculate crown parameters for a community at a specific age
#'
#' @param dat The results of run_scm_collect
#' @param tr A table of input traits
#' @param age Years since disturbance
#' @param lat Latitude (degrees)
#' @param map Mean annual precipitation (mm)
#' @param mat Mean annual temperature (degC)
#' 
#' @importFrom dplyr
#' @export
#'

shape_forest <- function(dat, tr, age, lat = -35, map = 1000, mat = 20) {
  
  result <- dat %>%
    plant::tidy_patch() %>% 
    plant:::expand_state()
  
  # Run plantLitter model and add field
  ######################
  
  #Create function indiCrowns
  indiCrowns <- function(H, eta, lat = -35, map = 1000, mat = 20) {
    Hc <- H*0.91*(1-exp(-0.11*eta))
    He <- H*0.96*(1-exp(-0.21*eta))
    Ht <- H*0.96*(1-exp(-0.45*eta))
    W <- 0.2057*H+0.248*sqrt(abs(lat))+0.0005964*map+0.01202*mat-1.977
    Af <- W*(Ht-He)+W*0.5*(H-Ht)+W*0.5*(He-Hc)
    
    res <- data.frame(matrix(ncol = 6, nrow = 1))
    x <- c("top", "base", "he", "ht", "W", "frontalArea")
    colnames(res) <- x
    res[1,]<-c(H, Hc, He, Ht, W, Af)
    
    return(res)
  }
  
  # Interpolate to set time
  full_data <- result$species
  study_age <- interpolate_to_times(full_data, age)
  study_age <- study_age[order(study_age$height),]%>%
    dplyr::filter(time == age, !is.na(height))
  
  # Add individual crown descriptors
  steps <- nrow(study_age)
  sp <- data.frame(matrix(ncol = 6, nrow = steps))
  colnames(sp) <- c("height", "base", "he", "ht", "w", "frontalArea")
  e <- as.vector(tr$eta) # Crown shape
  
  for (r in 1:steps) {
    sp[r,] <- indiCrowns(H=study_age$height[r], eta=e[as.numeric(study_age$species[r])],lat = lat, map = map, mat = mat)
  }
  all <- dplyr::left_join(study_age, sp, by = "height")
  
  # Update species names if tr table used
  #  if(!missing(tr)) {
  #    all <- left_join(all,tr, by = c("species" = "Species")) %>%
  #      mutate(species = name)
  #  }
  
  return(all)
}


#' Divide a plant community into strata
#' 
#' Creates a series of pseudo-transects to randomly sample the output data,
#' then divides these into 2-4 strata using k-means clustering, and
#' choosing the maximum number of significantly divided strata
#'
#' @param dat The output from shape_forest
#' @param sample Proportion of cohorts to test (0-1)
#' @param transects Number of repeats
#' 
#' @importFrom dplyr frame
#' @export
#'

stratify_community <- function(dat, tr, age, lat = -35, map = 1000, mat = 20, sample = 0.5, transects = 10){
  
  # Randomly sample results
  vegStrat <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(vegStrat) <- c("Species", "height", "base", "he", "ht")
  crowns <- shape_forest(dat, tr, age, lat, map, mat)
  steps <- nrow(crowns)
  samples <- ceiling(steps * sample)
  rand <- data.frame(matrix(ncol = transects, nrow = steps))
  for (t in 1:transects) {
    for (r in 1:steps) {
      # Randomise sampling
      rand[r,t] <- runif(n = 1, min = min(crowns$density), max = max(crowns$density))
    }
    # Order plants by the most likely to be sampled (cover x random)
    crowns$weight<-crowns$density*rand[,t]
    samp <- crowns[order(-crowns$weight),] %>%
      select(species, height, base, he, ht)
    # Select the chosen number of samples
    samp <- samp[1:samples, ]
    vegStrat <- rbind(vegStrat,samp)
  }
  
  # Calculate strata from random sample data
  vegStrat$Point <- 1
  vegStrat$top <- vegStrat$height
  vegStrat <- as.data.frame(vegStrat)
  vegStrat <- frame::frameStratify(veg = vegStrat, spName = "species") %>%
    group_by(Stratum, species) %>%
    summarise_if(is.numeric, mean)%>%
    select(Stratum, species, height)
  
  # Interpolate to the stratum division heights, then add these to original data
  breaks <- nrow(vegStrat)-1
  heights <- list()
  for (b in 1:breaks) {
    heights[b]<-vegStrat$height[b]
  }
  heights[breaks+1] <- max(crowns$height, na.rm = TRUE)
  heights <- as.numeric(heights)
  tidy_species_new <- interpolate_to_heights(crowns, heights) %>%
    select(!step)%>%
    mutate(cohort = NA)
  f<- rbind(tidy_species_new, crowns) %>% 
    drop_na() %>%
    mutate(Stratum = NA,
           d = NA)
  
  # Measure density
  out <- f[0,]
  outA <- f[0,]
  for (sp in unique(vegStrat$species)) {
    s <- filter(vegStrat, species == sp)
    outA <- filter(f, species == sp)
    s$height[nrow(s)] <- 1000
    # Identify the stratum of each plant
    outA$Stratum <- cut(as.numeric(outA$height), breaks = c(0,s$height),
                        labels = s$Stratum,
                        include.lowest = TRUE)
    # Integrate density
    outA <- outA[order(outA$height),]
    for (h in 1:nrow(outA)) {
      if (h == 1) {
        outA$d[h] <- max(((outA$density[h]/2) * outA$height[h]),0)
      } else {
        outA$d[h] <- max((mean(outA$density[h], outA$density[h-1]) * (outA$height[h] - outA$height[h-1])),0)
      }
    }
    
    out <- rbind(out, outA) 
  }
  return(out)
  
}


#' Wrapper to grow a forest in plant ready for fire modelling
#' 
#' Grows the forest from an imported table of traits
#' Runs tidy_patch and expand_state
#'
#' @param dat An input table listing species traits
#' @param B_lf1 Potential CO2 photosynthesis at average leaf nitrogen
#' @param dist Mean disturbance interval for investigating trends
#' 
#' @importFrom dplyr
#' @export
#'


grow_forest <- function(dat, B_lf1 = 0.8273474, dist = 1000) {
  
  plant::plant_log_console()
  
  # Collect plant traits
  params <- scm_base_parameters("FF16")
#  param$birth_rate <- 
    params$disturbance_mean_interval <- dist
  eta <- as.vector(dat$eta) # Crown shape
  lma <- as.vector(dat$lma) # kgm−2
  rho <- as.vector(dat$rho) # Wood density - kgm−3
  theta <- as.vector(dat$theta)  # Sapwood per unit leaf area
  a_l1 <- as.vector(dat$a_l1) # Height of plant with leaf area of 1m2
  a_l2 <- as.vector(dat$a_l2) # Exponent of relationship between height and leaf area
  a_b1 <- as.vector(dat$a_b1) # Ratio of bark area to sapwood area
  hmat <- as.vector(dat$hmat) # Height at maturation (m)
  a_f1 <- as.vector(dat$a_f1) # Maximum allocation to reproduction
  a_f2 <- as.vector(dat$a_f2) # Parameter determining rate of change in r(x,ml) around Hmat
  FF16_hyperpar1 <- make_FF16_hyperpar(B_lf1 = B_lf1) # Potential CO2 photosynthesis at average leaf Nitrogen - mold−1 m−2
  
  patch <- expand_parameters(trait_matrix(x = c(eta, lma, rho, theta, a_b1, hmat, a_l1, a_l2, a_f1, a_f2), 
                                          trait_name = c("eta", "lma", "rho", "theta", "a_b1", "hmat", "a_l1", "a_l2", "a_f1", "a_f2")), 
                             params, 
                             FF16_hyperpar1,
                             mutant = FALSE)
  
  patch <- build_schedule(patch)
  dat <- run_scm_collect(patch)
  
  return(dat)
}

#' Updates species names from an optional table
#'
#' @param comm The output from stratify_community
#' @param tr An optional table of input traits
#'  

updateSpecies <- function(comm, tr){
  
  # Update species names if tr table used
  if(!missing(tr)) {
    comm <- left_join(comm,tr, by = c("species" = "Species")) %>%
      mutate(species = name)
  } 
}