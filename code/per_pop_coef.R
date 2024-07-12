# create a base map to use in plotting per-population parameters
# The map goes from 23 W to 337 E and 55 S to 75 N, which captures all samples.
world_map <- map_data("world", wrap=c(-23, 337), ylim=c(-55, 75))
base_world <- ggplot() +
  xlab("") + ylab("") +
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group),
               colour="#ffffff00", fill="#00002020") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = '#ffffff00', colour = '#ffffff00'), 
        axis.line = element_line(colour = "#ffffff00"), legend.position="bottom",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

# for each individual location:
#  evaluate the contributions to phonemic distance to neighbors (w/in X km) 
#  based on its genetic distance, x-chr distance, and geographic distance

per_pop_nuclear_coef_fxn <- function(
    focal_pop = 1,
    spat_thresh = 2500,
    omit_pops = NA) {
  
  # only use on pairs of populations that include the focal pop
  focal_dists <- as.matrix(nuc_list$spat) * 0
  focal_dists[, focal_pop] <- focal_dists[focal_pop, ] <- 1
  focal_dists[focal_pop, focal_pop] <- 0
  
  if (all(! is.na(omit_pops))) {
    focal_dists[, omit_pops] <- focal_dists[omit_pops, ] <- 0
  }
  
  # select population pairs within the spatial threshold
  thresh <- which((nuc_list$spat <= spat_thresh
                   ) & as.dist(focal_dists) == 1 &
                    (nuc_list$spat > 0))
  
  # skip pops with fewer than 7 neighboring populations
  if (length(thresh) <= 5) {
    
    to_ret <- rep(NA, 4)
    names(to_ret) <- c(
      "spat_coef", "gene_coef",
      "gene_x_coef", "a")
    return(list(coefficients=to_ret))
  }
  
  # number of phonemes per population pair (unique, shared)
  pop_phon <- nuc_list$pop_phon[thresh, ]

  # population pair spatial distance matrix
  pop_dist_geo <- nuc_list$spat[thresh]
  pop_dist_geo <- (pop_dist_geo / sd(pop_dist_geo))
  # population pair autosomal genetic distance matrix
  pop_dist_gene <- nuc_list$gene[thresh]
  pop_dist_gene <- lm(pop_dist_gene ~ pop_dist_geo)$resid
  pop_dist_gene <- (pop_dist_gene / sd(pop_dist_gene))
  # population pair X-chrom. genetic distance matrix
  pop_dist_gene_x <- nuc_list$gene_x[thresh]
  pop_dist_gene_x <- lm(pop_dist_gene_x ~ pop_dist_gene + pop_dist_geo)$resid
  pop_dist_gene_x <- (pop_dist_gene_x / sd(pop_dist_gene_x))
  
  # initialize STAN model
  stan_model <- quap(
    alist(
      phon_shared ~ dbinom(phon_total, p),
      logit(p) <- a + spat_coef * spat +
        gene_coef * gene + gene_x_coef * gene_x,
      a ~ dnorm(0, 3),  # intercept
      spat_coef ~ dnorm(0, 0.5),  # spatial coefficient
      gene_coef ~ dnorm(0, 0.5),  # autosomal genetic coefficient
      gene_x_coef ~ dnorm(0, 0.5)  # X chrom genetic coefficient
    ),
    data = list(
      phon_shared = pop_phon[, 1],
      phon_total = pop_phon[, 2],
      spat = pop_dist_geo,
      gene = pop_dist_gene,
      gene_x = pop_dist_gene_x
    ),
    start = list(
      spat_coef = 0,
      gene_coef = 0,
      gene_x_coef = 0
    )
  )
  mu <- link(stan_model)
  mu_mean <- apply(mu, 2, mean)
  posterior <- extract.samples(stan_model)
  stan_resid <- pop_phon[, 1] / pop_phon[, 2] - mu_mean
  
  df <- data.frame(
    phon_shared = pop_phon[, 1],
    phon_total = pop_phon[, 2],
    spat = pop_dist_geo,
    gene = pop_dist_gene,
    gene_x = pop_dist_gene_x
  )
  

  # if ethnographic traits vary:
  p_val <- summary(glm(
    # binomial glm 
    cbind(phon_shared, phon_total - phon_shared) ~
      spat + gene + gene_x,
    family = binomial,
    data = df))$coefficient
  p_val <- p_val[2:4, 4]
  
  return_list <- list(posterior = posterior, resid = stan_resid,
                      thresh = thresh, p_values = p_val)
  return(return_list)
};
# per_pop_coef_results <- lapply(1:nrow(nuc_list$meta), per_pop_nuclear_coef_fxn)
# names(per_pop_coef_results) <- nuc_list$pop_order


##
run_per_pop_nuclear_coef <- function() {
  
  # run for all populations with data, then restrict to populations with
  #   endogamy/community marriage organization data available:
  for (ethn_trait in c(NA, "endo")) {
    thresholds_for_per_pop <- c(1500, 2500, 5000, 10000)
    # run the per-pop coefficient results for each spatial distance threshold
    
    omit_pops <- NA
    
    # select pops with community marriage organization data available
    endo_pops <- (
      nuc_list$meta[, paste0("Code..EA015.Community.marriage.organization.",
                        ".1.Exogamous..2.Agamous..3.Endogamous.")] == 3)
    endo_pops <- ifelse(endo_pops == TRUE, "endo", "nonendo")
    
    # if using ethnographic data, select populations with data available
    if (! is.na(ethn_trait)) {
      omit_pops <- which(is.na(endo_pops))
    }
    per_thresh_per_pop_coef_results <- lapply(
      thresholds_for_per_pop,
      function(n) {
        to_ret <- lapply(
          1:nrow(nuc_list$meta),
          per_pop_nuclear_coef_fxn,
          n,  # threshold
          omit_pops)
        names(to_ret) <- nuc_list$meta$pop_name
        to_ret
    })
    names(per_thresh_per_pop_coef_results) <- as.character(thresholds_for_per_pop)
    
    # plotting 3 maps per distance threshold
    threshold_maps <- list()
    for (i in 1:length(thresholds_for_per_pop)) {
      # one plot per coefficient:
      for (coef_for_plotting in c("spat_coef", "gene_coef", "gene_x_coef")) {
        coef_for_plotting_full <- switch(
          coef_for_plotting,
          "spat_coef" = "Spatial",
          "gene_coef" = "Autosomal",
          "gene_x_coef" = "X chrom."
        )
        plotting_table <- data.frame(
          long = nuc_list$meta$longitude,
          lat = nuc_list$meta$latitude,
          ethn = endo_pops,
          colors = sapply(per_thresh_per_pop_coef_results[[i]], function(n) {
            if (all(is.null(n$posterior))) {
              return(NA)
            } else {
              # return 0 for the coefficient value if the credibility interval
              #  intersects with 0 (i.e. no good evidence that the coefficient
              #  value is non-zero).
              if (prod(PCI(n$posterior[[coef_for_plotting]])) > 0) {
                return(
                  # find the median value of the posterior for each coefficient
                  (function (n) ifelse(is.null(n), NA, n))(
                    median(n$posterior[[coef_for_plotting]])))
              } else {
                return(0)
              }}
          }))
        # update latitudes for proper wrapping around the world
        plotting_table[, "long"] <- plotting_table[, "long"] %>% {
            ifelse(. < -25, . + 360, .)}
        # reorder so that the points with the greatest values (positive 
        #  or negative) are plotted on top
        plotting_table <- plotting_table[
          order(abs(plotting_table[, "colors"]), na.last = FALSE), ]
        if (is.na(ethn_trait)) {
          # if using all populations, initialize the world map:
          colored_world <-
            base_world +
            geom_point(data=plotting_table,
                       aes(x = long, y = lat, fill = colors),
                       size=3, color = "#000000d0", shape = 21, show.legend = T)
        } else {
          # if using only using populations with ethnographic data,
          #  initialize the world map and include a shape-scale to label
          #  ethnographic data for each population
          plotting_table <- plotting_table[which(!is.na(
            plotting_table$ethn
          )) ,]
          colored_world <- 
            base_world +
            geom_point(data=plotting_table,
                       aes(x = long, y = lat, fill = colors, shape = ethn),
                       size=3, color = "#000000d0", show.legend = T) + 
            scale_shape_manual(
              breaks = c("nonendo", "endo"),
              values = c(21, 22),
              labels = c("Non-Endogamous", "Endogamous")
            )
        }
        colored_world <- colored_world +
          # add a color scale to label coefficient value for each population
          scale_fill_steps2(
            midpoint = 0, low = '#1b7837', mid = '#f8f8f8', high = '#762a83',
            na.value = '#ffa0a090', nice.breaks = T, name = paste(
              "Effect of", coef_for_plotting_full,
              "var. on \nlanguage,", thresholds_for_per_pop[i], "km threshold"
            ), limits = c(-0.4, 0.4),
            labels = function(breaks) {
              legend_labels = c()
              for (b in breaks) {
                if (((b * 10) %% 2) < 1e-10) {
                  legend_labels <- c(legend_labels, as.character(b))
                } else {
                  legend_labels <- c(legend_labels, "")
                }}
              legend_labels},
            breaks = c(
              seq(-0.5, -0.05, 0.05),
              seq(0.05, 0.4, 0.05)
            )) +
          guides(fill = guide_colorbar(
            title.position = "bottom", ticks = F
          ))
        
        threshold_maps[[length(threshold_maps) + 1]] <- colored_world
        if (is.na(ethn_trait)) {
          plot_file_name <- paste0(
            "../figures/nuc_plots/nuc_map_", thresholds_for_per_pop[i],
            "_", coef_for_plotting, ".pdf"
          )
        } else {
          plot_file_name <- paste0(
            "../figures/nuc_plots/nuc_map_endo_", thresholds_for_per_pop[i],
            "_", coef_for_plotting, ".pdf"
          )
        }
        ggsave(plot_file_name, colored_world, device = "pdf",
               width = 10, height = 6, units = "in", dpi = 300)
      }
    }
  }
}; run_per_pop_nuclear_coef()


# For mitochondrial data, find the coefficients for each population of 
#  Prob(syllable sharing) as a function of mtDNA and geographic distance b/w
#  populations, and how their coefficients vary based on the ethnographic traits
#  found in these populations (post-marital residence, descent patterns, and 
#  community marriage organization)
mt_per_pop <- function(
    focal_pop = 1,
    spat_thresh = 2500,
    ethn_trait = NA) {
  
  focal_dists <- mt_data$min_samp * 0
  focal_dists[, focal_pop] <- focal_dists[focal_pop, ] <- 1
  focal_dists[focal_pop, focal_pop] <- 0
  
  # if subsetting by ethnographic trait, only include indivviduals with a non-NA
  #  value for that trait.
  if (! is.na(ethn_trait)) {
    for (i in which(is.na(switch(
      ethn_trait,
      endoexo = mt_data$meta$Code..EA015.Community.marriage.organization..1.Exogamous..2.Agamous..3.Endogamous.,
      residence = mt_data$meta$Residence,
      descent = mt_data$meta$Descent
    )))) {
      focal_dists[i, ] <- 0
      focal_dists[, i] <- 0
    }
  }
  
  # find pairs of populations within the spatial threshold distance and where
  #  both populations have at least 3 individuals sampled
  thresh <- which(
    mt_data$geo_dist <= spat_thresh &
      as.dist(focal_dists) == 1 &
      as.dist(mt_data$min_samp) >=3)
  
  # if there are fewer than 6 population pairs, return NA (not enough samples)
  if (length(thresh) <= 5) {
    to_ret <- rep(NA, 3)
    names(to_ret) <- c(
      "(Intercept)",
      "pop_dist_geo", "pop_dist_gene")
    return(list(coefficients=to_ret))
  }

  # phonemes shared for pairs of populations  
  pop_phon <- mt_data$mt_phon[thresh, ]
  
  # spatial distances for pairs of populations
  pop_dist_geo <- mt_data$geo_dist[thresh]
  pop_dist_geo <- (pop_dist_geo / sd(pop_dist_geo))
  
  # mtDNA distances for pairs of populations
  pop_dist_gene <- mt_data$mt_dist[thresh]
  pop_dist_gene <- lm(pop_dist_gene ~ pop_dist_geo)$resid
  pop_dist_gene <- (pop_dist_gene / sd(pop_dist_gene))
  
  stan_model <- quap(
    alist(
      phon_shared ~ dbinom(phon_total, p),
      logit(p) <- a +
        spat_coef * spat +
        gene_coef * gene,
      a ~ dnorm(0, 3),
      spat_coef ~ dnorm(0, 0.5),
      gene_coef ~ dnorm(0, 0.5)
    ),
    data = list(
      phon_shared = pop_phon[, 1],
      phon_total = pop_phon[, 2],
      spat = pop_dist_geo,
      gene = pop_dist_gene
    )
  )
  # run the model
  mu <- link(stan_model, n = 1e3)
  # get the model posterior
  posterior <- extract.samples(stan_model)
  return_list <- list(posterior = apply(posterior, 2, median),
                      post_range = apply(posterior, 2, PCI))
  return(return_list)
};

# run and draw plots of the mtDNA model run on a per-population basis
run_and_plot_mt_per_pop <- function() {
  
  ##
  plot_list_mt_residence_descent <- list(); plot_num <- 0
  # for each trait being analyzed:
  for (ethn_trait in c("residence", "descent", "endoexo")) {
    # run the model at a number of spatial thresholds:
    mt_per_pop_res_1000 <- lapply((1:nrow(mt_data$meta)), mt_per_pop, 1000, ethn_trait)
    mt_per_pop_res_1500 <- lapply((1:nrow(mt_data$meta)), mt_per_pop, 1500, ethn_trait)
    mt_per_pop_res_2000 <- lapply((1:nrow(mt_data$meta)), mt_per_pop, 2000, ethn_trait)
    mt_per_pop_res_2500 <- lapply((1:nrow(mt_data$meta)), mt_per_pop, 2500, ethn_trait)
    for (thresh in c(1000, 1500, 2000, 2500)) {
      # at each threshold, get the coefficient values for each pop
      plot_num <- plot_num + 1
      plotting_table <- cbind(
        mt_data$pop_longlat,
        value = sapply(get(paste0("mt_per_pop_res_", thresh)), function(d) {
          if (is.null(d$posterior)) {
            return(NA)
          }
          if (prod(d$post_range) < 0) {
            return(0)
          }
          return(d$posterior[3])
        }),
        as.factor(mt_data$meta$Descent %>% sapply(function(d) ifelse(is.na(d), "None", d))),
        as.factor(as.vector(
          mt_data$meta$Residence %>% sapply(function(d) ifelse(is.na(d), "None", d)))),
        as.factor(as.vector(mt_data$meta[
          , paste0("Code..EA015.Community.marriage.organization",
                   "..1.Exogamous..2.Agamous..3.Endogamous.")] %>% as.character %>%
            sapply(function(d) switch(
              d,
              "1" = "Exogamous",
              "2" = "Agamous",
              "3" = "Endogamous",
              "None"
            ))))
      ); colnames(plotting_table) <-
        c("lat", "long", "value", "descent", "residence", "endoexo"); plotting_table <-
        as.data.frame(plotting_table);plotting_table$long <- plotting_table$long + 360 * (
          plotting_table$long < -25); plotting_table <- 
        plotting_table[order(abs(plotting_table$value), na.last = F), ]
      if (ethn_trait == "descent") {
        colored_world <-
          base_world +
          geom_point(data=plotting_table,
                     aes(x = long, y = lat, fill = value,
                         shape = descent),
                     size=3, color = "#000000d0", show.legend = T) +
          scale_shape_manual(
            "Mode of descent",
            breaks = c("Matrilineal", "NonMatrilineal", "None"),
            labels = c("Female-biased", "Non-Female-biased", "No data"),
            values = c(22, 24, 23)) +
          scale_fill_steps2(
            midpoint = 0,
            low = '#1b7837', mid = '#f8f8f8', high = '#762a83', na.value = '#ffa0a090',
            nice.breaks = T, name = paste(
              "Effect of mtDNA on language,\n", thresh, "km threshold"),
            limits = c(-0.4, 0.4), labels = function(breaks) {
              legend_labels = c()
              for (b in breaks) {
                if (((b * 10) %% 2) < 1e-10) {
                  legend_labels <- c(legend_labels, as.character(b))
                } else {
                  legend_labels <- c(legend_labels, "")
                }
              }
              legend_labels},
            breaks = c(
              seq(-0.5, -0.05, 0.05),
              seq(0.05, 0.4, 0.05)
            )) +
          guides(color = guide_colorbar(
            title.position = "bottom", ticks = F, barwidth = 7
          ))
      }
      if (ethn_trait == "residence") {
        colored_world <-
          base_world +
          geom_point(data=plotting_table,
                     aes(x = long, y = lat, fill = value,
                         shape = residence),
                     size=3, color = "#000000d0", show.legend = T) +
          scale_shape_manual(
            "Marital residence with kin",
            breaks = c("FemaleKin Residence", "NonFemaleKin Residence", "None"),
            labels = c("with female kin", "not with female kin", "No data"),
            values = c(22, 24, 23)) +
          scale_fill_steps2(
            midpoint = 0,
            low = '#1b7837', mid = '#f8f8f8', high = '#762a83', na.value = '#ffa0a090',
            nice.breaks = T, name = paste(
              "Effect of mtDNA on language,\n", thresh, "km threshold"),
            limits = c(-0.4, 0.4), labels = function(breaks) {
              legend_labels = c()
              for (b in breaks) {
                if (((b * 10) %% 2) < 1e-10) {
                  legend_labels <- c(legend_labels, as.character(b))
                } else {
                  legend_labels <- c(legend_labels, "")
                }
              }
              legend_labels},
            breaks = c(
              seq(-0.5, -0.05, 0.05),
              seq(0.05, 0.4, 0.05)
            )) + 
          guides(color = guide_colorbar(
            title.position = "bottom", ticks = F, barwidth = 7
          ))
      
      }
      if (ethn_trait == "endoexo") {
        colored_world <-
          base_world +
          geom_point(data=plotting_table,
                     aes(x = long, y = lat, fill = value,
                         shape = endoexo),
                     size=3, color = "#000000d0", show.legend = T) +
          scale_shape_manual(
            "Community marriage organization",
            breaks = c("Exogamous","Agamous", "Endogamous", "None"),
            labels = c("Exogamous","Agamous", "Endogamous", "None"),
            values = c(21, 22, 24, 23)) +
          scale_fill_steps2(
            midpoint = 0,
            low = '#1b7837', mid = '#f8f8f8', high = '#762a83', na.value = '#ffa0a090',
            nice.breaks = T, name = paste(
              "Effect of mtDNA on language,\n", thresh, "km threshold"),
            limits = c(-0.4, 0.4), labels = function(breaks) {
              legend_labels = c()
              for (b in breaks) {
                if (((b * 10) %% 2) < 1e-10) {
                  legend_labels <- c(legend_labels, as.character(b))
                } else {
                  legend_labels <- c(legend_labels, "")
                }
              }
              legend_labels},
            breaks = c(
              seq(-0.5, -0.05, 0.05),
              seq(0.05, 0.4, 0.05)
            )) + 
          guides(color = guide_colorbar(
            title.position = "bottom", ticks = F, barwidth = 7
          ))
        }
      plot_list_mt_residence_descent[[plot_num]] <- colored_world
    }
  }
  
  plot_num <- 0; for (ethn_trait in c("residence", "descent", "endoexo")) {
    for (thresh in c(1000, 1500, 2000, 2500)) {
      plot_num <- plot_num + 1
      ggsave(paste0("../figures/map_mt_", ethn_trait, "_", thresh, ".pdf"),
             plot_list_mt_residence_descent[[plot_num]], device = "pdf",
             height = 6.5, width = 12, dpi = 300, units = "in")
    }}
}; run_and_plot_mt_per_pop()
