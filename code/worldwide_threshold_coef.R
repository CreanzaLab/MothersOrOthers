if (! dir.exists('../figures')) dir.create('../figures')
if (! dir.exists('../figures/nuc_plots/')) dir.create('../figures/nuc_plots/')

# function that finds the best-fit coefficients for pairs of populations across
#  the world using nuclear genetic data
nuclear_coefficient_fxn <- function(
    spat_thresh = 10000, nuclear_data = nuc_list, continents = NA) {
  
  # select population pairs within the spatial threshold
  thresh <- which((nuclear_data$spat <= spat_thresh) &
                    (nuclear_data$spat > 0))
  
  # select population pairs on the continents provided
  if (!all(is.na(continents))) {
    thresh <- thresh[nuclear_data$continent[thresh] %in% continents]
  }
  
  # continue analysis only if 6 or more population pairs exist, otherwise
  #  return NA
  if (length(thresh) <= 6) {
    to_ret <- rep(NA, 4)
    names(to_ret) <- c(
      "(Intercept)", "pop_dist_geo",
      "pop_dist_gene", "pop_dist_gene_x")
    return(list(coefficients=to_ret))
  }
  
  # unique/total shared phonemes across the population pairs
  pop_phon <- nuclear_data$pop_phon[thresh, ]
  
  # spatial distances across population pairs
  pop_dist_geo <- log1p(nuclear_data$spat[thresh])
  pop_dist_geo <- (pop_dist_geo / sd(pop_dist_geo))
  # genetic (autosomal) distances across pop pairs
  pop_dist_gene <- nuclear_data$gene[thresh]
  pop_dist_gene <- lm(pop_dist_gene ~ pop_dist_geo)$resid
  pop_dist_gene <- (pop_dist_gene / sd(pop_dist_gene))
  # genetic (x chr.) distances across pop pairs
  pop_dist_gene_x <- nuclear_data$gene_x[thresh]
  pop_dist_gene_x <- lm(pop_dist_gene_x ~ pop_dist_gene + pop_dist_geo)$resid
  pop_dist_gene_x <- (pop_dist_gene_x / sd(pop_dist_gene_x))
  
  # initialize a STAN model
  stan_model <- quap(
    # model structure
    alist(
      # binomial model, predicting the number of shared phonemes based on
      #  possible phonemes (in either language)
      phon_uniq ~ dbinom(phon_total, p),
      logit(p) <- a + spat_coef * spat +
        gene_coef * gene + gene_x_coef * gene_x,
      a ~ dnorm(0, 3),
      spat_coef ~ dnorm(0, 0.5),
      gene_coef ~ dnorm(0, 0.5),
      gene_x_coef ~ dnorm(0, 0.5)
    ),
    # inputs
    data = list(
      phon_uniq = pop_phon[, 1],
      phon_total = pop_phon[, 2],
      spat = pop_dist_geo,
      gene = pop_dist_gene,
      gene_x = pop_dist_gene_x
    ),
    # initial parameters
    start = list(
      spat_coef = 0,
      gene_coef = 0,
      gene_x_coef = 0
    )
  )
  # compute values for the STAN model
  mu <- link(stan_model)
  mu_mean <- apply(mu, 2, mean)
  posterior <- extract.samples(stan_model)
  stan_resid <- pop_phon[, 1] / pop_phon[, 2] - mu_mean
  
  # freq approach for p-values:
  p_val <- summary(glm(
    # binomial glm with unique features V1 and shared features (V2 - V1)
    cbind(V1, V2 - V1) ~ pop_dist_geo + pop_dist_gene + pop_dist_gene_x,
    # make a temporary data frame with phonemes (unique/total), spatial dists.,
    #  genetic dists, and x-chromosome genetic dists.
    data = as.data.frame(cbind(
      pop_phon, pop_dist_geo, pop_dist_gene, pop_dist_gene_x)),
    family = binomial))$coefficient[2:4, ]
  coef_val <- p_val[, 1]
  p_val <- p_val[, 4]
  return_list <- list(posterior = posterior, resid = stan_resid,
                      thresh = thresh, p_values = p_val, coef_values = coef_val)
  return(return_list)
};

###

# run and plot the nuclear genetic model using pairs of
#  populations from around the world
run_worldwide_nuclear_coef_mdl <- function() {
  # the geographic/spatial thresholds at which to calculate model values
  model_distances <- (4:40) * 250
  # output of the model for each threshold
  nuclear_coef_out <- lapply(model_distances, function(distance)
    nuclear_coefficient_fxn(distance))
  
  # take credible intervals of each coefficient (spatial, gene, and gene_x)
  #  i.e. the middle 89% of the posterior distribution for each coefficient; 
  #  this posterior describes likely values of each coefficient in the model.
  coefs_ci <- lapply(nuclear_coef_out, function(d) lapply(
    d$posterior, PI))
  coefs_ci <- sapply(coefs_ci, function(d) 
    d[c("spat_coef", "gene_coef", "gene_x_coef")])
  # find the coefficient values favored by the model (median of the posterior)
  coefs <- apply((sapply(nuclear_coef_out, function(d)
    lapply(d$posterior, median))), 1, as.numeric)
  # create a table of p-values made from equivalent non-bayesian models 
  p_val_table <- sapply(nuclear_coef_out, function(d) d$p_values)
  coef_val_table <- sapply(nuclear_coef_out, function(d) d$coef_values)
  
  # data frame for plotting
  coef_plotting_df <- data.frame(
    coef_distance = unlist(lapply(model_distances, rep, 3)),
    Coefficient = rep(c("spat_coef", "gene_coef", "gene_x_coef"), nrow(coefs)),
    coef_value = as.vector(t(coefs[, 1:3])),
    coef_ci_min = sapply(coefs_ci, function(d) d[1]),
    coef_ci_max = sapply(coefs_ci, function(d) d[2]),
    p_values = as.vector(p_val_table),
    coef_values = as.vector(coef_val_table)
  )
  
  # create a line plot with the posterior distributions of each coefficient
  nuc_ww_plot <- ggplot(data = coef_plotting_df, aes(
    # for this plot: x-axis is threshold distance, y-axis is the median value
    #  of the coefficient posterior at that threshold.
    x = coef_distance, y = coef_value, group = Coefficient)) +
    xlab("Threshold for spatial distance (km)") +
    ylab("Relationship to linguistic distances") +
    guides(fill = element_blank(),
           colour = NULL) +
    geom_line(aes(color = Coefficient)) +
    scale_x_continuous(
      lim = range(model_distances),
      breaks = (1:10) * 1000,
      expand = expansion(mult = c(0, 0))) +
    scale_color_manual(
      breaks = c("spat_coef", "gene_coef", "gene_x_coef"),
      values = c("black", "blue", "red"), labels = c(
        "Spatial distance",
        "Add'l effect of\nAutosomal distance",
        "Add'l effect of\nX-chrom. distance")) +
    # Add a visualization of the middle 89% percentiles of the posterior for
    #  each coefficient.
    geom_ribbon(aes(
      x = coef_distance,
      ymin = coef_ci_min, ymax = coef_ci_max,
      group = Coefficient, fill = Coefficient)) +
    scale_fill_manual(
      breaks = c("spat_coef", "gene_coef", "gene_x_coef"),
      values = c('#00000030', '#0000e030', '#e0000030'),
      labels = c(
        "Spatial distance",
        "Add'l effect of\nAutosomal distance",
        "Add'l effect of\nX-chrom. distance")) +
    geom_hline(yintercept = 0, color = '#80808030') +
    theme_classic()
  
  ggsave("../figures/nuc_plots/thresh_worldwide_log.pdf", nuc_ww_plot,
         width=7, height=5, units='in')
  # 
  # # plot p-values
  # nuc_ww_pval_plot <- ggplot(data = coef_plotting_df, aes(
  #   x = coef_distance, y = coef_values, group = Coefficient)) +
  #   xlab("Threshold for spatial distance (km)") +
  #   # ylab("log10(p-value)") +
  #   ylab("coefficient value") +
  #   guides(fill = element_blank(),
  #          colour = NULL) +
  #   geom_line(aes(color = Coefficient)) +
  #   scale_x_continuous(
  #     lim = range(model_distances),
  #     breaks = (1:10) * 1000,
  #     expand = expansion(mult = c(0, 0))) +
  #   # scale_y_continuous(trans = "log10") +
  #   scale_color_manual(
  #     breaks = c("spat_coef", "gene_coef", "gene_x_coef"),
  #     values = c("black", "blue", "red"), labels = c(
  #       "Spatial distance",
  #       "Add'l effect of\nAutosomal distance",
  #       "Add'l effect of\nX-chrom. distance")) +
  #   geom_hline(yintercept = 0.0, color = '#80808030') +
  #   theme_classic()
  # 
  # ggsave("../figures/nuc_plots/thresh_worldwide_coef.pdf", nuc_ww_pval_plot,
  #        width=7, height=5, units='in')
}; run_worldwide_nuclear_coef_mdl()


# function that finds the best fit coefficients for pairs of populations across
#  the world using mtDNA data
mt_coef_fxn <- function(
    spat_thresh = 10000, continents = NA,
    trait="matrilocality") {
  
  # which pairs of populations to use
  thresh <- which(mt_data$geo_dist <= spat_thresh &
                    as.dist(mt_data$min_samp) >=3)
  
  # subset continents
  if (!all(is.na(continents))) {
    thresh <- thresh[mt_data$continent[thresh] %in% continents]
  }
  
  # select relevant ethnographic trait
  if (trait == "matriliny") {
    pop_ethnographic <- ceiling(as.dist(mt_data$matrilineal)[thresh])
  }
  if (trait == "matrilocality") {
    pop_ethnographic <- ceiling(as.dist(mt_data$matrilocal)[thresh])
  }
  thresh <- thresh[! is.na(pop_ethnographic)]
  
  pop_ethnographic <- as.factor(pop_ethnographic[! is.na(pop_ethnographic)])
  
  # make sure there are more than six neighboring populations
  if (length(thresh) <= 6) {
    to_ret <- list(
      a = rep(NA, 2),
      spat_coef = matrix(NA, 5, 2),
      gene_coef = matrix(NA, 5, 2)
    )
    return(list(posterior = to_ret, resid = NA, thresh = NA))
  }
  
  # phoneme counts for the population pairs
  pop_phon <- mt_data$mt_phon[thresh, ]
  
  # spatial distances for the population pairs
  pop_dist_geo <- mt_data$geo_dist[thresh]
  pop_dist_geo <- (pop_dist_geo / sd(pop_dist_geo))
  # genetic (mtDNA) distances for the population pairs
  pop_dist_gene <- mt_data$mt_dist[thresh]
  pop_dist_gene <- lm(pop_dist_gene ~ pop_dist_geo)$resid
  pop_dist_gene <- (pop_dist_gene / sd(pop_dist_gene))
  
  # initialize the STAN model
  stan_model <- quap(
    alist(
      # Binomial(shared phonemes, total phonemes,p=
      #  logit(intercept + effect of geography + effect of genetics))
      #  for population pairs where at least one has the ethnographic trait
      #  in question.
      phon_shared ~ dbinom(phon_total, p),
      logit(p) <- a +
        spat_coef[ethn] * spat +
        gene_coef[ethn] * gene,
      a ~ dnorm(0, 3),  # intercept
      # spatial coefficient for ethn. trait = {0, 1}
      spat_coef[ethn] ~ dnorm(0, 0.5),
      # genetic coefficient for ethn. trait = {0, 1}
      gene_coef[ethn] ~ dnorm(0, 0.5)
    ),
    data = list(
      phon_shared = pop_phon[, 1],
      phon_total = pop_phon[, 2],
      spat = pop_dist_geo,
      gene = pop_dist_gene,
      ethn = pop_ethnographic
    )
  )
  # run the model and get a posterior
  mu <- link(stan_model, n = 1e4)
  mu_mean <- apply(mu, 2, mean)
  posterior <- extract.samples(stan_model)
  stan_resid <- pop_phon[, 1] / pop_phon[, 2] - mu_mean
  
  # freq approach for p-values:
  df <- data.frame(
    phon_shared = pop_phon[, 1],
    phon_total = pop_phon[, 2],
    spat = pop_dist_geo,
    gene = pop_dist_gene,
    ethn = as.factor(pop_ethnographic)
  )
  
  if (length(unique(pop_ethnographic)) > 1) {
    # if ethnographic traits vary:
    p_val <- summary(glm(
      # binomial glm 
      cbind(phon_shared, phon_total - phon_shared) ~
        spat * ethn + gene * ethn - ethn,
      family = binomial,
      data = df))$coefficient
    p_val <- p_val[2:nrow(p_val), ]
  } else {
    # if all ethnographic traits are identical:
    p_val <- summary(glm(
      # binomial glm 
      cbind(phon_shared, phon_total - phon_shared) ~
        spat + gene,
      family = binomial,
      data = df))$coefficient[2:3, ]
  }
  
  coef_val <- p_val[, 1]
  p_val <- p_val[, 4]
  
  return_list <- list(posterior = posterior, resid = stan_resid,
                      thresh = thresh, p_values = p_val, coef_values = coef_val)
  return(return_list)
};

# run and plot the mtDNA genetic model using pairs of
#  populations from around the world, including ethnographic traits
#  and regional subsetting
run_worldwide_mt_coef_mdl <- function() {
  
  # get regional subsets
  Americas <- (function(x) x[grepl("America", x, )])(unique(mt_data$continent))
  Africa <- (function(x) x[grepl("Africa", x, )])(unique(mt_data$continent))
  Central_East_Asia_and_Pacific <- (function(x) x[
    grepl("Asia|Pacific", x, )])(unique(mt_data$continent))
  West_Eurasia <- (function(x) x[grepl("Europe|Middle", x, )])(
    unique(mt_data$continent))
  
  # the function will generate a plot for each ethnographic trait and
  #  regional subset
  mt_plot_list <- list()
  mt_pval_plot_list <- list()
  plot_num <- 0
  for (ethn_trait in c("matrilocality", "matriliny")) {
    for (region_subset in c(
      "All",
      "Central & East Asia and Pacific",
      "Americas", "Africa", "West Eurasia"
      )) {
      print(region_subset)
      plot_num <- plot_num + 1
      if (region_subset == "All") {
        continents <- NA
        continents_name <- "Worldwide"
      } else {
        continents <- get(gsub("( |-|,|&)+", "_", region_subset))
        continents_name <- paste("in", region_subset)
      }
      
      # run the model for the specified region and ethnographic trait
      model_distances <- (4:40) * 250
      mt_coefmdl_stan <- lapply(model_distances, function(d) mt_coef_fxn(
        d, continents=continents, trait=ethn_trait))
      
      mt_coefs <- list()
      mt_model_ci <- list()
      # for both spat_coef & gene_coef for each corresponding ethnographic
      #  trait, find the median value and the credibility interval from
      #  their posteriors.
      for (cf in c("spat_coef", "gene_coef")) {
        if (region_subset != "West Eurasia") {
          mt_coefs[[cf]] <- sapply(
            mt_coefmdl_stan,
            function(d) {
              median(d$posterior[[cf]][, 1])
            }
          )
          mt_coefs[[paste0(cf, "_m")]] <- sapply(
            mt_coefmdl_stan,
            function(d) {
              median(d$posterior[[cf]][, 2] - d$posterior[[cf]][, 1])
            }
          )
          mt_model_ci[[cf]] <- sapply(
            mt_coefmdl_stan,
            function(d) {
              if (all(is.na(d$posterior[[cf]][, 1]))) {
                return(as.numeric(c(NA, NA)))
              } else {
                PI(d$posterior[[cf]][, 1])
              }
            }
          )
          mt_model_ci[[paste0(cf, "_m")]] <- sapply(
            mt_coefmdl_stan,
            function(d) {
              if (all(is.na(d$posterior[[cf]][, 1]))) {
                return(as.numeric(c(NA, NA)))
              } else {
                PI(d$posterior[[cf]][, 2] - d$posterior[[cf]][, 1])
              }
            }
          )
        } else {
          # For west Eurasia, only extract a single vector instead of two (
          #  there are only non-matrilocal and non-matrilineal populations here,
          #  so the function returns a slightly different output)
          mt_coefs[[cf]] <- sapply(
            mt_coefmdl_stan,
            function(d) {
              median(d$posterior[[cf]])
            }
          )
          mt_model_ci[[cf]] <- sapply(
            mt_coefmdl_stan,
            function(d) {
              PI(d$posterior[[cf]])
            }
          )
        }
      }
      
      # create tables for plotting data
      if (region_subset != "West Eurasia") {
        # extract coefficient posteriors and p-values
        p_val_table <- sapply(mt_coefmdl_stan, function(d) {
          n <- rep(NA, 4)
          n[1:length(d$p_values)] <- d$p_values
          n
          })[c(1, 3, 2, 4), ]
        coef_val_table <- sapply(mt_coefmdl_stan, function(d) {
          n <- rep(NA, 4)
          n[1:length(d$coef_values)] <- d$coef_values
          n
        })[c(1, 3, 2, 4), ]
        
        # collect data into table for ggplot to work with:
        mt_coef_plotting_df <- data.frame(
          # spatial thresholds
          coef_distance = unlist(lapply(model_distances, rep, 4)),
          # coefficient names
          Coefficient = rep(c("spat_coef", "spat_coef_m", "gene_coef", "gene_coef_m"),
                            length(mt_coefs$spat_coef)),
          # coefficient values
          coef_value = as.vector(t(as.matrix(as.data.frame(mt_coefs))[, 1:4])),
          # coefficient credibility intervals
          coef_ci_min = as.vector(t(sapply(mt_model_ci, function(d) d[1, ])[, 1:4])),
          coef_ci_max = as.vector(t(sapply(mt_model_ci, function(d) d[2, ])[, 1:4])),
          # p-values
          p_vals = as.vector(p_val_table),
          coef_vals = as.vector(coef_val_table)
        )
      } else {
        p_val_table <- sapply(mt_coefmdl_stan, function(d) d$p_values[1:2])
        coef_val_table <- sapply(mt_coefmdl_stan, function(d) d$coef_values[1:2])
        mt_coef_plotting_df <- data.frame(
          # spatial thresholds
          coef_distance = unlist(lapply(model_distances, rep, 2)),
          # coefficient names
          Coefficient = rep(c("spat_coef", "gene_coef"),
                            length(mt_coefs$spat_coef)),
          # coefficient values
          coef_value = as.vector(t(as.matrix(as.data.frame(mt_coefs))[, 1:2])),
          # coefficient credibility intervals
          coef_ci_min = as.vector(t(sapply(mt_model_ci, function(d) d[1, ])[, 1:2])),
          coef_ci_max = as.vector(t(sapply(mt_model_ci, function(d) d[2, ])[, 1:2])),
          # p-values
          p_vals = as.vector(p_val_table),
          coef_vals = as.vector(coef_val_table)
        )
      }
      
      ##
      # a line plot with threshold distances on the x-axis and coefficient
      #  values on the y-axis for mtDNA data
      mt_plot_list[[plot_num]] <- ggplot(data = mt_coef_plotting_df, aes(
        x = coef_distance, y = coef_value, group = Coefficient)) +
        xlab("Spatial distance threshold for pairs of populations (km)") +
        ylab(paste0("Contribution to linguistic variation\n(", continents_name, ")")) +
        guides(fill = element_blank(),
               colour = NULL) +
        geom_line(aes(color = Coefficient, linetype = Coefficient)) +
        scale_x_continuous(
          lim = range(model_distances),
          breaks = (0:10) * 1000,
          expand = expansion(mult = c(0, 0)))
      if (region_subset != "West Eurasia") {
        mt_plot_list[[plot_num]] <- mt_plot_list[[plot_num]] + scale_linetype_manual(
          values = c("solid", "solid", "dashed", "dashed"),
          breaks = c("spat_coef", "gene_coef", "spat_coef_m", "gene_coef_m"),
          labels = c(
            "Spatial distance",
            "mtDNA distance\n(non-spatial component)",
            paste0("Spatial distance\n(add'l effect of\n", ethn_trait, ")"),
            paste0("mtDNA distance\n(add'l effect of\n", ethn_trait, ")"))) +
          scale_color_manual(
            values = c("black", "red", "black", "red"),
            breaks = c("spat_coef", "gene_coef", "spat_coef_m", "gene_coef_m"),
            labels = c(
              "Spatial distance",
              "mtDNA distance\n(non-spatial component)",
              paste0("Spatial distance\n(add'l effect of\n", ethn_trait, ")"),
              paste0("mtDNA distance\n(add'l effect of\n", ethn_trait, ")"))) +
          geom_ribbon(aes(
            x = coef_distance,
            ymin = coef_ci_min, ymax = coef_ci_max,
            group = Coefficient, fill = Coefficient)) +
          scale_fill_manual(
            values = c('#00000030', '#e0000030', '#00000030', '#e0000030'),
            breaks = c("spat_coef", "gene_coef", "spat_coef_m", "gene_coef_m"),
            labels = c(
              "Spatial distance",
              "mtDNA distance\n(non-spatial component)",
              paste0("Spatial distance\n(add'l effect of\n", ethn_trait, ")"),
              paste0("mtDNA distance\n(add'l effect of\n", ethn_trait, ")"))) +
          geom_hline(yintercept = 0, color = '#80808030') +
          theme_classic()
      } else {
        mt_plot_list[[plot_num]] <- mt_plot_list[[plot_num]] + scale_linetype_manual(
          values = c("solid", "solid"),
          breaks = c("spat_coef", "gene_coef"),
          labels = c(
            "Spatial distance",
            "mtDNA distance\n(non-spatial component)")) +
          scale_color_manual(
            values = c("black", "red"),
            breaks = c("spat_coef", "gene_coef"),
            labels = c(
              "Spatial distance",
              "mtDNA distance\n(non-spatial component)")) +
          geom_ribbon(aes(
            x = coef_distance,
            ymin = coef_ci_min, ymax = coef_ci_max,
            group = Coefficient, fill = Coefficient)) +
          scale_fill_manual(
            values = c('#00000030', '#e0000030'),
            breaks = c("spat_coef", "gene_coef"),
            labels = c(
              "Spatial distance",
              "mtDNA distance\n(non-spatial component)")) +
          geom_hline(yintercept = 0, color = '#80808030') +
          theme_classic()
      }
      
      mt_pval_plot_list[[plot_num]] <- ggplot(data = mt_coef_plotting_df, aes(
        x = coef_distance, y = coef_vals, group = Coefficient)) +
        xlab("Spatial distance threshold for pairs of populations (km)") +
        ylab(paste0("coefficient values \n(", continents_name, ")")) +
        guides(fill = element_blank(),
               colour = NULL) +
        geom_line(aes(color = Coefficient, linetype = Coefficient)) +
        scale_x_continuous(
          lim = range(model_distances),
          breaks = 2 * (0:5) * 1000,
          expand = expansion(mult = c(0, 0))) # +
        # scale_y_continuous(
        #   transform = "log",
        #   breaks = c(10 ** c(-3, -5, (-1 - 10 * (0:20))))
        #   ) + ggforce::facet_zoom(ylim = c(0.2, 0.00001),
        #                           zoom.size = 1,
        #                           horizontal = F)
      if (region_subset != "West Eurasia") {
        mt_pval_plot_list[[plot_num]] <- mt_pval_plot_list[[plot_num]] +
          scale_linetype_manual(
            values = c("solid", "solid", "dashed", "dashed"),
            breaks = c("spat_coef", "gene_coef", "spat_coef_m", "gene_coef_m"),
            labels = c(
              "Spatial distance",
              "mtDNA distance\n(non-spatial component)",
              paste0("Spatial distance\n(add'l effect of\n", ethn_trait, ")"),
              paste0("mtDNA distance\n(add'l effect of\n", ethn_trait, ")"))) +
          scale_color_manual(
            values = c("black", "red", "black", "red"),
            breaks = c("spat_coef", "gene_coef", "spat_coef_m", "gene_coef_m"),
            labels = c(
              "Spatial distance",
              "mtDNA distance\n(non-spatial component)",
              paste0("Spatial distance\n(add'l effect of\n", ethn_trait, ")"),
              paste0("mtDNA distance\n(add'l effect of\n", ethn_trait, ")"))) +
          geom_hline(yintercept = 0.0, color = '#80808030') +
          theme_classic()
      } else {
        mt_pval_plot_list[[plot_num]] <- mt_pval_plot_list[[plot_num]] +
          scale_linetype_manual(
            values = c("solid", "solid"),
            breaks = c("spat_coef", "gene_coef"),
            labels = c(
              "Spatial distance",
              "mtDNA distance\n(non-spatial component)")) +
          scale_color_manual(
            values = c("black", "red"),
            breaks = c("spat_coef", "gene_coef"),
            labels = c(
              "Spatial distance",
              "mtDNA distance\n(non-spatial component)")) +
          geom_hline(yintercept = 0.0, color = '#80808030') +
          theme_classic()}
    } # end region_subset
  } # end ethn_trait
    
      
  # # save the plots made with mtDNA data
  # plot_num <- 0; for (
  #   ethn_trait in c("matrilocality", "matriliny")) {
  #   for (region_subset in c("world", "ceasiapacific",
  #                           "americas", "africa", "weurasia")) {
  #     plot_num <- plot_num + 1
  #     ggsave(filename = paste0("../figures/plot_mt_", region_subset, "_", ethn_trait, ".pdf"),
  #            plot = mt_plot_list[[plot_num]], device = "pdf", dpi = 300,
  #            height = 4, width = 8)
  #   }
  # }
  #     
  # save the p-value plots made with mtDNA data
  plot_num <- 0; for (
    ethn_trait in c("matrilocality", "matriliny")) {
    for (region_subset in c(
      "world" # , "ceasiapacific",
                  # "americas", "africa", "weurasia"
      )) {
      plot_num <- plot_num + 1
      ggsave(filename = paste0("../figures/plot_mt_", region_subset, "_", ethn_trait, "_coef.pdf"),
             plot = mt_pval_plot_list[[plot_num]], device = "pdf", dpi = 300,
             height = 4, width = 8)
    }
  }
}; run_worldwide_mt_coef_mdl()
  