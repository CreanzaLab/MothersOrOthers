# function that finds the best fit coefficients for pairs of populations across
#  the world using nuclear genetic data
nuclear_coefficient_ethn_fxn <- function(
    spat_thresh = 10000, nuclear_data = nuc_list, continents = NA) {
  
  # select population pairs within the spatial threshold
  thresh <- which((nuclear_data$spat <= spat_thresh) &
                    (nuclear_data$spat > 0))
  thresh <- thresh[
    ! (is.na(nuclear_data$endogamous[thresh]) +
         is.na(nuclear_data$exogamous[thresh]))]
  
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
  pop_dist_geo <- nuclear_data$spat[thresh]
  pop_dist_geo <- (pop_dist_geo / sd(pop_dist_geo))
  # genetic (autosomal) distances across pop pairs
  pop_dist_gene <- nuclear_data$gene[thresh]
  pop_dist_gene <- lm(pop_dist_gene ~ pop_dist_geo)$resid
  pop_dist_gene <- (pop_dist_gene / sd(pop_dist_gene))
  # genetic (x chr.) distances across pop pairs
  pop_dist_gene_x <- nuclear_data$gene_x[thresh]
  pop_dist_gene_x <- lm(pop_dist_gene_x ~ pop_dist_gene + pop_dist_geo)$resid
  pop_dist_gene_x <- (pop_dist_gene_x / sd(pop_dist_gene_x))
  
  # whether the population pair is exogamous or endogamous
  ethno_trait <- nuclear_data$endogamous # nuclear_data$exogamous + 2 *
  ethno_trait <- as.factor(ethno_trait[thresh])
  
  # initialize a STAN model
  stan_model <- quap(
    # model structure
    alist(
      # binomial model, predicting the number of shared phonemes based on
      #  possible phonemes (in either language)
      phon_uniq ~ dbinom(phon_total, p),
      logit(p) <- a + spat_coef[ethn] * spat +
        gene_coef[ethn] * gene + gene_x_coef[ethn] * gene_x,
      a ~ dnorm(0, 3),
      spat_coef[ethn] ~ dnorm(0, 0.5),
      gene_coef[ethn] ~ dnorm(0, 0.5),
      gene_x_coef[ethn] ~ dnorm(0, 0.5)
    ),
    # inputs
    data = list(
      phon_uniq = pop_phon[, 1],
      phon_total = pop_phon[, 2],
      spat = pop_dist_geo,
      gene = pop_dist_gene,
      gene_x = pop_dist_gene_x,
      ethn = ethno_trait
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
  return_list <- list(posterior = posterior, resid = stan_resid, thresh = thresh)
  return(return_list)
};

###

# run and plot the nuclear genetic model using pairs of
#  populations from around the world
run_worldwide_nuclear_coef_ethn_mdl <- function() {
  # the geographic/spatial thresholds at which to calculate model values
  model_distances <- (4:40) * 250
  # output of the model for each threshold
  nuclear_coef_out <- lapply(model_distances, function(distance)
    nuclear_coefficient_ethn_fxn(distance))
  
  # take credible intervals of each coef (spatial, gene, and gene_x)
  coefs <- list()
  coefs_ci <- list()
  # spat_coef & gene_coef
  for (cf in c("spat_coef", "gene_coef", "gene_x_coef")) {
    coefs[[cf]] <- sapply(
      nuclear_coef_out,
      function(d) {
        median(d$posterior[[cf]][, 1])
      }
    )
    coefs[[paste0(cf, "_endo")]] <- sapply(
      nuclear_coef_out,
      function(d) {
        median(d$posterior[[cf]][, 2] - d$posterior[[cf]][, 1])
      }
    )
    # coefs[[paste0(cf, "_exo")]] <- sapply(
    #   nuclear_coef_out,
    #   function(d) {
    #     median(d$posterior[[cf]][, 3] - d$posterior[[cf]][, 1])
    #   }
    # )
    coefs_ci[[cf]] <- sapply(
      nuclear_coef_out,
      function(d) {
        if (all(is.na(d$posterior[[cf]][, 1]))) {
          return(as.numeric(c(NA, NA)))
        } else {
          PI(d$posterior[[cf]][, 1])
        }
      }
    )
    coefs_ci[[paste0(cf, "_endo")]] <- sapply(
      nuclear_coef_out,
      function(d) {
        if (all(is.na(d$posterior[[cf]][, 1]))) {
          return(as.numeric(c(NA, NA)))
        } else {
          PI(d$posterior[[cf]][, 2] - d$posterior[[cf]][, 1])
        }})
    # coefs_ci[[paste0(cf, "_exo")]] <- sapply(
    #   nuclear_coef_out,
    #   function(d) {
    #     if (all(is.na(d$posterior[[cf]][, 1]))) {
    #       return(as.numeric(c(NA, NA)))
    #     } else {
    #       PI(d$posterior[[cf]][, 3] - d$posterior[[cf]][, 1])
    #     }})
  }
  
  # data frame for plotting
  coef_plotting_df <- data.frame(
    coef_distance = unlist(lapply(model_distances, rep, 6)), # 9
    Coefficient = rep(c("spat_coef", "spat_coef_endo", # "spat_coef_exo",
                        "gene_coef", "gene_coef_endo", # "gene_coef_exo",
                        "gene_x_coef", "gene_x_coef_endo"), #, "gene_x_coef_exo"),
                      length(coefs$spat_coef)),
    coef_value = as.vector(t(as.matrix(as.data.frame(coefs))[, 1:6])),  # 9
    coef_ci_min = as.vector(t(sapply(coefs_ci, function(d) d[1, ])[, 1:6])),  # 9
    coef_ci_max = as.vector(t(sapply(coefs_ci, function(d) d[2, ])[, 1:6]))  # 9
  )
  
  nuc_ww_plot <- ggplot(data = coef_plotting_df, aes(
    x = coef_distance, y = coef_value, group = Coefficient)) +
    xlab("Threshold for spatial distance (km)") +
    ylab("Relationship to linguistic distances") +
    guides(fill = element_blank(),
           colour = NULL) +
    geom_line(aes(color = Coefficient, linetype = Coefficient)) +
    scale_x_continuous(
      lim = range(model_distances),
      breaks = (1:10) * 1000,
      expand = expansion(mult = c(0, 0))) +
    scale_color_manual(
      breaks = c("spat_coef", "spat_coef_endo", # "spat_coef_exo",
                 "gene_coef", "gene_coef_endo", # "gene_coef_exo",
                 "gene_x_coef", "gene_x_coef_endo"), #, "gene_x_coef_exo"),
      values = c("#000000", "#000000", # "#c0c0c0",
                 "#0000ff", "#0000ff", # "#00ffff",
                 "#ff0000", "#ff0000"),  #"#ff00a0"),
      labels = c(
        "Spatial distance", "Spat. (endo)", # "Spat. (exo)",
        "Autosomal distance", "Auto. (endo)", # "Auto. (exo)",
        "X-chrom. distance", "X-chr. (endo)")) + #, "X-chr. (exo)")) +
    geom_ribbon(aes(
      x = coef_distance,
      ymin = coef_ci_min, ymax = coef_ci_max,
      group = Coefficient, fill = Coefficient)) +
    scale_fill_manual(
      breaks = c("spat_coef", "spat_coef_endo", # "spat_coef_exo",
                 "gene_coef", "gene_coef_endo", # "gene_coef_exo",
                 "gene_x_coef", "gene_x_coef_endo"), #, "gene_x_coef_exo"),
      values = c("#00000030", "#00000030", # "#c0c0c0",
                 "#0000ff30", "#0000ff30", # "#00ffff",
                 "#ff000030", "#ff000030"), # "#ff00a030"),
      labels = c(
        "Spatial distance", "Spat. (endo)", # "Spat. (exo)",
        "Autosomal distance", "Auto. (endo)", # "Auto. (exo)",
        "X-chrom. distance", "X-chr. (endo)" # "X-chr. (exo)"
        )) +
    scale_linetype_manual(
      breaks = c("spat_coef", "spat_coef_endo", # "spat_coef_exo",
                 "gene_coef", "gene_coef_endo", # "gene_coef_exo",
                 "gene_x_coef", "gene_x_coef_endo"), #, "gene_x_coef_exo"),
      values = c("solid", "dashed",
                 "solid", "dashed",
                 "solid", "dashed"), 
      labels = c(
        "Spatial distance", "Spat. (endo)", # "Spat. (exo)",
        "Autosomal distance", "Auto. (endo)", # "Auto. (exo)",
        "X-chrom. distance", "X-chr. (endo)" # "X-chr. (exo)"
      )) +
    geom_hline(yintercept = 0, color = '#80808030') +
    theme_classic()
  
  ggsave("../figures/nuc_plots/thresh_worldwide_endogamy.pdf", nuc_ww_plot,
         width=7, height=5, units='in')
}; run_worldwide_nuclear_coef_ethn_mdl()

###

# function that finds the best fit coefficients for pairs of populations across
#  the world using mtDNA data
mt_coef_cmo_fxn <- function(
    spat_thresh = 1000, continents = NA) {
  
  # which pairs of populations to use
  thresh <- which(mt_data$geo_dist <= spat_thresh &
                  as.dist(mt_data$min_samp) >=3)
  
  # subset continents
  if (!all(is.na(continents))) {
    thresh <- thresh[mt_data$continent[thresh] %in% continents]
  }
  
  # select relevant ethnographic trait
  pop_ethnographic <- as.dist(mt_data$endogamous)[thresh] #  * 2 + mt_data$exogamous
  thresh <- thresh[! is.na(pop_ethnographic)]

  pop_ethnographic <- as.factor(pop_ethnographic[! is.na(pop_ethnographic)])

  # make sure there are at least 7 neighboring populations
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
      phon_shared ~ dbinom(phon_total, p),
      logit(p) <- a +
        spat_coef[ethn] * spat +
        gene_coef[ethn] * gene,
      a ~ dnorm(0, 3),
      spat_coef[ethn] ~ dnorm(0, 0.5),
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

  p_val <- summary(glm(
    # binomial glm 
    cbind(phon_shared, phon_total - phon_shared) ~
      spat * ethn + gene * ethn - ethn,
    family = binomial,
    data = df))$coefficient
  p_val <- p_val[2:nrow(p_val), ]
  coef_val <- p_val[, 1]
  p_val <- p_val[, 4]
  
  return_list <- list(posterior = posterior, resid = stan_resid,
                      thresh = thresh, p_values = p_val, coef_values = coef_val)
  return(return_list)
};

#
run_worldwide_mt_coef_cmo_mdl <- function() {
  
  # get possible regional subsets
  Americas <- (function(x) x[grepl("America", x, )])(unique(mt_data$continent))
  Africa <- (function(x) x[grepl("Africa", x, )])(unique(mt_data$continent))
  Central_East_Asia_and_Pacific <- (function(x) x[
    grepl("Asia|Pacific", x, )])(unique(mt_data$continent))
  West_Eurasia <- (function(x) x[grepl("Europe|Middle", x, )])(
    unique(mt_data$continent))

  plot_num <- 0
  mt_plot_list <- list()
  mt_pval_plot_list <- list()
  for (region_subset in c("All", "Central & East Asia and Pacific",
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
      mt_coefmdl_stan <- lapply(model_distances, function(d) mt_coef_cmo_fxn(
        d, continents=continents))
      
      mt_coefs <- list()
      mt_coefs_ci <- list()
      for (cf in c("spat_coef", "gene_coef")) {
        mt_coefs[[cf]] <- sapply(
          mt_coefmdl_stan,
          function(d) {
            median(d$posterior[[cf]][, 1])
          }
        )
        # mt_coefs[[paste0(cf, "_exo")]] <- sapply(
        #   mt_coefmdl_stan,
        #   function(d) {
        #     median(d$posterior[[cf]][, 2] - d$posterior[[cf]][, 1])
        #   }
        # )
        mt_coefs[[paste0(cf, "_endo")]] <- sapply(
          mt_coefmdl_stan,
          function(d) {
            median(d$posterior[[cf]][, 2] - d$posterior[[cf]][, 1])
          }
        )
        mt_coefs_ci[[cf]] <- sapply(
          mt_coefmdl_stan,
          function(d) {
            if (all(is.na(d$posterior[[cf]][, 1]))) {
              return(as.numeric(c(NA, NA)))
            } else {
              PI(d$posterior[[cf]][, 1])
            }
          }
        )
        # mt_coefs_ci[[paste0(cf, "_exo")]] <- sapply(
        #   mt_coefmdl_stan,
        #   function(d) {
        #     if (all(is.na(d$posterior[[cf]][, 1]))) {
        #       return(as.numeric(c(NA, NA)))
        #     } else {
        #       PI(d$posterior[[cf]][, 2] - d$posterior[[cf]][, 1])
        #     }})
        mt_coefs_ci[[paste0(cf, "_endo")]] <- sapply(
          mt_coefmdl_stan,
          function(d) {
            if (all(is.na(d$posterior[[cf]][, 1]))) {
              return(as.numeric(c(NA, NA)))
            } else {
              PI(d$posterior[[cf]][, 2] - d$posterior[[cf]][, 1])
            }})
      }
      
      p_val_table <- sapply(mt_coefmdl_stan, function(d) {
        n <- rep(NA, 4)
        n[1:length(d$p_values)] <- d$p_values
        n
      })[c(1:4), ]
      coef_val_table <- sapply(mt_coefmdl_stan, function(d) {
        n <- rep(NA, 4)
        n[1:length(d$coef_values)] <- d$coef_values
        n
      })[c(1:4), ]
      mt_coef_plotting_df <- data.frame(
        # spatial thresholds
        coef_distance = unlist(lapply(model_distances, rep, 4)),  # 6
        # coefficient names
        Coefficient = rep(c("spat_coef", "spat_coef_endo",
                            # "spat_coef_exo", "gene_coef_exo",
                            "gene_coef", "gene_coef_endo"),
                          length(mt_coefs$spat_coef)),
        # coefficient values
        coef_value = as.vector(t(as.matrix(as.data.frame(mt_coefs))[, 1:4])),  # 6
        # coefficient credibility intervals
        coef_ci_min = as.vector(t(sapply(mt_coefs_ci, function(d) d[1, ])[, 1:4])),
        coef_ci_max = as.vector(t(sapply(mt_coefs_ci, function(d) d[2, ])[, 1:4])),
        p_vals = as.vector(p_val_table),
        coef_vals = as.vector(coef_val_table)
      )
      
      ##
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
          expand = expansion(mult = c(0, 0))) +
        # mt_plot_list[[plot_num]] 
        scale_linetype_manual(
          values = c("solid", "solid",
                     # "dashed", "dashed",
                     "dashed", "dashed"),
          breaks = c("spat_coef", "gene_coef",
                     # "spat_coef_exo", "gene_coef_exo",
                     "spat_coef_endo", "gene_coef_endo"),
          labels = c(
            "Spatial distance",
            "mtDNA distance\n(non-spatial component)",
            paste0("Spatial distance\n(add'l effect of\nendogamy)"),
            paste0("mtDNA distance\n(add'l effect of\nendogamy)"))) +
        scale_color_manual(
            values = c("black", "red",
                       # "#000000", "#ff0000",
                       "#000000", "#ff0000"),
            breaks = c("spat_coef", "gene_coef",
                       # "spat_coef_exo", "gene_coef_exo",
                       "spat_coef_endo", "gene_coef_endo"),
            labels = c(
              "Spatial distance",
              "mtDNA distance\n(non-spatial component)",
              paste0("Spatial distance\n(add'l effect of\nendogamy)"),
              paste0("mtDNA distance\n(add'l effect of\nendogamy)"))) +
          geom_ribbon(aes(
            x = coef_distance,
            ymin = coef_ci_min, ymax = coef_ci_max,
            group = Coefficient, fill = Coefficient)) +
          scale_fill_manual(
            values = c('#00000030', '#e0000030',
                       '#00000030', '#e0000030'),
                       # "#b0b0b030", "#ff800030"),
            breaks = c("spat_coef", "gene_coef",
                       # "spat_coef_exo", "gene_coef_exo",
                       "spat_coef_endo", "gene_coef_endo"),
            labels = c(
              "Spatial distance",
              "mtDNA distance\n(non-spatial component)",
              paste0("Spatial distance\n(add'l effect of\nendogamy)"),
              paste0("mtDNA distance\n(add'l effect of\nendogamy)"))) +
          geom_hline(yintercept = 0, color = '#80808030') +
          theme_classic()
      
      # mt_pval_plot_list[[plot_num]] <- ggplot(data = mt_coef_plotting_df, aes(
      #   x = coef_distance, y = p_vals, group = Coefficient)) +
      #   xlab("Spatial distance threshold for pairs of populations (km)") +
      #   ylab(paste0("p-value \n(", continents_name, ")")) +
      #   guides(fill = element_blank(),
      #          colour = NULL) +
      #   geom_line(aes(color = Coefficient, linetype = Coefficient)) +
      #   scale_x_continuous(
      #     lim = range(model_distances),
      #     breaks = 2 * (0:5) * 1000,
      #     expand = expansion(mult = c(0, 0))) +
      #   scale_y_continuous(
      #     transform = "log",
      #     breaks = c(10 ** c(-3, -5, (-1 - 10 * (0:20))))
      #   ) + ggforce::facet_zoom(ylim = c(0.2, 0.00001),
      #                           zoom.size = 1,
      #                           horizontal = F) +
      #   scale_linetype_manual(
      #     values = c("solid", "solid", "dashed", "dashed"),
      #     breaks = c("spat_coef", "gene_coef", "spat_coef_endo", "gene_coef_endo"),
      #     labels = c(
      #       "Spatial distance",
      #       "mtDNA distance\n(non-spatial component)",
      #       paste0("Spatial distance\n(add'l effect of\nendogamy)"),
      #       paste0("mtDNA distance\n(add'l effect of\nendogamy)"))) +
      #   scale_color_manual(
      #     values = c("black", "red", "black", "red"),
      #     breaks = c("spat_coef", "gene_coef", "spat_coef_endo", "gene_coef_endo"),
      #     labels = c(
      #       "Spatial distance",
      #       "mtDNA distance\n(non-spatial component)",
      #       paste0("Spatial distance\n(add'l effect of\nendogamy)"),
      #       paste0("mtDNA distance\n(add'l effect of\nendogamy)"))) +
      #   geom_hline(yintercept = 0.05, color = '#80808030') +
      #   theme_classic()
    }

  
  plot_num <- 0;
  for (region_subset in c("world" , "ceasiapacific",
                          "americas", "africa", "weurasia"
                          )) {
    plot_num <- plot_num + 1
    ggsave(filename = paste0("../figures/plot_mt_", region_subset, "_exoendo.pdf"),
           plot = mt_plot_list[[plot_num]], device = "pdf", dpi = 300,
           height = 4, width = 8)
    # ggsave(filename = paste0("../figures/plot_mt_", region_subset, "_exoendo_coef.pdf"),
    #        plot = mt_pval_plot_list[[plot_num]], device = "pdf", dpi = 300,
    #        height = 4, width = 8)
  }
}; run_worldwide_mt_coef_cmo_mdl()
