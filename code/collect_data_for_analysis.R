# generate a list with all of the info necessary to analyze nuclear genetic data
nuclear_data_fxn <- function() {
  
  # get metadata for pops with autosomal data
  metadata <- read.csv("../data/nuc_populations_withDPLACEmatches.csv")
  metadata$ID <- sub(" $", "", metadata$ID)  # remove trailing white-space
  
  # get pops with more 2 or more individuals were sampled
  pops_sampled <- unlist(read.table("../data/nuclear_pops_2orMoreIndividuals.txt"))
  metadata <- metadata[(metadata$ID %in% pops_sampled), ]
  
  # get the linguistic Ruhlen IDs associated with each population
  ruhlen_nums <- lapply(
    pops_sampled,
    function(pop) as.numeric(unlist(strsplit(
      metadata$ruhlen_id[(metadata$ID == pop)], ";")))
  )
  
  # geographic distances
  ruhlen_geo_dist_matrix <- read.table(
    "../data/GeographicDistanceRuhlen2082.csv",
    row.names=1, header=T, sep=",")
  
  get_dist <- function(pop1_ids, pop2_ids, dist_matrix) {
    d <- c()
    for (i in pop1_ids) {
      for (j in pop2_ids) {
        i <- which(rownames(dist_matrix) == as.character(i))
        j <- which(rownames(dist_matrix) == as.character(j))
        if (!identical(i, integer(0)) && !identical(j, integer(0))) {
          d <- c(d, dist_matrix[i, j])
        }
      }
    }
    return(mean(d))
  }
  
  geo_dist_matrix <- matrix(0, length(ruhlen_nums), length(ruhlen_nums))
  for (i in 1:(length(ruhlen_nums) - 1)) {
    for (j in (i + 1):length(ruhlen_nums)) {
      # print(c(ruhlen_nums[i], ruhlen_nums[j]))
      geo_dist_matrix[i, j] <- geo_dist_matrix[j, i] <- get_dist(
        ruhlen_nums[[i]], ruhlen_nums[[j]], ruhlen_geo_dist_matrix)
    }
  }
  colnames(geo_dist_matrix) <- rownames(geo_dist_matrix) <-
    pops_sampled
  rm(ruhlen_geo_dist_matrix)
  
  ## get ruhlen data
  
  ruhlen_phonemes = read.table("../data/merritt_phoneme_matrix.txt",
                               comment.char = "#", sep="\t", row.names = 1,
                               header = F )
  phoneme_data = read.table("../data/merritt_phoneme_info.txt",
                            comment.char = "#", sep="\t", row.names = 1,
                            header = T)
  
  pop_continent <- sapply(
    ruhlen_nums,
    function(r) {
      if (length(r) > 1) r <- r[1];
      return(ruhlen_phonemes[which(rownames(ruhlen_phonemes) == r), 6])
    }
  )
  
  phon = matrix(0, length(ruhlen_nums), ncol(ruhlen_phonemes) - 8)
  rownames(phon) <- pops_sampled
  colnames(phon) <- ruhlen_phonemes[1, 9:ncol(ruhlen_phonemes)]
  for (i in 1:length(ruhlen_nums)) {
    nums <- ruhlen_nums[[i]]
    for (j in 1:ncol(phon)) {
      phon[i, j] <- mean(as.numeric(
        ruhlen_phonemes[rownames(ruhlen_phonemes) %in% as.character(nums), 8 + j],
        na.rm = T))
    }
  }
  
  auto_fst <- 0 * as.matrix(geo_dist_matrix)
  X_fst <- auto_fst
  auto_fst_table <- read.table("../data/nuclear.fst.summary")
  X_fst_table <- read.table("../data/nuclear.x.fst.summary")
  
  for (i in rownames(geo_dist_matrix)) {
    for (j in rownames(geo_dist_matrix)) {
      if (i != j) {
        pos <- which(auto_fst_table[, 1] == i & auto_fst_table[, 2] == j)
        if (length(pos) == 1) {
          auto_fst[i, j] <- auto_fst[j, i] <- auto_fst_table[pos, 3]
          X_fst[i, j] <- X_fst[j, i] <- X_fst_table[pos, 3]
        }
      }
    }
  }
  
  rm(auto_fst_table, X_fst_table)
  
  # for each pair of populations, find the number of phonemes that are present
  #  in only one language (first column) and the number of phonemes that are
  #  present in either of the two languages (second column)
  pop_phon <- cbind(
    as.dist(round(sapply(
      1:nrow(phon),
      function(n) sapply(
        1:nrow(phon),
        function(m) sum(round(abs(phon[n, ] - phon[m, ])))
      )))),
    as.dist(round(sapply(
      1:nrow(phon),
      function(n) sapply(
        1:nrow(phon),
        function(m) sum(round(phon[n, ]) | round(phon[m, ]))
      ))))
  )
  
  gene_x <- as.dist(X_fst)
  gene <- as.dist(auto_fst)
  spat <- as.dist(geo_dist_matrix)

  continent <- as.matrix(sapply(
    pop_continent, function(i)
      sapply(pop_continent, function(j) {
        if (i == j) {
          return(i)
        } else {
          return(paste0(i, "_", j))
        }
      })))
  continent <- continent[lower.tri(continent)]
  
  pop_exogamous <- metadata[, paste0(
    "Code..EA015.Community.marriage.organization..1.Exogamous.",
    ".2.Agamous..3.Endogamous.")] == 1
  pop_endogamous <- metadata[, paste0(
    "Code..EA015.Community.marriage.organization..1.Exogamous.",
    ".2.Agamous..3.Endogamous.")] == 3
  exogamous <- endogamous <- nrow(metadata) %>% matrix(NA, ., .)
  for (i in 1:(length(pop_exogamous)) - 1) {
    for (j in (i + 1):length(pop_exogamous))
    {
      endogamous[i, j] <- endogamous[j, i] <-
        max(pop_endogamous[i], pop_endogamous[j])
      exogamous[i, j] <- exogamous[j, i] <- 
        max(pop_exogamous[i], pop_exogamous[j]) * (1 - endogamous[i, j])
    }
  }
  
  dat_list <- list(
    meta = metadata,
    pop_order = pops_sampled,
    spat = spat,
    gene = gene,
    gene_x = gene_x,
    pop_phon = pop_phon,
    endogamous = endogamous,
    exogamous = exogamous,
    continent = continent)
  
  dat_list
}; nuc_list <- nuclear_data_fxn()

# mtDNA data processing
make_mt_dists <- function() {
  # load FASTA data
  fasta <- ape::read.FASTA("../data/mtDNA_5141.fasta")
  fasta_pops <- read.table("../data/mtDNApopNames_5141.txt",
                           header = T, sep = "\t")$"PopNames5141"
  # filter for populations 
  pops_used <- rle(as.character(sort(fasta_pops))) %>% { .$values[.$lengths > 2] }
  use_fasta_pops <- which(
    as.character(fasta_pops) %in% pops_used[
      pops_used %in% read.csv("../data/mt_populations_withDPLACEmatches.csv",
                              header = T)$Pop_Names_Ethno])
  
  fasta_dist <- ape::dist.dna(fasta[use_fasta_pops])
  
  fasta_pops <- fasta_pops[use_fasta_pops]
  
  mt_dists <- matrix(0, nrow = length(pops_used), ncol = length(pops_used))
  rownames(mt_dists) <- colnames(mt_dists) <- pops_used
  for (i in 2:ncol(mt_dists)) {
    pops_i <- which(as.character(fasta_pops) == colnames(mt_dists)[i])
    for (j in 1:(i - 1)) {
      pops_j <- which(as.character(fasta_pops) == colnames(mt_dists)[j])
      
      mt_dists[i, j] <- mt_dists[j, i] <- mean(fasta_dist[
        pops_i, pops_j
      ])
    }
  }
  
  write.table(x = mt_dists, file = "../big_mt/mtDNA_pop_dists.txt")
  
}

mt_data_fxn <- function() {
  # get population data
  mt_meta <- read.csv("../data/mt_populations_withDPLACEmatches.csv",
                      header = T)
  mt_meta <- mt_meta[order(mt_meta$Pop_Names_Ethno), ]
  mt_meta[mt_meta == ""] <- NA
  
  mt_meta <- mt_meta[mt_meta$Num_Samples >= 3, ]
  
  # load mtDNA distances
  mt_dist <- read.table(
    "../data/mtDNA_pop_dists.txt")
  rownames(mt_dist) <- colnames(mt_dist) <- gsub(" $", "", rownames(mt_dist))
  
  # remove and reorder
  mt_dist <- rownames(mt_dist) %in% mt_meta$Pop_Names_Ethno %>% {
    mt_dist[., .]}
  mt_dist <- order(rownames(mt_dist)) %>% {mt_dist[., .]} %>% as.dist
  
  # get phonemes for these populations
  ruhlen_phonemes <- read.table(
    "../data/merritt_phoneme_matrix.txt",
    comment.char = "#", sep="\t", header = T)
  ruhlen_rows <- sapply(
    mt_meta$Ruhlen,
    function(n) which(ruhlen_phonemes$Merritt_NUMBER == n))
  ruhlen_phonemes <- ruhlen_phonemes[
    ruhlen_rows,
  ]
  
  # get data for each population from the ruhlen database
  pop_continent <- ruhlen_phonemes$GeographicRegion
  pop_longlat <- ruhlen_phonemes[, c("Longitude", "Latitude")]
  pop_matriliny <- 1 * (mt_meta$Descent == "Matrilineal")
  pop_matrilocality <- 1 * (mt_meta$Residence == "FemaleKin Residence")
  pop_exogamous <- mt_meta[, paste0(
    "Code..EA015.Community.marriage.organization..1.Exogamous.",
    ".2.Agamous..3.Endogamous.")] == 1
  pop_endogamous <- mt_meta[, paste0(
    "Code..EA015.Community.marriage.organization..1.Exogamous.",
    ".2.Agamous..3.Endogamous.")] == 3
  
  phon_mt <- ruhlen_phonemes[, 10:ncol(ruhlen_phonemes)]
  phon_mt <- phon_mt[, colSums(phon_mt) != 0]
  
  pairwise_nonshared <- matrix(0, nrow(phon_mt), nrow(phon_mt))
  pairwise_total <- matrix(0, nrow(phon_mt), nrow(phon_mt))
  
  geo_dist <- as.dist(read.table(
    "../data/GeographicDistanceRuhlen2082.csv",
    row.names=1, header=T, sep=",")[
      ruhlen_rows %>% sapply(function(d) ifelse(d == 1875, 1909, d)),
      ruhlen_rows %>% sapply(function(d) ifelse(d == 1875, 1909, d))])
  
  continent <- as.matrix(sapply(
    pop_continent, function(i) {
      conts <- (pop_continent == i)
      diff_cont <- which(! conts)
      conts[conts] <- i
      conts[diff_cont] <- paste0(i, "_", pop_continent[diff_cont])
      return(conts)
    }))
  continent <- continent[lower.tri(continent)]
  
  # for each population pair, identify whether the populations have female-
  #   -baised descent and female-biased postmarital residence
  matrilineal <- sapply(
    pop_matriliny, function(i)
      sapply(pop_matriliny, function(j) {
        if (is.na(i) | is.na(j)) {
          return(NA)
        } else {
          return(mean(i, j))
        }
      }))
  matrilocal <- sapply(
    pop_matrilocality, function(i)
      sapply(pop_matrilocality, function(j) {
        if (is.na(i) | is.na(j)) {
          return(NA)
        } else {
          return(mean(i, j))
        }
      }))
  
  exogamous <- endogamous <- matrilineal
  for (i in 1:(length(pop_exogamous)) - 1) {
    for (j in (i + 1):length(pop_exogamous))
    {
      endogamous[i, j] <- endogamous[j, i] <-
        max(pop_endogamous[i], pop_endogamous[j])
      exogamous[i, j] <- exogamous[j, i] <- 
        max(pop_exogamous[i], pop_exogamous[j]) * (1 - endogamous[i, j])
    }
  }
  
  # find the minimum number of samples for each population pair (for quality
  #  control)
  min_samp <- sapply(
    mt_meta$Num_Samples, function(i)
      sapply(mt_meta$Num_Samples, function(j) {
        min(i, j)
      }))
  
  # for each population pair, find the number of syllables unique to either
  #  population's language, and the number of syllables found in either
  #  language.
  for (n in 2:nrow(phon_mt)) {
    for (m in 1:(n - 1)) {
      pairwise_nonshared[m, n] <- pairwise_nonshared[n, m] <- 
        sum(round(abs(phon_mt[n, ] - phon_mt[m, ])))
      pairwise_total[m, n] <- pairwise_total[n, m] <- 
        sum(round(phon_mt[n, ]) | round(phon_mt[m, ]))
    }
  }
  
  mt_pairwise_lang <- cbind(
    as.dist(pairwise_nonshared),
    as.dist(pairwise_total))
  
  # return a list of data associated with the populations and their 
  #  characteristics for use in analysis and plotting
  mt_data <- list(
    meta = mt_meta,
    pop_longlat = pop_longlat,
    pop_continent = pop_continent,
    continent = continent,
    geo_dist = geo_dist,
    mt_dist = mt_dist,
    mt_phon = mt_pairwise_lang,
    matrilineal = matrilineal,
    matrilocal = matrilocal,
    exogamous = exogamous,
    endogamous = endogamous,
    min_samp = min_samp
  )
  
  mt_data
}; mt_data <- mt_data_fxn()
