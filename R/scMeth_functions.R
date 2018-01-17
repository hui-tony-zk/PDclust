cores_to_use <- 30

#' Divide cores
#'
#' @description Utility function to help parallelize jobs in R with less overhead.
#'     Likely wouldn't ever be used by the user.
#' @param total Total number of jobs to parallelize
#' @param ncores CPUs to parallelize over (must be >1)
#'
#' @return Returns a data.frame with two columns (from and to) and rows equal to the number of cores
#'
#' @examples divide_cores(total = 100, ncores = 2)
#'
divide_cores <- function(total, ncores = cores_to_use) {
  if (cores_to_use <= 1) stop("ncores <=1 - make sure to use 2 or more cores for parallelization")
  if (total <= ncores){
    return(data.frame(from=1:(total),to=1:(total)))
  } else if (total <= ncores*2) {
    from=c(1:(total-1))
    from=c(from)[seq(from = 1,to=length(from), by=2)]
    to=c(1:total)
    to=c(to)[seq(from = 2,to=length(to), by=2)]
    df = data.frame(from, to)
    if (total %% 2 == 1) {
      df[nrow(df),2]<- total
    }
    return(df)
  } else {
    batch = total %/% ncores
    from = seq(from = 1, to = total-batch, by = batch)[1:ncores-1]
    to = c(seq(from = batch, to = total, by = batch), total)[1:length(from)]
    total_batch <- data.frame(cbind(from, to))
    # Add on the remainder
    total_batch <- rbind(total_batch, data.frame(
      from=max(total_batch$to)+1,
      to=total)
    )
  }
  return(total_batch)
}

#' Create pairwise comparisons
#'
#' @description Create pairwise comparisons between single-cells
#' @param cpg A list of named data frames containing CpG calls. See details for required format of dataframes. Required.
#' @param digital Whether or not to discard non-binary CpG calls. Useful in single-cells as it's very unlikely that a single-cell contains a heterozygous methylation call. Defaults to TRUE.
#' @param ncores Number of cores to parallelize over. Defaults to 1
#' @param calcdiff Whether or not to directly calculate the average difference (if TRUE), or to return a list of dataframes containing pairwise common CpGs (if FALSE). Defaults to TRUE. Usually you don't want false unless you wish to do something else with all the pairwise data
#'
#' @return A list of dataframes if \code{calcdiff} is FALSE. Otherwise, a dataframe containing the pairwise dissimiarlties if \code{calcdiff} is TRUE.
#' @export
#' @import foreach doMC dplyr
#' @importFrom stats setNames
#'
#' @details
#' Each dataframe containing CpG calls must have the following four columns:
#' 1. Chromsome column, named "chr"
#' 2. Start/Position column, named "start"
#' 3. Percentage or fractional methylation column, named "meth" (between 0-100 or 0-1)
#'
create_pairwise_master <- function(cpg, digital = TRUE, cores_to_use = 2, calcdiff = TRUE){
  doMC::registerDoMC(cores_to_use)
  # Generate combinations
  comb.names <- utils::combn(names(cpg),2)
  total_batch <- divide_cores(ncol(comb.names), cores_to_use)
  # Combine
  pairwise <- foreach(i=1:nrow(total_batch), .combine = c) %dopar% {
    start <- total_batch[i,1]
    end <- total_batch[i,2]
    merge_bind=vector("list", end-start+1)
    for (f in start:end) {
      name1 <- comb.names[1,f]
      name2 <- comb.names[2,f]
      one <- cpg[[name1]] %>% select(1,2,4)
      two <- cpg[[name2]] %>% select(1,2,4)
      if (digital) {
        one <- filter(one, meth == 0 | meth == max(meth))
        two <- filter(two, meth == 0 | meth == max(meth))
      }
      merge_temp <- inner_join(one, two, by=c("chr", "start")) %>% setNames(c("chr","pos",paste0(name1), paste0(name2)))
      if (calcdiff) {
        merge_temp <- get_diff_df(merge_temp)
      }
      merge_bind[[f-start+1]] <- merge_temp
      names(merge_bind)[f-start+1] <- paste0(name1,"_",name2)
    }
    merge_bind
  }
  if (calcdiff) {
    pairwise <- suppressWarnings(bind_rows(pairwise))
    return(tbl_df(pairwise))
  } else {
    return(pairwise)
  }
}

#' Get Dissimilarity from a pairwise common data frame
#'
#' @param df A dataframe containing the pairwise common CpGs. Required.
#'
#' @return A 1-row dataframe containing two samples compared, the total number of pairwise CpGs, and the pairwise dissimilarity
#' @import dplyr
#' @importFrom stats setNames
#'
#' @details
#' Typically this function is not used by the user. It is called by \code{create_pairwise_master} when \code{calc_diff} is TRUE.
#'
#' The input dataframe requires four columns in the specific order:
#' 1. Chromosome
#' 2. Start
#' 3. Name of the 1st comparitor (e.g. cell1)
#' 4. Name of the 2nd comparitor (e.g. cell2)
#'

get_diff_df <- function(df) {
  diff.temp <- data.frame(
    x = names(df)[3],
    y = names(df)[4],
    total = nrow(df),
    pear_corr = c(cor(df[[3]], df[[4]], method = "pearson")),
    pairwise_dissimilarity_total = sum(abs(df[[3]]-df[[4]]))
  ) %>%
    mutate(pairwise_dissimilarity = pairwise_dissimilarity_total/total)
  diff.temp <- diff.temp %>%
    select(x, y, num_cpgs = total, pairwise_dissimilarity)
  return(diff.temp)
}

#' Convert pairwise dissimilarity data.frame to a dissimilarity matrix
#'
#' @param master_diff The data.frame containing pairwise dissimilarites. Typically the output of \code{create_pairwise_master}. Required.
#' @param measure The name of the column to use as the values for the pairwise dissimilarity matrix. Defaults to the pairwise_dissimilarity. Can be changed if user wants another measure (e.g. correlation, etc).
#' @param diag The value to put in the diagonal. Defaults to NA.
#' @param sample_subset A character vector of sample names to subset. Optional.
#'
#' @return A data.matrix representing the pairwise dissimilarites between every cell
#' @export
#' @import dplyr
#'
convert_to_dissimilarity_matrix <- function(master_diff, measure = "pairwise_dissimilarity", diag = NA, sample_subset = NULL) {
  if (! is.null(sample_subset)) {
    master_diff <- master_diff %>% filter(x %in% sample_subset, y %in% sample_subset)
  }
  test <- data.frame(master_diff)[,c("x","y",measure)]
  colnames(test)[3] <- "value"
  list_samples <- unique(c(as.character(test$x), as.character(test$y)))
  test.rev <- test %>% transform(x=y,y=x,value=value)
  test.same <- data.frame(x=list_samples,y=list_samples,value=diag)
  test <- unique(rbind(test.same, test, test.rev))

  x <- data.matrix(tidyr::spread(test, key = y, value = value)[,-1])
  rownames(x) <- colnames(x)
  return(x)
}

#' Clustering dissimilarities
#'
#' @param dissimilarity_matrix A matrix outputted by `convert_to_dissimilarity_matrix()`
#' @param num_clusters Number of clusters to divide into
#'
#' @return an list containing: 1) an hclust object and 2) cluster assignments of samples
#' @export
#' @import dplyr
#' @importFrom stats setNames
#'
cluster_dissimilarity <- function(dissimilarity_matrix, num_clusters) {
  if (nrow(dissimilarity_matrix) < num_clusters) stop("You have selected more clusters than samples!")
  hclust_obj <- hclust(dist(dissimilarity_matrix, method = "euclidean"), method = "ward.D2")
  cluster_assignments <- cutree(hclust_obj, k = num_clusters) %>% as.data.frame() %>% setNames("cluster")
  cluster_assignments$cluster <- as.factor(cluster_assignments$cluster)
  return(list(hclust_obj = hclust_obj, cluster_assignments = cluster_assignments))
}

#' Visualize clusters
#'
#' @param dissimilarity_matrix A matrix outputted by `convert_to_dissimilarity_matrix()`
#' @param cluster_results The result from `cluster_dissimilarity()`. If left empty, no clustering results will be visualized.
#'
#' @return A ggplot-ready data.frame
#' @export
#'
visualize_clusters <- function(dissimilarity_matrix, cluster_labels) {
  x <- as.dist(dissimilarity_matrix) %>%
    cmdscale() %>% as.data.frame()
  if (!missing(cluster_labels)) {
    x <- x %>%
      merge(cluster_labels, by="row.names", all.x=TRUE)
  } else {
    warning("Warning: no cluster_labels specified!")
  }
  return(tbl_df(x))
}

#' Create in silico merged bulk profiles from single-cell files
#'
#' @param cpg_all The list containing all CpG calls in data.frame format. Required.
#' @param cluster_assignments The cluster_assignments outputed from cluster_dissimilarity(). Required
#' @param desired_cluster The number desired from the input cluster assignments. Required
#'
#' @return a BSmooth object ready for smoothing
#' @import dplyr
#' @importFrom stats setNames
#' @export
#'
#' @details Uses the bsseq package to perform in silico merging of single-cell CpG calls. Requires the bsseq R package to be installed
#'
merge_cpgs <- function(cpg_all, cluster_assignments, desired_cluster) {
  if (!require("bsseq")) {
    stop("'bsseq' package needed for this function to work. Please install it at http://bioconductor.org/packages/release/bioc/html/bsseq.html.",
         call. = FALSE)
  }
  # get group members
  cluster_members <- cluster_assignments %>%
    subset(cluster == desired_cluster) %>%
    row.names()
  # get group CpGs
  tmp <- cpg_all[cluster_members]
  tmp <- lapply(seq_along(tmp), function(i) {
    x = tmp[[i]]
    y = BSseq(M = as.matrix(x$meth/100), Cov = as.matrix(rep(1, nrow(x))), chr = x$chr, pos = x$start, sampleNames = names(tmp)[[i]])
    return(y)
  })
  tmp2 = tmp[[1]]
  for (i in 2:length(tmp)) {
    tmp2 = combineList(list(tmp2, tmp[[i]]))
    samples = sampleNames(tmp2)
    tmp2 = collapseBSseq(tmp2, columns = rep(paste0("group", desired_cluster), length(samples)) %>% setNames(samples))
    gc()
  }
  return(tmp2)
}

"cpg_files"
