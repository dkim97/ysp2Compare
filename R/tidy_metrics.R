#' @title Tidy data from gds file
#' @description Tidy data to long format with column sample, variant, and the
#'   metric name. internal use only
#'
#' @param data_mat a data matrix in sample x variant format
#' @param values_cname name of the metric -- eg RealDepth
#' @param rnames vector to name the columns, eg sample.ids
#' @param cnames vector to name the rows, eg variant_ids
#'
#' @return a tidy dataframe in long format with columns sample, variant, metric
#'
#' @importFrom dplyr as_tibble
#' @importFrom tidyr pivot_longer
tidy_metrics = function(data_mat, values_cname, rnames, cnames){

  errors = list(
    dimension_error = paste0("vcf_to_qtlseqr_table error: rownames do not ",
                             "match data matrix. see tidy_metrics() function")
  )

  if(length(rnames) != dim(data_mat)[1] | length(cnames) != dim(data_mat)[2]){
    stop(errors$dimension_error)
  }

  data_mat %>%
    `rownames<-`(rnames) %>%
    `colnames<-`(cnames) %>%
    dplyr::as_tibble(rownames = "sample") %>%
    tidyr::pivot_longer(-sample, names_to = "variant", values_to = values_cname)
}
