# snpeff_df is one of the dataframes fread in from snpeff genes.txt
# sample_name will be appended in a column called "sample"
pivot_snpeff = function(snpeff_df, sample_name){
  # snpeff_df is input
  snpeff_df %>%
    # this drops the first column, which starts with a # and is redundant with
    # GeneId
    dplyr::select(-`#GeneName`) %>%
    pivot_longer(-c(GeneId, TranscriptId, BioType),
                 names_to = 'variant_type',
                 values_to = 'count') %>%
    mutate(sample = sample_name)
}
