library(tidyverse)
library(SeqArray)
library(SeqVarTools)
library(here)

# read in functions from other files
source(here("R/pivot_snpeff.R"))

# Read in data. the data is stored in the local project directory data

dk8_snpeff_genes = read_tsv(
  here("data/snpeff/DK8_dusted.recode_sorted_sorted_snpEff.genes.txt"),
  skip=1
)

dk15_snpeff_genes = read_tsv(
  here('data/snpeff/DK15_dusted.recode_sorted_sorted_snpEff.genes.txt'),
  skip=1
)


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

# pivot the snpeff data to long format where the variant types are stored in
# column variant_type and number of those variants with a given
# gene/transcript/biotype is stored in count
dk8_snpeff_genes_long = pivot_snpeff(dk8_snpeff_genes, "DK8")
dk15_snpeff_genes_long = pivot_snpeff(dk15_snpeff_genes, "DK15")

# full join the data, replace nas appropriately
dk8_dk15_compare = dk8_snpeff_genes_long %>%
  full_join(dk15_snpeff_genes_long,
            by=c('GeneId', 'TranscriptId', 'BioType', 'variant_type')) %>%
  replace_na(list(
    count.x = 0,
    sample.x = "DK8",
    count.y = 0,
    sample.y = "DK15"
  )) %>%
  dplyr::rename(
    dk8_count = count.x,
    dk15_count = count.y
  ) %>%
  dplyr::select(-c('sample.x', 'sample.y'))

dk8_dk15_compare %>%
  filter(dk8_count != 0, dk15_count != 0)

dk8_dk15_compare %>%
  filter(dk8_count == 0, dk15_count != 0)

dk8_not_dk15 = dk8_dk15_compare %>%
  filter(dk8_count != 0, dk15_count == 0)

SeqArray::seqVCF2GDS(
  here("data/snpeff/DK8_dusted.recode_sorted_sorted_snpEff.ann.vcf"),
  here("data/snpeff/DK8_dusted.recode_sorted_sorted_snpEff.ann.gds")
)

dk8_gds = SeqArray::seqOpen(
  here("data/snpeff/DK8_dusted.recode_sorted_sorted_snpEff.ann.gds"))

seqSetFilter(dk8_gds,
             variant.sel = cumsum(x$length) ==
               which(str_detect(x$data, 'CKF44_00476')))

seqResetFilter(dk8_gds)



