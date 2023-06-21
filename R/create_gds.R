library(tidyverse)
#a metapackage - makes data easier to handle (e.g. reading in CSV files)
library(SeqArray)
library(SeqVarTools)
#tools for reading vcf files
library(here)
#a convenience package - makes pathfinding easier

dir.create(here("data/cnvpytor"))
dir.create(here("data/snpeff"))

#read in data
#data is stored in the local project directory data, but more can be found
#in Z:\Active\lab_members\data_from_brentlab\variant_pipeline_results\run_6453
dk8_snpeff_genes = read_tsv(
  here("data/snpeff/DK8_dusted.recode_sorted_sorted_snpEff.genes.txt"),
  skip=1)
dk15_snpeff_genes = read_tsv(
  here("data/snpeff/DK15_dusted.recode_sorted_sorted_snpEff.genes.txt"),
  skip=1)

#snpeff_df is one of the dataframes read in from snpeff_genes.txt
#sample_name will be appended in a colum called "sample"
pivot_snpeff = function(snpeff_df, sample_name){
#  tmp_df = dplyr::select(snpeff_df, -"#GeneName")
#  long_df = pivot_longer(tmp_df, ...)
#  return(long_df)
#the above three lines are essentially what the below five lines do
  snpeff_df %>%
    #drops the first column, which starts with a # and is redundant with
    #GeneId
    dplyr::select(-"#GeneName") %>%
    pivot_longer(-c(GeneId, TranscriptId, BioType),
                 names_to = "variant_type",
                 values_to = "count") %>%
    mutate(sample = sample_name)
#adding a column
}
#the below code is defined by the above function
#dk8_snpeff_genes %>%
#  dplyr::select(-"#GeneName") %>%
#  pivot_longer(-c(GeneId, TranscriptId, BioType),
#               names_to = "variant_type",
#               values_to = "count") %>% view()
#it is easier bioinformatically to work in long format instead of wide formats
#first argument is always a dataframe - in this case, select
#the "pipe" (%>%) directs the argument towards

dk8_snpeff_genes_long = pivot_snpeff(dk8_snpeff_genes, "DK8")
dk15_snpeff_genes_long = pivot_snpeff(dk15_snpeff_genes, "DK15")

dk8_dk15_comp = dk8_snpeff_genes_long %>%
  full_join(dk15_snpeff_genes_long,
            by = c("GeneId", "TranscriptId", "BioType", "variant_type")) %>%
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
  dplyr::select(-c("sample.x", "sample.y"))
#joins the two tables by the defined columns

dk8_dk15_comp %>%
  filter(dk8_count != 0, dk15_count != 0)

dk8_not_dk15 = dk8_dk15_comp %>%
  filter(dk8_count == 0, dk15_count != 0)

not_dk8_dk15 = dk8_dk15_comp %>%
  filter(dk8_count != 0, dk15_count == 0)

SeqArray::seqVCF2GDS(
  here("data/snpeff/DK8_dusted.recode_sorted_sorted_snpEff.ann.vcf"),
  here("data/snpeff/DK8_dusted.recode_sorted_sorted_snpEff.ann.gds")
)

dk8_gds = SeqArray::seqOpen(
  here("data/snpeff/DK8_dusted.recode_sorted_sorted_snpEff.ann.gds")
)

