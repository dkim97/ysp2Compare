# Give your self a nice description of what your script does -- this will be
# helpful to you-in-the-future
#
# This script serves as an example of how to extract data from a VCF file into
# a more usable format. In this case, it turns a VCF file into a GDS file,
# and then uses SeqArray and SeqVarTools to extract the variant info (chr, pos,
# depth, annotations, ...) and parse that data into a dataframe. The final
# product of this script is a dataframe that can be used to filter variants
# in the DK samples based on how their genotypes compare to the YSP2-8 genotype
# at the same location.

library(tidyverse)
library(SeqArray)
library(SeqVarTools)
library(here)

# read in functions from other files
source(here("R/pivot_snpeff.R"))
source(here("R/tidy_metrics.R"))

# Read in data. Note that this data may not exist on your machine -- you
# should copy the data into the project data directory, or update the path,
# or both, as appropriate
group_1_genes = read_tsv(
  here("data/snpeff/group_1_merged_dusted.recode_sorted_sorted_snpEff.genes.txt"),
  skip=1
)

# pivot the snpeff data to long format where the variant types are stored in
# column variant_type and number of those variants with a given
# gene/transcript/biotype is stored in count
group_1_long = pivot_snpeff(group_1_genes, "group_1")

# This translates the vcf files to a gds file. A GDS file is a format which
# the bioconductor packages SeqArray (this is one of the more commonly used
# softwares for handling VCF data) uses to manipulate VCF file data
#
# I wrap this in a try catch b/c if the gds file exists, an error is thrown
# and execution halts. You can always manually delete the old gds file if you
# want to overwrite it (i think there is also an overwrite argument --
# check ?seqArray). But, if the file does exist, really I just want to ignore
# the error and move on. With the tryCatch block, the error will still be
# printed, so we can see it, but it won't halt execution
tryCatch(
  SeqArray::seqVCF2GDS(
    here("data/snpeff/group_1_merged_dusted.recode_sorted_sorted_snpEff.ann.vcf"),
    here("data/snpeff/group_1_merged_dusted.recode_sorted_sorted_snpEff.ann.gds")),
  error = function(e){
    message(e)
  })


# After creating the GDS file, we now 'open' it. Note that this does not
# open the entire file, but rather creates a file handle. The difference is
# that we are not reading the data into RAM. This means that we can handle
# enormous VCF files, which is a common use case with this type of data
group_1_gds = SeqArray::seqOpen(
  here("data/snpeff/group_1_merged_dusted.recode_sorted_sorted_snpEff.ann.gds"))

# set a filter such that we only look at loci with at most 1 alternative allele
# In the group_1 data, there is only a single position with more than 1 alt
message("Setting single allele loci filter...")
single_allele_fltr =
  SeqArray::seqGetData(group_1_gds,
                       'variant.id')[SeqVarTools::nAlleles(group_1_gds) == 2]
message(sprintf("...total variants: %s",
                length(SeqVarTools::nAlleles(group_1_gds))))
message("...FALSE refers to the number of multi allelic loci")
print(table(SeqVarTools::nAlleles(group_1_gds) == 2))
SeqArray::seqSetFilter(group_1_gds,
                       variant.id=single_allele_fltr)

# This is the set of 'fields' we're going to extract from the GDS file. Note
# that the return type of each of these 'fields' will not be the same -- eg,
# genotype returns a 3 dimensional matrix, chromosome returns a vector, ...
data_to_extract = c("chromosome", "position", "$ref", "$alt",
                    "annotation/qual", "genotype", "annotation/info/ANN",
                    "annotation/info/DP", "annotation/format/DP",
                    "annotation/format/RO", "annotation/format/AO",
                    "sample.id")

# `seqGetData` is one of the methods we can use to extract data from the
# GDS file. Using `seqGetData` extracts the data from the file and returns it
# as a list in memory. In the example below, we now have the annotations data
# in memory and can start to parse it
group_1_annotations = seqGetData(group_1_gds, data_to_extract)

# create unique IDs for the variants
variant_ids = paste0("var_",
                     seq(1,dim(group_1_annotations$`annotation/format/DP`)[2]))


# The annotations are returned in a list of two items, `length` and `data`.
# The `length` list looks something like 3, 1, 6, 2, ... where the number
# represents the number of times the variant at that index is annotated to a
# different spnEFF annotation. These annotations are stored in the list `data`.
# So, variant 1 has 3 different annotations in `data`, variant 2 has 1
# annotation, and so on.
#
# The goal of the code below is to extract a "key" column from the
# annotations to the variants. If we consider the variant index to be the key,
# what we want to know is that the first 3 entries of the `data` map to the
# variant at index 1
variant_index = map2(
  seq(length(group_1_annotations$`annotation/info/ANN`$length)),
  group_1_annotations$`annotation/info/ANN`$length,
  ~rep(.x,.y)) %>%
  unlist()

# this regex will be used to extract, from a annotation string which might
# look like:
# "C|synonymous_variant|LOW|CKF44_00014|CKF44_00014|transcript|gnl|..."
# the following:
# "C", "synonymous_variant", "LOW", "CKF44_00014"
snpeff_annotation_regex = "^([^|]+)\\|([^|]+)\\|([^|]+)\\|([^|]+)"

# create a dataframe where the first column is the variant name, the second
# column is the full annotation string, and the `alt`, `effect`, `impact`,
# and `gene` columns are created by parsing out the relevant information from
# the annotation string
variant_annotation_df = tibble(
  variant = paste('var', variant_index, sep="_"),
  annotation = group_1_annotations$`annotation/info/ANN`$data) %>%
  mutate(alt = str_match(annotation, snpeff_annotation_regex)[, 2],
         effect = str_match(annotation, snpeff_annotation_regex)[, 3],
         impact = str_match(annotation, snpeff_annotation_regex)[, 4],
         gene = str_match(annotation, snpeff_annotation_regex)[, 5])

# this is part of a function from the BSA package which was set up as a
# convenience to parse some of this data. It is reused here -- see
# R/tidy_metrics.R
data_to_transform = list(
  genotype     = group_1_annotations$genotype[1,,],
  RealDepth    = group_1_annotations$`annotation/format/DP`,
  Reference    = group_1_annotations$`annotation/format/RO`,
  Alternative1 = group_1_annotations$`annotation/format/AO`$data)
# result is a long data frame with columns sample, variant, RealDepth,
# Reference (which is the reference allele depth) and Alternative1, which is
# the Alternative depth. Note that if single_allele_loci_only is not set,
# then this is not actually Alternative1, but just Alternative
depth_metrics_df = suppressMessages(purrr::map2(data_to_transform,
                                                names(data_to_transform),
                                                tidy_metrics,
                                                group_1_annotations$sample.id,
                                                variant_ids) %>%
                                      purrr::reduce(dplyr::left_join))

# Now we're going to create our analysis data -- in this case, we're going to
# join the depth df to a temporary dataframe with the chromosome, position, etc
# information. Then we join the annotations dataframe and filter for positions
# which have RealDepth (actual reads aligning over a given position,
# not counting reads which aligned with a gap over that variant position) of
# at least 10
group_1_annotated_df = depth_metrics_df %>%
  left_join(
  tibble(variant = variant_ids,
         chr = group_1_annotations$chromosome,
         pos = group_1_annotations$position,
         ref = group_1_annotations$`$ref`,
         alt = group_1_annotations$`$alt`)) %>%
  left_join(variant_annotation_df %>%
              select(-annotation, -alt),
            by = c('variant'),
            relationship = 'many-to-many') %>%
  filter(!is.na(RealDepth) & RealDepth >= 10) %>%
  # exclude any marker sequences
  filter(!chr %in% c("NAT", "G418")) %>%
  # factor the `impact` column -- note that I ran everything up to the
  # `filter` step above, and then looked at the unique values in `impact`
  # to figure out the levels. The order of the labels determines how
  # this column will sort, with "HIGH" being highest impact and "MODIFIER"
  # the lowest
  mutate(impact = factor(impact,levels = c('HIGH', 'MODERATE',
                                           'LOW', 'MODIFIER'))) %>%
  # `HIGH` on top
  arrange(impact)

# Now we're going to extract just the ysp2_8 data and reform it into a
# dataframe with columns `variant` and `ysp2_8_genotype`
ysp2_8_genotype_df = group_1_annotated_df %>%
  filter(sample == 'ysp2_8') %>%
  select(variant, genotype) %>%
  dplyr::rename(ysp2_8_genotype = genotype)

# And finally, we remove the ysp2_8 data from the group_1 data, and then join
# the ysp2_genotype_df so that we now have an additional column called
# `ysp_2_8_genotype`
#
# We filter for "complete cases" (or, only keep rows with no NAs in any column)
# so that we're only looking at variants for which we have reliable data over
# both a given sample, and ysp2_8
dk_samples_analysis_df = group_1_annotated_df %>%
  filter(sample != 'ysp2_8') %>%
  left_join(ysp2_8_genotype_df,
            by = 'variant',
            relationship = 'many-to-many') %>%
  filter(complete.cases(.))

# At this point, we can use the fields `ysp2_8_genotype` and `genotype`
# (which refers to the DK sample genotype) to look at variants which are
# called as either REF or ALT in the DK sample and not in the ysp2_8 sample,
# or vice versa

# for example, this would select variants which are ALT in the DK samples
# but REF in the ysp2_8 strain
dk_samples_analysis_df %>%
  filter(genotype == 1, ysp2_8_genotype == 0) %>%
  View()

# you might wish to only look at those which have predicted impact HIGH
dk_samples_analysis_df %>%
  filter(genotype == 1, ysp2_8_genotype == 0, impact=="HIGH") %>%
  View()

# and now you might wish to save this to a variable
high_impact_dk_alt_ysp2_ref = dk_samples_analysis_df %>%
  filter(genotype == 1, ysp2_8_genotype == 0, impact=="HIGH")

# or even write it to file
# NOTE: uncomment to run this, but comment it out to prevent yourself from
# accidently overwriting data
#
# write_csv(high_impact_dk_alt_ysp2_ref,
#           here("output/high_impact_dk_alt_ysp2_ref.csv"))

# finally -- you want to try to keep your script(s) clean. Lines 171 to ~191
# above represent what I might do as I go along, but I would likely clean it
# up and only save this. If I wanted to keep the variable
# `high_impact_dk_alt_ysp2_ref` for later, I'd do this in two steps rather
# than one.
#
# NOTE!! I typically comment out my "write" statements, so that I don't
# accidently overwrite them
#
# dk_samples_analysis_df %>%
#   filter(genotype == 1, ysp2_8_genotype == 0) %>%
#   write_csv(here("output/dk_alt_ysp2_ref.csv"))
#
# dk_samples_analysis_df %>%
#   filter(genotype == 0, ysp2_8_genotype == 1) %>%
#   write_csv(here("output/dk_ref_ysp2_alt.csv"))

# There comes a point where there are too many lines in a single script for it
# to be managable. Try to think about this as you are exploring the data --
# have you reached a natural breakpoint? Maybe it is time to move over to a
# new script that is focused on a single analytic taks, with a descriptive
# name that matches it. Or, possibly you want to create a notebook so that you
# can mix text, code and plots?
