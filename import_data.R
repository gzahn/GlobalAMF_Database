library(tidyverse)
library(ShortRead)
library(phyloseq) # for last bit

# find files
sh_files <- 
list.files("/home/gzahn/Desktop/GIT_REPOSITORIES/GlobalAMF_Database/SH",
           full.names = TRUE, pattern = "sh_list.txt$")
fasta_files <- 
list.files("/home/gzahn/Desktop/GIT_REPOSITORIES/GlobalAMF_Database/Fasta",
           full.names = TRUE, pattern = ".fas$")
metadata_files <- 
list.files("/home/gzahn/Desktop/GIT_REPOSITORIES/GlobalAMF_Database/Metadata",
           full.names = TRUE, pattern = "_sample_list.txt$")


# tests to make sure same number of files, in same order
sh_names <- sh_files %>% basename() %>% str_remove("_sh_list.txt")
fa_names <- fasta_files %>% basename() %>% str_remove("species_") %>% str_remove("\\.fas")
md_names <- metadata_files %>% basename() %>% str_remove("_sample_list.txt")
length_same <- length(metadata_files) == length(sh_files) & length(metadata_files) == length(fasta_files)
order_same <- identical(sh_names,fa_names) & identical(sh_names, md_names)

# make new list to hold results
species_data_list <- list()

# run tests on file order
if(length_same & order_same){

  # for-loop to process all 3 files for each species, if test passes
for(i in seq_along(metadata_files)){
  metadata <- read_delim(metadata_files[i])
  fasta <- readFasta(fasta_files[i])
  SH <- read_delim(sh_files[i])
  sample_ids <- fasta@id %>% as.character() %>% str_split("\\|") %>% map_chr(2) %>% str_remove("SampleID_") %>% as.numeric()
  markers <- fasta@id %>% as.character() %>% str_split("\\|") %>% map_chr(4)
  sequences <- fasta@sread %>% as.character()

# build new vars into dataframe
  seq_data <- 
    data.frame(id=sample_ids,
               marker=markers,
               sequence=sequences) %>% 
    mutate(marker = marker %>% str_remove("marker_")) %>% 
    full_join(SH) # join SH (species hypothesis ID) and taxonomy with metadata and DNA sequence

  df <- metadata %>% 
    full_join(seq_data, by="id", multiple='all') %>% 
    mutate(pH = as.numeric(pH), # take care to force data types for columns in acse of odd entries
           abundances = as.numeric(abundances),
           year_of_sampling = as.numeric(year_of_sampling),
           MAT = as.numeric(MAT),
           MAP = as.numeric(MAP),
           latitude = as.numeric(latitude),
           longitude = as.numeric(longitude),
           primers = as.character(primers),
           id = as.character(id))
  # assign data frame to list
  species_data_list[[i]] <- df %>% as.data.frame()
  
  }
}

#inspect for sanity
species_data_list %>% map("Species") %>% map_chr(unique) # should be all species

#save intermediate file (it's rather large and unwieldy)
saveRDS(species_data_list,"/home/gzahn/Desktop/GIT_REPOSITORIES/GlobalAMF_Database/species_data_list.RDS")

# pull them all into one big data frame
df <- purrr::reduce(species_data_list,full_join)

#inspect for sanity...do biomes match locations?
glimpse(df)
df %>% 
  dplyr::select(longitude,latitude,Biome) %>% 
  unique.data.frame() %>% 
  ggplot(aes(longitude,latitude,color=Biome)) +
  geom_point() +
  scale_color_viridis_d()

# save that file
saveRDS(df,"/home/gzahn/Desktop/GIT_REPOSITORIES/GlobalAMF_Database/full_database_df.RDS")
#reload point
df <- readRDS("/home/gzahn/Desktop/GIT_REPOSITORIES/GlobalAMF_Database/full_database_df.RDS")

# build phyloseq object for downstream microbiome work

# build taxonomy table
tax <- 
df %>% 
  dplyr::select(Kingdom,Phylum,Class,Order,Family,Genus,Species,accession) %>% 
  unique.data.frame() %>% 
  tax_table()
# add accessions as taxa names
taxa_names(tax) <- tax[,8]
# rename taxa ranks
colnames(tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Accession")

# build metadata table
met <- 
df %>% 
  select(all_of(c("id","paper_id","primers","longitude","latitude","continent",
                  "sample_type","marker_size","Biome","MAT","MAP","pH",
                  "year_of_sampling","plants_dominant","target_gene",
                  "sequencing_platform","marker"))) %>% 
  unique.data.frame() %>% 
  sample_data()

# number of unique samples
df$id %>% unique %>% length
# add sample names to slot
sample_names(met) <- met$id

# build "otu" table 
df %>% glimpse
x <- df %>% 
  dplyr::select(id,accession) %>% 
  unique.data.frame()
x <- table(x$id,x$accession) %>% as("matrix") %>% as.data.frame()
# x %>% View
# place samples (rows) in same order as metadata,
# and taxa (cols) in same order as tax_table
otu <- x[sample_names(met),taxa_names(tax)] %>% otu_table(taxa_are_rows = FALSE)

# sanity check
identical(sample_names(met),sample_names(otu))
identical(taxa_names(tax),taxa_names(otu))

# make sure sequences are in same order as the accessions that they belong to
# there are going to be multiple sequences for each accession
# find consensus sequence for each accession
x <- df %>% 
  dplyr::select(sequence,accession) %>% 
  unique.data.frame()

#inspect
x %>% 
  group_by(accession) %>% 
  summarize(l=length(sequence)) %>% 
  pluck("l") %>% plot() # number of unique seqs in each accession

# can consensus seq be built without aligning?
x %>% filter(accession == "FN547496") %>% 
  pluck("sequence") %>% 
  DNAStringSet() %>% 
  DECIPHER::ConsensusSequence() %>% 
  as.character() # yeah, but there are some accessions that have LOTS of seqs


consensus_df <- 
x %>% 
  group_by(accession) %>% 
  summarize(consensus = sequence %>% 
              rep(2) %>% # double each sequence list to account for accessions only represented by 1 sequence 
              DNAStringSet() %>% 
              # DECIPHER::AlignSeqs() %>% # skip alignment step
              DECIPHER::ConsensusSequence() %>% 
              as.character())
row.names(consensus_df) <- consensus_df$accession
# reorder to match otu table
consensus_df <- consensus_df[taxa_names(otu),]

beepr::beep()
# build phyloseq object
ps <- phyloseq(otu_table(otu),
               tax_table(tax),
               sample_data(met))
# add reference sequence slot (from consensus seqs)
ps@refseq <- consensus_df$consensus %>% DNAStringSet()

ps@otu_table
ps@sam_data
ps@tax_table

# save object as RDS
saveRDS(ps,"/home/gzahn/Desktop/GIT_REPOSITORIES/GlobalAMF_Database/ps_object_not_cleaned.RDS")
