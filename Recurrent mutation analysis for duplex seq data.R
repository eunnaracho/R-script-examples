#Adapted from R script written by Annette Dodge and written as an extension to Annette's mutation analysis
#Uses .mut files processed by Annette's script which list mutations by nucleotide
#Identifies single nucleotide variants that are present in more than one sample ("recurrent mutations")
#Outputs graphs showing recurrent mutations by sample, by chromosome, and by frequency.

library(ggplot2)
library(ggpp)
library(ggpattern)
library(tidyverse)
library(GenomicRanges)
library(fuzzyjoin)
library(writexl)
library(stringr)
library(ggthemes)
library(gtools)
library(flextable)
library(grid)
library(cowplot)

#Chemical name and concentration unit (for labels in figures and file names)
chemical <- "hiPSC liver KBrO3"
conc_unit <- "\U03bcM"

#Combine .mut files in one data frame
path <- "C:/Users/EUCHO/Desktop/NIPH DS data/hiPSC_liver"

filenames_list <- list.files(path= path, full.names=TRUE)

All <- lapply(filenames_list,function(filename){
  print(paste("Merging",filename,sep = " "))
  read.delim(filename)
})

dat <- do.call(rbind.data.frame, All)

#Add a dose column
sampledat <- read.delim("C:/Users/EUCHO/Desktop/NIPH DS data/sample_list_liver.txt")
sampledat$dose <- factor(sampledat$dose)
dat <- dplyr::left_join(dat, sampledat)

##############################
#After uploading and compiling .mut files (raw list of nucleotides and mutations) in the dataframe "dat", Annette's script is used identify analyzable nucleotides
#This script uses the "dat" dataframe processed by Annette's script
##############################

##########Recurrent mutation analysis##########

#Filter all mutations observed in the same position in 2 or more samples at Variant allele frequency (VAF) <= 0.01
recurrent_removed <- dat %>%
  filter(filter != "EndRepairFillInArtifact" & filter != "EndRepairFillInArtifact,v2" & filter != "EndRepairFillInArtifact,NM8.0,v2" & filter != "EndRepairFillInArtifact,NM8.0") %>%
  filter(VAF <= 0.01) %>%
  group_by(start, variation_type, alt) %>%
  add_tally(name = "recurrence") %>%
  filter((recurrence == 1 & variation_type != "no_variant")|variation_type == "no_variant") %>%
  subset(select = -c(recurrence))

#Isolate recurrent snv, indels, mnv, and sv observed in 2 or more samples at VAF 0.01
recurrent_mut_all <- dat %>%
  #filter(dose != 0) %>%
  filter(filter != "EndRepairFillInArtifact" & filter != "EndRepairFillInArtifact,v2" & filter != "EndRepairFillInArtifact,NM8.0,v2" & filter != "EndRepairFillInArtifact,NM8.0") %>%
  filter(variation_type != "no_variant" & VAF <= 0.01) %>%
  group_by(start, variation_type, alt) %>% 
  add_tally(name = "recurrence") %>%
  filter(recurrence >= 2) %>%
  ungroup()

recurrent_mut_hifreq <- dat %>%
  group_by(start, variation_type, alt) %>%
  filter(VAF > 0.01) %>%
  add_tally(name = "recurrence") %>%
  filter(recurrence >= 2) %>%
  ungroup()

write_xlsx(recurrent_mut_all, paste0(chemical,"_Recurrent mutations_all.xlsx"))

#Adds 'chemical' column for use in the "Overlapping recurrent mutations" code
recurrent_mut_all2 <- select(recurrent_mut_all, contig:variation_type, alt, alt_depth, depth, subtype, recurrence) %>%
  mutate(chemical = chemical)

write_xlsx(recurrent_mut_all2, paste0(chemical,"_Recurrent mutations_all_2.xlsx"))

#Calculate the number of recurring mutations by subtype
recurrent_mut_spec <- recurrent_mut_all['subtype'] 

#Function that converts purine bases to pyrimidine bases
convert_AG_subtype <- function(x){
  x['subtype'][x['subtype']== "G>T"] <- "C>A"
  x['subtype'][x['subtype']== "G>A"] <- "C>T"
  x['subtype'][x['subtype']== "G>C"] <- "C>G" 
  x['subtype'][x['subtype']== "A>T"] <- "T>A"
  x['subtype'][x['subtype']== "A>G"] <- "T>C" 
  x['subtype'][x['subtype']== "A>C"] <- "T>G" 
  x['subtype'][x['subtype']== "."] <- "Indel/mnv/sv"
  return(x)
}

recurrent_mut_spec <- convert_AG_subtype(recurrent_mut_spec)

recurrent_mut_spec <- recurrent_mut_spec %>%
  group_by(subtype) %>%
  summarise(number_mut = n())

recurrent_mut_spec$subtype <- factor(recurrent_mut_spec$subtype, levels=c('C>A', 'C>G','C>T', 'T>A', 'T>C', 'T>G', 'Indel/mnv/sv'))

#write.table(recurrent_mut_spec, "Recurrent mutation spectra.txt", col.names=NA)

#Plot the number of recurrent mutations by type
ggsave(filename = paste0(chemical,"_Recurrent by subtype.pdf"),
       ggplot(data = recurrent_mut_spec, mapping = aes(x = subtype, y = number_mut)) + 
         geom_bar(stat = "identity") +
         geom_text(aes(label = number_mut), 
                   position = position_stacknudge(y = 50)),
       width = 20, height = 15, units = "cm")


#Counts the number of recurring snv and indels in each sample
recur_mut_spec_per_samp <- select(recurrent_mut_all, subtype, sample, dose)

recur_mut_spec_per_samp <- convert_AG_subtype(recur_mut_spec_per_samp)  

recur_mut_spec_per_samp <- recur_mut_spec_per_samp %>%
  group_by(sample, subtype) %>% 
  mutate(num_recurring_mut = n()) %>%
  ungroup() %>%
  select (subtype:num_recurring_mut, dose)

recur_mut_spec_per_samp <- recur_mut_spec_per_samp[!duplicated(recur_mut_spec_per_samp), ]
recur_mut_spec_per_samp$subtype <- factor(recur_mut_spec_per_samp$subtype, levels=c('C>A', 'C>G','C>T', 'T>A', 'T>C', 'T>G', 'Indel/mnv/sv'))

factor_reordered <- factor(recur_mut_spec_per_samp[["sample"]],
                           levels = unique(recur_mut_spec_per_samp[["sample"]][mixedorder(recur_mut_spec_per_samp[["dose"]])]),
                           ordered = FALSE)

recur_mut_spec_per_samp[["sample"]] <- factor_reordered 

#Plot the spectrum of recurrent mutations by the number of occurrences
ggsave(filename = paste0(chemical,"_Recurrent mut spec by sample_VC.pdf"),
       ggplot(data = recur_mut_spec_per_samp, mapping = aes(x = subtype, y = num_recurring_mut, fill = subtype)) + 
         geom_bar(stat = "identity") + 
         labs(title = "Recurrent mutation spectra by sample",
              x = "Subtype", 
              y = "Number of recurrent mutations",
              fill = "Subtype") +
         theme(text = element_text(size = 10),
               axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1)) +
         facet_grid(~sample),
       width = 25, height = 15, units = "cm")

# Total number of recurrent mutations in each sample
recur_mut_per_samp1 <- recur_mut_spec_per_samp %>%
  group_by(sample) %>%
  mutate(mut_total = sum(num_recurring_mut)) %>%
  select(sample, dose, mut_total)

recur_mut_per_samp1 <- recur_mut_per_samp1[!duplicated(recur_mut_per_samp1), ]

write.table(recur_mut_per_samp1, paste0(chemical, "_Recurrent mutations in treated samples.txt"))

recur_mut_summary <- cbind(recur_mut_per_samp1, Depth["depth"])

#Calculate the percentage of each subtype in each sample
recur_mut_per_samp2 <- recur_mut_spec_per_samp %>%
  group_by(sample) %>%
  mutate(mut_total = sum(num_recurring_mut)) %>%
  select(sample, dose, subtype, num_recurring_mut, mut_total) %>%
  mutate(percentage = num_recurring_mut/mut_total*100)

write.csv(recur_mut_per_samp2, paste0(chemical, "_Recurrent mutation by subtype in each samples.csv"))

ggsave(filename = paste0(chemical,"_Recurrent mut spec by sample_VC.pdf"),
       ggplot(data = recur_mut_per_samp2, mapping = aes(x = subtype, y = percentage, fill = subtype)) + 
         geom_bar(stat = "identity") + 
         labs(title = paste0(chemical, " Recurrent mutation spectra by sample_VC"),
              x = "Subtype", 
              y = "Percentage",
              fill = "Subtype") +
         theme(text = element_text(size = 10),
               axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1)) +
         facet_grid(~sample),
       width = 25, height = 15, units = "cm")

#Counts the number of recurring mutations in each target region
recur_per_target <- recurrent_mut_all %>%
  group_by(contig) %>%
  summarise(num_recurring_mut1 = n())

recur_per_target <- left_join(recur_per_target, genic_regions) %>% 
  select("contig", "location_relative_to_genes", "num_recurring_mut1")

recur_per_target <- recur_per_target[order(-recur_per_target$num_recurring_mut1), ]

order <- recur_per_target$contig

#Plot the number of recurrent mutations by target site
ggsave(filename = paste0(chemical,"_Recurrent mut per target.pdf"),
       ggplot(data = recur_per_target, 
              mapping = aes(x = factor(contig, levels = order),
                            y = num_recurring_mut1,
                            fill = location_relative_to_genes)) + 
         geom_bar(stat="identity") +
         labs(x = "Target Site", 
              y = "Number of recurrent mutations",
              fill = "Target site type") +
         geom_text(aes(label = num_recurring_mut1), 
                   position = ggpp::position_stacknudge(y = 15),
                   size = 3),
       width = 25, height = 15, units = "cm")

#Isolate unique snv 
unique_snv_only <- dat %>%
  group_by (start, variation_type) %>%
  add_tally(name = "recurrence") %>%
  filter (recurrence == 1 & VAF <= 0.01 & variation_type == "snv")

unique_snv_only <- convert_AG_subtype(unique_snv_only)

unique_snv_only <- unique_snv_only %>%
  group_by(subtype) %>%
  summarise(number_mut = n())

unique_snv_only$subtype <- factor(unique_snv_only$subtype, levels=c('C>A', 'C>G','C>T', 'T>A', 'T>C', 'T>G', 'Indel/mnv/sv'))

ggsave(filename = paste0(chemical,"_Unique snv by subtype.pdf"),
       ggplot(data = unique_snv_only, mapping = aes (x = subtype, y = number_mut)) + 
         geom_bar(stat = "identity") +
         geom_text(aes(label = number_mut), 
                   position = position_stacknudge(y = 50)),
       width = 20, height = 15, units = "cm")


#Mutations by the number of recurrence
#Group mutations by the number of recurrence in the dataset 
recur_distribution <- dplyr::select(recurrent_mut_all, subtype, sample, dose, recurrence)

recur_distribution <- convert_AG_subtype(recur_distribution)

recur_distribution <- recur_distribution %>%
  group_by(recurrence) %>%
  add_tally(name = "Number_mutation") 

simple_distribution <- group_by (recur_distribution, recurrence) %>%
  summarise(Number_mutation = n())

colnames(simple_distribution) <- c("Recurrence", "Number of mutations")

recur_distribution_type <- group_by(recur_distribution, recurrence, subtype) %>%
  add_tally(name = "Number_mutation") %>%
  dplyr::select(subtype, recurrence, Number_mutation)

recur_distribution_type <- recur_distribution_type[!duplicated(recur_distribution_type), ]

recur_distribution_type$subtype <- factor(recur_distribution_type$subtype, levels=c('C>A', 'C>G','C>T', 'T>A', 'T>C', 'T>G', 'Indel/mnv/sv'))

plot_recurrence <- ggplot(data = recur_distribution_type, mapping = aes(x = subtype, y = Number_mutation, fill = subtype)) + 
  geom_bar(stat = "identity") + 
  labs(title = "Mutations by recurrence",
       x = "Subtype", 
       y = "Number of mutations",
       fill = "Subtype") +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1)) +
  facet_grid(~recurrence)

recur_table <- simple_distribution %>% 
  flextable::flextable() %>% 
  as_raster()

table <- ggplot() + 
  theme_void() + 
  annotation_custom(rasterGrob(recur_table), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

plot_w_table <- plot_grid(plot_recurrence, table, nrow = 1, ncol = 2, rel_widths = c(4, 1) )

ggsave(filename = paste0(chemical,"_Mutations by recurrence.pdf"), plot_w_table, width = 30, height = 15, units = "cm") 
