rm(list=ls()) 
library(magrittr)
library(stringi)
library(dplyr)
library(Biostrings) 
library(table1)
library(lubridate) 
library(tidyr)
library(parallel) 
library(readr) 
library(microbiome)
library(decontam)
library(MASS)
library(GGally)
library(FactoMineR)
library(vegan) 
library(fpc) 
library(Nonpareil)
library(gridExtra) 
library(cowplot)
library(ape) 
library(logistf) 
library(forcats) 
library(tibble) 
library(treeclimbR)
library(TreeSummarizedExperiment)
library(tidytree)
library(pROC)

n_clust=3

#Rename group for publication graphs
pub_rename <- function(group){
  group[which(group=="Febrile")] <- "Febrile control"
  group
}
#Maximum from each column
colMax <- function(df){
  apply(df, 2, max) %>% unlist()
}
#Number formatting
fnum <- function(x, d=3){
  formatC(signif(x,digits=d), digits=d,format="fg", flag="#")
}

#Excluded on basis of quality
exclude_qual<-c("479-1510181")
#Two swabs from same case (second is later and excluded from primary analyses)
duplicates<-c("69-4968-Fb", "71-4968")

read_groups_ordered <- c("Human", "Bacteria", "Fungus", "Virus", "Archaea")
major_OTUs_ordered <- c("Homo sapiens", "Bacteria", "Fungi", "Viruses", "Archaea")

palette <- c("#E69F00", "#56B4E9", "#009E73")

#PREPARATORY 1
#Read in kraken and bracken data, preprocess

#List all kraken (GTDB) reports and cut down to only have latest version per sample 
kraken_GTDB_files <- list.files("read_binning/kraken/conf", pattern="*_kraken_gtdb_report.txt", full.names=TRUE)

reports <- lapply(kraken_GTDB_files, function(f){
  r<-read.table(f, sep="\t", header=FALSE, fill=TRUE, quote="")
  r$level <- stri_extract_first_regex(r$V6, "[ ]+") %>% nchar() / 2
  r$level[is.na(r$level)] <- 0
  r$file <- stri_extract_first_regex(f, "[^/]+(?=_kraken_gtdb_report[.]txt)")
  r
})
reports %<>% bind_rows()
reports$database="GTDB_conf"

kraken_GTDB_files_noconf <- list.files("read_binning/kraken/noconf", pattern="*_kraken_gtdb_report.txt", full.names=TRUE)

reports_noconf <- lapply(kraken_GTDB_files_noconf, function(f){
  r<-read.table(f, sep="\t", header=FALSE, fill=TRUE, quote="")
  r$level <- stri_extract_first_regex(r$V6, "[ ]+") %>% nchar() / 2
  r$level[is.na(r$level)] <- 0
  r$file <- stri_extract_first_regex(f, "[^/]+(?=_kraken_gtdb_report[.]txt)")
  r
})
reports_noconf %<>% bind_rows()
reports_noconf$database="GTDB_noconf"

reports %<>% bind_rows(reports_noconf)

#Extract sample IDs
reports$sample <- stri_extract_first_regex(reports$file, "^[^.]+")

#Update column names
colnames(reports)[1:6] <- c("cum_percentage", "cum_reads", "reads", "rank", "OTU_ID", "OTU_name")

#Attach sample groups
reports %<>% mutate(group=case_when(
  stri_count_regex(sample, "[-]3[0-9]{3}$") > 0 ~ "KD",
  nchar(stri_extract_first_regex(sample, "(?<=[0-9][-])[0-9]+")) == 4 ~ "Febrile",
  sample %in% c("507", "508") ~ "Negative control",
  TRUE ~ "KD"))

#Clear up OTU names
reports$OTU_name %<>% trimws()
reports$OTU_name %<>% stri_replace_first_regex("^[a-z]__", "")

#Same for Bracken
bracken_GTDB_files  <- list.files("read_binning/bracken/", pattern="*_bracken_gtdb_output.txt", full.names=TRUE)

bracken_reports <- lapply(bracken_GTDB_files, function(f){
  r<-read.table(f, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE)
  r$file <- stri_extract_first_regex(f, "[^/]+(?=_bracken_gtdb_output[.]txt)")
  r
})

bracken_reports %<>% bind_rows()
bracken_reports$database="GTDB"

#Add samples and groups
bracken_reports$sample <- stri_extract_first_regex(bracken_reports$file, "^[^.]+")
bracken_reports$group <- case_when(
  stri_count_regex(bracken_reports$sample, "[-]3[0-9]{3}$") > 0 ~ "KD",
  nchar(stri_extract_first_regex(bracken_reports$sample, "(?<=[0-9][-])[0-9]+")) == 4 ~ "Febrile",
  bracken_reports$sample %in% c("507", "508") ~ "Negative control",
  TRUE ~ "KD")

#Reads per sample
filter(reports, database=="GTDB_conf") %>%
  group_by(file, sample) %>%
  summarise(reads=sum(reads)) %>%
  arrange(file) -> reads
quantile(reads[!reads$sample %in% c("507","508"),]$reads, c(0.25,0.5,0.75))
mean(reads[!reads$sample %in% c("507","508"),]$reads)
sum(reads$reads)

#PREPARATORY 2
#Create metadata table

#Create sample lookup
filter(reports, OTU_ID==9606) %>%
  select(sample, group) %>% unique() -> sample_lookup
sample_lookup$basename <- stri_replace_first_regex(sample_lookup$sample, "[-]?Fb", "")
library_info <- read.csv("metadata/nonclinical/library_info.csv")
library_info$ID %<>% stri_extract_first_regex("^[^ ]+")
library_info$ID[library_info$ID=="219-163006.2"] <- "219-163006-2"
library_info$sample_name <- sample_lookup$sample[match(library_info$ID, sample_lookup$basename)]

sample_lookup$concentration <- as.character(library_info$Conc.ng.ul[match(sample_lookup$sample, library_info$sample_name)])
sample_lookup$concentration[sample_lookup$concentration == "Too Low"] <- 2
sample_lookup$concentration %<>% as.numeric()
sample_lookup$PCR_cycles <- library_info$PCR.Cycles[match(sample_lookup$sample, library_info$sample_name)]
sample_lookup$batch <- library_info$Batch[match(sample_lookup$sample, library_info$sample_name)]
sample_lookup$pool <- library_info$Pool[match(sample_lookup$sample, library_info$sample_name)]
sample_lookup$position <- library_info$Pos[match(sample_lookup$sample, library_info$sample_name)]
sample_lookup$plate <- stri_extract_first_regex(sample_lookup$position, "[0-9]+$") %>% as.numeric()
sample_lookup$col <- stri_extract_first_regex(sample_lookup$position, "[A-Z]")
sample_lookup$row <- stri_extract_first_regex(sample_lookup$position, "[0-9]+") %>% as.numeric()

sample_sheet_files <- list.files("metadata/nonclinical/sequencing_runs/", recursive=TRUE, pattern="H*/*SampleSheet.csv", full.names=TRUE)
sample_indices <- lapply(sample_sheet_files, function(x){
  r<-read.csv(x, stringsAsFactors=FALSE, skip=20)
  r$flowcell <- stri_extract_first_regex(x, "H[A-Z]+")
  r
}) %>% bind_rows() %>% select(Sample_Name, index, index2) %>% unique()

sample_lookup$P5_index <- as.character(reverseComplement(DNAStringSet(sample_indices$index2[match(sample_lookup$sample, sample_indices$Sample_Name)])))
sample_lookup$P7_index <- sample_indices$index[match(sample_lookup$sample, sample_indices$Sample_Name)]

metadata_febrile <- read.csv("metadata/clinical/metadata_febrile.csv", stringsAsFactors=FALSE)
metadata_febrile$doo %<>% as.Date(format="%d/%m/%Y")
metadata_febrile$pldate %<>% as.Date(format="%d/%m/%Y")
metadata_febrile$lab_ill_day <- metadata_febrile$pldate - metadata_febrile$doo + 1
metadata_febrile$lab_ill_day[which(metadata_febrile$lab_ill_day<0)] <- metadata_febrile$lab_ill_day[which(metadata_febrile$lab_ill_day<0)]+365
metadata_febrile$lab_ill_day %<>% as.numeric()
metadata_KD <- read.csv("metadata/clinical/metadata_KD_USA.csv", stringsAsFactors=FALSE)
#Sample ID was mistranscribed and persists in filenames
metadata_KD %<>% filter(ID!=153037) %>% mutate(ID=case_when(
  ID==153035 ~ 153037,
  .default=ID
))

metadata_KD$doo %<>% as.Date(format="%d/%m/%Y")
metadata_KD$ID %<>% stri_replace_all_fixed(".", "-")
metadata_SMH_KD <- read.csv("metadata/clinical/metadata_KD_UK.csv", stringsAsFactors=FALSE)
for(echo in 1:5){
  tmp<-metadata_SMH_KD %>% select(contains(paste0("echo_", echo))) %>% apply(1,function(x){
    res<-names(x)[which(x==1)]
    if(length(res)==0){
      res<-NA
    }
    res
  }) %>% lapply(paste, collapse=";") %>% stri_replace_all_regex("ip_echo_[0-9]_results___","")
  metadata_SMH_KD[[paste0("echo_", echo)]] <- tmp
}
metadata_SMH_KD %<>% mutate(
  aneurysm=case_when(
    grepl("aneurysm", paste(echo_1,echo_2,echo_3,echo_4,echo_5)) ~ TRUE,
    echo_1=="NA" & echo_2=="NA" & echo_3=="NA" & echo_4=="NA" & echo_5 == "NA" ~ NA,
    TRUE ~ FALSE
  )
)
metadata_SMH_KD$date_of_birth %<>% as.Date(format="%d/%m/%Y")
metadata_SMH_KD$ip_date_present_loc_hosp %<>% as.Date(format="%d/%m/%Y")
metadata_SMH_KD$ip_date_transf_ref_hosp %<>% as.Date(format="%d/%m/%Y")
metadata_SMH_KD$ip_date_fever_onset %<>% as.Date(format="%d/%m/%Y")
metadata_SMH_KD$date_adm <- if_else(is.na(metadata_SMH_KD$ip_date_present_loc_hosp), metadata_SMH_KD$ip_date_transf_ref_hosp, metadata_SMH_KD$ip_date_present_loc_hosp)
metadata_SMH_KD$doo <- metadata_SMH_KD$ip_date_fever_onset
metadata_SMH_KD$age <- as.numeric((metadata_SMH_KD$date_adm - metadata_SMH_KD$date_of_birth) / 365.25) %>% round(1)
metadata_SMH_KD$gender <- ifelse(metadata_SMH_KD$gender=="male",1,2) #check correct coding!
metadata_SMH_KD$group <- "KD"
metadata_SMH_KD %<>% mutate(eth=case_when(
  final_ethnicity=="asian_other" ~ 1,
  final_ethnicity=="black_african" ~ 2,
  final_ethnicity=="black_carribean" ~ 2,
  final_ethnicity=="mixed_other" ~ 6,
  final_ethnicity=="mixed_white_african" ~ 6,
  final_ethnicity=="other" ~ 9,
  final_ethnicity=="white_british" ~ 3,
  final_ethnicity=="white_other" ~ 3,
  is.na(final_ethnicity) ~ NA_real_
))
metadata_SMH_KD %<>% rename(ID=subject_id) %>% mutate(ID=as.character(ID))
metadata_febrile$pcrp %<>% as.numeric()
metadata_febrile$pcrp <- metadata_febrile$pcrp * 10
metadata_KD$pcrp %<>% as.numeric()
metadata_KD$pcrp <- metadata_KD$pcrp * 10
metadata_febrile$pwbc %<>% as.numeric()
metadata_febrile$ppolys %<>% as.numeric()
metadata_KD$pwbc %<>% as.numeric()
metadata_KD$ppolys %<>% as.numeric()
metadata_KD$grouped_syndrome<-"KD"

metadata_febrile$country <- "USA"
metadata_KD$country <- "USA"
metadata_SMH_KD$country <- "UK"
metadata_SMH_KD$grouped_syndrome <- "KD"
metadata_SMH_KD$pcrp <- metadata_SMH_KD$ip_lab_1_crp
metadata_SMH_KD$ppolys <- metadata_SMH_KD$ip_lab_1_neutro
metadata_SMH_KD$pwbc <- metadata_SMH_KD$ip_lab_1_wcc
metadata_SMH_KD$lab_ill_day <- as.numeric(as.Date(metadata_SMH_KD$ip_npa_ts_mbiol_date, format="%d/%m/%Y") - metadata_SMH_KD$doo + 1)
metadata_febrile$group <- "Febrile"
metadata_febrile$eth %<>% as.numeric()
metadata_KD$group <- "KD"
metadata_febrile$age %<>% as.numeric()
metadata_KD$age %<>% as.numeric()
metadata_KD$gender %<>% as.numeric()
metadata_febrile$gender %<>% as.numeric()
metadata_febrile$ID %<>% as.character()
metadata_KD$ID %<>% as.character()
metadata_KD$aneurysm <- as.numeric(metadata_KD$zworstever) >= 2.5
metadata_febrile$aneurysm <- NA
metadata <- bind_rows(select(metadata_febrile, ID, gender, age, eth, country, pcrp, ppolys, pwbc, grouped_syndrome, doo, lab_ill_day, aneurysm),
                      select(metadata_KD,ID, gender, age, eth, country, pcrp, ppolys, pwbc, grouped_syndrome, doo, lab_ill_day, aneurysm),
                      select(metadata_SMH_KD, ID, gender, age, eth, country, pcrp, ppolys, pwbc, grouped_syndrome, doo, lab_ill_day, aneurysm))

numeric_ID <- stri_extract_first_regex(sample_lookup$sample, "[0-9]{4}[0-9]*([-][0-9]+)?")
sample_lookup %<>% cbind(metadata[match(numeric_ID, metadata$ID), -1])

sample_lookup %<>% mutate(season=case_when(
  between(month(doo), 3,5) ~ "Spring",
  between(month(doo), 6,8) ~ "Summer",
  between(month(doo), 9,11) ~ "Autumn",
  !is.na(doo) ~ "Winter",
  TRUE ~ NA_character_
))

antibiotic <- read.csv("metadata/clinical/USA_antibiotics.csv", stringsAsFactors=FALSE)
colnames(antibiotic) <- c("ID", "Index", "Group", "Antibiotic", "Amoxicillin", "Coamoxiclav", "Azithromycin", "Cefazolin", "Cephalexin", "Clindamycin", "Ceftriaxone", "Other_antibiotic", "Other_antibiotic_name", "X")
numeric_ID <- as.numeric(stri_extract_first_regex(sample_lookup$sample, "[0-9]{4}[0-9]*"))
sample_lookup %<>% cbind(antibiotic[match(numeric_ID, antibiotic$Index), c("Antibiotic"), drop=FALSE])

swab_dates <- read.csv("metadata/clinical/USA_swab_dates.csv", stringsAsFactor=FALSE) %>%
  mutate(swab_day=as.Date(swab_day, format="%d/%m/%Y"))

sample_lookup %<>% mutate(ID=stri_extract_first_regex(sample, "(?<=[-])[0-9]+"),
                          swab_date=swab_dates$swab_day[match(ID, as.character(floor(swab_dates$ID)))],
                          swab_date=case_when(
                            is.na(swab_date) ~ doo + lab_ill_day,
                            .default=swab_date
                          ),
                          swab_ill_day=as.numeric(swab_date-doo+1)) %>%
  select(-ID)

#Quality
quality<-read.table("metadata/nonclinical/fastq_quality.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(quality) <- c("file","NA", "format","type","num_seqs","sum_len","min_len","avg_len","max_len","Q1","Q2","Q3","sum_gap","N50","Q20","Q30")
quality$sample <- stri_extract_first_regex(quality$file, "^[^_]+")
quality$pair <- stri_extract_first_regex(quality$file, "R[12]")
quality %<>% pivot_wider(id_cols=sample, names_from=pair, values_from=Q20)
sample_lookup$R1_Q20 <- quality$R1[match(sample_lookup$sample, quality$sample)]
sample_lookup$R2_Q20 <- quality$R2[match(sample_lookup$sample, quality$sample)]

#Classify rows into groups
get_group_rows <- function(group){
  group_rows <- lapply(which(reports$OTU_name == group), function(x){
    domains<-which(reports$rank%in%c("D","U","R","R1","R2")) 
    next_domain <- domains[which(domains > x)[1]]
    if(is.na(next_domain))
      next_domain <- nrow(reports)+1
    x:(next_domain-1)
  }) %>% unlist()
}

reports$bacteria <- FALSE
reports$bacteria[get_group_rows("Bacteria")] <- TRUE  

reports$fungus <- FALSE
reports$fungus[get_group_rows("Fungi")] <- TRUE  

reports$virus <- FALSE
reports$virus[get_group_rows("Viruses")] <- TRUE

reports$archaea <- FALSE
reports$archaea[get_group_rows("Archaea")] <- TRUE

reports[ ,c("OTU_ID", "bacteria", "virus", "fungus", "archaea")] %>% unique() ->
  class_OTUs

bracken_reports %<>% merge(class_OTUs, by.x="taxonomy_id", by.y="OTU_ID", all.x=TRUE)

bracken_reports %<>% mutate(read_group = case_when(
  bacteria ~ "Bacteria",
  fungus ~ "Fungus",
  virus ~ "Virus",
  archaea ~ "Archaea",
  name == "Homo sapiens" ~ "Human",
  TRUE ~ "Other"
))

reports %<>% mutate(read_group = case_when(
  bacteria ~ "Bacteria",
  fungus ~ "Fungus",
  virus ~ "Virus",
  archaea ~ "Archaea",
  OTU_name == "Homo sapiens" ~ "Human",
  TRUE ~ "Other"
))

groups=c("bacteria", "fungus", "virus", "archaea")
bracken_reports %<>% group_by(sample, read_group) %>% mutate(new_est_read_prop=new_est_reads/sum(new_est_reads)) %>% ungroup()

#SECTION 1
#Demographics
filter(sample_lookup, group %in% c("KD","Febrile") &! sample %in% c(duplicates[2], exclude_qual)) %>%
  mutate(Sex=ifelse(gender==1, "Male", "Female")) %>%
  mutate(Ethnicity=case_when(
    eth==0 ~ "Unknown",
    eth==1 ~ "Asian",
    eth==2 ~ "Black/African American",
    eth==3 ~ "Caucasian",
    eth==4 ~ "Hispanic",
    eth==6 ~ "Multiple",
    eth==7 ~ "American Indian/Alaska Native",
    eth==8 ~ "Native Hawaiian or Other Pacific Islander",
    eth==9 ~ "Other",
    is.na(eth) ~ "Unknown",
    TRUE ~ "Mistake"),
    Season=season,
    `CRP (mg/dL)`=pcrp,
    `PMNs (10<sup>9</sup>/L)`=ifelse(country=="USA", ppolys*pwbc/100, ppolys)
  ) %>%
  rename(`Age (y)`=age, Country=country, Group=group, `WBC (10<sup>9</sup>/L)`=pwbc, `Day of illness`=swab_ill_day, `Antibiotics before sampling`=Antibiotic, `Coronary artery aneurysm`=aneurysm) %>%
  mutate(`Antibiotics before sampling`=as.logical(`Antibiotics before sampling`)) %>%
  table1(~ Sex + `Age (y)` + Country + Season + Ethnicity + `Day of illness` + `CRP (mg/dL)` + `WBC (10<sup>9</sup>/L)` + `PMNs (10<sup>9</sup>/L)` + `Antibiotics before sampling` + `Coronary artery aneurysm` | Group, data=., overall=FALSE, render.continuous="Median [Q1, Q3]") ->
  this_table1_PUB

as.character(this_table1_PUB) %>% cat(file="output/PUB_table1.html")

#SECTION 2
#Contaminants (both reference library and sequencing)
#Mislabelled human
reports %>%
  group_by(database, sample) %>%
  mutate(cum_percentage=cum_reads/sum(reads))%>%
  ungroup() %>%
  filter(rank=="S") %>%
  pivot_wider(id_cols=c(sample, group), names_from=c(OTU_name, database), names_sep="+", values_from=cum_percentage, values_fill=list(cum_percentage=0)) ->
  total_abundance_matrix

species <- colnames(total_abundance_matrix)[-grep("Homo sapiens", colnames(total_abundance_matrix), fixed=TRUE)][-c(1:2)]

total_abundance_matrix[,-c(1,2)] %<>% apply(2, function(x){
  x[x==0] <- min(x[x!=0]) / 2
  x
})

cl<-makeCluster(n_clust)
clusterExport(cl, "total_abundance_matrix")
clusterEvalQ(cl, {library(magrittr)})

parLapply(cl, species, function(s){
  s_parts <- strsplit(s, "+", fixed=TRUE)[[1]]
  model<-as.formula(paste0("`", s, "`~ `Homo sapiens+", s_parts[2], "`")) %>%
    glm(data=total_abundance_matrix) %>% summary()
  list(species=s_parts[1], database=s_parts[2], intercept=model$coefficients[1,1], intercept_p=model$coefficients[1,4],  coef=model$coefficients[2,1], p=model$coefficients[2,4], r2=1-(model$deviance/model$null.deviance))
}) -> human_assoc

human_assoc %<>% bind_rows()

stopCluster(cl)

human_assoc %<>% mutate(p.log10=-log10(p)*sign(coef)) %>%
  pivot_wider(id_cols=species, names_from=database, values_from=c(intercept, intercept_p, coef, p, r2, p.log10))

(read_file("resources/ref_genome_contamination/taxids_with_human_HG.txt") %>%
    stri_split_fixed("\n"))[[1]] %>%
  as.numeric() -> taxids_human_HG
taxids_human_HG <- taxids_human_HG[!is.na(taxids_human_HG)]

(read_file("resources/ref_genome_contamination/taxids_with_human_noHG.txt") %>%
    stri_split_fixed("\n"))[[1]] %>%
  as.numeric() -> taxids_human_noHG
taxids_human_noHG <- taxids_human_noHG[!is.na(taxids_human_noHG)]

OTU_known_human <- data.frame(OTU_ID=c(taxids_human_HG, taxids_human_noHG), HG_match=c(rep(TRUE, length(taxids_human_HG)), rep(FALSE, length(taxids_human_noHG)))) %>%
  mutate(OTU_name=reports$OTU_name[match(OTU_ID, reports$OTU_ID)]) %>% filter(!duplicated(OTU_ID))

human_assoc %<>% merge(OTU_known_human, by.x="species", by.y="OTU_name", all.x=TRUE, all.y=FALSE)

filter(human_assoc, (r2_GTDB_conf>=0.1 & coef_GTDB_conf>0) | (r2_GTDB_noconf>=0.1 & coef_GTDB_noconf>0)) %>%
  arrange(p_GTDB_noconf) ->
  potential_human_spurious

table(contaminated=!is.na(potential_human_spurious$OTU_ID), HG_match=potential_human_spurious$HG_match)

bracken_reports %>% group_by(sample) %>%
  summarise(bacteria=sum(new_est_reads[bacteria]), human=new_est_reads[name=="Homo sapiens"], total=sum(new_est_reads),
            human_prop=human/total, bacteria_prop=bacteria/total,
            human_assoc=sum(new_est_reads[name %in% potential_human_spurious$species]),
            human_assoc_bacteria=sum(new_est_reads[name %in% potential_human_spurious$species & bacteria])) %>%
  mutate(human_assoc_prop=human_assoc / total,
         human_assoc_human_prop=human_assoc / human,
         human_assoc_bacteria_prop = human_assoc_bacteria / bacteria) ->
  human_assoc_abundance

#Reallocate human-associated kraken reads
reports %>% filter(database=="GTDB_conf") %>%
  group_by(sample) %>%
  summarise(human_assoc=sum(reads[OTU_name %in% potential_human_spurious$species])) ->
  kraken_human_assoc_abundance

reports %>% filter(database=="GTDB_conf" & !OTU_ID %in% potential_human_spurious$OTU_ID) %>%
  mutate(reads=reads+ifelse(OTU_ID==9606, kraken_human_assoc_abundance$human_assoc[match(sample, kraken_human_assoc_abundance$sample)], 0)) ->
  kraken_corrected_conf

#Reallocate human associated reads
bracken_reports %<>% filter(!name %in% potential_human_spurious$species) %>%
  mutate(new_est_reads=new_est_reads+ifelse(taxonomy_id==9606, human_assoc_abundance$human_assoc[match(sample, human_assoc_abundance$sample)], 0))

#Count reads in aggregate
filter(reports, OTU_name %in% c("Homo sapiens", "Bacteria", "Fungi", "Viruses", "Archaea", "unclassified") &! sample %in% c("507","508")) %>%
  group_by(OTU_name) %>%
  summarise(reads=sum(cum_reads)) %>%
  arrange(OTU_name) %>%
  mutate(read_prop = reads/sum(reads))

filter(reports, database=="GTDB_conf" & OTU_name %in% c("Homo sapiens", "Bacteria", "Fungi", "Viruses", "Archaea") &! sample %in% c("507","508")) %>%
  left_join(select(ungroup(reads) %>% mutate(total_reads=reads), sample, total_reads), by="sample") %>%
  group_by(OTU_name, database, read_group, group) %>%
  summarise(p25=quantile(cum_reads*100/total_reads,0.25), median=median(cum_reads*100/total_reads), p75=quantile(cum_reads*100/total_reads,0.75)) %>%
  mutate(median=paste0(fnum(median), "% (", fnum(p25), "-", fnum(p75), ")")) %>%
  select(-p25,-p75) %>%
  pivot_wider(id_cols=OTU_name, names_from=group, values_from=median) %>%
  mutate(OTU_name=ordered(OTU_name, levels=major_OTUs_ordered))%>%
  arrange(OTU_name) %>%
  rename(`Organism group` = OTU_name) ->
  org_prop_table

#Sequencing contaminants
reads_per_sample <- group_by(bracken_reports, sample) %>% summarise(reads=sum(new_est_reads))
group_reads_per_sample <- group_by(bracken_reports, sample, read_group) %>%
  summarise(reads=sum(new_est_reads), reads_prop=reads/reads_per_sample$reads[reads_per_sample$sample==sample[1]]) %>%
  mutate(read_group=paste0(tolower(read_group), "_prop")) %>%
  pivot_wider(id_cols=sample, names_from=read_group, values_from=reads_prop)
sample_lookup %<>% merge(group_reads_per_sample, by="sample", all.x=TRUE)
sample_lookup$effective_human_conc <- sample_lookup$human_prop * sample_lookup$concentration
sample_lookup$effective_bacteria_conc <- sample_lookup$bacteria_prop * sample_lookup$concentration
sample_lookup$effective_archaea_conc <- sample_lookup$archaea_prop * sample_lookup$concentration
sample_lookup$effective_fungus_conc <- sample_lookup$fungus_prop * sample_lookup$concentration
sample_lookup$effective_virus_conc <- sample_lookup$virus_prop * sample_lookup$concentration

#Classify contaminants
bracken_reports %>% 
  group_by(sample) %>%
  mutate(prop=new_est_reads/sum(new_est_reads), human_prop=new_est_reads[name=="Homo sapiens"]/sum(new_est_reads)) %>%
  ungroup() %>%
  filter(bacteria) %>%
  select(sample, human_prop, name, prop) %>%
  pivot_wider(id_cols=c(sample, human_prop), names_from=name, values_from=prop, values_fill=0) ->
  contam_matrix_new
contam_matrix_new$plate <- sample_lookup$plate[match(contam_matrix_new$sample, sample_lookup$sample)]
contam_matrix_new$row <- sample_lookup$row[match(contam_matrix_new$sample, sample_lookup$sample)]
contam_matrix_new$col <- sample_lookup$col[match(contam_matrix_new$sample, sample_lookup$sample)]
contam_matrix_new$plate_row <- paste(contam_matrix_new$plate, contam_matrix_new$row, sep="-")
contam_matrix_new$P5_index <- sample_lookup$P5_index[match(contam_matrix_new$sample, sample_lookup$sample)]
contam_matrix_new$P7_index <- sample_lookup$P7_index[match(contam_matrix_new$sample, sample_lookup$sample)]

cl<-makeCluster(n_clust)
clusterExport(cl, "contam_matrix_new")

parLapply(cl, colnames(contam_matrix_new)[-c(1,2,(ncol(contam_matrix_new)-5):ncol(contam_matrix_new))], function(s){
  f <- as.formula(paste0("`", s, "` ~ P5_index + P7_index"))
  f_basic <- as.formula(paste0("`", s, "` ~ 1"))
  glm(f, data=contam_matrix_new[!is.na(contam_matrix_new$P7_index) & !is.na(contam_matrix_new$P5_index),]) -> model
  glm(f_basic, data=contam_matrix_new[!is.na(contam_matrix_new$P7_index) & !is.na(contam_matrix_new$P5_index),]) -> model_basic
  if(length(model$coefficients)>1){
    list(species=s, p=anova(model, model_basic, test="F")$`Pr(>F)`[2], r2=1-(model$deviance/model$null.deviance))
  }else{
    list(species=s, p=NA, r2=NA)
  }
}) -> contam_coef

stopCluster(cl)

contam_coef %<>% bind_rows()
contam_coef$q <- p.adjust(contam_coef$p, method="BH")

contam_coef %<>% arrange(q)

filter(contam_coef, q <= 0.0005) %>% arrange(-q) %>% as.data.frame() %>% pull(species) -> bacteria_tag_sp

#Potential missing contaminants
contam_proteobacteria_genera<-("Afipia
Aquabacterium
Asticcacaulis
Aurantimonas
Beijerinckia
Bosea
Bradyrhizobium
Brevundimonas
Caulobacter
Craurococcus
Devosia
Hoeflea
Mesorhizobium
Methylobacterium
Novosphingobium
Ochrobactrum
Paracoccus
Pedomicrobium
Phyllobacterium
Rhizobium
Roseomonas
Sphingobium
Sphingomonas
Sphingopyxis
Acidovorax
Azoarcus
Azospira
Burkholderia
Comamonas
Cupriavidus
Curvibacter
Comamonas
Duganella
Herbaspirillum
Janthinobacterium
Leptothrix
Limnobacter
Massilia
Methylophilus
Methyloversatilis
Oxalobacter
Pelomonas
Polaromonas
Ralstonia
Schlegelella
Sulfuritalea
Undibacterium
Variovorax
Enhydrobacter
Enterobacter
Escherichia
Nevskia
Pseudoxanthomonas
Psychrobacter
Stenotrophomonas
Xanthomonas" %>% strsplit("\n"))[[1]]

contam_genus_grep_str=paste0("^(", paste(contam_proteobacteria_genera, collapse="|"), ")")

classify_contaminants <- function(this_group, mode="conc", zeroes="replace", additional=character(0), annotate=TRUE, output_proportion=TRUE, p=0.05, exclude=character(0)){
  filter(bracken_reports, read_group==this_group &! name %in% exclude) %>%
    pivot_wider(id_cols=sample, names_from=name, values_from=new_est_read_prop, values_fill=c(new_est_read_prop=0)) ->
    group_matrix
  if(mode=="conc"){
    conc <- sample_lookup[match(group_matrix$sample, sample_lookup$sample), paste0("effective_", tolower(this_group),"_conc")]
    conc <- conc * rowSums(group_matrix[,-1])
  }else{
    conc <- sample_lookup[match(group_matrix$sample, sample_lookup$sample), paste0(tolower(this_group),"_prop")]
  }
  contam_matrix <- as.matrix(group_matrix[!is.na(conc),-1])
  rownames(contam_matrix) <- group_matrix$sample[!is.na(conc)]
  contam_matrix <- contam_matrix[,colSums(contam_matrix)>0]
  contam_matrix <- contam_matrix / rowSums(contam_matrix)
  if(zeroes=="replace")
    contam_matrix %<>% apply(2, function(x){
      x[x==0] <- min(x[x!=0]) / 2
      x
    })
  isContaminant(contam_matrix, conc[!is.na(conc)]) %>% arrange(p) -> contam
  
  #Annotate contaminants
  contam %<>% add_rownames(var="species")
  contam %<>%
    mutate(mean_abundance=colMeans(contam_matrix[,match(species, colnames(contam_matrix))]),
           median_abundance=apply(contam_matrix[,match(species, colnames(contam_matrix))], 2, function(x){
             result<-median(x)
             if(result==0){
               result <- 0.5 * min(x[x!=0]) / 2
             }
             result
           }),
           max_abundance=colMax(contam_matrix[,match(species, colnames(contam_matrix))])) 
  
  contam$p[contam$species %in% additional] <- -1
  
  contam$mode<-paste(mode, zeroes, sep="_")
  
  if(annotate){
    #Annotate bracken and kraken reports
    bracken_reports <<- mutate(bracken_reports, contaminant = case_when(
      !read_group == this_group ~ contaminant,
      database != "GTDB" ~ NA,
      name %in% contam$species[which(contam$p<=p)] ~ TRUE,
      TRUE ~ FALSE
    ))
    reports <<- mutate(reports, contaminant = case_when(
      !read_group == this_group ~ contaminant,
      database != "GTDB" ~ NA,
      OTU_name %in% contam$species[which(contam$p<=p)] ~ TRUE,
      TRUE ~ FALSE
    ))
  }
  if(output_proportion){
    #Contaminant proportion versus estimated bacterial DNA biomass
    rowSums(group_matrix[,contam$species[which(contam$p <= p)]]) -> total_contam_prop
    contam_prop <<- bind_rows(contam_prop, data.frame(read_group=this_group, sample=group_matrix$sample, total_contam_prop=total_contam_prop, conc=conc))
  }
  contam
}

bacteria_contam_z_prop <- classify_contaminants("Bacteria",
                                                mode="prop",
                                                zeroes="noreplace",
                                                annotate=FALSE,
                                                output_proportion=FALSE)
bacteria_contam_z_conc <- classify_contaminants("Bacteria",
                                                mode="conc",
                                                zeroes="noreplace",
                                                annotate=FALSE,
                                                output_proportion=FALSE)
bacteria_contam_z_conc_exclude <- classify_contaminants("Bacteria",
                                                        mode="conc",
                                                        zeroes="noreplace",
                                                        annotate=FALSE,
                                                        output_proportion=FALSE,
                                                        exclude=unique(bacteria_contam_z_prop$species[bacteria_contam_z_prop$p<=0.05], bacteria_tag_sp))
bacteria_contam_z_conc_exclude$mode <- paste0(bacteria_contam_z_conc_exclude$mode, "_exclude")

bind_rows(bacteria_contam_z_conc, bacteria_contam_z_conc_exclude, bacteria_contam_z_prop) %>% mutate(contam=p<=0.05) %>%
  mutate(contaminant=!(p>=0.05 | is.na(p))) %>%
  pivot_wider(id_cols=c(species), names_from=mode, values_from=c(p, contaminant), values_fill=list(p=1, contaminant=FALSE)) ->
  contam_bacteria

contam_bacteria %<>% merge(bacteria_contam_z_prop[,c("species", "mean_abundance", "median_abundance")], by="species")
contam_bacteria %<>% mutate(contaminant_tag=ifelse(species %in% bacteria_tag_sp, TRUE, FALSE),
                                 contaminant_genus=grepl(contam_genus_grep_str, species))

#Manual review of abundant bacteria
filter(bracken_reports, group != "Negative control" & bacteria &! name %in% contam_bacteria$species[contam_bacteria$any_contam]) %>%
  group_by(sample) %>%
  mutate(RA=new_est_reads/sum(new_est_reads)) %>%
  filter(RA >= 0.1) %>%
  group_by(name) %>%
  summarise(n=n()) %>%
  arrange(-n) ->
  suspicious_abundant_bacteria

contam_bacteria %<>% mutate(contaminant_abund=ifelse(species %in% c("Syntrophomonas methylbutyratica", "Sphingobacterium multivorum", "Bifidobacterium breve"), TRUE, FALSE))
contam_bacteria$any_contam_at_all <- select(contam_bacteria, starts_with("contaminant")) %>% rowSums() > 0

table(contam_bacteria$any_contam_at_all)

reports$contaminant <- NA
bracken_reports$contaminant <- NA

bracken_reports <- mutate(bracken_reports, contaminant = case_when(
  !read_group == "Bacteria" ~ contaminant,
  database != "GTDB" ~ NA,
  name %in% contam_bacteria$species[contam_bacteria$any_contam_at_all] ~ TRUE,
  TRUE ~ FALSE
))
reports <- mutate(reports, contaminant = case_when(
  !read_group == "Bacteria" ~ contaminant,
  !database %in% c("GTDB_conf", "GTDB_noconf") ~ NA,
  OTU_name %in% contam_bacteria$species[contam_bacteria$any_contam_at_all] ~ TRUE,
  TRUE ~ FALSE
))

#Contam proportion per sample
bracken_reports %>% group_by(sample, group) %>%
  summarise(total_reads=sum(new_est_reads), human_prop=new_est_reads[name=="Homo sapiens"]/total_reads, bacteria_contam_reads=sum(new_est_reads[contaminant & bacteria])) %>%
  mutate(conc=sample_lookup$concentration[match(sample, sample_lookup$sample)]) %>%
  mutate(contam_prop=bacteria_contam_reads/total_reads) -> contam_prop_total

quantile(contam_prop_total$contam_prop, c(0.25,0.5,0.75))
tapply(contam_prop_total$contam_prop, contam_prop_total$group, quantile, c(0.25,0.5,0.75))
wilcox.test(contam_prop_total$contam_prop[contam_prop_total$group=="KD"], contam_prop_total$contam_prop[contam_prop_total$group=="Febrile"])

#Recalculate estimated read proportion within taxonomic groups
bracken_reports %<>% group_by(sample, read_group) %>%
  mutate(new_est_read_prop=new_est_read_prop/sum(new_est_read_prop)) %>%
  ungroup()

filter(bracken_reports, read_group=="Bacteria") %>%
  pivot_wider(id_cols=sample, names_from=name, values_from=new_est_read_prop, values_fill=c(new_est_read_prop=0)) ->
  group_matrix
rowSums(group_matrix[,contam_bacteria$species[contam_bacteria$any_contam_at_all]]) -> total_contam_prop
conc <- sample_lookup[match(group_matrix$sample, sample_lookup$sample), "effective_bacteria_conc"]
prop <- sample_lookup[match(group_matrix$sample, sample_lookup$sample), "bacteria_prop"]
sample_lookup$bacteria_contam_prop <- contam_prop$total_contam_prop[contam_prop$read_group=="Bacteria"][match(sample_lookup$sample, contam_prop$sample[contam_prop$read_group=="Bacteria"])]

bracken_reports %>% filter(bacteria & !contaminant) %>%
  group_by(sample) %>%
  summarise(reads=sum(new_est_reads)) ->
  bacteria_reads_nocontam

sample_lookup$bacteria_reads_nocontam <- bacteria_reads_nocontam$reads[match(sample_lookup$sample, bacteria_reads_nocontam$sample)]

sample_lookup$total_reads <- reads$reads[match(sample_lookup$sample, reads$sample)]

#SECTION 3
#Descriptive analyses

#Nonpareil curves
npo<-list.files(path="intermediates/nonpareil/", pattern="*.npo", full.names=TRUE)
npo_samples<-stri_extract_first_regex(npo, "(?<=[.][.][/][.][.][/]output[/]nonhuman_trim[/]np_100[/])[^.]+")
npo_groups <- case_when(
  stri_count_regex(npo_samples, "[-]3[0-9]{3}$") > 0 ~ "KD",
  nchar(stri_extract_first_regex(npo_samples, "(?<=[0-9][-])[0-9]+")) == 4 ~ "Febrile",
  npo_samples %in% c("507", "508") ~ "Negative control",
  TRUE ~ "KD")

colours <- ifelse(npo_groups=="KD", "orange", "sky blue")

np_set<-Nonpareil.set(npo, labels=npo_samples, col=colours, plot.opts=list(plot.observed=FALSE, plot.model=TRUE, legend.opts=NULL, curve.alpha=0.3, model.alpha=0.5), correction.factor=TRUE)
legend(x = 1e10, y = 0.7,               
       legend = c("KD", "Febrile"),
       col = c("orange", "sky blue"),
       lty=1,
)

#Bacterial microbiome
sample_lookup$exclude_descriptive <- sample_lookup$sample %in% exclude_qual

bracken_reports %>% filter(bacteria &! sample%in%exclude_qual) %>%
  group_by(sample) %>%
  mutate(RA=new_est_reads / sum(new_est_reads)) %>%
  pivot_wider(id_cols=sample, names_from=name, values_from=RA, values_fill=list(RA=0)) ->
  all_bacteria_table

all_bacteria_matrix <- as.matrix(all_bacteria_table[,-1])
rownames(all_bacteria_matrix) <- all_bacteria_table$sample

bracken_reports %>% filter(bacteria & !contaminant &! sample %in% exclude_qual) %>%
  group_by(sample) %>%
  mutate(RA=new_est_reads / sum(new_est_reads)) %>%
  pivot_wider(id_cols=sample, names_from=name, values_from=RA, values_fill=list(RA=0)) ->
  bacteria_table

bacteria_matrix <- as.matrix(bacteria_table[,-1])
rownames(bacteria_matrix) <- bacteria_table$sample

table(colMax(bacteria_matrix[!rownames(bacteria_matrix) %in% c("507", "508"), ])>=0.001)

sample_lookup$bacteria_diversity <- diversity(bacteria_matrix)[match(sample_lookup$sample, rownames(bacteria_matrix))]
sample_lookup$bacteria_evenness <- evenness(t(bacteria_matrix), index="pielou", zeroes=FALSE)$pielou[match(sample_lookup$sample, rownames(bacteria_matrix))]

quantile(sample_lookup$bacteria_diversity[sample_lookup$group!="Negative control" &! sample_lookup$exclude_descriptive], c(0.25,0.5,0.75))
tapply(sample_lookup$bacteria_diversity[sample_lookup$group!="Negative control" &! sample_lookup$exclude_descriptive], sample_lookup$group[sample_lookup$group!="Negative control" &! sample_lookup$exclude_descriptive], quantile, c(0.25,0.5,0.75))
wilcox.test(sample_lookup$bacteria_diversity[sample_lookup$group=="KD" &! sample_lookup$exclude_descriptive], sample_lookup$bacteria_diversity[sample_lookup$group=="Febrile" &! sample_lookup$exclude_descriptive])

quantile(sample_lookup$bacteria_evenness[sample_lookup$group!="Negative control" &! sample_lookup$exclude_descriptive], c(0.25,0.5,0.75))
tapply(sample_lookup$bacteria_evenness[sample_lookup$group!="Negative control" &! sample_lookup$exclude_descriptive], sample_lookup$group[sample_lookup$group!="Negative control" &! sample_lookup$exclude_descriptive], quantile, c(0.25,0.5,0.75))
wilcox.test(sample_lookup$bacteria_evenness[sample_lookup$group=="KD" &! sample_lookup$exclude_descriptive], sample_lookup$bacteria_evenness[sample_lookup$group=="Febrile" &! sample_lookup$exclude_descriptive])

#Graph table 7
bracken_reports %>% filter(name!="Homo sapiens" & !contaminant &! sample%in%exclude_qual) %>%
  filter(read_group=="Bacteria") %>%
  mutate(Genus=stri_extract_first_regex(name, "[^ _]+")) %>%
  group_by(sample, Genus) %>%
  summarise(new_est_reads=sum(new_est_reads)) %>%
  mutate(`Relative abundance`=new_est_reads / sum(new_est_reads)) %>%
  pivot_wider(id_cols=sample, names_from=Genus, values_from=`Relative abundance`, values_fill=list(`Relative abundance`=0)) %>%
  pivot_longer(-sample, names_to="Genus", values_to="Relative abundance") %>%
  group_by(Genus) %>%
  mutate(median_RA=median(`Relative abundance`)) %>%
  arrange(-median_RA) %>%
  ungroup() %>%
  group_by(-median_RA, Genus) %>%
  mutate(rank=cur_group_id()) %>%
  filter(rank<=10) %>%
  mutate(Genus=ordered(Genus, levels=group_keys(.) %>% pull(Genus))) %>%
  group_by(sample) %>%
  mutate(Group=sample_lookup$group[sample_lookup$sample==sample[1]],
         Group=pub_rename(Group)) %>%
  ungroup() %>%
  filter(Group %in% c("KD", "Febrile control")) %>%
  ggplot(aes(x=Genus, y=`Relative abundance`, fill=Group))+
  geom_boxplot(position="dodge")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_fill_manual(values=palette)

ggsave("PUB_abundant_genera.svg", device="svg", path="output", height=4, width=8)

tree_text<-read_file("resources/GTDB/bac120_r95.sp_labels.tree") %>%
  stri_replace_all_regex("['][0-9]+[.][0-9]+[:]([^']+)[']", "$1") %>%
  stri_replace_all_fixed("; ", "|") %>%
  stri_replace_all_fixed(" ", "_")
all_tree <- read.tree(text=tree_text)
all_tree$tip.label %<>% substring(5, nchar(.) - 1)

colnames(all_bacteria_matrix) %<>% stri_replace_all_fixed(" ", "_")
m_data_all <- all_bacteria_matrix %>% as.data.frame()
contam_tree <- drop.tip(all_tree, which(!all_tree$tip.label %in% colnames(all_bacteria_matrix)), collapse.singles=FALSE)
contam_tree$node.label[grep("^[0-9]+$", contam_tree$node.label)] <- grep("^[0-9]+$", contam_tree$node.label) + length(contam_tree$tip.label)

m_contam <- m_data_all
m_contam[,contam_bacteria$species[!contam_bacteria$any_contam_at_all] %>% stri_replace_all_fixed(" ", "_")] <- 0
m_contam <- m_contam / rowSums(m_contam)

m_noncontam <- m_data_all
m_noncontam[,contam_bacteria$species[contam_bacteria$any_contam_at_all] %>% stri_replace_all_fixed(" ", "_")] <- 0
m_noncontam <- m_noncontam / rowSums(m_noncontam)

all_contam_node <- showNode(tree=contam_tree, only.leaf=FALSE)
sample_contam_lse <- TreeSummarizedExperiment(assays = list(t(m_contam),
                                                            t(m_noncontam)),
                                              rowTree = contam_tree)
sample_contam_tse <- aggValue(x=sample_contam_lse, rowLevel=all_contam_node, FUN=sum)
all_phyla <- grep("^p__", contam_tree$node.label) + length(contam_tree$tip.label)
tse_subset <- subsetByNode(sample_contam_tse, all_phyla)

contam_data <- tse_subset@assays@data@listData[[1]] %>% t()
noncontam_data <- tse_subset@assays@data@listData[[2]] %>% t()
colnames(contam_data) <- paste0("contam:", showNode(contam_tree)[all_phyla] %>% names %>% stri_extract_first_regex("(?<=p__)[^|]+"))
colnames(noncontam_data) <- paste0("noncontam:", showNode(contam_tree)[all_phyla] %>% names %>% stri_extract_first_regex("(?<=p__)[^|]+"))
bind_cols(as.data.frame(contam_data), as.data.frame(noncontam_data)) %>%
  rownames_to_column(var="sample") %>%
  pivot_longer(-sample, names_to=c("category", "phylum"), names_sep=":", values_to="RA") %>%
  mutate(phylum=stri_extract_first_regex(phylum, "^[^_]+")) %>%
  group_by(sample, category, phylum) %>%
  summarise(RA=sum(RA)) ->
  sample_contam

filter(sample_contam, category=="noncontam") %>% group_by(phylum) %>% summarise(median_RA=median(RA)) %>% arrange(-median_RA) %>% pull(phylum) %>% head(n=10) -> top_phyla
filter(sample_contam, category=="noncontam" & phylum==top_phyla[1]) %>% arrange(-RA) %>% pull(sample) -> sample_order
filter(sample_contam, category=="noncontam" & phylum==top_phyla[1]) %>% pull(RA) %>% quantile(c(0,0.25,0.5,0.75,1))

#PUB phylum
filter(sample_contam, category=="noncontam") %>%
  mutate(Group=m_metadata$group[match(sample, rownames(m_metadata))]) %>%
  filter(Group!="Negative control") %>%
  mutate(Group=pub_rename(Group)) %>%
  rename(Phylum=phylum, `Relative abundance` = RA) %>%
  group_by(Phylum) %>%
  mutate(median_RA=median(`Relative abundance`)) %>%
  arrange(-median_RA) %>%
  ungroup() %>%
  group_by(-median_RA, Phylum) %>%
  mutate(rank=cur_group_id()) %>%
  filter(rank<=10) %>%
  mutate(Phylum=ordered(Phylum, levels=group_keys(.) %>% pull(Phylum))) %>%
  ggplot(aes(x=Phylum, y=`Relative abundance`, fill=Group))+
  geom_boxplot(position="dodge")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_fill_manual(values=palette)

ggsave("PUB_abundant_phyla.svg", path="output", device="svg")    

#SECTION 4
#MDS and ANCOVA

#UniFrac distances

bacteria_phylo_otu <- bacteria_matrix %>% apply(2, function(x){
  x[x==0] <- min(x[x!=0]) / 2
  x
})
bacteria_phylo_otu <- bacteria_phylo_otu[,colMax(bacteria_phylo_otu)>=0.001] %>% log10()
colnames(bacteria_phylo_otu) %<>% stri_replace_all_fixed(" ", "_")
tree <- drop.tip(all_tree, which(!all_tree$tip.label %in% colnames(bacteria_phylo_otu)), collapse.singles=TRUE)

bacteria_phylo_otu %<>% otu_table(taxa_are_rows=FALSE)
rownames(sample_lookup) <- sample_lookup$sample
phylo_sample_lookup <- sample_data(sample_lookup)
bacteria_phylo<-phyloseq(bacteria_phylo_otu, phylo_sample_lookup, tree)
bacteria_phylo_distance <- phyloseq::distance(bacteria_phylo, method="wunifrac")

species_mds<-metaMDS(bacteria_phylo_distance, k=2, maxit=500)

#PUB MDS
as.data.frame(species_mds$points) %>%
  rownames_to_column(var="sample") %>%
  mutate(Group=sample_lookup$group[match(sample, sample_lookup$sample)],
         Group=pub_rename(Group)) %>%
  ggplot(aes(x=MDS1, y=MDS2, colour=Group))+
  geom_point()+
  stat_ellipse(type="norm")+
  theme_light()+
  scale_colour_manual(values=palette)

ggsave("PUB MDS.svg", device="svg", path="output", height=4, width=8)

vd2 <- phyloseq::distance(subset_samples(bacteria_phylo, group!="Negative control"), method="wunifrac")

(permanova_m <- adonis2(vd2 ~ country+gender+age+group, data=sample_lookup[match(labels(vd2), sample_lookup$sample), ]))
permanova_m$SumOfSqs[1:4]/permanova_m$SumOfSqs[length(permanova_m$SumOfSqs)]

#Export
m_metadata <- sample_lookup %>%
  filter(group != "Negative control" &! sample %in% exclude_qual) %>%
  select(sample, age, group, gender, country, eth, season, human_prop, concentration, effective_bacteria_conc, effective_human_conc, batch, doo, plate, row, col, pool, lab_ill_day, bacteria_prop, bacteria_reads_nocontam, total_reads, Antibiotic, pcrp, pwbc, ppolys) %>%
  rename(ID=sample)

m_data <- bacteria_matrix %>% as.data.frame()
m_data <- m_data[m_metadata$ID,]

row.names(m_metadata)<-m_metadata$ID
m_metadata %<>% select(-ID)

#SECTION 5
#Differential RA analysis
#Replace bacteria_matrix with all_bacteria_matrix to run analysis with contaminants included
colnames(bacteria_matrix) %<>% stri_replace_all_fixed(" ", "_")

m_metadata <- m_metadata[!rownames(m_metadata) == duplicates[2], ]
m_metadata$concentration[is.na(m_metadata$concentration)] <- median(m_metadata$concentration, na.rm=TRUE)

m_data <- bacteria_matrix %>% as.data.frame()
m_data <- m_data[rownames(m_metadata),]

tree <- drop.tip(all_tree, which(!all_tree$tip.label %in% colnames(bacteria_matrix)), collapse.singles=FALSE)

tree$node.label[grep("^[0-9]+$", tree$node.label)] <- grep("^[0-9]+$", tree$node.label) + length(tree$tip.label)

my_lse <- TreeSummarizedExperiment(assays = list(t(m_data)),
                                   colData = m_metadata,
                                   rowTree = tree)

all_node <- showNode(tree=tree, only.leaf=FALSE)
my_tse <- aggValue(x=my_lse, rowLevel=all_node, FUN=sum)

min_abundance=0.001
min_prevalence = 0.1

#Get all node data (summed up tree)
all_node_data <- my_tse@assays@data@listData[[1]] %>% t()

#Get node IDs to exclude
exclude_nodes <- colnames(all_node_data)[colMax(all_node_data) < min_abundance | colSums(all_node_data >0) < min_prevalence * nrow(all_node_data)] %>% substring(7) %>% as.numeric()

#Create tree table
tree_tbl <- as_tibble(tree)

#Remove rows with nodes to be excluded
tree_tbl_trim <- tree_tbl %>% filter(!parent %in% exclude_nodes &! node %in% exclude_nodes)
#Classify tips (occur only in nodes column, not parent)
tree_tbl_trim$tip <- !tree_tbl_trim$node %in% tree_tbl_trim$parent
#Create indices column
tree_tbl_trim$indices <- NA_integer_
#Give tips new numbers ascending from 1
tree_tbl_trim$indices[tree_tbl_trim$tip] <- 1:length(which(tree_tbl_trim$tip))
#Extract internal nodes and sort by original ID
node_indices <- unique(tree_tbl_trim$parent) %>% sort()
#Create indices for internal nodes by finding in ordered list
tree_tbl_trim$indices[!tree_tbl_trim$tip] <- match(tree_tbl_trim$node[!tree_tbl_trim$tip], node_indices) + length(which(tree_tbl_trim$tip))
indices <- tree_tbl_trim$node[order(tree_tbl_trim$indices)]
#Replace parent node IDs based on matching with nodes
tree_tbl_trim$parent <- tree_tbl_trim$indices[match(tree_tbl_trim$parent, tree_tbl_trim$node)]
#Replace descendent node IDs with new ones
tree_tbl_trim$node <- tree_tbl_trim$indices
tree_trim <- as.phylo(tree_tbl_trim[,1:4])
trim_node_data <- t(all_node_data[,colMax(all_node_data) >= min_abundance  & colSums(all_node_data >0) >= min_prevalence * nrow(all_node_data)])
rownames(trim_node_data) %<>% substring(7) %>%
  as.numeric() %>%
  match(indices) %>%
  paste0("alias_", .)
my_tse_trim <- TreeSummarizedExperiment(assays = list(trim_node_data),
                                        colData = m_metadata,
                                        rowTree = tree_trim)

all_node <- showNode(tree=tree_trim, only.leaf=FALSE)

this_data <- m_metadata

cl <- makeCluster(n_clust)
clusterEvalQ(cl, {
  library(dplyr)
  library(TreeSummarizedExperiment)
  library(treeclimbR)})
clusterExport(cl, c("my_tse_trim", "min_prevalence", "this_data", "min_abundance"))

parLapply(cl, all_node, function(x){
  cat(x)
  cat(" ")
  tse_subset <- subsetByNode(my_tse_trim, x)
  this_RA <- tse_subset@assays@data@listData[[1]] %>% t()
  if(length(which(this_RA!=0)) / nrow(this_RA) < min_prevalence | max(this_RA) < min_abundance){
    #Should never happen
    NULL
    stop("OTU has too low abundance or prevalence. Shouldn't happen!")
  }else{
    this_RA[this_RA==0] <- min(this_RA[this_RA!=0])/2
    this_data$RA <- this_RA[match(rownames(this_data), rownames(this_RA))]
    model<-glm(log(RA) ~ group + country + age + gender, data=this_data) %>% summary()
    model_sens<-glm(log(RA) ~ group + country + age + gender, data=this_data[this_data$bacteria_contam_prop < 0.1,]) %>% summary()
    c<-list(mode=rep("continuous", nrow(model$coefficients)-1), node=rep(x, nrow(model$coefficients)-1), covariate=row.names(model$coefficients)[-1], coef=model$coefficients[,"Estimate"][-1], p=model$coefficients[,"Pr(>|t|)"][-1], p_sens=model_sens$coefficients[1:(nrow(model$coefficients)),"Pr(>|t|)"][-1], note=NA_character_)
    c_sens<-list(mode=rep("continuous_sens", nrow(model_sens$coefficients)-1), node=rep(x, nrow(model_sens$coefficients)-1), covariate=row.names(model_sens$coefficients)[-1], coef=model_sens$coefficients[,"Estimate"][-1], p=model_sens$coefficients[,"Pr(>|t|)"][-1], note=NA_character_)
    bind_rows(c, c_sens)
  }
}) -> results

results %<>% bind_rows()
results$direction <- sign(results$coef)

select(results, mode, covariate) %>% unique() ->
  combs

clusterExport(cl, c("results", "combs", "tree_trim"))

parLapply(cl, 1:nrow(combs), function(i){
  these_results <- filter(results, mode==combs$mode[i] & covariate==combs$covariate[i])
  
  my_cand <- getCand(tree = rowTree(my_tse_trim), score_data = these_results, 
                     node_column = "node", p_column = "p",
                     threshold = 0.05,
                     sign_column = "direction", message = TRUE)
  
  my_candL <- my_cand$candidate_list
  
  # evaluate candidates
  my_best <- evalCand(tree = tree_trim, levels = my_cand$candidate_list, 
                      score_data = these_results, node_column = "node",
                      p_column = "p", sign_column = "direction")
  infoCand(object = my_best)
  
  my_best
}) -> all_best

stopCluster(cl)

all_outB <- lapply(all_best, topNodes, n = Inf, p_value = 1) %>% bind_rows()

all_outB$name <- convertNode(tree_trim, all_outB$node, use.alias=FALSE)

treeMat=matTree(tree)
treeMatNames<-matrix((c(tree$tip.label, tree$node.label)[treeMat]), ncol=ncol(treeMat))

treeMat_trim=matTree(tree_trim)
treeMatNames_trim<-matrix((c(tree_trim$tip.label, tree_trim$node.label)[treeMat_trim]), ncol=ncol(treeMat_trim))

all_outB$descendants <- lapply(all_outB$node, function(n){
  n2<-indices[n]
  matches<-which(treeMat==n2, arr.ind=TRUE)
  apply(matches, 1, function(x){
    names<-treeMatNames[x[1], x[2]:1]
    names <- names[!grepl("^[0-9]+$", names)]
    names[1]
  }) %>% unique() -> full
  
  matches<-which(treeMat_trim==n, arr.ind=TRUE)
  apply(matches, 1, function(x){
    names<-treeMatNames_trim[x[1], x[2]:1]
    names <- names[!grepl("^[0-9]+$", names)]
    names[1]
  }) %>% unique() -> partial
  
  extended <- full[!full %in% partial]
  
  result <- paste(partial, collapse=",")
  if(length(extended)>0)
    result %<>% paste0(" (", paste(extended, collapse=","), ")")
  result
}) %>% unlist()


all_outB$n_descendants <- lapply(all_outB$node, function(n){
  nrow(which(treeMat_trim==n, arr.ind=TRUE))
}) %>% unlist()

all_outB$descendants %<>% stri_replace_all_regex("[|][^,]+", "")

tse_subset <- subsetByNode(my_tse_trim, all_outB$node)
mean_RA <- rowMeans(tse_subset@assays@data@listData[[1]])
mean_RA_KD <- rowMeans(tse_subset@assays@data@listData[[1]][,rownames(this_data)[this_data$group=="KD"]])
mean_RA_Febrile <- rowMeans(tse_subset@assays@data@listData[[1]][,rownames(this_data)[this_data$group=="Febrile"]])
median_RA <- rowMeans(tse_subset@assays@data@listData[[1]])
all_outB$mean_RA <- mean_RA[match(paste0("alias_", all_outB$node), names(mean_RA))]
all_outB$mean_RA_KD <- mean_RA_KD[match(paste0("alias_", all_outB$node), names(mean_RA_KD))]
all_outB$mean_RA_Febrile <- mean_RA_Febrile[match(paste0("alias_", all_outB$node), names(mean_RA_Febrile))]
all_outB$median_RA <- median_RA[match(paste0("alias_", all_outB$node), names(median_RA))]

all_outB %<>% mutate(level=stri_extract_first_regex(descendants, "[a-z](?=__)"), level=case_when(
  level=="p" ~ "Phylum",
  level=="g" ~ "Genus",
  level=="c" ~ "Class",
  level=="o" ~ "Order",
  level=="f" ~ "Family",
  is.na(level) ~ "Species"
), descendants=stri_replace_all_regex(descendants,"[a-z]__", "") %>% stri_replace_all_regex("_(?=[a-z])", " ") %>% stri_replace_all_fixed(",", ", "))

all_outB %<>% group_by(mode, covariate) %>%
  mutate(adj.p_sens=ifelse(mode%in%c("continuous", "dichotomous"), p.adjust(p_sens, method="BH"), NA))

table(all_outB$covariate[all_outB$adj.p <= 0.2], sign(all_outB$coef[all_outB$adj.p <= 0.2]), all_outB$mode[all_outB$adj.p <= 0.2])
table(all_outB$covariate[all_outB$adj.p <= 0.05], sign(all_outB$coef[all_outB$adj.p <= 0.05]), all_outB$mode[all_outB$adj.p <= 0.05])

table(all_outB$covariate, all_outB$level)

mutate(all_outB, row_num=1:n()) %>%
  ungroup() %>%
  filter(direction==1, mode=="continuous", covariate=="groupKD", adj.p<=0.2) %>%
  arrange(-coef) %>%
  select(row_num, descendants, level, coef, adj.p, mean_RA_KD, mean_RA_Febrile, adj.p_sens_2) %>%
  mutate(coef=fnum(coef, 2),
         adj.p=fnum(adj.p,2),
         mean_RA_KD=fnum(mean_RA_KD*100,2),
         mean_RA_Febrile=fnum(mean_RA_Febrile*100,2),
         adj.p_sens_2=fnum(adj.p_sens_2,2)) %>%
  relocate(level, .after=descendants) %>%
  rename(Organisms=descendants, Coefficient=coef, `Q value`=adj.p, `Mean RA Febrile`=mean_RA_Febrile, `Mean RA KD`=mean_RA_KD, `Taxonomic level`=level, `Q value sensitivity`=adj.p_sens_2) %>%
  select(-row_num) %>%
  write.table("", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mutate(all_outB, row_num=1:n()) %>%
  ungroup() %>%
  filter(direction==-1, mode=="continuous", covariate=="groupKD", adj.p<=0.2) %>%
  arrange(coef) %>%
  select(row_num, descendants, level, coef, adj.p, mean_RA_KD, mean_RA_Febrile, adj.p_sens_2) %>%
  mutate(coef=fnum(coef, 2),
         adj.p=fnum(adj.p,2),
         mean_RA_KD=fnum(mean_RA_KD*100,2),
         mean_RA_Febrile=fnum(mean_RA_Febrile*100,2),
         adj.p_sens_2=fnum(adj.p_sens_2,2)) %>%
  relocate(level, .after=descendants) %>%
  rename(Organisms=descendants, Coefficient=coef, `Q value`=adj.p, `Mean RA Febrile`=mean_RA_Febrile, `Mean RA KD`=mean_RA_KD, `Taxonomic level`=level, `Q value sensitivity`=adj.p_sens_2) %>%
  select(-row_num) %>%
  write.table("", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

mutate(all_outB, row_num=1:n()) %>%
  filter(direction==1, mode=="continuous", covariate=="groupKD", adj.p<=0.2) %>%
  arrange(-coef) %>% ungroup() %>%
  select(node, descendants) -> tmp

nodes <- tmp$node
names(nodes) <- tmp$descendants

tse_subset <- subsetByNode(my_tse_trim, nodes)
these_RA <- tse_subset@assays@data@listData[[1]] %>% t()
these_RA %<>% as.data.frame()
colnames(these_RA) <- c(tree_trim$tip.label,tree_trim$node.label)[as.numeric(substring(colnames(these_RA),7,100))]
colnames(these_RA) %<>% stri_replace_all_fixed("_", " ")
these_RA$Abiotrophia <- these_RA$`Abiotrophia defectiva` + these_RA$`Abiotrophia sp001815865`
these_RA %<>% mutate(across(everything(), ~ifelse(.==0, min(.[.!=0])/2, .)))
these_RA$age <- this_data$age[match(rownames(these_RA), rownames(this_data))]
these_RA$sample=rownames(these_RA)
these_RA$Group=m_metadata$group[match(these_RA$sample, rownames(m_metadata))]

glm(log(`Abiotrophia`) ~ Group + country + age + gender, data=these_RA) %>% summary()
glm(log(`Abiotrophia`) ~ Group + country + ifelse(age<1.5,age,1.5) + gender, data=these_RA) %>% summary()

glm(log(`Neisseria sp000090875`) ~ Group + country + age + gender, data=these_RA) %>% summary()
glm(log(`Neisseria sp000090875`) ~ Group + country + ifelse(age<1.5,age,1.5) + gender, data=these_RA) %>% summary()

glm(log(`Prevotella oris`) ~ Group + country + age + gender, data=these_RA) %>% summary()
glm(log(`Prevotella oris`) ~ Group + country + gender, data=these_RA) %>% summary()
glm(log(`Rothia dentocariosa`) ~ Group + country + age + gender, data=these_RA) %>% summary()
glm(log(`Rothia dentocariosa`) ~ Group + country + gender, data=these_RA) %>% summary()
glm(log(`Streptococcus mitis AD`) ~ Group + country + age + gender, data=these_RA) %>% summary()
glm(log(`Streptococcus mitis AD`) ~ Group + country + gender, data=these_RA) %>% summary()
glm(log(`Streptococcus pseudopneumoniae O`) ~ Group + country + age + gender, data=these_RA) %>% summary()
glm(log(`Streptococcus pseudopneumoniae O`) ~ Group + country + gender, data=these_RA) %>% summary()
glm(log(`Streptococcus sanguinis G`) ~ Group + country + age + gender, data=these_RA) %>% summary()
glm(log(`Streptococcus sanguinis G`) ~ Group + country + gender, data=these_RA) %>% summary()
glm(log(`Streptococcus sp000831085`) ~ Group + country + age + gender, data=these_RA) %>% summary()
glm(log(`Streptococcus sp000831085`) ~ Group + country + gender, data=these_RA) %>% summary()

these_RA %<>% pivot_longer(-c(sample, age, Group, country, gender), names_to="Species", values_to="RA")
mutate(these_RA, Group=pub_rename(Group)) %>%
  ggplot(aes(x=age, y=RA, colour=Group))+
  geom_point()+
  geom_smooth(se=FALSE)+
  scale_y_log10()+
  theme_light()+
  xlab("")+
  ylab("Relative abundance")+
  facet_wrap(~Species, ncol=3)+
  coord_cartesian(xlim=c(0,8))+
  scale_x_continuous(breaks=1:15)+
  xlab("Age")+
  scale_colour_manual(values=palette)

ggsave("PUB_candidate_microbe_age_associations.svg", device="svg", path="output", height=8, width=9)

filter(these_RA, Species %in% c("Abiotrophia", "Lautropia mirabilis", "Rothia dentocariosa")) %>%
  mutate(Group=pub_rename(Group)) %>%
  ggplot(aes(x=age, y=RA, colour=Group))+
  geom_point()+
  geom_smooth(se=FALSE)+
  scale_y_log10()+
  theme_light()+
  xlab("")+
  ylab("Relative abundance")+
  facet_wrap(~Species, ncol=3)+
  coord_cartesian(xlim=c(0,8))+
  scale_x_continuous(breaks=1:15)+
  xlab("Age")+
  scale_colour_manual(values=palette)

ggsave("PUB_selected_microbe_age_associations.svg", device="svg", path="output", height=6, width=8)

assoc_nodes<-c(Abiotrophia=1872, `Lautropia mirabilis`=41, `Rothia dentocariosa`=514)

tse_subset <- subsetByNode(my_tse_trim, assoc_nodes)
these_RA <- tse_subset@assays@data@listData[[1]] %>% t()
colnames(these_RA) <- names(assoc_nodes)[match(as.numeric(substring(colnames(these_RA),7,100)), assoc_nodes)]
these_RA %<>% apply(2, function(x){
  x[x==0] <- min(x[x!=0])/2
  x
})
these_RA %<>% as.data.frame()
these_RA$age <- this_data$age[match(rownames(these_RA), rownames(this_data))]
these_RA$sample=rownames(these_RA)
these_RA$Group=m_metadata$group[match(these_RA$sample, rownames(m_metadata))]

lm(log10(`Abiotrophia`) ~ log10(`Lautropia mirabilis`), data=these_RA) %>% summary()
lm(log10(`Abiotrophia`) ~ log10(`Rothia dentocariosa`), data=these_RA) %>% summary()
lm(log10(`Lautropia mirabilis`) ~ log10(`Rothia dentocariosa`), data=these_RA) %>% summary()

these_RA$A_residual <- lm(log10(`Abiotrophia`) ~ ifelse(age<1.5,age,1.5), data=these_RA)$residuals
these_RA$Lm_residual <- lm(log10(`Lautropia mirabilis`) ~ ifelse(age<1.5,age,1.5), data=these_RA)$residuals
these_RA$Rd_residual <- lm(log10(`Rothia dentocariosa`) ~ age, data=these_RA)$residuals
these_RA$All_residuals <- these_RA$A_residual + these_RA$Lm_residual + these_RA$Rd_residual

glm(Group=="KD" ~ A_residual + Lm_residual + Rd_residual, data=these_RA, family=binomial) %>% summary()
glm(Group=="KD" ~ All_residuals, data=these_RA, family=binomial) %>% summary()

these_RA %>%
  mutate(Group=pub_rename(Group)) %>%
  ggplot(aes(y=All_residuals, x=Group, colour=Group, fill=Group))+
  geom_violin()+
  labs(y="Summed residual relative abundances")+
  scale_fill_manual(values=palette)

ggsave("PUB_summed_residuals_by_group.svg", device="svg", path="output", height=4, width=4)

roc <- roc(these_RA$Group=="KD", these_RA$All_residuals)
ggroc(roc) +
  theme_minimal()+
  labs(x="Specificity", y="Sensitivity")

ggsave("PUB_ROC.svg", device="svg", path="output", height=4, width=4)

#SECTION 6
#PanPhlAn

organism <- "Lautropia"
#or
organism <- "Abiotrophia"
#or
organism <- "Rothia_dentocariosa"
panphlan_file <- paste0("intermediates/panphlan/", organism, "_profile_stringent.tsv")

panphlan_org <- organism
panphlan <- read.table(panphlan_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
colnames(panphlan) %<>% stri_extract_first_regex("^[^.]+")
row.names(panphlan) <- panphlan[,1]
panphlan <- panphlan[,-1] %>% t()
panphlan <- panphlan[!rownames(panphlan) %in% c(exclude_qual, duplicates[2]), ]

panphlan_df <- as.data.frame(panphlan) %>%
  rownames_to_column(var="sample") %>%
  mutate(
    reference=grepl("^REF", sample),
    group=case_when(
      grepl("bin$", sample) ~ "MAG",
      grepl("^REF", sample) ~ "Reference",
      TRUE ~ sample_lookup$group[match(sample, sample_lookup$sample)]),
    country=sample_lookup$country[match(sample, sample_lookup$sample)])

panphlan_df %<>% select(-starts_with("PC"))
panphlan_df %<>% bind_cols(pca_df[match(panphlan_df$sample, rownames(pca_df)),])

filter(panphlan_df, group %in% c("KD", "Febrile")) %>%
  rename(Group=group) %>%
  mutate(Group=pub_rename(Group)) %>%
  ggplot(aes(x=PC1, y=PC2, colour=Group))+
  geom_point()+
  scale_colour_manual(values=palette)

ggsave(paste0("PUB_panphlan_", organism, ".svg"), device="svg", path="output", height=4, width=6)

filter(panphlan_df, group=="KD") %>% select(contains("UniRef")) %>% rowSums() %>% median() -> median_KD
filter(panphlan_df, group=="Febrile") %>% select(contains("UniRef")) %>% rowSums() %>% median() -> median_Febrile
ratio <- median_KD / median_Febrile

panphlan_df %>% select(sample, group, country, starts_with("UniRef")) %>%
  pivot_longer(-c(sample,group, country), names_to="COG", values_to="presence") %>%
  group_by(COG) %>%
  summarise(n_KD=round(length(which(group=="KD" & presence==1))/ratio,0), d_KD=length(which(group=="KD")), n_Febrile=length(which(group=="Febrile" & presence==1)), d_Febrile=length(which(group=="Febrile"))) %>%
  mutate(n_KD=ifelse(n_KD>d_KD, d_KD, n_KD)) %>%
  filter(n_KD+n_Febrile >= 10) %>%
  group_by(COG) %>%
  mutate(p=fisher.test(matrix(c(n_KD, d_KD-n_KD, n_Febrile, d_Febrile-n_Febrile), ncol=2))$p.value,
         KD_vs_febrile=n_KD*d_Febrile/(d_KD*n_Febrile)) %>%
  arrange(p) %>%
  ungroup() %>%
  mutate(q=p.adjust(p, method="BH")) ->
  panphlan_diff

table(KD_up=panphlan_diff$KD_vs_febrile>1, sig=panphlan_diff$p<=0.05)
table(KD_up=panphlan_diff$KD_vs_febrile>1, sig=panphlan_diff$q<=0.05)
table(KD_up=panphlan_diff$KD_vs_febrile>1, sig=panphlan_diff$q<=0.2)

panphlan_genomic <- read.table(paste0("intermediates/panphlan/", panphlan_org, "_pangenome.tsv"), sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(panphlan_genomic) <- c("COG", "unknown", "genome", "contig", "start", "end")

filter(panphlan_diff, p<=0.05) %>%
  pull(COG) -> sig_COGs

filter(panphlan_genomic) %>% group_by(contig) %>%
  mutate(sig=COG %in% sig_COGs) %>%
  mutate(n_genes=length(which(sig))) %>%
  arrange(-n_genes, contig, start) %>%
  mutate(gene_n=1:n()) %>%
  ungroup() %>%
  mutate(contig_n=as.integer(fct_inorder(contig))) -> panphlan_gene_pos

gene_runs <- rle(paste(panphlan_gene_pos$genome, panphlan_gene_pos$contig, panphlan_gene_pos$COG))
gene_run_starts <- c(1, cumsum(gene_runs$lengths))
gene_run_starts <- gene_run_starts[-length(gene_run_starts)]
mapply(function(start, length){
  seq(from=start+1, to=start+length-1)
}, gene_run_starts[gene_runs$lengths>1], gene_runs$lengths[gene_runs$lengths>1], SIMPLIFY=FALSE) %>% unlist() ->
  run_censor

panphlan_gene_pos <- panphlan_gene_pos[-run_censor, ]

rle(panphlan_gene_pos$sig) -> panphlan_sig_gene_rle

table(panphlan_sig_gene_rle$lengths[panphlan_sig_gene_rle$value])

panphlan_df %>% select(sample, group, starts_with("UniRef")) %>%
  filter(group %in% c("KD", "Febrile")) %>%
  pivot_longer(-c(sample, group), names_to="COG", values_to="present") %>%
  group_by(COG) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(total_genes=sum(present)) -> panphlan_long

panphlan_gene_pos %<>% merge(panphlan_long, by="COG", all.x=TRUE, sort=FALSE)

#filter(panphlan_gene_pos) -> panphlan_test_data

panphlan_gene_pos %>%
  select(genome, contig, COG, gene_n) %>%
  distinct() %>%
  filter(!duplicated(COG)) %>%
  mutate(ID=paste(genome, contig, COG, gene_n)) ->
  distinct_COGs

mutate(panphlan_gene_pos, ID=paste(genome, contig, COG, gene_n)) %>%
  filter(ID %in% distinct_COGs$ID) ->
  panphlan_counts

panphlan_gene_pos %<>% select(-sample, -group, -present, -total_genes)

d_KD=length(which(panphlan_df$group=="KD"))
d_Febrile=length(which(panphlan_df$group=="Febrile"))

cl<-makeCluster(n_clust)
clusterEvalQ(cl, {library(dplyr)
  library(magrittr)
  library(forcats)})
clusterExport(cl, c("panphlan_gene_pos", "distinct_COGs", "panphlan_counts", "sample_lookup", "d_KD", "d_Febrile", "ratio"))

g<-"all"
clusterExport(cl, c("g"))
parLapply(cl, 1:500, function(x){
  cat(x)
  if(x > 1){
    panphlan_counts$sample %<>% as.factor() %>% fct_shuffle()
    levels(panphlan_counts$sample) %<>% sort()
    panphlan_counts$sample %<>% as.character()
    panphlan_counts$group <- sample_lookup$group[match(panphlan_counts$sample, sample_lookup$sample)]
  }
  if(g == "all")
    g<-unique(panphlan_gene_pos$genome)
  
  panphlan_counts %>% filter(genome %in% g) %>%
    group_by(genome, contig, COG, gene_n, contig_n) %>%
    summarise(mean_diff=mean(present[group=="KD"])-mean(present[group=="Febrile"]),
              n_KD=round(length(which(present[group=="KD"]==1))/ratio, 0),
              n_KD=ifelse(n_KD>d_KD, d_KD, n_KD),
              n_Febrile=length(which(present[group=="Febrile"]==1)),
              p=fisher.test(matrix(c(n_KD, d_KD-n_KD, n_Febrile, d_Febrile-n_Febrile), ncol=2))$p.value) %>%
    mutate(log_p=-log10(p) * sign(mean_diff),
           cluster=sign(mean_diff)*(abs(log_p)>=-log10(0.05))*contig_n) %>%
    ungroup() %>%
    select(-genome, -contig, -gene_n, -contig_n) ->
    panphlan_stats
  
  filter(panphlan_gene_pos, genome%in%g) %>%
    group_by(genome, contig, COG, gene_n, contig_n) %>%
    summarise() %>%
    merge(panphlan_stats, by=c("COG"), all.x=TRUE, sort=FALSE) %>%
    arrange(genome, gene_n) -> tmp
  
  rle(tmp$cluster) -> rle_cluster
  
  cluster_starts<-c(1,cumsum(rle_cluster$lengths)+1)
  cluster_starts <- cluster_starts[-length(cluster_starts)]
  
  mapply(function(start,length){
    sum(abs(tmp$log_p[seq(from=start, to=start+length-1)]))
  }, cluster_starts[rle_cluster$values!=0], rle_cluster$lengths[rle_cluster$values!=0], SIMPLIFY=FALSE) %>% unlist() ->
    cluster_masses
  
  if(x==1){
    list(tmp, cluster_starts[rle_cluster$values!=0], rle_cluster$lengths[rle_cluster$values!=0], cluster_masses, max(cluster_masses))
  }else{
    max(cluster_masses)
  }
}) -> cluster_mass

panphlan_genome_results <- cluster_mass[[1]][[1]]
cluster_starts <- cluster_mass[[1]][[2]]
cluster_lengths <- cluster_mass[[1]][[3]]
cluster_masses <- cluster_mass[[1]][[4]]
cluster_mass[[1]] <- cluster_mass[[1]][[5]]
cluster_mass %<>% unlist()
cluster_mass_threshold <- quantile(cluster_mass[-1],0.95)

cluster_masses[which(cluster_masses >= cluster_mass_threshold)]
cluster_starts[cluster_masses >= cluster_mass_threshold]
cluster_lengths[cluster_masses >= cluster_mass_threshold]
range(panphlan_genome_results$gene_n)
mapply(function(s, l){
  seq(from=s, to=s+l-1)
},cluster_starts[cluster_masses >= cluster_mass_threshold], cluster_lengths[cluster_masses >= cluster_mass_threshold], SIMPLIFY=FALSE) %>% unlist() ->
  cluster_rows
panphlan_genome_results[cluster_rows,]

stopCluster(cl)

panphlan_genome_results %<>% ungroup() %>% mutate(same_set=ifelse(c(1, diff(gene_n))==1, 0, 1),
                                                  cluster_index=cumsum(same_set)+1) %>%
  select(-same_set, -cluster)

split(panphlan_genome_results$COG, panphlan_genome_results$cluster_index) ->
  contiguous_gene_sets

remove <- FALSE
remain <- 1:length(contiguous_gene_sets)
for(i in 1:length(contiguous_gene_sets)){
  others <- remain[remain!=i]
  for(j in others){
    if(all(contiguous_gene_sets[[i]] %in% contiguous_gene_sets[[j]])){
      remove<-TRUE
      break
    }
  }
  if(remove){
    remain <- remain[remain!=i]
    remove <- FALSE
  }
}

contiguous_gene_sets[remain] %>% unlist() %>% unique() ->
  sig_COGs

panphlan_genome_results %<>% filter(cluster_index %in% remain)

COG_groups <- data.frame(COG=contiguous_gene_sets[remain] %>% unlist(), cluster=rep(1:length(remain), lapply(contiguous_gene_sets[remain], length) %>% unlist()), stringsAsFactors=FALSE)

reduced_contiguous_gene_sets <- contiguous_gene_sets[remain]

#SECTION 7
#StrainPhlAn

plots <- list()

for(organism in c("Abiotrophia_defectiva", "Lautropia_mirabilis", "Rothia_dentocariosa")){
  distmat <- read.table(paste0("intermediates/strainphlan/s__", organism, ".StrainPhlAn3_concatenated.dist"), skip=8, fill=TRUE, sep="\t", strip.white=TRUE)
  # remove the first column of the data as it is blank
  distmat[1] <- NULL
  # get the header as the last column of the data as a character vector
  header <- lapply(distmat[,ncol(distmat)], as.character) %>% stri_extract_first_regex("[^ ]+")
  # remove the last column from the data as it has been stored as a header
  distmat[ncol(distmat)] <- NULL
  # remove the current last column from the data as it is blank
  distmat[ncol(distmat)] <- NULL
  # add the sample names to the columns and rows of the data matrix
  rownames(distmat) <- unlist(header)
  colnames(distmat) <- unlist(header)
  # make symmetric, add lower triangle to upper triangle
  distmat[lower.tri(distmat)] <- t(distmat)[lower.tri(distmat)]
  which(is.nan(as.matrix(distmat)), arr.ind=TRUE) %>% as.vector() %>% table() ->
    nan_counts
  which(nan_counts == sum(nan_counts)/2) %>% names() %>% as.numeric() -> remove
  if(length(remove) > 0){
    distmat <- distmat[-remove, -remove]
  }
  distmat <- distmat[which(rownames(distmat)!=duplicates[2]), which(rownames(distmat)!=duplicates[2])]
  distmat <- distmat[rownames(distmat) %in% rownames(m_metadata) | grepl("^G", rownames(distmat)), colnames(distmat) %in% rownames(m_metadata) | grepl("^G", colnames(distmat))]
  distmat[!grepl("^G", rownames(distmat)),!grepl("^G", colnames(distmat))] %>% as.dist() -> dist_proper
  adonis2(dist_proper ~ group, data=m_metadata[match(labels(dist_proper), rownames(m_metadata)), ])
  e.sir.pcoa <- cmdscale(distmat, eig = T )
  # variance explained 
  variance <- head(eigenvals(e.sir.pcoa)/sum(eigenvals(e.sir.pcoa)))
  x_variance <- as.integer(variance[1]*100)
  y_variance <- as.integer(variance[2]*100)
  sum(variance[1:2])
  # get scores for plotting
  e.sir.scores <- as.data.frame( e.sir.pcoa$points )
    e.sir.scores.meta <- e.sir.scores
  # change colnames
  colnames(e.sir.scores.meta)[c(1,2)] <- c( "PCo1", "PCo2")
  e.sir.scores.meta$sample <- rownames(e.sir.scores.meta)
  e.sir.scores.meta$Group <- m_metadata$group[match(rownames(e.sir.scores.meta), rownames(m_metadata))]
  e.sir.scores.meta$Group[grep("^G", e.sir.scores.meta$sample)] <- "Reference"
  # plot ordination
  plots[[length(plots)+1]] <- mutate(e.sir.scores.meta, Group=pub_rename(Group)) %>% ggplot(aes(PCo1, PCo2, color=Group)) + 
    geom_point(size = 4, alpha = 0.75) + theme_classic() + 
    theme(axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
          axis.ticks = element_blank(), axis.text = element_blank()
    )+
    xlab(paste0("PCo1 (",x_variance,"% variance)")) + ylab(paste0("PCo2 (",y_variance,"% variance)"))+
    scale_colour_manual(values=palette)
  apply(distmat[-grep("^G", rownames(distmat)),grep("^G", colnames(distmat))], 1, min) %>% quantile(c(0, 0.05,0.25,0.5,0.75,0.95,1))
  for(i in 1:ncol(distmat)){
    distmat[i,i]<-NA
  }
  unlist(distmat[-grep("^G", rownames(distmat)),-grep("^G", colnames(distmat))]) %>% quantile(c(0, 0.05,0.25,0.5,0.75,0.95,1), na.rm=TRUE)
}
leg <- get_legend(plots[[1]])
plots %<>% lapply(function(x){
  x+theme(legend.position="none")
})
plot_grid(plots[[1]], plots[[2]], plots[[3]], nrow=1, ncol=3) %>%
  plot_grid(leg, ncol=2, rel_widths=c(1,0.2))
ggsave("PUB_strainphlan_pca.svg", device="svg", path="output", width=10, height=5)

#SECTION 8
#Assembly based analayis of representative gene clusters

qa_files <- list.files(path="intermediates/assembly/checkm/", full.names=TRUE, recursive=FALSE)
qa <- lapply(qa_files, function(x){
  t<-read.table(x, sep="\t", header=TRUE, stringsAsFactors=FALSE, comment.char="")
  t$sample <- stri_extract_first_regex(x, "[0-9]+[-][0-9]+([-][0-9])?(Fb|[-]Fb)?")
  t
})
qa %<>% bind_rows()

qa$complete_90 <- qa$Completeness >= 90
qa$low_contam <- qa$Contamination <= 10
table(complete_50=qa$Completeness >= 50, contam_10=qa$Contamination < 10, !qa$sample == exclude_qual)
table(complete_90=qa$Completeness >= 90, contam_10=qa$Contamination < 5, !qa$sample== exclude_qual)

contig_taxonomy_kraken<-read.table("intermediates/assembly/contig_kraken.txt", sep="\t", stringsAsFactors=FALSE, header=FALSE)
colnames(contig_taxonomy_kraken) <- c("sample_contig", "name_taxid", "length", "confidence", "RTL_confidence")
contig_taxonomy_kraken %<>%
  mutate(sample=stri_extract_first_regex(sample_contig, "^[^.]+"),
         contig=stri_extract_first_regex(sample_contig, "(?<=[.]).+"),
         species=stri_extract_first_regex(name_taxid, "^.+(?=[ ][(])"),
         taxid=as.numeric(stri_extract_first_regex(name_taxid, "(?<=taxid[ ])[0-9]+")),
         genus=stri_extract_first_regex(species, "^[^_ ]+"))

contig_taxonomy_kraken %<>% filter(!sample == exclude_qual)

class_OTUs %<>% mutate(read_group = case_when(
  bacteria ~ "Bacteria",
  fungus ~ "Fungus",
  virus ~ "Virus",
  archaea ~ "Archaea",
  OTU_ID==9606 ~ "Human",
  OTU_ID==0 ~ "Unclassified",
  TRUE ~ "Other"
))
contig_taxonomy_kraken$read_group <- class_OTUs$read_group[match(contig_taxonomy_kraken$taxid, class_OTUs$OTU_ID)]
table(contig_taxonomy_kraken$read_group, useNA="ifany")

select(reports, OTU_ID, OTU_name, level, rank) %>% unique() ->
  levels_ranks

contig_taxonomy_kraken$level <- levels_ranks$level[match(contig_taxonomy_kraken$taxid, levels_ranks$OTU_ID)]
contig_taxonomy_kraken$rank <- levels_ranks$rank[match(contig_taxonomy_kraken$taxid, levels_ranks$OTU_ID)]

contig_bins <- read.table("intermediates/assembly/contig_bins.txt", sep="\t" , header=FALSE, stringsAsFactors=FALSE)
contig_taxonomy_kraken$bin <- contig_bins$V1[match(contig_taxonomy_kraken$sample_contig, contig_bins$V2)]

contam_bacteria$OTU_ID <- levels_ranks$OTU_ID[match(contam_bacteria$species, levels_ranks$OTU_name)]

contig_taxonomy_kraken$species_contaminant <- contig_taxonomy_kraken$taxid %in% contam_bacteria$OTU_ID[contam_bacteria$any_contam_at_all]

filter(bracken_reports, group != "Negative control" & bacteria) %>%
  mutate(genus=stri_extract_first_regex(name, "^[^ ]+"), contaminant=name %in% contam_bacteria$species[contam_bacteria$any_contam_at_all]) %>%
  group_by(group, sample, genus) %>%
  summarise(new_est_reads_contam=sum(new_est_reads[contaminant]), new_est_reads=sum(new_est_reads)) %>%
  group_by(sample) %>%
  mutate(contam_prop=new_est_reads_contam/new_est_reads) %>%
  group_by(genus) %>%
  summarise(mean=mean(contam_prop), median=(median(contam_prop))) %>%
  mutate(OTU_ID = levels_ranks$OTU_ID[match(genus, levels_ranks$OTU_name)]) %>%
  filter(mean>=0.5) %>%
  arrange(-mean) ->
  genera_mostly_contam

contig_taxonomy_kraken$genus_contaminant <- contig_taxonomy_kraken$taxid %in% genera_mostly_contam$OTU_ID

MAG_taxonomy_sourmashk31 <- read.table("intermediates/assembly/sourmash_lca_rs202_k31.txt", sep=",", header=TRUE, stringsAsFactors=FALSE) %>%
  mutate(sample=stri_extract_first_regex(ID, "^[^.]+"),
         Bin.Id=stri_extract_first_regex(ID, "bin[.][0-9]+"),
         OTU_name=case_when(
           species!="" ~ species,
           genus != "" ~ genus,
           family != "" ~ family,
           order != "" ~ order,
           class != "" ~ class,
           phylum != "" ~ phylum,
           superkingdom != "" ~ superkingdom,
           TRUE ~ NA_character_
         ) %>% substring(4),
         level=case_when(
           species!="" ~ "species",
           genus != "" ~ "genus",
           family != "" ~ "family",
           order != "" ~ "order",
           class != "" ~ "class",
           phylum != "" ~ "phylum",
           superkingdom != "" ~ "superkingdom",
           TRUE ~ "unidentified"))

MAG_taxonomy_sourmashk31 %<>% filter(!sample == exclude_qual)

contig_taxonomy_kraken$bin_OTU_name <- MAG_taxonomy_sourmashk31$OTU_name[match(stri_extract_first_regex(contig_taxonomy_kraken$sample_contig, "^.+(?=_k127)")%>%paste0("_hg_bin.", contig_taxonomy_kraken$bin), MAG_taxonomy_sourmashk31$ID)]

contig_taxonomy_kraken$bin_contaminant <- contig_taxonomy_kraken$bin_OTU_name %in% genera_mostly_contam$genus | contig_taxonomy_kraken$bin_OTU_name %in% contam_bacteria$species[contam_bacteria$any_contam_at_all]

contig_taxonomy_kraken$contaminant <- contig_taxonomy_kraken$species_contaminant | contig_taxonomy_kraken$genus_contaminant | contig_taxonomy_kraken$bin_contaminant

contig_taxonomy_kraken %>% filter(!is.na(bin)) %>%
  group_by(sample, bin, species, read_group) %>%
  summarise(n=n(), contaminant=any(contaminant)) %>%
  group_by(sample, bin) %>%
  arrange(sample, bin, -n) %>%
  summarise(OTU_names=paste(species, n, sep=" ", collapse="\n"), contaminant=any(contaminant),
            read_groups=paste(sort(unique(read_group)), collapse=" ")) ->
  MAG_taxonomy_kraken

MAG_taxonomy_kraken %<>% filter(!sample == exclude_qual)

qa %<>% merge(mutate(MAG_taxonomy_kraken, Bin.Id=paste0("bin.", bin)), by=c("sample", "Bin.Id"), all.x=TRUE, sort=FALSE)
qa %<>% merge(select(MAG_taxonomy_sourmashk31, sample, Bin.Id, OTU_name), by=c("sample", "Bin.Id"), all.x=TRUE, sort=FALSE)

genes <- read.table("intermediates/all_genes_nt_cluster_results_95.tsv", sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(genes) <- c("cluster_rep_gene", "member_gene")
tmp<-rle(genes$cluster_rep_gene)
tmp$values <- 1:length(tmp$values)
genes$cluster_index <- inverse.rle(tmp)
genes %<>%
  mutate(sample_contig=stri_extract_first_regex(member_gene, "^[0-9]+[-][0-9]+(-?Fb)?([-][0-9])?[.][0-9]{4}[-][0-9]{2}[-][0-9]{2}_[0-9]+(_hg)?_k127_[0-9]+"),
         OTU_name=contig_taxonomy_kraken$species[match(sample_contig, contig_taxonomy_kraken$sample_contig)],
         contaminant=contig_taxonomy_kraken$contaminant[match(sample_contig, contig_taxonomy_kraken$sample_contig)])

group_by(genes, cluster_index, OTU_name) %>%
  summarise(n=n(), n_contam=length(which(contaminant))) %>%
  arrange(-n) %>%
  group_by(cluster_index) %>%
  summarise(size=sum(n), n_contam=sum(n_contam), contam_prop=n_contam/size, OTU_names=paste(OTU_name, n, sep=" ", collapse="\n")) %>%
  arrange(-size) ->
  gene_cluster_summary

group_by(genes, cluster_index) %>%
  mutate(sample=stri_extract_first_regex(member_gene, "^[^.]+")) %>%
  summarise(n_samples=n_distinct(sample)) ->
  gene_cluster_samples

gene_cluster_summary$n_samples <- gene_cluster_samples$n_samples[match(gene_cluster_summary$cluster_index, gene_cluster_samples$cluster_index)]

filter(gene_cluster_summary, grepl("Homo[ ]sapiens", OTU_names)) %>%
  mutate(human_status=stri_replace_all_regex(OTU_names,
                                             "Homo[ ]sapiens|cellular[ ]organisms|NA|unclassified|Opisthokonta|root|[0-9]+",
                                             "") %>% trimws() %>% equals("") %>% ifelse("Human","Mixed")) %>%
  select(cluster_index, human_status)->
  human_cluster_status

gene_cluster_summary$human_status <- human_cluster_status$human_status[match(gene_cluster_summary$cluster_index, human_cluster_status$cluster_index)]
gene_cluster_summary$human_status[is.na(gene_cluster_summary$human_status)] <- "No human"

gene_cluster_summary %<>% mutate(contamination=case_when(
  contam_prop == 0 ~ "No contaminants", 
  contam_prop < 0.1 ~ "Low contamination", 
  contam_prop < 0.9 ~ "Middling contamination",
  TRUE ~ "High contamination"))

table(gene_cluster_summary$human_status)
table(gene_cluster_summary$contamination[gene_cluster_summary$human_status=="No human"])

quantile(gene_cluster_summary$size, c(0,0.05,0.25,0.5,0.75,0.95,1))
table(gene_cluster_summary$size >=5)/nrow(gene_cluster_summary)

kmer_path="intermediates/jellyfish_count"
kmer_files <- list.files(path=kmer_path, pattern=".*counts[.]tsv")

kmer_files <- kmer_files[!grepl(exclude_qual, kmer_files, fixed=TRUE)]

first <- TRUE
colnames <- NULL
kmer_data <- lapply(kmer_files, function(f){
  d<-read.table(file.path(kmer_path, f), colClasses=c(ifelse(first, "character", "NULL"), "numeric"), sep="\t", header=FALSE, stringsAsFactors=FALSE)
  if(first){
    colnames(d) <- c("rep_gene", stri_extract_first_regex(f, "^[^.]+"))
    colnames <<- d$rep_gene
  }
  else
    colnames(d) <- stri_extract_first_regex(f, "^[^.]+")
  first <<-FALSE
  cat(".")
  values <- d[,ncol(d)]
  values
}) %>% do.call(rbind,.)
colnames(kmer_data) <- colnames
rownames(kmer_data)<-stri_extract_first_regex(kmer_files, "^[^.]+")

gene_cluster_summary$rep_gene <- genes$cluster_rep_gene[match(gene_cluster_summary$cluster_index, genes$cluster_index)]

kmer_data <- kmer_data[,-which(colnames(kmer_data) %in%  gene_cluster_summary$rep_gene[gene_cluster_summary$human_status != "No human"])]

all_rep_genes <- colnames(kmer_data)

bacterial_depth <- data.frame(sample=sample_lookup$sample, group=sample_lookup$group, log_cum_reads=log10(sample_lookup$bacteria_reads_nocontam), stringsAsFactors=FALSE) %>%
  filter(group != "Negative control") %>% filter(sample != exclude_qual)

bacterial_depth %<>% cbind(kmer_data[match(bacterial_depth$sample, rownames(kmer_data)),])
covariate_cols<-which(colnames(bacterial_depth)[1:17] %in% c("group", "log_cum_reads"))

cluster<-makeCluster(n_clust)
clusterExport(cluster, c("bacterial_depth", "covariate_cols"))
parLapply(cluster, 4:ncol(bacterial_depth), function(col){
  col_name<-colnames(bacterial_depth)[col]
  formula <- as.formula(paste0("`", col_name, "` ~ group + log_cum_reads"))
  model<-summary(glm(formula, data=bacterial_depth[,c(covariate_cols, col)]))
  model$coefficients[2,"Pr(>|t|)"] * sign(model$coefficients[2,"Estimate"])
}) %>% unlist() ->
  p_vals

stopCluster(cluster)

names(p_vals) <- colnames(bacterial_depth)[4:ncol(bacterial_depth)]
directions <- sign(p_vals)
p_vals %<>% abs()
p_vals[which(directions==0)] <- 1

p_vals_reduced <- ifelse(colSums(bacterial_depth[,4:ncol(bacterial_depth)]) >= 20, p_vals, NA)

p_vals.bh <- p.adjust(p_vals)
p_vals_reduced.bh <- p.adjust(p_vals_reduced)

table(p_vals.bh <= 0.5, directions)
table(p_vals_reduced.bh <= 0.5, directions)

bacterial_depth$total_genes <- rowSums(bacterial_depth[,-(1:3)]>0)
bacterial_depth$total_kmers <- rowSums(bacterial_depth[,-(1:3)])
bacterial_depth %<>% relocate(c(total_genes, total_kmers), .after=log_cum_reads)
covariate_cols<-which(colnames(bacterial_depth)[1:17] %in% c("group", "log_cum_reads", "total_genes", "total_kmers"))

cluster<-makeCluster(n_clust)
clusterEvalQ(cluster, {library(logistf)})
clusterExport(cluster, c("bacterial_depth", "covariate_cols"))

cols_to_test<-which(colSums(bacterial_depth[,6:ncol(bacterial_depth)]>0) >= 10)+5
names(cols_to_test) <- colnames(bacterial_depth)[cols_to_test]
Sys.time()
parLapply(cluster, cols_to_test, function(col){
  col_name<-colnames(bacterial_depth)[col]
  formula <- as.formula(paste0("`", col_name, "` > 0 ~ group + log(total_genes)"))
  model<-logistf(formula, bacterial_depth[, c(covariate_cols, col)])
  model$prob[2] * sign(model$coefficients[2])
}) %>% unlist() ->
  p_vals_lr
Sys.time()

stopCluster(cluster)

directions_lr <- sign(p_vals_lr)
p_vals_lr %<>% abs()
p_vals_lr[which(directions_lr==0)] <- 1

p_vals_lr.bh <- p.adjust(p_vals_lr)

table(p_vals_lr <= 0.05, directions_lr)
table(p_vals_lr.bh <= 0.05, directions_lr)
table(p_vals_lr.bh <= 0.2, directions_lr)

KD_up <- names(p_vals_lr.bh[which(p_vals_lr.bh<=0.05 & directions_lr==1)])


#Representative gene detection in samples
set.seed(1984)
(lapply(sample(colnames(kmer_data), 10000), function(col){
  sample <- stri_extract_first_regex(col, "^[^.]+")
  if(sample %in% rownames(kmer_data))
    kmer_data[sample,col] > 0
  else
    NULL
}) %>% unlist() %>% table() -> tmp)
tmp/sum(tmp)

con <- file("../../output/megahit/all_genes_nt_cluster_results_95.fasta", "rt")
next_header<-character(0)
kmer_presences <- vector(mode="list", length=0)
for(i in 1:100){
  next_cluster<-FALSE
  lines <- next_header
  last_line_header<-FALSE
  while(!next_cluster){
    new_line <- scan(con, what=character(), nmax=1, sep="\n", nlines=1)
    if(last_line_header & substring(new_line, 1,1)==">"&length(lines)>2){
      next_cluster <- TRUE
      next_header<-c(lines[length(lines)], new_line)
      lines <- lines[1:(length(lines)-1)]
    }else{
      lines %<>% c(new_line)
    }
    if(substring(new_line, 1,1)==">"){
      last_line_header<-TRUE
    }else{
      last_line_header<-FALSE
    }
  }
  
  lines[-1] %>% paste(collapse="\n") ->
    fasta
  fasta_conn <- textConnection(fasta, open = "r")
  cluster_sequences<-read.FASTA(fasta_conn) %>% as.character() %>% lapply(paste0, collapse="")
  names(cluster_sequences) %<>% stri_extract_first_regex("^[^ ]+")
  close(fasta_conn)
  cluster_rep<-substring(lines[1],2)
  kmers<-system(paste0('grep -F -A 30 "', cluster_rep,'" ../../output/megahit/all_genes_nt_cluster_results_95_rep_jf_unique_kmer_reduced.fasta'), intern=TRUE, show.output.on.console = FALSE)
  kmers <-kmers[2:(grep("^>", kmers)[2]-1)]
  kmer_presence <- lapply(tolower(kmers), grep, cluster_sequences, fixed=TRUE)
  kmer_absence <- lapply(tolower(kmers), grep, cluster_sequences, fixed=TRUE, invert=TRUE)
  new_index<-length(kmer_presences)+1
  data.frame(cluster_rep_seq=cluster_rep,
             gene=names(cluster_sequences)[c(unlist(kmer_presence), unlist(kmer_absence))],
             kmers=c(rep(kmers, lapply(kmer_presence, length) %>% unlist()), rep(kmers, lapply(kmer_absence, length) %>% unlist())),
             presence=c(rep(TRUE, length(unlist(kmer_presence))), rep(FALSE, length(unlist(kmer_absence)))),
             stringsAsFactors=FALSE) %>%
    pivot_wider(id_cols=c(cluster_rep_seq, gene), names_from=kmers, values_from=presence, values_fill=list(kmers=FALSE)) ->
    kmer_presences[[new_index]]
  kmer_presences[[new_index]]$length <- lapply(cluster_sequences[kmer_presences[[new_index]]$gene], nchar) %>% unlist()
  kmer_presences[[new_index]]$any_kmer <- rowSums(kmer_presences[[new_index]][,kmers])>0
}
close(con)

lapply(kmer_presences, select, cluster_rep_seq, any_kmer, length) %>%
  bind_rows() %>%
  group_by(cluster_rep_seq) %>%
  summarise(members=n(), has_any_kmer=length(which(any_kmer)), has_any_kmer_prop=has_any_kmer/members, min_length=min(length), max_length=max(length), median_length=median(length), median_length_no_kmer=median(length[!any_kmer], na.rm=TRUE), median_length_kmer=median(length[any_kmer], na.rm=TRUE))->
  kmer_presence_summary

quantile(kmer_presence_summary$has_any_kmer_prop, c(0,0.25,0.5,0.75,1))
quantile(kmer_presence_summary$median_length_no_kmer, c(0,0.25,0.5,0.75,1), na.rm=TRUE)
quantile(kmer_presence_summary$median_length_kmer, c(0,0.25,0.5,0.75,1))