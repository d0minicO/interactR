## this script is for processing the published interactors databases
## to allow searching a protein of itnerests interactome from many databases at once

## so far covering

## STRING (0.4 med conf cutoff, experimental evidence only)
## Human cell map (BFDR threshold=0.1)
## BioGRID
## BioPlex Huttlin AP-MS interactomes
## OpenCell https://opencell.czbiohub.org/download
## Human reference interactome Y2H (Luck 2020 Nature)
## HIPPIE v2.3 http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php

library(tidyverse)
library(magrittr)
library(HGNChelper)
library(data.table)



################
#### inputs ####
################

base = "C:/Users/dowens/OneDrive/Postdoc/Projects/interactR/"

dat_dir = paste0(base,"data/")
dir.create(dat_dir,showWarnings = F)


######## STRING-DB ##########

## string db download
string_file = paste0(dat_dir,"9606.protein.links.full.v11.5.txt")
prot_info_file = paste0(dat_dir,"9606.protein.info.v11.5.txt")
prot_seqs_file = paste0(dat_dir, "9606.protein.sequences.v11.5.fa")




######## HUTTLIN ##########

## BioPlex Huttlin AP-MS
huttlin_file1 = paste0(dat_dir,"BioPlex_interactionList_v2.tsv")
huttlin_file2 = paste0(dat_dir,"BioPlex_interactionList_v4a.tsv")
huttlin_file3 = paste0(dat_dir,"BioPlex_HCT116_Network_5.5K_Dec_2019.tsv")
huttlin_file4 = paste0(dat_dir,"BioPlex_293T_Network_10K_Dec_2019.tsv")


######## BIOGRID ##########

## load the full bioGRID database
grid_file = paste0(dat_dir,"BIOGRID-ALL-4.4.213.tab3.txt")

# load the preanalyzed data as analyzing takes ages and has already been done
## data saved as Rds
grid_dat_file = paste0(dat_dir,"BioGrid_table_hgnc_compliant.Rds")



######## HUMAN CELL MAP ##########

#GO et al bio id human cell map
cellmap_file = paste0(dat_dir,"saint-latest_humancellmap.txt")

BFDR_thresh = 0.1



######## OPEN CELL ##########

opencell_file = paste0(dat_dir,"opencell-protein-interactions.csv")



######## HURI ##########

huri_file = paste0(dat_dir,"HI-union.tsv")



######## HIPPIE ##########

hip_file = paste0(dat_dir,"hippie_current.txt")


######## HGNC naming ##########

## load new HGNC map
# downloaded April 2022 https://github.com/waldronlab/HGNChelper/blob/master/data/hgnc.table.rda
load(file=paste0(dat_dir,"hgnc.table.rda"))
#"https://github.com/d0minicO/unimachR/blob/main/hgnc.table.rda?raw=TRUE"
hgnc.table %<>%
  dplyr::select(1,2)


#################
#### outputs ####
#################


out_dir = paste0(base,"interactR_results/")
dir.create(out_dir,showWarnings = F)



#################################################################################
########################### data load and prep ##################################
#################################################################################




###################################
#### HUTTLIN INTERACTOMES PREP ####
###################################

## prepare the data

# bioplex 1
a = read_tsv(huttlin_file1) %>%
  dplyr::select(`Symbol A`,`Symbol B`) %>%
  dplyr::rename(A=`Symbol A`,B=`Symbol B`)

# bioplex 2
b = read_tsv(huttlin_file2) %>%
  dplyr::select(SymbolA,SymbolB) %>%
  dplyr::rename(A=SymbolA,B=SymbolB)

# bioplex 3
d = read_tsv(huttlin_file3) %>%
  dplyr::select(SymbolA,SymbolB) %>%
  dplyr::rename(A=SymbolA,B=SymbolB)

# bioplex 3
e = read_tsv(huttlin_file4) %>%
  dplyr::select(SymbolA,SymbolB) %>%
  dplyr::rename(A=SymbolA,B=SymbolB)


hut_dat = rbind.data.frame(a,b,d,e)

# make non redundant
hut_dat %<>%
  distinct()

## correct hgnc_symbols
hut_dat$A =
  HGNChelper::checkGeneSymbols(hut_dat$A,map=hgnc.table) %>%
  dplyr::pull(Suggested.Symbol)

## correct hgnc_symbols
hut_dat$B =
  HGNChelper::checkGeneSymbols(hut_dat$B,map=hgnc.table) %>%
  dplyr::pull(Suggested.Symbol)


## drop NAs
hut_dat %<>%
  drop_na()

## for ambiguosly named (that look like eg KAT14 /// PET117)
## just keep first protein mapped
hut_dat %<>%
  separate(A, into="A_1",sep=" ///",remove=F) %>%
  separate(B, into="B_1",sep=" ///",remove=F) %>%
  #filter(grepl("///",A)|grepl("///",B)) %>%
  dplyr::select(A_1,B_1) %>%
  rename(A=A_1,B=B_1)

# make non redundant
hut_dat %<>%
  distinct()


########################
#### STRING-DB PREP ####
########################

## could use string package but for now just using downloaded string database
#https://string-db.org/cgi/download?sessionId=b6lgRJ5JMFCs&species_text=Homo+sapiens October 13th 2022

string_dat = fread(string_file)
prot_info = fread(prot_info_file)


## different confidences in combined_score and meaning, from string-db.org
# highest confidence > 900
# high confidence > 700
# med confidence is > 400
# low confidence is > 150

# just keep the interactions with experimental evidence
hi =
  string_dat %>%
  filter(experiments>400)


## join the two tables using dt

# set keys to help join faster
setkey(hi, "protein1")
setkey(prot_info,"#string_protein_id")

## join to get protein1 name, but keep other column names too
joined = 
  hi[prot_info, nomatch = 0] %>%
  #dplyr::select(protein1,protein2,combined_score,preferred_name) %>%
  dplyr::rename(protein1_name = preferred_name)

# join again to get protein2 name, dt not working so dplyr join
prot_info %<>%
  dplyr::rename(protein2 =`#string_protein_id`)

joined %<>%
  left_join(prot_info, by = "protein2") %>%
  dplyr::select(protein1,protein1_name,protein2,preferred_name,combined_score) %>%
  dplyr::rename(protein2_name = preferred_name)


histring = joined

# tidy
histring %<>%
  rename(A=protein1_name,B=protein2_name) %>%
  dplyr::select(A,B)

## correct hgnc_symbols
histring$A =
  HGNChelper::checkGeneSymbols(histring$A,map=hgnc.table) %>%
  dplyr::pull(Suggested.Symbol)

## correct hgnc_symbols
histring$B =
  HGNChelper::checkGeneSymbols(histring$B,map=hgnc.table) %>%
  dplyr::pull(Suggested.Symbol)


## drop NAs
histring %<>%
  drop_na()

## for ambiguosly named (that look like eg KAT14 /// PET117)
## just keep first protein mapped
histring %<>%
  separate(A, into="A_1",sep=" ///",remove=F) %>%
  separate(B, into="B_1",sep=" ///",remove=F) %>%
  #filter(grepl("///",A)|grepl("///",B)) %>%
  dplyr::select(A_1,B_1) %>%
  rename(A=A_1,B=B_1)


## only keep unique
## actually they are all unique already!
histring %<>%
  distinct()


######################
#### BIOGRID PREP ####
######################

# load data
#grid_dat = read_tsv(grid_file)

## correct names
## takes ages!!!
#grid_dat$`Official Symbol Interactor A`  =
#  HGNChelper::checkGeneSymbols(grid_dat$`Official Symbol Interactor A`,map=hgnc.table) %>%
#  dplyr::pull(Suggested.Symbol)

#grid_dat$`Official Symbol Interactor B`  =
#  HGNChelper::checkGeneSymbols(grid_dat$`Official Symbol Interactor B`,map=hgnc.table) %>%
#  dplyr::pull(Suggested.Symbol)


## clean up
#grid_dat %<>%
#  dplyr::rename(A= `Official Symbol Interactor A`) %>%
#  dplyr::rename(B = `Official Symbol Interactor B`) 

## load preanalyzed data
grid_dat = readRDS(grid_dat_file)


## keep just the physical evidence interactions (interactors)
grid_dat %<>%
  filter(`Experimental System Type`=="physical")


## tidy
grid_dat %<>%
  dplyr::select(A,B)


## drop NAs
grid_dat %<>%
  drop_na()

## for ambiguosly named (that look like eg KAT14 /// PET117)
## just keep first protein mapped
grid_dat %<>%
  separate(A, into="A_1",sep=" ///",remove=F) %>%
  separate(B, into="B_1",sep=" ///",remove=F) %>%
  #filter(grepl("///",A)|grepl("///",B)) %>%
  dplyr::select(A_1,B_1) %>%
  rename(A=A_1,B=B_1)

## only keep unique
grid_dat %<>%
  distinct()



#############################
#### HUMAN CELL MAP PREP ####
#############################

# read data
cellmap = read_tsv(cellmap_file)

# filter on given BFDR
cellmap %<>%
  filter(BFDR<BFDR_thresh)

## correct names
cellmap$Bait =
  HGNChelper::checkGeneSymbols(cellmap$Bait,map=hgnc.table) %>%
  dplyr::pull(Suggested.Symbol)

cellmap$PreyGene =
  HGNChelper::checkGeneSymbols(cellmap$PreyGene,map=hgnc.table) %>%
  dplyr::pull(Suggested.Symbol)


## drop NAs
cellmap %<>%
  drop_na()

## tidy
cellmap %<>%
  rename(A=Bait,B=PreyGene) %>%
  dplyr::select(A,B)

## for ambiguosly named (that look like eg KAT14 /// PET117)
## just keep first protein mapped
cellmap %<>%
  separate(A, into="A_1",sep=" ///",remove=F) %>%
  separate(B, into="B_1",sep=" ///",remove=F) %>%
  #filter(grepl("///",A)|grepl("///",B)) %>%
  dplyr::select(A_1,B_1) %>%
  rename(A=A_1,B=B_1)

## only keep unique
cellmap %<>%
  distinct()



########################
#### OPEN CELL PREP ####
########################

# read data
opencell = read_csv(opencell_file)

## tidy
opencell %<>%
  rename(A=target_gene_name,B=interactor_gene_name) %>%
  dplyr::select(A,B)


## correct names
opencell$A =
  HGNChelper::checkGeneSymbols(opencell$A,map=hgnc.table) %>%
  dplyr::pull(Suggested.Symbol)

opencell$B =
  HGNChelper::checkGeneSymbols(opencell$B,map=hgnc.table) %>%
  dplyr::pull(Suggested.Symbol)


## drop NAs
opencell %<>%
  drop_na()


## for ambiguosly named (that look like eg KAT14 /// PET117)
## just keep first protein mapped
opencell %<>%
  separate(A, into="A_1",sep=" ///",remove=F) %>%
  separate(B, into="B_1",sep=" ///",remove=F) %>%
  #filter(grepl("///",A)|grepl("///",B)) %>%
  dplyr::select(A_1,B_1) %>%
  rename(A=A_1,B=B_1)

## only keep unique
opencell %<>%
  distinct()




###################
#### HURI PREP ####
###################

# read data
huri = read_tsv(huri_file,col_names=c("A","B"))

nrow(huri)

huri %>%
  distinct() %>%
  nrow()



# get a non redundant mapping of ENSID to Symbol
pcgenes =
  grch38 %>%
  filter(biotype=="protein_coding") %>%
  dplyr::select(ensgene,symbol) %>%
  distinct()


## correct names
pcgenes$symbol =
  HGNChelper::checkGeneSymbols(pcgenes$symbol,map=hgnc.table) %>%
  dplyr::pull(Suggested.Symbol)


pcgenes %<>%
  drop_na()

# some ENSIDs match to multiple gene symbols
# checked a few on ensembl and they just look like alternate sequence genes
pcgenes %>%
  group_by(symbol) %>%
  dplyr::count() %>%
  arrange(desc(n))


## convert ENSID to symbols
## lose 28 interactions mapping A
## lose 307 interactions mapping B
## very little so happy to proceed
## checked some lost and was not-well annotated read-through gene
huri %<>%
  inner_join(pcgenes,by=c("A"="ensgene")) %>%
  dplyr::select(symbol,B) %>%
  rename(A=symbol) %>%
  inner_join(pcgenes,by=c("B"="ensgene")) %>%
  dplyr::select(A,symbol) %>%
  rename(B=symbol)

## drop NAs
huri %<>%
  drop_na()


## for ambiguosly named (that look like eg KAT14 /// PET117)
## just keep first protein mapped
huri %<>%
  separate(A, into="A_1",sep=" ///",remove=F) %>%
  separate(B, into="B_1",sep=" ///",remove=F) %>%
  #filter(grepl("///",A)|grepl("///",B)) %>%
  dplyr::select(A_1,B_1) %>%
  rename(A=A_1,B=B_1)

## only keep unique
huri %<>%
  distinct()



nrow(huri)
# 63665 final interactions



#####################
#### HIPPIE PREP ####
#####################

# read data
hipdat = 
  read_tsv(hip_file,col_names=c(NA,"entrezA",NA,"entrezB")) %>%
  dplyr::select("entrezA","entrezB") %>%
  filter(entrezA!=entrezB)


nrow(hipdat)

hipdat %>%
  distinct() %>%
  nrow()


# get a non redundant mapping of entrezID to Symbol
pcgenes_entrez =
  grch38 %>%
  filter(biotype=="protein_coding") %>%
  dplyr::select(entrez,symbol) %>%
  distinct()


## correct names
pcgenes_entrez$symbol =
  HGNChelper::checkGeneSymbols(pcgenes_entrez$symbol,map=hgnc.table) %>%
  dplyr::pull(Suggested.Symbol)


pcgenes_entrez %<>%
  drop_na()

# some entrez ids match to multiple gene symbols
# checked a few on ensembl and they just look like alternate sequence genes
pcgenes_entrez %>%
  group_by(symbol) %>%
  dplyr::count() %>%
  arrange(desc(n))


## convert entrez to symbols
hipdat %<>%
  inner_join(pcgenes_entrez,by=c("entrezA"="entrez")) %>%
  dplyr::select(symbol,entrezB) %>%
  rename(A=symbol) %>%
  inner_join(pcgenes_entrez,by=c("entrezB"="entrez")) %>%
  dplyr::select(A,symbol) %>%
  rename(B=symbol)

## drop NAs
hipdat %<>%
  drop_na()


## for ambiguosly named (that look like eg KAT14 /// PET117)
## just keep first protein mapped
hipdat %<>%
  separate(A, into="A_1",sep=" ///",remove=F) %>%
  separate(B, into="B_1",sep=" ///",remove=F) %>%
  #filter(grepl("///",A)|grepl("///",B)) %>%
  dplyr::select(A_1,B_1) %>%
  rename(A=A_1,B=B_1)

## only keep unique
hipdat %<>%
  distinct()

nrow(hipdat)
# 801698 final interactions



###################################################
#### COMBINE INTERACTOMES AND ONLY KEEP UNIQUE ####
###################################################


## databases loaded

cellmap
grid_dat
histring
hut_dat
opencell
huri
hipdat

# create a list and loop through to process
# only keep unique non-redundant interactions

dat_list=list(cellmap,grid_dat,histring,hut_dat,opencell,huri,hipdat)
names(dat_list) = c("cellmap","biogrid","string","bioplex","opencell","huri","hippie")

int_nums = tibble()
i=1
for(i in 1:length(dat_list)){
  
  dat=dat_list[[i]]
  this_dat = names(dat_list)[i]
  
  cat(this_dat,"\n")
  
  # nrows of the redundant original list
  prefiltnum = nrow(dat)
  
  
  # only keep unique non-redundant interactions
  # sort each row alphabetically
  dat_sort <- t(apply(dat, 1, sort))
  
  # remove duplicates
  dat_unique <- unique(dat_sort)
  
  dat = 
    as.data.frame(dat_unique) %>%
    tibble() %>%
    rename(A=V1,B=V2)
  
  postfiltnum = nrow(dat)
  
  temp = tibble(prefiltnum,postfiltnum,this_dat)
  int_nums %<>% rbind.data.frame(temp)
  
  
}


## combine tables into one and make non-redundant
ints =
  rbind.data.frame(cellmap,
                 grid_dat,
                 histring,
                 hut_dat,
                 opencell,
                 huri,
                 hipdat)


## make unique
ints_sort <- t(apply(ints, 1, sort))

# remove duplicates
ints_unique <- unique(ints_sort)

ints = 
  as.data.frame(ints_unique) %>%
  tibble() %>%
  rename(A=V1,B=V2)



### save data object of the complete interactome
saveRDS(ints,file=paste0(dat_dir,"InteractR_interactome.Rds"))
