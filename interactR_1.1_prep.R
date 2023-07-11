## this script is for data cleaning, harmonizing, and processing the published interactome databases
## to allow searching for a protein of interest in the interactomrs from many databases at once
## within the interactR shiny app
## https://dominico.shinyapps.io/interactR/

## database versions included in this version 1.1 (June 28 2023)

## STRING v12 physical links network, experimental evidence only (not text mining) https://stringdb-downloads.org/download/protein.physical.links.full.v12.0/9606.protein.physical.links.full.v12.0.txt.gz (11.8 Mb). STRING now allows you to view the physical-only interaction network. The physical network will only display edges between the proteins for which we have evidence of their binding or forming a physical complex.
## BioGRID 4.4.222 https://thebiogrid.org

## Human cell map v1 (BFDR threshold=0.1) https://humancellmap.org/
## BioPlex 1.0, 2.0, & 3.0 https://bioplex.hms.harvard.edu/
## OpenCell https://opencell.czbiohub.org/download
## Human Reference interactome http://www.interactome-atlas.org/
## HIPPIE v2.3 http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php


library(tidyverse)
library(magrittr)
library(HGNChelper)
library(data.table)
library(annotables)
library(Biostrings)


################
#### inputs ####
################

base = "C:/Users/dowens/OneDrive/Postdoc/Projects/interactR/"

dat_dir = paste0(base,"data/")
dir.create(dat_dir,showWarnings = F)


# uniprot reviewed human proteome (downloaded within FragPipe)
prot_file = paste0(dat_dir,"2023-04-19-decoys-reviewed-contam-UP000005640.fa")


######## STRING-DB ##########

## string db download
string_file = paste0(dat_dir,"9606.protein.physical.links.full.v12.0.txt")
prot_info_file = paste0(dat_dir,"9606.protein.info.v12.0.txt")




######## HUTTLIN ##########

## BioPlex Huttlin AP-MS
huttlin_file1 = paste0(dat_dir,"BioPlex_interactionList_v2.tsv")
huttlin_file2 = paste0(dat_dir,"BioPlex_interactionList_v4a.tsv")
huttlin_file3 = paste0(dat_dir,"BioPlex_HCT116_Network_5.5K_Dec_2019.tsv")
huttlin_file4 = paste0(dat_dir,"BioPlex_293T_Network_10K_Dec_2019.tsv")


######## BIOGRID ##########

## load the full bioGRID database
#grid_file = paste0(dat_dir,"BIOGRID-ALL-4.4.222.tab3.txt")

# load the preanalyzed data as analyzing takes ages and has already been done
## data saved as Rds
grid_dat_file = paste0(dat_dir,"Biogrid_4.4.222_precleaned.Rds")



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
  dplyr::rename(A=`Symbol A`,B=`Symbol B`) %>%
  mutate(ref="PMID:26186194")

# bioplex 2
b = read_tsv(huttlin_file2) %>%
  dplyr::select(SymbolA,SymbolB) %>%
  dplyr::rename(A=SymbolA,B=SymbolB) %>%
  mutate(ref="PMID:28514442")

# bioplex 3
d = read_tsv(huttlin_file3) %>%
  dplyr::select(SymbolA,SymbolB) %>%
  dplyr::rename(A=SymbolA,B=SymbolB) %>%
  mutate(ref="PMID:33961781")

# bioplex 3
e = read_tsv(huttlin_file4) %>%
  dplyr::select(SymbolA,SymbolB) %>%
  dplyr::rename(A=SymbolA,B=SymbolB) %>%
  mutate(ref="PMID:33961781")


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
  dplyr::select(A_1,B_1,ref) %>%
  dplyr::rename(A=A_1,B=B_1)

# make non redundant
hut_dat %<>%
  distinct()

## add database info and experiment_type info
hut_dat %<>%
  mutate(ref=ref,
         experiment_type = "IP-MS",
         db="bioplex.hms.harvard.edu")


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
  filter(experiments>700)


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
  dplyr::rename(A=protein1_name,B=protein2_name) %>%
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
  dplyr::rename(A=A_1,B=B_1)


## only keep unique
## actually they are all unique already!
histring %<>%
  distinct()


## add database and ref info
histring %<>%
  mutate(ref="PMID:30476243",
         experiment_type = "unknown",
         db="string-db.org") %>%
  as_tibble()



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

## save preanalyzed data
#saveRDS(grid_dat,paste0(dat_dir,"Biogrid_4.4.222_precleaned.Rds"))




## load preanalyzed data
grid_dat = readRDS(grid_dat_file)

## keep just the physical evidence interactions (interactors)
grid_dat %<>%
  filter(`Experimental System Type`=="physical")

## check the different physical experiment type options
unique(grid_dat$`Experimental System`)

# rank order them
## most are IP-MS and Y2H
test =
  grid_dat %>%
  group_by(`Experimental System`) %>%
  dplyr::count()

## tidy
grid_dat %<>%
  dplyr::select(A,B,`Experimental System`,`Publication Source`) %>%
  dplyr::rename(experiment_type=`Experimental System`,ref=`Publication Source`)

## rename the experiment_types to match with the other conventions
grid_dat %<>%
  mutate(experiment_type=if_else(
    experiment_type=="Affinity Capture-MS",
    "IP-MS",
    experiment_type
))


## drop NAs
grid_dat %<>%
  drop_na()

## for ambiguosly named (that look like eg KAT14 /// PET117)
## just keep first protein mapped
grid_dat %<>%
  separate(A, into="A_1",sep=" ///",remove=F) %>%
  separate(B, into="B_1",sep=" ///",remove=F) %>%
  #filter(grepl("///",A)|grepl("///",B)) %>%
  dplyr::select(A_1,B_1,experiment_type,ref) %>%
  dplyr::rename(A=A_1,B=B_1)

## only keep unique
grid_dat %<>%
  distinct()

## add database and ref info
grid_dat %<>%
  mutate(ref=gsub("PUBMED","PMID",ref),
       experiment_type = experiment_type,
       db="thebiogrid.org") %>%
  as_tibble()




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
  dplyr::rename(A=Bait,B=PreyGene) %>%
  dplyr::select(A,B)

## for ambiguosly named (that look like eg KAT14 /// PET117)
## just keep first protein mapped
cellmap %<>%
  separate(A, into="A_1",sep=" ///",remove=F) %>%
  separate(B, into="B_1",sep=" ///",remove=F) %>%
  #filter(grepl("///",A)|grepl("///",B)) %>%
  dplyr::select(A_1,B_1) %>%
  dplyr::rename(A=A_1,B=B_1)

## only keep unique
cellmap %<>%
  distinct()

## add database and ref info
cellmap %<>%
  mutate(ref="PMID:34079125",
         experiment_type = "Proximity Label-MS",
         db="humancellmap.org") %>%
  as_tibble()


########################
#### OPEN CELL PREP ####
########################

# read data
opencell = read_csv(opencell_file)

## tidy
opencell %<>%
  dplyr::rename(A=target_gene_name,B=interactor_gene_name) %>%
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
  dplyr::rename(A=A_1,B=B_1)

## only keep unique
opencell %<>%
  distinct()

## add database and ref info
opencell %<>%
  mutate(ref="PMID:35271311",
         experiment_type = "IP-MS",
         db="opencell.czbiohub.org") %>%
  as_tibble()



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
  dplyr::rename(A=symbol) %>%
  inner_join(pcgenes,by=c("B"="ensgene")) %>%
  dplyr::select(A,symbol) %>%
  dplyr::rename(B=symbol)

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
  dplyr::rename(A=A_1,B=B_1)

## only keep unique
huri %<>%
  distinct()



nrow(huri)
# 63243 final interactions

## add database and ref info
huri %<>%
  mutate(ref="PMID:32296183",
         experiment_type = "Two-hybrid",
         db="interactome-atlas.org") %>%
  as_tibble()


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
  dplyr::rename(A=symbol) %>%
  inner_join(pcgenes_entrez,by=c("entrezB"="entrez")) %>%
  dplyr::select(A,symbol) %>%
  dplyr::rename(B=symbol)

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
  dplyr::rename(A=A_1,B=B_1)

## only keep unique
hipdat %<>%
  distinct()

nrow(hipdat)
# 799808 final interactions

## add database and ref info
hipdat %<>%
  mutate(ref="PMID:27794551",
         experiment_type = "unknown",
         db="cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/") %>%
  as_tibble()


###################################################
#### COMBINE INTERACTOMES AND ONLY KEEP UNIQUE ####
###################################################


## databases loaded

## save a complete list of interactors where A and B does mean bait and prey


## now actually do the combining of the tables into one and double check to make non-redundant
ints =
  rbind.data.frame(cellmap,
                   grid_dat,
                   histring,
                   hut_dat,
                   opencell,
                   huri,
                   hipdat)


## remove the self interactors as they are meaningless
ints %<>%
  filter(A!=B)





###### UNIPROT ID #######

## add the uniprot IDs by using uniprot reference proteome

## load uniprot human proteome
sequences <- readAAStringSet(prot_file)

# Convert to a data frame and tidy up metadata
prot =
  as.data.frame(sequences) %>%
  rownames_to_column("meta") %>%
  as_tibble() %>%
  separate(meta,into=c( "prefix",
                        "ID",
                        "rest"),
           sep="\\|") %>%
  separate(rest, into=c(NA,"species"),sep="OX=",remove=F) %>%
  separate(rest, into=c(NA,"Gene"),sep="GN=",remove=T) %>%
  mutate(species=sub(" .*", "", species)) %>%
  mutate(Gene=sub(" .*", "", Gene)) %>%
  dplyr::rename(seq=x)


## got this reference proteome from within fragpipe, so contains contaminants and reversered peptides
## remove these
prot %<>%
  filter(prefix!="rev_sp") %>%
  filter(species=="9606") %>%
  dplyr::select(Gene,ID)


## join to the main ints table
ints %<>%
  left_join(prot,by=c("A"="Gene")) %>%
  dplyr::rename(A_Gene=A,A_ID=ID) %>%
  left_join(prot,by=c("B"="Gene")) %>%
  dplyr::rename(B_Gene=B,B_ID=ID) %>%
  dplyr::select(A_Gene,A_ID,B_Gene,B_ID,ref,experiment_type,db)


### save data object of the complete interactome
## this is the file for interactR shiny app
saveRDS(ints,file=paste0(dat_dir,"InteractR_interactome_full.Rds"))




###########################################
## NOT NECESSARY PAST HERE FOR INTERACTR ##
###########################################







##### now interested in a NON-REDUNDANT list of interactors

## where protein1-protein2 will only be considered once, and protein2-protein1 will be ignored

## add dummy info to later columns to allow sorting to work to remove duplicate entries

# function to add dummy column info to allow sorting
dummycols = function(data){
  data %<>%
    mutate(ref=paste0("zzzzzzzzzz_",ref),
           experiment_type = paste0("zzzzzzzzzzz_",experiment_type),
           db=paste0("zzzzzzzzzzzz_",db))
  return(data)
}


cellmap2 = dummycols(cellmap)
grid_dat2 = dummycols(grid_dat)
histring2 = dummycols(histring)
hut_dat2 = dummycols(hut_dat)
opencell2 = dummycols(opencell)
huri2 = dummycols(huri)
hipdat2 = dummycols(hipdat)


# Calculate the numbers of unique non-redundant interactions in each database
## exclude the self interactions
# create a list of the tables and loop through to process
## note this is just to calculate, and doesn't actually do any joining

dat_list=list(cellmap2,grid_dat2,histring2,hut_dat2,opencell2,huri2,hipdat2)
names(dat_list) = c("cellmap","biogrid","string","bioplex","opencell","huri","hippie")

int_nums = tibble()
i=1
for(i in 1:length(dat_list)){
  
  dat=dat_list[[i]]
  this_dat = names(dat_list)[i]
  
  cat(this_dat,"\n")
  
  # nrows of the redundant original list
  prefiltnum = nrow(dat)
  
  ## remove the self interactions as that is meaningless
  dat %<>%
    filter(A!=B)
  
  # only keep unique non-redundant interactions
  # sort each row alphabetically
  dat_sort <- t(apply(dat, 1, sort))
  
  # remove duplicates
  dat_unique <- unique(dat_sort)
  
  dat = 
    as.data.frame(dat_unique) %>%
    tibble() %>%
    dplyr::rename(A=V1,B=V2)
  
  postfiltnum = nrow(dat)
  
  temp = tibble(prefiltnum,postfiltnum,this_dat)
  int_nums %<>% rbind.data.frame(temp)
  
  
}



### now create a non-redundant version
## remove ref, experiment_type and database info for just non-redundant interactome count
ints2 =
  ints %>%
  dplyr::select(A,B) %>%
  distinct()


# only keep unique non-redundant interactions
# sort each row alphabetically
dat_sort <- t(apply(ints2, 1, sort))

# remove duplicates
dat_unique <- unique(dat_sort)

ints2 = 
  as.data.frame(dat_unique) %>%
  tibble() %>%
  dplyr::rename(A=V1,B=V2)


### save data object of the complete interactome
saveRDS(ints2,file=paste0(dat_dir,"InteractR_interactome.Rds"))


