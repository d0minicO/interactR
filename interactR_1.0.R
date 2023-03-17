## this script is for looking for published interactors from many databases at once

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
library(annotables)


#################
#### OPTIONS ####
#################

options(max.print=100)



###################
#### FUNCTIONS ####
###################

interactR <- function(A_search,hut_dat,histring,grid_dat,cellmap,export=T,geneGrep=NULL){
  
  ################################
  #### DO INTERACTOMES SEARCH ####
  ################################
  
  ## Huttlin
  A_hut =
    hut_dat %>%
    filter((SymbolA %in% A_search | SymbolB %in% A_search)) %>%
    drop_na()
  
  ######################
  #### SAVE OUTPUTS ####
  ######################
  
  ## save counts info
  ##  checkif we have entries for this protein
  ## if so then create the directory and save the files
  c_hut = nrow(A_hut)
  c_string = nrow(A_string)
  c_grid = nrow(A_grid)
  c_map = nrow(A_map)
  
  total = sum(c_hut,c_string,c_grid,c_map)

  # df of numbers of interactors
  counts = tibble(
    prot=A_search,
    Hut=c_hut,
    String=c_string,
    Grid=c_grid,
    CellMap=c_map
  )
  
  ## combined interactor table, grepping on a certain gene name / family root
  p = A_hut %>%
    filter(grepl(geneGrep,SymbolA)|grepl(geneGrep,SymbolB)) %>%
    dplyr::select(SymbolA,SymbolB) %>%
    distinct() %>%
    mutate(type="Huttlin_PMID:33961781") %>%
    rename(A=SymbolA,B=SymbolB)
  
  q = A_string %>% 
    filter(grepl(geneGrep,protein1_name)|grepl(geneGrep,protein2_name)) %>%
    dplyr::select(protein1_name,protein2_name) %>%
    distinct() %>%
    mutate(type="STRING") %>%
    rename(A=protein1_name,B=protein2_name)
  
  r = A_grid %>%
    filter(grepl(geneGrep,A)|grepl(geneGrep,B)) %>%
    dplyr::select(A,B,`Publication Source`) %>%
    rename(type=`Publication Source`) %>%
    distinct()
  
  s = A_map %>%
    #filter(grepl(geneGrep,Bait)|grepl(geneGrep,PreyGene)) %>%
    dplyr::select(Bait,PreyGene) %>%
    rename(A=Bait,B=PreyGene) %>%
    distinct() %>%
    mutate(type="HumanCellMap_PMID:34079125")
  
  out=rbind.data.frame(p,q,r,s)
  
  if(total>0&export){
    
    ## create the output directory for this protein
    prot_dir = paste0(out_dir,A_search,"/")
    dir.create(prot_dir,showWarnings = F)
    
    ## save the tables if entries
    write.table(A_hut,
                file = paste0(prot_dir,"Huttlin_coIPs_",A_search,".tsv"),
                sep="\t",
                quote=F,
                row.names = F,
                col.names = T)
    
    write.table(A_string,
                file = paste0(prot_dir,"STRING_hi.7",A_search,".tsv"),
                sep="\t",
                quote=F,
                row.names = F,
                col.names = T)
    
    write.table(A_grid,
                file = paste0(prot_dir,"BIOGRID_",A_search,".tsv"),
                sep="\t",
                quote=F,
                row.names = F,
                col.names = T)
  
    write.table(A_map,
                file = paste0(prot_dir,"HumanCellMap_",A_search,".tsv"),
                sep="\t",
                quote=F,
                row.names = F,
                col.names = T)
    
    ## combined table grepped
    
    write.table(out,
                file = paste0(prot_dir,"filtered_",geneGrep,"_", A_search,".tsv"),
                sep="\t",
                quote=F,
                row.names = F,
                col.names = T)
    
  }
  return(counts)
}



################
#### inputs ####
################


base = "C:/Users/dowens/OneDrive/Postdoc/Projects/interactR/"

dat_dir = paste0(base,"data/")
dir.create(dat_dir,showWarnings = F)


# load the saved session image containing data objects
load(paste0(dat_dir,"InteractR_data.Rdata"))



## load new HGNC map

# downloaded April 2022 https://github.com/waldronlab/HGNChelper/blob/master/data/hgnc.table.rda
load(file=paste0(dat_dir,"hgnc.table.rda"))

hgnc.table %<>%
  dplyr::select(1,2)


## tissue info data
tiss_file = paste0(dat_dir, "sample_info.csv")

tiss = fread(tiss_file)



##########################
#### PREANALYZED DATA ####
##########################

## protein
prot = readRDS(paste0(dat_dir,"protCounts_all_long_uniprotIDs_normed_geneSymbols_rmDup.Rds"))


## RNA
rna = readRDS(paste0(dat_dir,"RNACounts_v2_220613.Rds"))

# table of interactors ran March 2023
main_table = readRDS(paste0(dat_dir,"main_table.Rds"))




######## proteins to search ##########

## main A protein to get interactors for
#A_search = c("USP5")

## optional B proteins that are hypothesised to interact with A
## this helps to narrow down the interactome of A to just these proteins of interest 
## not implemented yet
#B_search = c("SMARCA2","SMARCA4")


## run on all protein coding genes (aka proteins) in ENSEMBL grch38
A_search =
  grch38 %>%
  filter(biotype=="protein_coding") %>%
  pull(symbol) %>%
  unique()

## correct names
A_search =
  HGNChelper::checkGeneSymbols(A_search,map=hgnc.table) %>%
  dplyr::select(Suggested.Symbol) %>%
  drop_na() %>%
  distinct() %>%
  separate(Suggested.Symbol,into="Suggested.Symbol2",sep=" ///") %>% # HGNC when amibguous names detected uses  /// delemiter which breaks file name so just sep and skip the second one (quick fix)
  filter(!grepl("\\:",Suggested.Symbol2)) %>% # also have some like HGNC:18790 which kills the file naming
  dplyr::pull(Suggested.Symbol2)


# filtering : didnt work so remove another way

# Find the indices of elements containing ":"
#idx = grep("\\:", A_search)

# Subset the vector to remove elements containing ":"
#A_search = A_search[-idx]
  


#################
#### outputs ####
#################


out_dir = paste0(base,"interactR_results2/")
dir.create(out_dir,showWarnings = F)


  

############################
#### InteractR Searches ####
############################


## single protein
target = "FTSJ3"
# use geneGrep to pull interactors of a certain type
geneGrep = "WDR"

out = interactR(target,hut_dat,histring,grid_dat,cellmap,export=T,geneGrep)




## run just on a subset that are not yet done
## by listing the files in the output folder
## and excluding them as these are already processed
#done = list.files(out_dir,include.dirs = T)
#A_search=A_search[!A_search %in% done]

#set up empty table
main_table = tibble()
i=100
prog = 0 # initialize progress as 0
for(i in 1:length(A_search)){
  
  prot=A_search[i]
  
  ## display percentage progress
  prog2 = round(i/length(A_search),3)*100
  
  if(prog2>prog){
    cat(prog2,"%","\n")
  }
  #set prog to prog2 to stop overwrite
  prog=prog2

  temp_out = interactR(prot,hut_dat,histring,grid_dat,cellmap,export=F)
  
  # add to main table
  main_table %<>% rbind.data.frame(temp_out)
}

#saveRDS(main_table,file=paste0(dat_dir,"main_table.Rds"))

#main_table = readRDS(paste0(dat_dir,"main_table.Rds"))



## add a totals column
main_table %<>%
  mutate(total=Hut+Grid+CellMap)

### exploring the data
counts =
  main_table %>%
  group_by(total) %>%
  dplyr::count()

counts %>%
  mutate(type=if_else(total>0,"light","dark")) %>%
  group_by(type) %>%
  summarise(total_type=sum(n))

dark_prots =
  main_table %>%
  filter(total==0)


write.table(dark_prots,
            file=paste0(base,"Dark_interactome_proteins.txt"),
            quote=F,
            sep="\t",
            col.names = T,
            row.names = F)


counts2 =
  counts %>%
  mutate(group = case_when(
    total<10~as.character(total),
    total>9&total<20~("10 to 20"),
    total>19&total<100~("20 to 100"),
    total>99~("100+")
  )) %>%
  group_by(group) %>%
  summarise(count_in_group = sum(n))


counts2$group=factor(counts2$group,levels=c(as.character(0:10),"10 to 20","20 to 100","100+"))


ggplot(counts2,aes(group,count_in_group,fill=group))+
  geom_bar(stat="identity",colour="black",size=0.1)+
  theme_bw()+
  labs(y="Number of proteins",x="Reported PPIs",fill=NULL)+
  #scale_fill_brewer(palette = "Greys")+
  scale_fill_grey(start=0,end=0.95)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave(paste0(base,"interactor_plots.pdf"),
       width=6,
       height=4)



############################################
#### QUANTIFY PROTEIN LEVEL IN LC-MS/MS ####
############################################

## uniprot code for filtering proteomics data on
target_id =   "P18206"
#"Q8NEZ5" # FBXO22
#"P60709" # ACTIN
#"Q8IVV7" # GID4
#"P18206" # Vinculin

#  hgnc_symbol
target_name = "FBXO22"

#"GID4"

## get total number of samples and number that detected our protein
s_nums =
  prot %>%
  dplyr::select(stripped_cell_line_name,Study.x) %>%
  distinct() %>%
  group_by(Study.x) %>%
  dplyr::count() %>%
  rename(total=n)

## filter either on uniprot or on hgnc_symbol
this_prot =
  prot %>%
  filter(hgnc_symbol=="GID4")
#filter(uniprot_id==target_id)

## quantify the numbers of samples
prot_quant =
  this_prot %>%
  dplyr::select(stripped_cell_line_name,Study.x) %>%
  distinct() %>%
  group_by(Study.x) %>%
  dplyr::count() %>%
  rename(this_prot=n) %>%
  left_join(s_nums) %>%
  mutate(perc=this_prot*100/total)


this_rna =
  rna %>%
  filter(Gene==target_name)


## for protein, just pick the main normal tissue study
this_prot %<>%
  filter(grepl("Goncalves|Nusinow|Jiang",Study.x))


