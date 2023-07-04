## interactR

## perform searches for interactions of a protein of interest

## so far covering

## STRING v11.5 physical links network, experimental evidence only (not text mining) https://string-db.org/
## Human cell map v1 (BFDR threshold=0.1) https://humancellmap.org/
## BioGRID 4.4.219 https://thebiogrid.org
## BioPlex 1.0, 2.0, & 3.0 https://bioplex.hms.harvard.edu/
## OpenCell https://opencell.czbiohub.org/download
## Human Reference interactome http://www.interactome-atlas.org/
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

# full list of interactors including metadata ie publication source
interactR_full = function(protein,ints_db,experiment_type="all",out_path=NA){
  
  ## function that takes a character vector of HGNC-compliant gene names
  ## and the database file
  ## default will return data as a tibble
  ## or if out_dir is set to a path, then will save a .tsv file
  
  for(target in protein){
    
    this_ints =
      ints_db %>%
      filter((A == target | B == target)) 
    
    # save this table if entries and if not returning data
    if(nrow(this_ints)>0&!is.na(out_path)){
      
      write.table(this_ints,
                  file = paste0(out_path,"InteractR_",target,"_full.tsv"),
                  sep="\t",
                  quote=F,
                  row.names = F,
                  col.names = T)
    }
  }
  
  return(this_ints)
  
}

# simple non-redundant list of just interactors
interactR_simple = function(protein,ints,experiment_type="all",out_path=NA){
  
  ## function that takes a character vector of HGNC-compliant gene names
  ## and the database file
  ## default will return data as a tibble
  ## or if out_dir is set to a path, then will save a .tsv file
  
  for(target in protein){
    
    this_ints =
      ints %>%
      filter((A == target | B == target))
    
    ## clean up to just have the interactors
    this_ints = unique(c(this_ints$A,this_ints$B))
    # remove itself from the list
    this_ints=sort(this_ints[this_ints!=target])
    
    # make a tibble and give nice colname
    this_ints = as_tibble(this_ints)
    colnames(this_ints) = paste0(target,"_interactors")
    
    # save this table if entries and if not returning data
    if(nrow(this_ints)>0&!is.na(out_path)){
      
      write.table(this_ints,
                  file = paste0(out_path,"InteractR_",target,"_simple.tsv"),
                  sep="\t",
                  quote=F,
                  row.names = F,
                  col.names = T)
    }
  }
  
  return(this_ints)
  
}


interactr.list.types = function(){
  ## function to list valid experiment_type options to search for in interactR
  ## mostly comes from biogrid classification
  
  types=c(
    "Proximity Label-MS","Two-hybrid","IP-MS","Affinity Capture-Luminescence","Affinity Capture-Western","Biochemical Activity",
    "Co-crystal Structure","Far Western","FRET",
    "Reconstituted Complex","Co-localization","Protein-peptide",
    "PCA","Affinity Capture-RNA","Protein-RNA",
    "Co-purification","Co-fractionation","unknown"
  )
  return(types)
}

################
#### inputs ####
################


base = "C:/Users/imnot/OneDrive/Postdoc/Projects/interactR/"

dat_dir = paste0(base,"data/")
dir.create(dat_dir,showWarnings = F)


# load the saved combined interactomes data WITH DB INFO
ints_db = readRDS(file=paste0(dat_dir,"InteractR_interactome_full.Rds"))

#write.table(ints_db,
#            file = paste0(base,"Interactome_full.tsv"),
#            quote=F,
#            row.names = F,
#            col.names = T,
#            sep="\t")

### load the non-redundant combined interactomes data
ints = readRDS(file=paste0(dat_dir,"InteractR_interactome.Rds"))


# table of interactors ran March 2023
int_nums = readRDS(paste0(dat_dir,"int_nums.Rds"))


## load new HGNC map
# downloaded April 2022 https://github.com/waldronlab/HGNChelper/blob/master/data/hgnc.table.rda
load(file=paste0(dat_dir,"hgnc.table.rda"))

hgnc.table %<>%
  dplyr::select(1,2)


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


out_dir = paste0(base,"interactR_2.0_results/")
dir.create(out_dir,showWarnings = F)



############################
#### InteractR Searches ####
############################

#######
## 1 ##
#######

## search for a single protein and save the output

## save output text file for each protein

target = A_search[1]

## on the full table
interactR_full("USP5",
               ints_db,
               out_path=out_dir)


## or on the simplified table
interactR_simple("USP5",
                 ints,
                 out_path=out_dir)



## don't save the output
interactR_full("STARD5",
               ints_db)




#######
## 2 ##
#######

## search and just keep the number of interactors for each protein

## return a df

#set up empty table to fill
int_nums = tibble()
i=1
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
  
  this_ints =
    ints %>%
    filter((A == prot | B == prot)) %>%
    nrow()
  
  temp_out = tibble(protein = prot,
                    int_num = this_ints)
  
  # add to main table
  int_nums %<>% rbind.data.frame(temp_out)
}

#saveRDS(int_nums,file=paste0(dat_dir,"int_nums.Rds"))

## load the precalculated interactor numbers
int_nums = readRDS(paste0(dat_dir,"int_nums.Rds"))


## export the list of dark interactome proteins to look at later
## defining as 0 or 1 unique interactors
dark_prots =
  int_nums %>%
  filter(int_num<2) %>%
  pull(protein)

#saveRDS(dark_prots,file=paste0(dat_dir,"dark_prots.Rds"))


write.table(dark_prots,
            file=paste0(dat_dir,"Dark_interactome_proteins.txt"),
            quote=F,
            sep="\t",
            col.names = T,
            row.names = F)



###################################

# EXPLORATORY PLOTTING #

## get the counts in each bin num of interactors
counts =
  int_nums %>%
  group_by(int_num) %>%
  dplyr::count()


counts %>%
  mutate(type=if_else(total>0,"light","dark")) %>%
  group_by(type) %>%
  summarise(total_type=sum(n))



counts2 =
  counts %>%
  mutate(group = case_when(
    int_num<10~as.character(int_num),
    int_num>9&int_num<20~("10 to 20"),
    int_num>19&int_num<100~("20 to 100"),
    int_num>99~("100+")
  )) %>%
  group_by(group) %>%
  summarise(count_in_group = sum(n))


counts2$group=factor(counts2$group,levels=c(as.character(0:10),"10 to 20","20 to 100","100+"))


ggplot(counts2,aes(group,count_in_group,fill=group))+
  geom_bar(stat="identity",colour="black",size=0.1)+
  theme_bw()+
  labs(y="Number of proteins",x="Reported PPIs",fill=NULL)+
  #scale_fill_brewer(palette = "Greys")+
  scale_fill_grey(start=0,end=1)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        text=element_text(size=6),
        panel.border=element_rect(size=0.1),
        axis.ticks = element_line(size=0.1))

ggsave(paste0(base,"interactor_plots.pdf"),
       width=2.5,
       height=2)


