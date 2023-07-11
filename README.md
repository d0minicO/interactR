# interactR v1.1
Tool for exploring human protein-protein interaction (interactome) databases

Code is provided as a supplement to the shiny app https://dominico.shinyapps.io/interactR/

Different experiment types can easily be excluded to help refine interactomes

#### Updated June 28 2023

## Data from the following sources:

- STRING v12.0 physical links network, experimental evidence only (not text mining) https://string-db.org/

- Human cell map (BFDR threshold=0.1) https://humancellmap.org/

- BioGRID v4.4.222 physical interactions https://thebiogrid.org

- BioPlex (IP-MS) https://bioplex.hms.harvard.edu/

- OpenCell (IP-MS) https://opencell.czbiohub.org/download

- Human Reference interactome (Y2H) http://www.interactome-atlas.org/

- HIPPIE v2.3 http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php)

- Interpro human domains annotation https://www.ebi.ac.uk/interpro/


---

[HGNChelper 0.8.1](https://github.com/waldronlab/HGNChelper) is used to correct gene names to allow different datasets to be harmonized

[annotables 0.2.0](https://github.com/stephenturner/annotables) is used to map gene IDs to gene symbols