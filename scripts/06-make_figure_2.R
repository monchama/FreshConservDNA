# code for figure 2, Monchamp et al.  FreshConservDNA project
# vincent fugere 2022

rm(list=ls())

# packages and session information

library(tidyverse)
library(ggtree)

# print(sessionInfo())
# 
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggtree_3.2.1    forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7     purrr_0.3.4     readr_2.0.1    
# [7] tidyr_1.1.3     tibble_3.1.4    ggplot2_3.3.5   tidyverse_1.3.1
# 
# loaded via a namespace (and not attached):
#   [1] treeio_1.18.1      tidyselect_1.1.1   haven_2.4.3        lattice_0.20-44    ggfun_0.0.6       
# [6] colorspace_2.0-2   vctrs_0.3.8        generics_0.1.3     gridGraphics_0.5-1 utf8_1.2.2        
# [11] rlang_0.4.11       pillar_1.6.2       glue_1.4.2         withr_2.4.2        DBI_1.1.1         
# [16] dbplyr_2.1.1       modelr_0.1.8       readxl_1.3.1       lifecycle_1.0.0    munsell_0.5.0     
# [21] gtable_0.3.0       cellranger_1.1.0   rvest_1.0.1        tzdb_0.1.2         parallel_4.1.1    
# [26] fansi_0.5.0        broom_0.7.9        Rcpp_1.0.9         scales_1.1.1       backports_1.2.1   
# [31] jsonlite_1.7.2     fs_1.5.0           hms_1.1.0          aplot_0.1.4        stringi_1.7.4     
# [36] grid_4.1.1         yulab.utils_0.0.4  cli_3.0.1          tools_4.1.1        magrittr_2.0.1    
# [41] lazyeval_0.2.2     patchwork_1.1.1    crayon_1.4.1       ape_5.6-2          pkgconfig_2.0.3   
# [46] tidytree_0.3.9     ellipsis_0.3.2     ggplotify_0.1.0    xml2_1.3.2         reprex_2.0.1      
# [51] lubridate_1.7.10   assertthat_0.2.1   httr_1.4.2         rstudioapi_0.13    R6_2.5.1          
# [56] nlme_3.1-152       compiler_4.1.1    

## data

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sp.list <- read_tsv('./output/2022-06-18_species_list.tsv')
dat.full <- read_tsv('./output/2022-06-18_bold_records.tsv')

colors <- ochRe::ochre_palettes[['lorikeet']] #colors for plotting

## create species-level dataset

dat.s <- dat.full %>% filter(!is.na(recordID)) %>% 
  group_by(species) %>% summarize(records=n())

dat.s <- left_join(sp.list, dat.s, by = 'species')

#need a list of phyla to fill missing records in full species list, where species with no bold records have no phyla
phylum.list = dat.full %>% select(class,phylum_name) %>% distinct(class, .keep_all = T) %>% rename('phylum' = phylum_name)
phylum.list[c(3,12),'phylum'] <- c('Chordata','Chordata')

dat.s <- left_join(dat.s,phylum.list)
dat.s <- dat.s %>% select(taxon, phylum, class:records)
dat.s$records[is.na(dat.s$records)] = 0

#add a column indicating if species has at least 3 records in BOLD
dat.s$represented = 0
dat.s$represented[dat.s$records > 4] = 1

#add a column for binary conservation status (at risk or not)
dat.s$at.risk <- 1
dat.s$at.risk[dat.s$assessment == 'Data Deficient'] <- 0
dat.s$at.risk[dat.s$assessment == 'Not At Risk'] <- 0
dat.s$at.risk[dat.s$assessment == 'Not Available'] <- 0

#add a column combining previous two metrics, i.e. species is both a risk & has 3+ BOLD records 
dat.s$at.risk.with.seq = 0
dat.s$at.risk.with.seq[dat.s$at.risk == 1 & dat.s$represented == 1] = 1

#add columns for threat data status (data deficient or not) & genetic resources for DD species
dat.s$data.def <- 0
dat.s$data.def[dat.s$assessment == 'Data Deficient'] <- 1
dat.s$data.def.with.seq = 0
dat.s$data.def.with.seq[dat.s$data.def == 1 & dat.s$represented == 1] = 1

## compute various order-level metrics for plot

dat.o <- dat.s %>%
  group_by(phylum,class,order) %>% 
  summarize(species = n_distinct(species),
            orecords=sum(records),
            median.records = median(records),
            species.with.rec = sum(represented), 
            species.at.risk = sum(at.risk),
            at.risk.represented = sum(at.risk.with.seq),
            data.deficient.species = sum(data.def),
            data.deficient.represented = sum(data.def.with.seq)) %>%
  drop_na(order) %>%
  mutate(logspecies = log10(species+1), # panel a, bars
         records.per.species = orecords/species, 
         log.records.per.species = log10(records.per.species+1), # panel a, points
         prct.represented = 100*(species.with.rec/species), # panel b
         prct.at.risk = 100*(species.at.risk/species), # panel c, light bars
         prct.at.risk.with.seq = 100*(at.risk.represented/species), # panel c, dark bars
         prct.dd = 100*(data.deficient.species/species), #panel d, light bars
         prct.dd.with.seq = 100*(data.deficient.represented/species)) %>% #panel d, bark bars
  ungroup()

## make a taxonomy-based tree

taxo <- dat.o %>% select(phylum,class,order)
taxo$kingdom <- 'Animalia'
taxo$kingdom[str_detect(taxo$phylum, 'phyta')] <- 'Plantae'
taxo$domain <- 'Eukaryota'
taxo$class[taxo$class == 'Teleostei'] <- 'Actinopterygii'
taxo$class[taxo$class == 'Chondrostei'] <- 'Actinopterygii'
taxo$class[taxo$class == 'Holostei'] <- 'Actinopterygii'
taxo <- taxo %>% mutate_all(as.factor)

frm <- ~domain/kingdom/phylum/class/order
tree <- ape::as.phylo(frm, data = taxo, collapse=FALSE)
tree$edge.length <- rep(1, nrow(tree$edge))

dat.o <- dat.o[match(tree$tip.label,dat.o$order),]
dat.o <- dat.o %>% ungroup %>% select(-phylum, -class)

## figure (7 X 6 pdf)

t <- ggtree(tree,color = 'gray80',layout="rect") + geom_tiplab(color='gray60',size=2) + 
  geom_nodelab(size=3, geom='text') + theme_tree2() + xlim_tree(6)

# I could not find a way to do everything in ggplot/ggtree. 
# I combine the following two pdfs in Inkscape, and then add a bit of 
# decoration (e.g. highlighting priority taxa)

#plot 1, tree of the right scale

p <- t + geom_facet(panel='# sp. (log)', data=dat.o, geom=geom_segment, 
                    mapping=aes(x=0, xend=logspecies, y=y, yend=y), size=2, color='gray80') +
  geom_facet(panel='seq. per sp. (log[x+1])', data=dat.o, geom=geom_point, 
             mapping=aes(x=log.records.per.species, y=y), size=1, shape=16, color=colors[2]) +
  geom_facet(panel='% with seq.', data=dat.o, geom=geom_segment, 
             mapping=aes(x=0, xend=prct.represented, y=y, yend=y), size=2, color=colors[4]) +
  geom_facet(panel='% at risk', data=dat.o, geom=geom_segment, 
             mapping=aes(x=0, xend=prct.at.risk, y=y, yend=y), size=2, color=alpha(colors[1],0.5)) +
  geom_facet(panel='at risk with seq.', data=dat.o, geom=geom_segment, 
             mapping=aes(x=0, xend=prct.at.risk.with.seq, y=y, yend=y), size=2, color=colors[1])

facet_labeller(p, c(Tree = "Taxa")) %>% facet_widths(widths=c(0.2,0.1,0.1,0.1,0.1,0.1))

#plot 2, weird tree but need to scale xlim = c(0,100) for DD panel

p <- t + geom_facet(panel='# sp. (log)', data=dat.o, geom=geom_segment, 
                    mapping=aes(x=0, xend=logspecies, y=y, yend=y), size=2, color='gray80') +
  geom_facet(panel='seq. per sp. (log[x+1])', data=dat.o, geom=geom_point, 
             mapping=aes(x=log.records.per.species, y=y), size=1, shape=16, color=colors[2]) +
  geom_facet(panel='% with seq.', data=dat.o, geom=geom_segment, 
             mapping=aes(x=0, xend=prct.represented, y=y, yend=y), size=2, color=colors[4]) +
  geom_facet(panel='% dd', data=dat.o, geom=geom_segment, 
             mapping=aes(x=0, xend=prct.dd, y=y, yend=y), size=2, color=alpha(colors[3],0.5)) +
  geom_facet(panel='dd with seq.', data=dat.o, geom=geom_segment, 
             mapping=aes(x=0, xend=prct.dd.with.seq, y=y, yend=y), size=2, color=colors[3]) + xlim(c(0,100))

facet_labeller(p, c(Tree = "Taxa")) %>% facet_widths(widths=c(0.2,0.1,0.1,0.1,0.1,0.1))

# taxa to outline in figure

dat.o %>% mutate(deficit = prct.at.risk - prct.at.risk.with.seq) %>% filter(deficit > 25) %>% pull(order)
dat.o %>% mutate(deficit = prct.dd - prct.dd.with.seq) %>% filter(deficit > 25) %>% pull(order)
