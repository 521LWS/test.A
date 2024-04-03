library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(magrittr)
library(ggnewscale)
library(ggtreeExtra)
library(RColorBrewer)


data <- read_tsv("data.xls") %>% select(-Group) %>% 
  group_by(id,Phylum) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm=TRUE))) %>% 
  pivot_wider(names_from = "Phylum",values_from = "Abundance")

exp <- data %>% pivot_longer(-id)

tree <- data %>% column_to_rownames(var="id") %>% 
  dist() %>% ape::bionj()

group <- read_tsv("data.xls") %>% select(id,Group) %>% 
  mutate(group="group")
  
ggtree(tree,branch.length = "none", layout = "circular",
       linetype = 1,size = 0.5, ladderize = T)+
  layout_fan(angle =180)+
  theme(plot.margin=margin(0,1,-7,0,"cm"))+
  geom_tiplab(offset=10.5,show.legend=FALSE,size=2.8,
              color = "black",starstroke = 0)+
  geom_fruit(data=exp,geom=geom_tile,
             mapping=aes(y=id,x=name,fill=value),
             pwidth=0.6,offset=0.02,
             axis.params=list(axis="x",text.angle=-90,text.size=2,hjust=0))+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdBu")))+
  new_scale_fill()+
  geom_fruit(data=group,geom=geom_tile,
             mapping=aes(y=id,x=group,fill=Group),color="white",
             pwidth=1,offset=0.4)+
  scale_fill_manual(values = c("#EDB749","#3CB2EC","#9C8D58"))


