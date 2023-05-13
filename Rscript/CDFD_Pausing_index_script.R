library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
x<- CTRL_PTPN6_Pol_II_pausing

x <- Ctrl_shp1_pol2

df = x %>% select(contains("pausing") | contains("Transcript")) %>% 
  rename_at(vars(1),substr,1,10) %>% 
  pivot_longer(-Transcript,names_to = "condition") %>% 
  mutate(condition=case_when(
    grepl("Ctrl",condition) ~ "Control",
    grepl("shp1",condition) ~ "shp1"
  )) %>%
  group_by(Transcript,condition) %>% summarize_all(median) 

####################################################################################


df = x %>% select(contains("pausing") | contains("Transcript")) %>% 
  rename_at(vars(1),substr,1,10) %>% 
  pivot_longer(-Transcript,names_to = "condition") %>% 
  mutate(condition=case_when(
    grepl("Control",condition) ~ "Control",
    grepl("PTPN6",condition) ~ "PTPN6_KO"
  )) %>%
  group_by(Transcript,condition) %>% summarize_all(median) 


########################################################################################






pdf("Pausing_Index_new.pdf",20,15)
p <- ggboxplot(df, x = "condition", y = "value",
               color = "condition", palette =c("#c72e29","#016392","#be9c2e","#098154")
)+
stat_compare_means(comparisons=list( c("Control", "shp1")))+
  theme_bw(base_size=15) +
  
  
  #Modify various aspects of the plot text and legend
  theme(
    legend.position="none",
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=14, face="bold", vjust=1),
    
    axis.text.x=element_text(angle=45, size=35, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=35, face="bold", vjust=0.5),
    axis.title=element_text(size=24, face="bold"),
    
    #Legend
    legend.key=element_blank(),     #removes the border
    legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
    legend.text=element_text(size=2),  #Text size
    title=element_text(size=12)) +      #Title text size
  
  #Change the size of the icons/symbols in the legend
  guides(colour=guide_legend(override.aes=list(size=1))) +
  
  #Set x- and y-axes labels
  xlab("Group") +
  ylab("Pausing Index")

p
dev.off()


Paused_genes = df[df$value > 5,]
Active_genes = df[df$value < 0.5,]

A <- Paused_genes[c(1)]
Paused_genes_unique = unique(A)

B <- Active_genes[c(1)]
Active_genes_unique = unique(B)
write_tsv(Active_genes_unique,"Active_genes_unique.txt")
