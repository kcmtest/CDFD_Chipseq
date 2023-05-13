library(ggplot2)
library(ggpmisc)

df1 <- read.csv('Promoter_analysis.txt',sep = "\t")


df <- log2(df1 + 1)
my.formula <- y ~ x

head(df)

pdf("Genebody_scatter.pdf",20,15)

p <- ggplot(data = df, aes(x = Ctrl_Promoter, y = PTPN6_KO_Promoter)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +         
  geom_point()+
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
  xlab("CTRL") +
  ylab("PTPN6")

p