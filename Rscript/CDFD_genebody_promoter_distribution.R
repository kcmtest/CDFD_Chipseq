require(reshape2)
require(ggplot2)

###############################
#df1 <- df[,-1]
#df1 <- log2(df1+1)
#df2 <- df[c(1)] 
#df3 <- cbind.data.frame(df2,df1)
##############################
#d <- read.csv('box_plot_data.txt',row.names = 1)
#dim(d)
#head(d)
#library(reshape)
#df <- apply(d,2,log2)
#head(df)
#write.csv(df,"BOX_LOG2_BOX.txt",row.names = TRUE,quote = F)

df1 <- read.csv('Genebody.txt',sep = "\t")
df1 <- Promoter
df1 <- genebody
df <- df1
df <- log2(df1 + 1)
head(df)
names(df)[1] = "Ctrl"
names(df)[2] = "shp1_ko"

ex <- melt(df, id.vars=c(NULL))
head(ex)
#colnames(ex) <- c("gene","group","exprs")
colnames(ex) <- c("group","exprs")


head(ex)


p <- ggplot(data=ex, aes(x=group, y=exprs)) +
  
  geom_boxplot(position=position_dodge(width=0.5),
               width=.8,
               outlier.shape=17, outlier.colour="red", 
               outlier.size=0.5, aes(fill=exprs)) +
  
  scale_fill_manual(values=c("red", "royalblue")) + #for boxplot
  stat_boxplot(geom="errorbar", width=.1)+
  #stat_summary(fun.y=median, colour="red", geom="point", aes(group = 1))+
  
  
  
  #geom_jitter(position=position_jitter(width=0.001), size=0.01, colour="blue") +
  
  #Choose which colours to use; otherwise, ggplot2 choose automatically
  #scale_color_manual(values=c("red3", "white", "blue")) + #for scatter plot dots
  #scale_fill_brewer() + #for boxplot
  
  #Add the scatter points (treats outliers same as 'inliers')
#geom_jitter(position=position_jitter(width=0.1), size=.5, colour="blue") +

#Set the size of the plotting window
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
  xlab("Cell Type") +
  ylab("Expression")
p

p+geom_violin(scale = "width")

library("ggpubr")
pdf("Genebody.pdf",20,15)

head(ex)
ggboxplot(ex, x = "group", y = "exprs",
          color = "group") +
  stat_compare_means(method = "wilcox.test",
                     label.x = 1.2, 
                     label.y = 10)+
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
  xlab("") +
  ylab("rpkm")

dev.off()


##################################################################GYAN##############
#colors <- c("brown", "maroon1", "black","darkgreen") # Number of color according to the number of groups
require(reshape2)
require(ggplot2)
require(RColorBrewer)
#require(dplyr)
require(tidyr)
require(tibble)

df <- read.csv("Fig2a input.txt",sep = ",")
head(df)


#df <- as.data.frame(Disease_module_2_box_plot_input)

Delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}
df <- Delete.na(df)

class(df)

head(df)

#df_melt=melt(df,id.vars=NULL)
df_melt = melt(df,id.vars = NULL)
#df_melt1 = melt(df1,id.vars = "class")

#head(df_melt1)
head(df_melt)
tail(df_melt)
###########################################################
df_melt$Group <- gsub("Ctrl_",'',df_melt$variable)
head(df_melt)
df_melt$Group <- gsub("PTPN6_KO_",'',df_melt$Group)
#df_melt$Group <- gsub("Blast",'', df_melt$Group)
head(df_melt)
###################################################################





factor(df_melt$Group)

head(df_melt)

df_melt$Group <- factor(df_melt$Group,levels = "genebody","Promoter")


head(df_melt$variable)

pdf("Fig2.pdf",20,15)

ggplot(df_melt, aes(variable,value)) +
  stat_boxplot(geom="errorbar", width=.5)+
  geom_boxplot(aes(fill=Group), position = position_dodge2(preserve = "total"))+
  #scale_fill_manual(values=c("grey","dodgerblue4","red","darkgreen","brown","maroon1")) + #for boxplot

  #scale_fill_manual(values=c("grey","dodgerblue4","red","darkgreen","brown","maroon1","Purple",
   #                          "turquoise","sienna","slate gray","peru","deep pink","cyan")) + #for boxplot
  
  
  scale_fill_manual(values=c("grey","dodgerblue4","red","darkgreen","brown","maroon1","Purple",
                             "turquoise","sienna","peru")) + #for boxplot
  
  
    # scale_fill_manual(values=c("grey","dodgerblue4","red","darkgreen","brown")) + #for boxplot
  geom_jitter(position=position_jitter(width=0.3), size=0.01, colour="black") +
  
  #scale_fill_manual(values = c("black","brown","maroon1","darkgreen"))+
  #stat_summary(fun.y=mean, colour="red", geom="line", aes(group = 1))+
  #geom_jitter(position = position_jitter(0.0001),size=.00001)+
  #theme_bw(base_size=15)+
  theme_bw()+
  
  theme(axis.text.x=element_text(angle = 90, size=20, face="bold", hjust = 1), 
        axis.text.y=element_text(angle=0, size=35, face="bold", vjust=0.5),
        plot.title = element_text(size=40, face="bold"), 
        legend.title=element_blank(), 
        legend.key.size=unit(1, "cm"),      #Sets overall area/size of the legend
        legend.text=element_text(size=30)) +
  
  labs(x="Cell_Type", y=expression(paste("Expression"), title="Expression")) + 
  guides(colour = guide_legend(override.aes = list(size=5)))

dev.off()

#####################################################
