# Summarize compleasm logs and plot per-species completeness

file_list<-list.files(path = "../bat1k_compleasm/", pattern = "log")

all_spp<-list()
for(q in 1:length(file_list)){
  df<-read.table(file=paste0("../bat1k_compleasm/",file_list[q]), sep = "\t")
  df2<-as.data.frame(t(df[c(1,5:8),]))
  df3<-as.data.frame(t(gsub("*.:","",gsub("%.*","",df2))))
  colnames(df3)<-c("SPP","S","D","F","M")
  all_spp[[q]]<-df3
}
all_spp_df<-do.call(rbind, all_spp)
all_spp_df[, 2:5] <- lapply(all_spp_df[, 2:5], as.numeric)
all_spp_df$C<-all_spp_df$S + all_spp_df$D
all_spp_df <- all_spp_df[, c(1, 6, 2:5)]

write.table(file="../bat1k_compleasm/all_species_compleasm_stats.txt", all_spp_df,
            quote = F, row.names = F, sep = "\t")

# Format for plot:
library(ggpubr)
library(viridis)

lplot_df<-list()
for(s in 1:30){
  df<-all_spp_df[all_spp_df$SPP == levels(factor(all_spp_df$SPP))[s],]
  df2<-as.data.frame(t(df))
  df2$SPP<-df2[1,1]
  df2$SCORE<-rownames(df2)
  df3<-df2[c(3:6),]
  colnames(df3)[1]<-"NUM"
  lplot_df[[s]]<-df3
}

plot_df<-do.call(rbind, lplot_df)
plot_df$NUM<-as.numeric(plot_df$NUM)
plot_df$SCORE <- factor(plot_df$SCORE, levels = c("S", "D", "F", "M"))
plot_df_or<-plot_df %>%
  mutate(SPP = factor(SPP, levels = species_order)) %>%  # Set factor levels for ordering
  arrange(SPP)  # Order the dataframe by the new factor levels
plot_df_or<-plot_df_or[plot_df_or$SPP != "HLrhiSin2",]

bp1<-ggbarplot(plot_df_or, 
          x = "SPP", 
          y = "NUM", 
          fill = "SCORE",  # Fill by SCORE
          color = "black",  # Outline color of the bars
          position = position_stack(),  
          palette = c("#31688EFF","#35B779FF","#FDE725FF","#440154FF"),  # Using viridis with 4 categories
          lab.size = 10,  
          xlab = "Species",  
          ylab = "Score",
          orientation = "horiz")
bp1

ggsave(filename = "../bat1k_compleasm/all_species_compleasm.pdf", plot = bp1, width = 12, height = 14)
