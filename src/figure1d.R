##################################################################
#Figure 1D - Percentage mutation rate by cancer type             #
#Merge of SCC exome sequencing with TCGA data from cBioportal    #
#See materials and methods for download and preprocessing instr. #
##################################################################

library(ggplot2)
Mutation_percentages_df = readRDS('Mutation_percentages_df.RDS')

Mutation_percentages_df$Percent_Mutated[which((Mutation_percentages_df$Percent_Mutated == 0))] = NA
cols <- rev(c(rep('darkred', 5), rep('red', 10),rep('orange', 10), rep('yellow', 10), 'deepskyblue1','navyblue'))
p <- ggplot(Mutation_percentages_df, aes(Gene, Cancer_Type)) + geom_tile(aes(fill = Percent_Mutated, height = 1, width = 1), color = 'white') + scale_fill_gradientn(na.value = "navyblue", colours = cols, breaks=seq(0,100,5), limits = c(0,100)) + theme(axis.text.x = element_text(angle = 90, hjust = 1), strip.background = element_blank(), strip.text.y = element_blank()) +
  facet_grid(Group~ Group2, scales ="free", space="free") +theme(legend.key.height = unit(1.5,'cm')) + ylab('Cancer Type') + xlab('Gene')
p
jpeg(file = paste('~/', 'Heatmap1D_rainbow_2' , '.jpeg', sep = ''),width = 6000, height = 2000, res=300)
print(p)
dev.off()

pdf(file = paste('~/', 'Heatmap1D_rainbow_2' , '.pdf', sep = ''), width = 30, height = 10)
print(p)
dev.off()