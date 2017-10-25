library(gridExtra)

YRI_abyss=2876134;
TSI_abyss=1799234;
YRI_TSI=TSI_abyss-990847;
all=2658964;
all_abyss=YRI_abyss+TSI_abyss-YRI_TSI;
all_all_abyss=all-523549; #127051 k-mer
YRI_all=YRI_abyss-1361077;
TSI_all=TSI_abyss-696649;

tiff("venn_all_all_abyss_25_49.tiff", width = 4, height = 4, units = 'in', res = 300);#pointsize=12);
venn.plot<- draw.pairwise.venn(all, all_abyss, all_all_abyss,fill=NULL, lwd=rep(10,2),col= c("azure4","navy" ),cat.cex=c(.6,.6),cex=1);
#venn.plot<- draw.pairwise.venn(all, all_abyss, all_all_abyss, c("Genotypes", expression(YRI*union(TSI))),fill = c("azure4","navy" ),cat.cex=c(.6,.6),cex=1);
#venn.plot <- draw.pairwise.venn(100, 100, 30, c("YRI", "TSI"),fill = c("darkblue", "azure4"));
grid.draw(venn.plot);

print( grid.arrange(gTree(children=venn.plot))); 

venn.plot[[5]]$label  <- "42.0%";  #3684521-2135415/3684521
# in baa only
venn.plot[[6]]$label <- "19.7%";  #523549/2658964
# intesection
venn.plot[[7]]$label <- "";  

# plot  
grid.newpage()
grid.draw(venn.plot)

dev.off();