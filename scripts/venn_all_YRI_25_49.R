library(gridExtra)
library(grDevices)

YRI_abyss=2876134;
TSI_abyss=1799234;
YRI_TSI=TSI_abyss-990847;
all=2658964;
all_abyss=YRI_abyss+TSI_abyss-YRI_TSI;
all_all_abyss=all-523549; #128767 k-mer
YRI_all=YRI_abyss-1361077;
TSI_all=TSI_abyss-696649;
all_minus_YRI=1005422;
YRI_minus_all=1361077;

tiff("venn_all_YRI_abyss_25_49.tiff", width = 4, height = 4, units = 'in', res = 300);#pointsize=12);

venn.plot<- draw.pairwise.venn(all, YRI_abyss, YRI_all,fill =NULL, lwd=rep(10,2),col= c("azure4","forestgreen"),cat.cex=c(.6,.6),cex=1);

#venn.plot<- draw.pairwise.venn(all, YRI_abyss, YRI_all, c("Genotypes", "YRI"),fill = c("azure4","forestgreen"),cat.cex=c(.6,.6),cex=.6);
#venn.plot <- draw.pairwise.venn(100, 100, 30, c("YRI", "TSI"),fill = c("darkblue", "azure4"));
grid.draw(venn.plot);

print( grid.arrange(gTree(children=venn.plot))); 

venn.plot[[5]]$label  <- "47.3%";	#1361077/2876134  
# in baa only
venn.plot[[6]]$label <- "37.8%";  #1005422/2658964
# intesection
venn.plot[[7]]$label <- "";  

# plot  
grid.newpage()
grid.draw(venn.plot)

dev.off();