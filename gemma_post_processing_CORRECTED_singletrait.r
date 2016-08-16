#GEMMA POST-GWAS CORRECTED

######################################################################
######################################################################
#download output files from FARM
#set the FARM directory names 
fibrpath="/home/maktaylo/fibr_gwas_maf10/"
traitname='plasticity_among_site'
pathdash="/"

#make an output directory on my machine
Sys.Date()
dir.create(file.path('/Users/mrktaylor531/Desktop/Projects/FIBR/gemma_results', paste(Sys.Date(),'_',traitname,sep='')))
outputdir = paste('/Users/mrktaylor531/Desktop/Projects/FIBR/gemma_results/',Sys.Date(),'_',traitname,sep='')

outputdir=paste(fibrpath,traitname,'/','output', sep='')

#figure out how to automatically download stuff from R in my machine

######################################################################
######################################################################
#now read in the output files

lmm = read.table(list.files(outputdir,pattern='lmm.assoc.txt',full.names=T), header=T)

head(lmm)
#manhattan plot of EMMA results

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7") #set plotting color palette
ggplot(data=lmm, aes(x = ps, y = -log10(p_wald), color=as.factor(chr))) + 
  geom_point(size=0.75) +
  facet_grid(chr ~ .) +
  xlab('position') +
  ylab('-log10(Wald test p-score)') +
  ggtitle(paste('LMM',traitname, sep=' ')) +
  theme(text = element_text(size=17)) +
  theme(axis.text=element_text(size=14,color='black'),
        axis.title=element_text(size=14,color='black')) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(plot.title = element_text(size = 24)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(10)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(legend.position = 'none') +
  theme(axis.line = element_line(colour = 'black')) +
  scale_colour_manual(values=cbbPalette) +
  scale_fill_manual(values=cbbPalette) +
  theme(panel.background = element_rect(colour = "black"))

#histogram of p-values
xint = quantile(-log10(lmm$p_wald), 0.99)

pwald_dist = ggplot(data=lmm, aes(x=-log10(p_wald))) + 
  geom_histogram(data=lmm, aes(x=-log10(p_wald)), alpha=0.8, binwidth = 0.1) +
  geom_vline(aes(xintercept = xint), colour="red", size=1.2) +
  xlab('Wald p-values') +
  ylab('') +
  ggtitle(paste(traitname, 'LMM', 'Wald')) +
  theme(axis.text=element_text(size=14,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(10)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(panel.background = element_rect(colour = "black"))

snp_num = 100 #how many of the top-associated SNPs to take for each trait!
trait='DTBclean' #write down the trait i want

################################################################################################
#read in GEMMA BSLMM results
#this loop pulls out the top 100 most associated SNPs

bslmmres = read.table(list.files(outputdir,pattern='output.param.txt',full.names=T), header=T)
head(bslmmres)

manplot_bslmm = ggplot(data=bslmmres, aes(x = ps, y = abs(gamma), color=as.factor(chr))) + 
  geom_point(size=0.75) +
  facet_grid(chr ~ .) +
  xlab('position') +
  ylab('absolute value (beta)') +
  ggtitle(paste('BSLMM beta', traitname,sep=' ')) +
  theme(text = element_text(size=17)) +
  theme(axis.text=element_text(size=14,color='black'),
        axis.title=element_text(size=14,color='black')) +
  theme(plot.title = element_text(size = 24)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(10)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(legend.position = 'none') +
  theme(axis.line = element_line(colour = 'black')) +
  scale_colour_manual(values=cbbPalette) +
  scale_fill_manual(values=cbbPalette) +
  theme(panel.background = element_rect(colour = "black"))

xintb = quantile(abs(bslmmres$beta), 0.99)

beta_dist = ggplot(data=bslmmres, aes(x=abs(beta))) + 
  geom_histogram(data=bslmmres, aes(x=abs(beta)), alpha=0.8) +
  geom_vline(aes(xintercept = xintb), colour="red", size=1.2) +
  xlab(paste(poo[1,1], 'beta',sep=' ')) +
  ggtitle(paste(poo[1,1], 'BSLMM', 'beta')) +
  theme(axis.text=element_text(size=14,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  #scale_x_continuous(breaks=number_ticks(10)) +
  #scale_y_continuous(breaks=number_ticks(10)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black'))

###compare LMM vs BSLMM beta results

head(lmm)
lmm_beta = ggplot(data=lmm, aes(x = ps, y = abs(beta), color=as.factor(chr))) + 
  geom_point(size=0.75) +
  facet_grid(chr ~ .) +
  xlab('position') +
  ylab(expression(paste('|',beta,'|',sep=''))) +
  ggtitle('LMM') +
  theme(text = element_text(size=17)) +
  theme(axis.text=element_text(size=14,color='black'),
        axis.title=element_text(size=14,color='black')) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(plot.title = element_text(size = 24)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(10)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(legend.position = 'none') +
  theme(axis.line = element_line(colour = 'black')) +
  scale_colour_manual(values=cbbPalette) +
  scale_fill_manual(values=cbbPalette) +
  theme(panel.background = element_rect(colour = "black"))

head(bslmmres)
bslmm_beta = ggplot(data=bslmmres, aes(x = ps, y = abs(beta), color=as.factor(chr))) + 
  geom_point(size=0.75) +
  facet_grid(chr ~ .) +
  xlab('position') +
  ylab(expression(paste('|',beta,'|',sep=''))) +
  ggtitle('BSLMM') +
  theme(text = element_text(size=17)) +
  theme(axis.text=element_text(size=14,color='black'),
        axis.title=element_text(size=14,color='black')) +
  theme(strip.background = element_rect(fill = 'white', color='black', size=1)) +
  theme(plot.title = element_text(size = 24)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(10)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(legend.position = 'none') +
  theme(axis.line = element_line(colour = 'black')) +
  scale_colour_manual(values=cbbPalette) +
  scale_fill_manual(values=cbbPalette) +
  theme(panel.background = element_rect(colour = "black"))

multiplot(lmm_beta, bslmm_beta, cols=2)

###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
#read in EMMA LMM results
#this loop pulls out the top 100 most associated SNPs
lmmres = do.call(rbind, lapply(1:length(lmm), function(i) {
  #i=1
  tryCatch({
    #i=1
    emmares = read.table(file=lmm[i], header=T)
    #re-order from most to least significant
    emmares.1 = emmares[order(emmares[14]), ]  
    emmares.2 = emmares.1[1:snp_num,] #take only top 100
    poo = data.frame(strsplit(as.character(lmm_names[i]), trait))
    data.frame(emmares.2, rep(poo[1,1],snp_num))
  }, error=function(e){})
}))

head(lmmres)
names(lmmres)[ncol(lmmres)]='site'

write.table(lmmres, file=paste(outputdir,'top_SNPs_LMM.txt',sep=''), row.names=F, quote=F,sep='\t')
lmmres = read.table(file=paste(outputdir,'top_SNPs_LMM.txt',sep=''), header=T)

lmmreps = do.call(rbind, lapply(1:nrow(lmmres), function(i) {
  #i=2
  poo = as.data.frame(match(lmmres$rs, lmmres$rs[i]))
  matchs = data.frame(table(poo))
  matchs = matchs[1,2]
  if (matchs > 1) {data.frame(lmmres[i,],matchs)}
}))
lmmreps

lmmreps = lmmreps[order(-lmmreps$matchs),]
head(lmmreps,15)

(length(unique(lmmreps$rs)) / length(unique(lmmres$rs)))*100 # = % of SNPs that overlap for a trait

################################################################################################
#read in GEMMA BSLMM results
#this loop pulls out the top 100 most associated SNPs

bslmm = 

###############################################################################################
###############################################################################################
#searh for overlaps between sites 
lmtraits = unique(lmmres$site)
bstraits = unique(bsres$site)

if (length(lmtraits) > length(bstraits)) {countnum = length(bstraits)} else {countnum=length(lmtraits)} 

#search for overlapping SNPs b/w LMM and BSLMM
overlapping_SNPs = do.call(rbind, lapply(1:countnum, function(i) {
  #i=4
  lematch = bstraits[match(lmtraits[i], bstraits)]
  
  if(is.na(match(lmtraits[i], bstraits))==F) { #if there is a matching trait in the bs/lmm
    lmtops = lmmres[lmmres$site==lmtraits[i],] #combine top snps for a trait between from LMM & BSLMM
    bstops = bsres[bsres$site==lematch,]
    
    allmatchs = do.call(rbind, lapply(1:nrow(lmtops), function(k) {
      #k=1
      if(is.na(match(lmtops$rs[k], bstops$rs))==F) {
        data.frame(lmtops[k,], bstops[match(lmtops$rs[k], bstops$rs),])
      }
    }))
    data.frame(allmatchs, rep(lematch, nrow(allmatchs)))
  }}))

#now summarize the results
names(overlapping_SNPs)[ncol(overlapping_SNPs)] = 'letrait'
redtraits = unique(overlapping_SNPs$letrait)

overlap_summary = do.call(rbind, lapply(1:length(redtraits), function(i) {
  #i=1
  data.frame(redtraits[i],((nrow(overlapping_SNPs[overlapping_SNPs$letrait==redtraits[i],])/snp_num)*100))
}))
names(overlap_summary)[1] = 'site'
names(overlap_summary)[2] = 'percent_overlap' #percent overlap is the number of SNPs a site has in common with the others

overlap_summary

library(ggplot2)
ggplot(data=overlap_summary, aes(x=site, y=percent_overlap, fill=site)) + 
  geom_bar(stat = "identity") +
  xlab('') +
  ylab('percent overlap') +
  ggtitle('Percent of top 100 SNPs for a site found in other sites') +
  theme(axis.text=element_text(size=20,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  scale_y_continuous(breaks=number_ticks(10)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #rotate labels

################################################################################################
################################################################################################
#make Manhattan plots and genome-wise beta distribution plots
#i=1
library(ggplot2)
pdf(paste('/Users/mrktaylor531/Desktop/Projects/FIBR/gemma_results/output_graphs_6_29_15/',trait,'_Manhattan_and_BSLMMparams.pdf', sep=''))

for (i in 1:length(lmm)) {
  
  #tryCatch({try catch to the error files aren't working!!!!!
  #get emma file from phenos dataframe 
  #i=1
  emmares = read.table(file=paste(lmm[i]), header=T)
  poo = data.frame(strsplit(as.character(lmm_names[i]), trait))
  
  manhattanplot = ggplot(data=emmares, aes(x = ps, y = -log10(p_wald), color=as.factor(chr))) + 
    geom_point(size=0.75) +
    facet_grid(chr ~ .) +
    xlab('position') +
    ylab('-log10(Wald test p-score)') +
    ggtitle(paste(poo[1,1], 'LMM', sep=' ')) +
    theme(text = element_text(size=17)) +
    theme(axis.text=element_text(size=14,color='black'),
          axis.title=element_text(size=14,color='black')) +
    theme(plot.title = element_text(size = 24)) +
    theme(panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    scale_x_continuous(breaks=number_ticks(10)) +
    #scale_y_continuous(breaks=number_ticks(10)) +
    theme(axis.line = element_line(colour = 'black')) +
    theme(legend.position = 'none')
  
  xint = quantile(-log10(emmares$p_wald), 0.99)
  
  pwald_dist = ggplot(data=emmares, aes(x=-log10(p_wald))) + 
    geom_histogram(data=emmares, aes(x=-log10(p_wald)), alpha=0.8, binwidth = 0.1) +
    geom_vline(aes(xintercept = xint), colour="red", size=1.2) +
    xlab(paste(poo[1,1], 'Wald p-value',sep=' ')) +
    ylab('') +
    ggtitle(paste(poo[1,1], 'LMM', 'Wald')) +
    theme(axis.text=element_text(size=14,color='black')) +
    theme(plot.title = element_text(size = 28)) +
    theme(panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(legend.text=element_text(size=24)) +
    theme(text = element_text(size=24)) +
    #scale_x_continuous(breaks=number_ticks(10)) +
    #scale_y_continuous(breaks=number_ticks(10)) +
    theme(axis.line = element_line(colour = 'black')) +
    theme(axis.line = element_line(colour = 'black'))
  
  ################################################################################################
  #BSLMM beta distributions
  
  bslmmres = read.table(bslmm[i], header=T)
  poo = data.frame(strsplit(as.character(bslmm_names[i]), trait))
  
  manplot_bslmm = ggplot(data=bslmmres, aes(x = ps, y = abs(beta), color=as.factor(chr))) + 
    geom_point(size=0.75) +
    facet_grid(chr ~ .) +
    xlab('position') +
    ylab('absolute value (beta)') +
    ggtitle(paste(poo[1,1], 'BSLMM beta', sep=' ')) +
    theme(text = element_text(size=17)) +
    theme(axis.text=element_text(size=14,color='black'),
          axis.title=element_text(size=14,color='black')) +
    theme(plot.title = element_text(size = 24)) +
    theme(panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    scale_x_continuous(breaks=number_ticks(10)) +
    #scale_y_continuous(breaks=number_ticks(10)) +
    theme(axis.line = element_line(colour = 'black')) +
    theme(legend.position = 'none')
  
  xintb = quantile(abs(bslmmres$beta), 0.99)
  
  beta_dist = ggplot(data=bslmmres, aes(x=abs(beta))) + 
    geom_histogram(data=bslmmres, aes(x=abs(beta)), alpha=0.8) +
    geom_vline(aes(xintercept = xintb), colour="red", size=1.2) +
    xlab(paste(poo[1,1], 'beta',sep=' ')) +
    ggtitle(paste(poo[1,1], 'BSLMM', 'beta')) +
    theme(axis.text=element_text(size=14,color='black')) +
    theme(plot.title = element_text(size = 28)) +
    theme(panel.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(legend.text=element_text(size=24)) +
    theme(text = element_text(size=24)) +
    #scale_x_continuous(breaks=number_ticks(10)) +
    #scale_y_continuous(breaks=number_ticks(10)) +
    theme(axis.line = element_line(colour = 'black')) +
    theme(axis.line = element_line(colour = 'black'))
  
  print(manhattanplot)
  print(manplot_bslmm)
  multiplot(pwald_dist, beta_dist, cols=1)
  
}

dev.off() 

#############################################################################################
#############################################################################################
#plot overlapping points betweeen top 100 Manhattan and top 100 beta distribution points

#plot top 100 overlapping SNPs between BSLMM and LMM for all sites
lmmres$method = rep('LMM', nrow(lmmres))
bsres$method = rep('BSLMM', nrow(bsres))
avg_ld = rbind(0,10000)
col2 = rbind(10,10)
avg_ld = data.frame(cbind(avg_ld, col2))
names(avg_ld)[1] = 'dist'
names(avg_ld)[2] = 'value'
head(avg_ld)

ggplot() + 
  geom_point(data=bsres, aes(x = ps, y = abs(beta), color='darkblue', aes=0.5), color='darkblue', aes=0.5) + #LMM results
  geom_point(data=lmmres, aes(x = ps, y = -log10(p_wald), color='darkgreen', aes=0.5), color='darkgreen', aes=0.5) + #BSLMM results
  geom_line(data=avg_ld, aes(x=dist, y=value, color='red'), color='red', size=3) +
  facet_grid(chr ~ .) +
  xlab('position') +
  ylab('absolute value (beta) or -log10(p-Wald)') +
  ggtitle('Compare top SNPs between methods') +
  theme(text = element_text(size=17)) +
  theme(axis.text=element_text(size=14,color='black'),
        axis.title=element_text(size=14,color='black')) +
  theme(plot.title = element_text(size = 24)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  scale_x_continuous(breaks=number_ticks(10)) +
  #scale_y_continuous(breaks=number_ticks(10)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(legend.position = 'none')

#########################################################################################################
#########################################################################################################
#find nearest genomic subclass from gemma results
tair9 = read.table(file='/Users/mrktaylor531/Desktop/Projects/FIBR/tair9_transcript_location_CORRECTED.txt', header=T) #read
kinds = unique(tair9$chromosome)

#define function to find nearest genomic subclass and the distance to the CENTER of this subclass
find_closest <- function(lmmres){
  closest = do.call(rbind, lapply(1:nrow(lmmres), function(i) {
    #i=1800
    chrom = tair9[tair9$chr==lmmres$chr[i],]
    matchs = chrom[abs(chrom$center_pos-lmmres$ps[i])==min(abs(chrom$center_pos-lmmres$ps[i])),]
    distance = lmmres$ps[i] - matchs[1,8]
    data.frame(lmmres[i,], distance ,matchs[1,]) #output SNP results along with closest gene
  }))
  return(closest)
}

bsclosest = find_closest(bsres)
lmmclosest = find_closest(lmmres)

#write out the closest genomic elements
write.table(lmmclosest, file=paste('/Users/mrktaylor531/Desktop/Projects/FIBR/gemma_results/post_gwas_analysis_8_3_15/top_',snp_num,'LMM_.txt', sep=''), row.names=F, sep='\t')
write.table(bsclosest, file=paste('/Users/mrktaylor531/Desktop/Projects/FIBR/gemma_results/post_gwas_analysis_8_3_15/top_',snp_num,'BSLMM_.txt', sep=''), row.names=F, sep='\t')

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
##search for genome kind enrichment (exon, mRNA, etc)

#define function to search for nearest genomic subclass to a SNP and then ask if this SNP is inside the genomic subclass (T/F)
enrichment = function(lmmres) {genome_kinds_summary = do.call(rbind, lapply(1:length(kinds), function(j) {
  #j=1
  tryCatch({ #suppress error because the top100 SNPs might not represent all possible genomic classes
    
    exon = tair9[tair9$chromosome==kinds[j],] #exon is actually any subclass of the genome
    
    closest_genes = do.call(rbind, lapply(1:nrow(lmmres), function(i) {
      #i=1640
      chrom = exon[exon$chr==lmmres$chr[i],]
      matchs = chrom[abs(chrom$center_pos-lmmres$ps[i])==min(abs(chrom$center_pos-lmmres$ps[i])),]
      matchs = matchs[1,]
      lmmres$ps[i] < matchs$end_pos & lmmres$ps[i] > matchs$start_pos #True means that the SNP is within the kind of genome
    }))
    
    sums = data.frame(table(closest_genes))
    data.frame(kinds[j], (sums[2,2]/nrow(closest_genes)))
  }, error=function(e){})
}))
return(genome_kinds_summary)
}


lmkinds_summary = enrichment(lmmres)
names(lmkinds_summary)[1] = 'class'
names(lmkinds_summary)[2] = 'percent'

bskinds_summary = enrichment(bsres)
names(bskinds_summary)[1] = 'class'
names(bskinds_summary)[2] = 'percent'

#define function to replace NA's with zeros in second column of the kinds summaries
replace_zero = function(df) {
  todo = do.call(rbind, lapply(1:nrow(df), function(i) {
    if(is.na(df$percent[i])==T) {poo = data.frame(df[i,1], 0)} else {poo = data.frame(df[i,])}
    row.names(poo) = i
    names(poo)[1] = 'class'
    names(poo)[2] = 'percent'
    poo
  }))
  return(todo)
}

lmkinds_sum = replace_zero(lmkinds_summary)
bskinds_sum = replace_zero(bskinds_summary)
lmkinds_sum$method = rep('lmm', rep(nrow(lmkinds_sum)))
bskinds_sum$method = rep('bslmm', rep(nrow(bskinds_sum)))

kinds_sum = rbind(lmkinds_sum, bskinds_sum)

########################################################################
#calculate the percentage of SNPs in each class at per_val percentile value for the empirical distribution of SNP classes
#see if there is an excess of a genomic class
perms = read.table(file='/Users/mrktaylor531/Desktop/Projects/FIBR/gemma_results/post_gwas_analysis_6_29_15/snp_chip_CLOSEST_PERMUTATION_ALLPERMS.txt', header=T)

permutationres = perms
leskinds = unique(permutationres$kind)
per_val = 0.90
permutation_sum = do.call(rbind, lapply(1:length(leskinds), function(i) {
  #i=1
  toavg = permutationres[permutationres$kind==leskinds[i],]
  data.frame(leskinds[i],quantile(na.omit(toavg$percent), per_val))
}))

names(permutation_sum)[1] = 'class'
names(permutation_sum)[2] = 'percent'
permutation_sum$method = rep(paste(per_val, '_percentile', sep=''), nrow(permutation_sum))

#######################################################################
#plot percent of top SNPs in each genomic subclass
ggplot(data=kinds_sum, aes(x=class, y=percent, fill=method)) + 
  geom_bar(stat = "identity", position="dodge") +
  xlab('') +
  ggtitle('Percent of top 100 SNPs in each genomic class') +
  theme(axis.text=element_text(size=20,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #rotate labels
  scale_fill_hue(l=45) #rotate labels

#####################################################################
#plot perent of top SNPs in each clas vs permutation subclass average
perm_sum = read.table(file='/Users/mrktaylor531/Desktop/Projects/FIBR/gemma_results/post_gwas_analysis_6_29_15/snp_chip_SUMMARY_PERMUTATION.txt', header=T)
perm_sum$method = rep('SNP_permutation', nrow(perm_sum))

kinds_sum = rbind(kinds_sum, perm_sum)

ggplot(data=kinds_sum, aes(x=class, y=percent, fill=method)) + 
  geom_bar(stat = "identity", position="dodge") +
  xlab('') +
  ggtitle('Percent of top 100 SNPs in each genomic class') +
  theme(axis.text=element_text(size=20,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #rotate labels
  scale_fill_hue(l=45) #rotate labels

#####################################################################
#plot with the percentile value
kinds_sum = rbind(kinds_sum, permutation_sum)

ggplot(data=kinds_sum, aes(x=class, y=percent, fill=method)) + 
  geom_bar(stat = "identity", position="dodge") +
  xlab('') +
  ggtitle('Percent of top 100 SNPs in each genomic class') +
  theme(axis.text=element_text(size=20,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #rotate labels
  scale_fill_hue(l=45) #rotate labels

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#NON-SYNONYMOUS VS SYNONYMOUS ENRICHMENT
#only if it is actually within proteins

source("http://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Athaliana.TAIR.TAIR9") #install package
library(BSgenome.Athaliana.TAIR.TAIR9)
library(BSgenome)
#biocLite("BSgenome") #install package

lmmres = read.table(file=paste(outputdir,'top_SNPs_LMM.txt',sep=''), header=T)
bres = read.table(file=paste(outputdir,'top_SNPs_BSLMM.txt',sep=''), header=T)

cds = read.table(file='/Users/mrktaylor531/Desktop/Projects/FIBR/tair9_cds.txt', header=T)
mRNA = read.table(file='/Users/mrktaylor531/Desktop/Projects/FIBR/tair9_mRNA.txt', header=T)
protein = read.table(file='/Users/mrktaylor531/Desktop/Projects/FIBR/tair9_protein.txt', header=T)

nonsyn = do.call(rbind, lapply(1:nrow(lmmres), function(i) {
  #i=1640
  chrom = protein[protein$chr==lmmres$chr[i],]
  
  matchs = chrom[abs(chrom$center_pos-lmmres$ps[i])==min(abs(chrom$center_pos-lmmres$ps[i])),]
  matchs = matchs[1,]
  if((lmmres$ps[i] < matchs$end_pos & lmmres$ps[i] > matchs$start_pos)==T){
    seq = getSeq(Athaliana, paste('Chr',lmmres$chr[i], sep=''), matchs$start_pos, matchs$end_pos)
    snp_pos = lmmres$ps[i] - matchs$start_pos #this is where the SNP lies in the seq
    
    #get the protein sequence for both alleles
    seq[snp_pos] = as.character(lmmres$allele1[i])
    aaseq1 = as.character(translate(seq))
    seq[snp_pos] = as.character(lmmres$allele0[i])
    aaseq2 = as.character(translate(seq))
    data.frame(lmmres[i,],aaseq1==aaseq2) #the last column if TRUE=synonymous, FALSE=non-synonymous
  }
}))

head(nonsyn)
names(nonsyn)[ncol(nonsyn)]='synonymous'  #the last column if TRUE=synonymous, FALSE=non-synonymous
sums = data.frame(table(nonsyn$synonymous))

sums[1,2]/nrow(nonsyn) #this is the percent of NS SNPs in the top 100 SNPs
#the false designation changes based on other things!!

#permutation of SYN vs non-SYN across SNP chip
#transfer file from FARM
# scp farm:/home/maktaylo/fibr_gwas/post_gwas_analysis/SYN_VS_NONSYN_PERMUTATION.txt /Users/mrktaylor531/Desktop/Projects/FIBR/SYN_VS_NONSYN_PERMUTATION.txt
synperm = read.table(file='/Users/mrktaylor531/Desktop/Projects/FIBR/SYN_VS_NONSYN_PERMUTATION.txt', header=T)

head(synperm)
xint = quantile(synperm$percent_nonsyn,0.9)

ggplot(data=synperm, aes(x=percent_nonsyn)) + 
  geom_histogram(alpha=0.9, colour = "black", fill='grey') +
  geom_vline(xintercept = xint, color='black', size=2) +
  geom_vline(xintercept = 0.7425, color='red', size=2) +
  geom_vline(xintercept = 0.25, color='red', size=2) +
  xlab('proportion') +
  ylab('') +
  ggtitle('Proportion of non-syn SNPs in 100 randomly sampled SNPs') +
  theme(axis.text=element_text(size=18,color='black')) +
  theme(plot.title = element_text(size = 28)) +
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(legend.text=element_text(size=24)) +
  theme(text = element_text(size=24)) +
  #scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(6)) +
  theme(axis.line = element_line(colour = 'black')) +
  theme(axis.line = element_line(colour = 'black'))



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#examine the allele frequency spectrum ON FARM

#across the entire SNP chip then in my own subsamples
#start FARM session : srun -p serial --pty R

chr1 = read.table(file='/home/maktaylo/fibr_gwas/genotypes/formatted_genotype/mdp_genotype_chr1.hmp.txt', header=T)
chr1[1:20,1:20]

allelefreq = do.call(rbind, lapply(1:nrow(chr1), function(i) {
  #i=1
  poo = data.frame(table(t(chr1[i,c(12:ncol(chr1))])))
  if(poo$Freq[1] >= poo$Freq[2]) {poo$Freq[2] / sum(poo$Freq)} else if(poo$Freq[1] <= poo$Freq[2]) {poo$Freq[1] / sum(poo$Freq)}
}))

#make column of allele position within the chromosomel
rsnum = do.call(rbind, lapply(1:nrow(chr1), function(i) {
  #i=2
  poo = data.frame(strsplit(as.character(chr1$rs[i]), '_'))
  data.frame(poo[2,])
}))
chr1 = cbind(rsnum, chr1)








