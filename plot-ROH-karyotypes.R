library(karyoploteR)
library(dplyr)
library(stringr)
# A script to plot ROH from the bed file outputted from the vcf parser 'roh-bin-genome.py'
# It accepts a multisample (stacked) txt file and output all samples on one karyogram.
# It works for up to 3 samples
# It plots samples from the bottom up, from the ideogram,
# sample order is obtained from a seperate text file generated from the 'roh-bin-genome.py' script

#read data, change columns and establish sample details and count
args <- commandArgs(TRUE)
finaldata <- data.frame(read.delim2(args[1], header=FALSE, row.names=NULL)[-1,])
colnames(finaldata) <-c('chrom','start','end','sample','run','zygosity','count','bases_per_snv','chunk')
sample_names <- as.character(read.delim2(args[2], header=FALSE, row.names=NULL)[-4,])
sample_names_count <- length(unique(finaldata$sample))

sample_plotting_dict <- list('1_start_1' = 0,
                             '1_end_1' = 1,
                             '2_start_1' = 0,
                             '2_start_2' = 0.5,
                             '2_end_1' = 0.48,
                             '2_end_2' = 0.97,
                             '3_start_1' = 0,
                             '3_start_2' = 0.33,
                             '3_start_3' = 0.72,
                             '3_end_1' = 0.31,
                             '3_end_2' = 0.70,
                             '3_end_3' = 0.98)

pdf(file = args[3],
    width = 10, 
    height = 8)
# Create plot
kp <- plotKaryotype(genome="hg19")
#kpDataBackground(kp,r0=0.5,r1=1,color = "#FFEEEE",data.panel = 2)
for (i in 1:length(sample_names)) {
  correct_sample_order <- rev(sample_names) # need to rev sample order because Kplot plots from bottom to top
  colours <- c('#CCFFAA','#FFAACC','#031CFC')
  # set the values according to number of samples
  start <- as.numeric(sample_plotting_dict[paste0(sample_names_count,'_start_',i)])
  end <- as.numeric(sample_plotting_dict[paste0(sample_names_count,'_end_',i)])
  
  filter(finaldata,str_detect(sample,sample_names[i])) %>%
    filter(.,str_detect(zygosity, "hom")) %>% 
    select(.,chrom,start,end) %>%
    filter(.,(end-start > 300000)) %>% 
    toGRanges(.)  %>%
    kpPlotRegions(kp, .,r0=start,r1=end,num.layers = 1,col=colours[i])
  # Make a one-off additon of the ROI  
  kpRect(kp, chr="chr4", x0=146864389, x1=152368608, y0=0, y1=0.98,border="red")
  kpAddLabels(kp,labels=correct_sample_order[i],r0=(i-0.49),r1=(i-0.03),cex=.35)
}

dev.off()
