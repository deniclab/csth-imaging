#This program reads in FlowJo count files (.csv) and produces violin plots of the data
#Red:Green ratios were calculated in FlowJo prior to data export

#set working directory
#point at folder containing exported data (.csv files) from FlowJo
setwd("/Users/Chris/Desktop/experiments/")

#call ggplot2
require(ggplot2)

#read in all files in a folder
temp = list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

#plot violin of calculated red:green ratios
ggplot() +
  
#start with negative control 1
  geom_violin(aes("negative control 1", export_sample1.csv),draw_quantiles = .5, fill="dark gray") +
  
#proceed with each sample of interest
  geom_violin(aes("sample name 1", export_sample1.csv), draw_quantiles = .5, fill="dark gray") + 
  geom_violin(aes("sample name 2", export_sample1.csv), draw_quantiles = .5, fill = "dark gray") +
  geom_violin(aes("sample name 3", export_sample1.csv), draw_quantiles = .5,  fill="dark gray") +
  geom_violin(aes("sample name 4", export_sample1.csv), draw_quantiles = .5, fill = "dark gray") +
  
#...
#continue with pattern until all samples are imported
  
#end with negative control 2
  geom_violin(aes("negative control 2", export_sample1.csv), draw_quantiles = .5,  fill="dark gray") +
  
#scaling of y axis (chose one)
  scale_y_log10() + 
    #or
  #scale_y_log10(limits=c(.1,50)) + #use this for custom scaling
  
  theme_bw() + ggtitle("Enter Plot Title here") + 
  theme(axis.text.x=element_text(angle = (-90), vjust = 0.5))+
  xlab("Samples") + ylab("Red:Green ratio") 

