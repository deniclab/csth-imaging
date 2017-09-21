# lc3_p62_timelapse_plotting.R

# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# load data
input_df <- read_csv('~/Dropbox/code/csth-imaging/output_files/lc3_p62_timelapse_raw_analysis_output.csv')

# list the different unique treatment identifiers
unique(input_df$treatment)

# checking to see how many foci failed parent cell assignment
ddply(.data=input_df, .variables='parent_cell', summarize,
      n=length(parent_cell))

# 1798 were mis-assigned; not too bad.

# typo in some of the filenames led to two separate torin treatment names:
# 1hrTort and 3hrTort. Adding a new column, tmt_time, to make plotting easier.

input_df$tmt_time <- NA
input_df$tmt_time[input_df$treatment == 'NoTreat'] <- 0
input_df$tmt_time[input_df$treatment == '1hrTor'] <- 1
input_df$tmt_time[input_df$treatment == '1hrTort'] <- 1
input_df$tmt_time[input_df$treatment == '3hrTor'] <- 3
input_df$tmt_time[input_df$treatment == '3hrTort'] <- 3



# forgot to add z- and oof-flags to this part of the dataset.
# i'm going to pull it out of the summary file.

summ_file <- read_csv('~/Dropbox/code/csth-imaging/output_files/lc3_p62_timelapse_summary_analysis_output.csv')
z_flagged <- subset(summ_file, flagged_z == 1)
input_df$flagged_z <- 0
for (r in 1:length(z_flagged$flagged_z)){  # iterate over rows
  input_df$flagged_z[input_df$channel == z_flagged$channel[r] &
                       input_df$filename == z_flagged$filename[r] &
                       input_df$im_number == z_flagged$im_number[r] &
                       input_df$parent_cell == z_flagged$parent_cell[r]] <- 1
}

# get rid of rows where parent cell assignment failed and where cell is on image edge
for_plt <- subset(input_df, parent_cell != 65535 & flagged_z == 0)
# summarize data.
# this is already done in the summary csv, but I'm re-doing it here so that
# I can work with one source data file for all of the plotting
# (I find it easier)

summ_df <- ddply(.data=for_plt,
                 .variables=c('cell_line','tmt_time','parent_cell','channel',
                              'im_number','czi_id'),
                 summarize, vol_mean = mean(volume, na.rm=T),
                 vol_sd = sd(volume, na.rm=T), 
                 int_mean = mean(intensity, na.rm=T),
                 int_sd = sd(intensity, na.rm=T),
                 n = length(intensity),  # number of foci per cell
                 n_overlap = sum(overlap),
                 scaling_factor = mean(scaling_factor)  # shouldn't change it
                 )

ggplot(summ_df, aes(x=factor(tmt_time), y=n)) +
  facet_grid(channel ~ cell_line) + geom_violin() +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank()) +
  labs(x='Torin treatment time (hr)', y='Number of foci per cell')

ggsave('~/Dropbox/code/csth-imaging/r_scripts/p62_timelapse_plots/p62_timecourse_foci_cts_violin.pdf',
       device=cairo_pdf(width=6, height=4))