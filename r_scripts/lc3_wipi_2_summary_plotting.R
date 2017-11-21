# lc3_wipi_2_summary_plotting.R


# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)
require(tidyverse)
# load data
input_df <- read_csv('~/Dropbox/code/csth-imaging/output_files/lc3_wipi_2_summary_analysis_output.csv')

# load in the corrected "z overlap" flags
z_flags <- read_csv('~/Dropbox/code/csth-imaging/output_files/lc3_wipi_2_z_flags.csv')
split_fnames <- strsplit(z_flags$filename, '/')
split_fnames <- data.frame(matrix(unlist(split_fnames), nrow = length(split_fnames), byrow = TRUE))
z_flags$filename <- split_fnames$X7
z_flags$im_number <- NULL
names(z_flags)[names(z_flags) == 'image'] <- 'im_number'
names(input_df)[names(input_df) == 'flagged_z'] <- 'old_z_flags'
input_df <- left_join(input_df, z_flags)
# remove cells that border on the edge of the image and defective parent cell assignment
for_plt <- subset(input_df, parent_cell != 65535 & flagged_z != 1)

# add a 0 overlap_counts label to all cells with 0 foci
for_plt$overlap_count[for_plt$count == 0] <- 0
for_plt$treatment[for_plt$treatment %in% c('1hrTort', '1hrTor')] <- 1
for_plt$treatment[for_plt$treatment == '0hrTor'] <- 0
for_plt$treatment <- factor(for_plt$treatment)
for_plt$cell_line[for_plt$cell_line == 'dVPS'] <- 'dVPS37A'
for_plt$cell_line[for_plt$cell_line == 'dTMEM'] <- 'dTMEM41B'
for_plt$cell_line[for_plt$cell_line == 'dTMEM+dFIP200'] <- 'dTMEM41B + dFIP200'
for_plt$cell_line <- factor(for_plt$cell_line,
                            levels=c('WT','dVPS37A','dTMEM41B', 'dTMEM41B + dFIP200'))


# Add a separate column indicating which cells are boxplot outliers for jittering
for_plt <- for_plt %>% group_by(channel, treatment, cell_line) %>%
  dplyr::mutate(count_outlier = count > quantile(count, 0.75) + IQR(count)*1.5 |
                  count < quantile(count, 0.25) - IQR(count)*1.5) %>% ungroup
for_plt <- for_plt %>% group_by(channel, treatment, cell_line) %>%
  dplyr::mutate(overlap_outlier = overlap_count > quantile(overlap_count, 0.75) + IQR(overlap_count)*1.5 |
                  overlap_count < quantile(overlap_count, 0.25) - IQR(overlap_count)*1.5) %>% ungroup

# count the # of cells in each sample
n_cells <- for_plt %>% tbl_df %>% group_by(channel, cell_line, treatment) %>% dplyr::summarise(n_cells = n())
View(n_cells)

# plot treatment vs # of foci in each channel for each cell line
ggplot(subset(for_plt, cell_line != 'dVPS37A'), aes(x=cell_line, y=count)) +
  facet_grid(channel ~ .) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(data=subset(for_plt, count_outlier==TRUE & cell_line != 'dVPS37A'), width=0.2, height=0, size=1) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank(),
        plot.title = element_text(hjust=0.5)) +
  labs(x='Torin treatment time (hr)',
       y='Number of foci per cell',
       title='Cell line')
ggsave('~/Dropbox/code/csth-imaging/r_scripts/lc3_wipi_2_plots/lc3_wipi_2_foci_cts_boxplot_2.pdf',
       device=cairo_pdf(width=2, height=4), useDingbats=F)

ggplot(subset(for_plt, cell_line != 'dVPS37A'), aes(x=cell_line, y=overlap_count)) +
  facet_grid(channel ~ .) + 
  geom_boxplot(outlier.shape=NA, fill = 'white') + geom_jitter(data=subset(for_plt, overlap_outlier==TRUE & cell_line != 'dVPS37A'), width=0.2, height=0, size=1) +
  scale_y_continuous(breaks=c(0,50,100)) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank(),
        plot.title = element_text(hjust=0.5)) +
  labs(x='Torin treatment time (hr)',
       y='Number of foci that overlap\nwith the opposite channel',
       title='Cell line')
ggsave('~/Dropbox/code/csth-imaging/r_scripts/lc3_wipi_2_plots/lc3_wipi_2_overlap_cts_boxplot_2.pdf',
  device=cairo_pdf(width=2, height=4), useDingbats=F)
