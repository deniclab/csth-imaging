# lc3_p62_summary_plotting.R


# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)

# load data
input_df <- read_csv('~/Dropbox/code/csth-imaging/output_files/lc3_p62_summary_analysis_output.csv')

# remove cells that border on the edge of the image and defective parent cell assignment
for_plt <- subset(input_df, parent_cell != 65535 & flagged_z != 1)

# add a 0 overlap_counts label to all cells with 0 foci
for_plt$overlap_count[for_plt$count == 0] <- 0


# plot treatment vs # of foci in each channel for each cell line
ggplot(for_plt, aes(x=treatment, y=count)) +
  facet_grid(channel ~ cell_line) + geom_violin() +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank()) +
  labs(y='Number of foci per cell')
#ggsave('~/Dropbox/code/csth-imaging/r_scripts/lc3_p62_plots/lc3_p62_foci_cts_violin.pdf',
#       device=cairo_pdf(width=6, height=4))

ggplot(for_plt, aes(x=treatment, y=overlap_count)) +
  facet_grid(channel ~ cell_line) + geom_violin() +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank()) +
  labs(y='Number of foci per cell that overlap\nwith foci in the other channel')
#ggsave('~/Dropbox/code/csth-imaging/r_scripts/lc3_p62_plots/lc3_p62_overlap_cts_violin.pdf',
#  device=cairo_pdf(width=6, height=4))
