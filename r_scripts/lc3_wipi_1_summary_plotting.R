# lc3_wipi_1_summary_plotting.R


# import dependencies
require(readr)
require(ggplot2)
require(plyr)
require(dplyr)
require(gridExtra)
require(Cairo)
require(reshape2)

# load data
input_df <- read_csv('~/Dropbox/code/csth-imaging/output_files/lc3_wipi_1_summary_analysis_output.csv')

# remove cells that border on the edge of the image and defective parent cell assignment
for_plt <- subset(input_df, parent_cell != 65535 & flagged_z != 1)

# add a 0 overlap_counts label to all cells with 0 foci
for_plt$overlap_count[for_plt$count == 0] <- 0


# plot treatment vs # of foci in each channel for each cell line
ggplot(for_plt, aes(x=treatment, y=count)) +
  facet_grid(channel ~ cell_line) + geom_boxplot() +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank()) +
  labs(y='Number of foci per cell')
#ggsave('~/Dropbox/code/csth-imaging/r_scripts/lc3_wipi_1_plots/lc3_wipi_1_foci_cts_violin.pdf',
#       device=cairo_pdf(width=6, height=4))

ggplot(for_plt, aes(x=treatment, y=overlap_count)) +
  facet_grid(channel ~ cell_line) + geom_boxplot() +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        strip.background = element_blank()) +
  labs(y='Number of foci per cell that overlap\nwith foci in the other channel')
#ggsave('~/Dropbox/code/csth-imaging/r_scripts/lc3_wipi_1_plots/lc3_wipi_1_overlap_cts_violin.pdf',
#  device=cairo_pdf(width=6, height=4))

# add a column for non-overlapping foci
for_plt$non_overlap_count <- for_plt$count - for_plt$overlap_count

# make a new df with one column containing 488_only, 561_only, and 561_overlap as the different data
# points for each object. make a row for interaction(im_number, id), but remove the label. sort by total
# foci.

# re-organize data so that the overlap count and non-overlap-count are in the same column (value), with
# a different variable that indicates which of the two each row corresponds to (variable).
melted_df <- melt(for_plt, measure.vars = c('overlap_count','non_overlap_count'))
# make a new ID variable column, which will be used to determine color for plotting, which contains
# both the channel and whether the value is indicating overlap or not.
melted_df$id_var <- interaction(melted_df$channel, melted_df$variable)
# sort each by # of foci
stacked_df <- subset(melted_df, id_var %in% c('488.non_overlap_count', '561.overlap_count',
                                              '561.non_overlap_count'))
get_tot_df <- dcast(stacked_df, im_number + parent_cell + czi_id + cell_line + treatment ~ id_var,
                    value.var = 'value')
get_tot_df$total <- get_tot_df$`488.non_overlap_count` + get_tot_df$`561.overlap_count` + get_tot_df$`561.non_overlap_count`
get_tot_df <- get_tot_df[order(get_tot_df$total), ]
get_tot_df$rank <- 1:length(get_tot_df$total)
get_tot_df <- melt(get_tot_df, measure.vars=c('488.non_overlap_count', '561.overlap_count',
                                              '561.non_overlap_count'))
names(get_tot_df)[names(get_tot_df) == 'variable'] <- 'id_var'
names(get_tot_df)[names(get_tot_df) == 'value'] <- 'foci_ct'
stacked_df <- merge(stacked_df, get_tot_df, by=intersect(names(stacked_df), names(get_tot_df)))
stacked_df$id_var <- factor(stacked_df$id_var,
                            levels=c('488.non_overlap_count', '561.overlap_count',
                                     '561.non_overlap_count'))
ggplot(stacked_df, aes(x = factor(rank),
                       y = foci_ct, color = id_var, fill=id_var)) + geom_bar(stat='identity') +
  facet_grid(. ~ interaction(treatment, cell_line), scales='free_x', space='free_x') +
  scale_color_manual(values=c('green','yellow','red')) +
  scale_fill_manual(values=c('green','yellow','red')) +
  theme(panel.background = element_rect(fill = 'black', color = 'white'),
        plot.background = element_rect(fill = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'white'),
        axis.ticks = element_line(color = 'white'),
        strip.text = element_text(color='white', angle=45),
        strip.background = element_blank(),
        axis.title = element_text(color = 'white'),
        axis.text.x = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(color='white')) +
  labs(x='Cell Line', y='Number of foci')
ggsave('~/Dropbox/code/csth-imaging/r_scripts/lc3_wipi_1_plots/lc3_wipi_ordered_bars.pdf',
       device = cairo_pdf(width=10, height=4))