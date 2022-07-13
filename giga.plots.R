
######  Generating plots ####

if(!exists('load')) {load <- T}

if (load) { load('all.RData') }
library(ggsci)
palette <- c(pal_nejm()(8), pal_npg()(10))
try({ dir.create(fig.dir) })

## Figure 2. Ratio of sub-genomic / genomic transcripts
#### ####
plot.data <- all.data %>% group_by(sample, category) %>% summarise(count=n()) %>% spread(category, count, fill=0)
#colnames(plot.data)[-1] <- c('non-genomic', 'genomic')
plot.data$ratio <- plot.data$'sub-genomic' / plot.data$genomic
plot.data <- merge(plot.data, metafilt[,c("sample_name", "h", "hpi", "Time", 'rep')], by.x='sample', by.y='sample_name')

plot.sum  <- plot.data[plot.data$hpi != 'dRNA', ] %>% 
  group_by(h, hpi, Time) %>% summarise(mean=mean(ratio), sd=sd(ratio), ymin=mean-sd, ymax=mean+sd )

gg.clust <- ggplot(plot.data[plot.data$hpi != 'dRNA', ]) + 
  geom_point(aes(x=Time, y=ratio)) + 
  geom_smooth(aes(x=Time, y=ratio)) +
  #geom_pointrange(aes(x=Time, y=mean, ymin=ymin, ymax=ymax, colour='mean'), size=0.25, data = plot.sum) +
  scale_fill_manual(values = palette) +
  scale_x_continuous(name = 'Time (hours past infection)') +
  scale_y_continuous(name = 'ratio (sgRNA/gRNA)') +
  theme_bw() + theme(legend.position = 'none') 
gg.clust

ggsave('giga.Fig2.jpg', width = 20, height = 14)
#### ####
##

## SuppFig S1 Violinplot of sub-genomic and genomic RNA lengths
#### ####
tr.stats <- tr.sp[,1:14] #merge(tr.sp[,1:7], tr.uni[,1:4], by='TR_ID')
colnames(tr.stats)[11] <- "Mapped Length (nt)"
tr.stats$method <- 'cDNA'
tr.stats$method[tr.stats$sample == 'dRNA' ] <- 'dRNA'

tr.stats$category <- NA
tr.stats$category[ tr.stats$is.genomic == T & tr.stats$is.subgenomic == F ] <- 'genomic'
tr.stats$category[ tr.stats$is.genomic == F & tr.stats$is.subgenomic == T ] <- 'sub-genomic'
tr.stats$category[ tr.stats$is.genomic == T & tr.stats$is.subgenomic == T ] <- 'sub-genomic'
tr.stats$category[is.na(tr.stats$category)] <- 'unclassified'

tr.stats <- merge(tr.stats, metafilt[c('sample_name', 'h')], by.y='sample_name', by.x='sample')

ggv.cDNA <- ggviolin(tr.stats[tr.stats$method == 'cDNA', ], size = 0.25,
                     'h', "Mapped Length (nt)", color = 'h', #palette = 'npg', 
                     add=c('boxplot', 'median', 'mean_sd'), fill='lightgrey') + 
  scale_color_viridis_d() +
  coord_cartesian(ylim=c(0, 32000)) +
  theme_bw() + 
  theme(legend.position='none',
        #axis.text.x = element_blank(),
        axis.title.x = element_blank()
  ) + 
  facet_grid(rows = vars(category), 
             cols = vars(method), scales = 'free_x')

ggv.dRNA <- ggviolin(tr.stats[tr.stats$method != 'cDNA', ], size = 0.25,
                     'method', "Mapped Length (nt)", color = 'method', palette = 'npg', 
                     add=c('boxplot', 'median', 'mean_sd'), fill='lightgrey') + 
  #scale_color_viridis_d() +
  coord_cartesian(ylim=c(0, 32000)) +
  theme_bw() + 
  theme(legend.position='none',
        #axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
  ) + 
  facet_grid(rows = vars(category), 
             cols = vars(method), scales = 'free_x')
ggv <- cowplot::plot_grid(ggv.cDNA, ggv.dRNA, rel_widths = c(6,1), align = 'vh', axis = 'tblr')

ggsave('giga.SuppFig_S1.jpg', ggv, width = 24, height = 12)
#### ####
##

#### ####
##
