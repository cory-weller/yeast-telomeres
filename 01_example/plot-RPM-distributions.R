#!/usr/bin/env R

library(data.table)
library(ggplot2)
library(foreach)
library(ggthemes)
library(showtext)

font_add_google('Atkinson Hyperlegible', 'Atkinson Hyperlegible')
showtext_auto()

# dat <- foreach(file=list.files(pattern='*.telo-counts.txt'), .combine='rbind') %do% {
#     sra_id <- strsplit(file, '_')[[1]][1]
#     tmp <- fread(file)
#     setnames(tmp, c('id', 'nReads', 'nTeloReads'))
#     tmp[, SRA := sra_id][]
# }
# dat[, RPM := 1e6 * (nTeloReads / nReads) ]


# dat[SRA=='SRR9330808', facetLabel := 'YJM981 x CBS2888']
# dat[SRA=='SRR9330831', facetLabel := 'CBS2888 x CLIB219']
# dat[SRA=='SRR9330832', facetLabel := 'CLIB219 x M22']
# dat[SRA=='SRR9330811', facetLabel := 'M22 x BY']
# dat[SRA=='SRR9330837', facetLabel := 'BY x RM']
# dat[SRA=='SRR9330810', facetLabel := 'RM x YPS163']
# dat[SRA=='SRR9330836', facetLabel := 'YPS163 x YJM145']
# dat[SRA=='SRR9330813', facetLabel := 'YJM145 x CLIB413']
# dat[SRA=='SRR9330815', facetLabel := 'CLIB413 x YJM978']
# dat[SRA=='SRR9330812', facetLabel := 'YJM978 x YJM454']
# dat[SRA=='SRR9330833', facetLabel := 'YJM454 x YPS1009']
# dat[SRA=='SRR9330814', facetLabel := 'YPS1009 x I14']
# dat[SRA=='SRR9330817', facetLabel := 'I14 x Y10']
# dat[SRA=='SRR9330816', facetLabel := 'Y10 x PW5']
# dat[SRA=='SRR9330830', facetLabel := 'PW5 x 273614']
# dat[SRA=='SRR9330809', facetLabel := '273614 x YJM981']

# dat[, facetLabel := factor(facetLabel, levels=c('YJM981 x CBS2888','CBS2888 x CLIB219',
#                             'CLIB219 x M22','M22 x BY','BY x RM','RM x YPS163','YPS163 x YJM145',
#                             'YJM145 x CLIB413','CLIB413 x YJM978','YJM978 x YJM454',
#                             'YJM454 x YPS1009','YPS1009 x I14','I14 x Y10',
#                             'Y10 x PW5','PW5 x 273614','273614 x YJM981'))]


# g <- ggplot(data=dat, aes(x=RPM)) + geom_histogram(bins=80) + theme_few(12) +
# labs(x='Telomeric Reads per 1M reads', y='N') +
# facet_wrap(~facetLabel, nrow=4) +
# theme(text=element_text(family="Atkinson Hyperlegible", size = 20))

# ggsave(g, file='telomeric-min8x.png', width=10, height=8, units='cm', dpi=300, scale=2)


SRA <- 'SRR9330808'

files <- list.files(pattern=paste0(SRA, '.*telo-counts.txt'))
teloCounts <- foreach(file=files, .combine='rbind') %do% {
    split_fn <- strsplit(file, split="_")[[1]]


    tmp <- fread(file, header=F)
    tmp[, SRA := split_fn[1]]
    tmp[, minRepeats := split_fn[2]]
    tmp[]
}

teloCounts[, c('sample','teloReads') := tstrsplit(V1, split=':')]
teloCounts[, V1 := NULL]

lineCounts <- fread(paste0(SRA, '_linecounts.txt'), header=F)
setnames(lineCounts, c('lines','sample'))
lineCounts <- lineCounts[sample != 'total']
lineCounts[, totalReads := lines/4]
lineCounts[, lines := NULL]
setkey(lineCounts, sample)
setkey(teloCounts, sample)

dat <- merge(lineCounts, teloCounts)
dat[, teloReads := as.numeric(teloReads)]
dat[, totalReads := as.numeric(totalReads)]
dat[, RPM := 1e6 * teloReads/totalReads]


dat[SRA=='SRR9330808', facetLabel := 'YJM981 x CBS2888']
dat[SRA=='SRR9330831', facetLabel := 'CBS2888 x CLIB219']
dat[SRA=='SRR9330832', facetLabel := 'CLIB219 x M22']
dat[SRA=='SRR9330811', facetLabel := 'M22 x BY']
dat[SRA=='SRR9330837', facetLabel := 'BY x RM']
dat[SRA=='SRR9330810', facetLabel := 'RM x YPS163']
dat[SRA=='SRR9330836', facetLabel := 'YPS163 x YJM145']
dat[SRA=='SRR9330813', facetLabel := 'YJM145 x CLIB413']
dat[SRA=='SRR9330815', facetLabel := 'CLIB413 x YJM978']
dat[SRA=='SRR9330812', facetLabel := 'YJM978 x YJM454']
dat[SRA=='SRR9330833', facetLabel := 'YJM454 x YPS1009']
dat[SRA=='SRR9330814', facetLabel := 'YPS1009 x I14']
dat[SRA=='SRR9330817', facetLabel := 'I14 x Y10']
dat[SRA=='SRR9330816', facetLabel := 'Y10 x PW5']
dat[SRA=='SRR9330830', facetLabel := 'PW5 x 273614']
dat[SRA=='SRR9330809', facetLabel := '273614 x YJM981']

dat[, facetLabel := factor(facetLabel, levels=c('YJM981 x CBS2888','CBS2888 x CLIB219',
                            'CLIB219 x M22','M22 x BY','BY x RM','RM x YPS163','YPS163 x YJM145',
                            'YJM145 x CLIB413','CLIB413 x YJM978','YJM978 x YJM454',
                            'YJM454 x YPS1009','YPS1009 x I14','I14 x Y10',
                            'Y10 x PW5','PW5 x 273614','273614 x YJM981'))]

dat[, minRepeats := factor(minRepeats, levels=c(5,10,20,30,35))]
                            
g <- ggplot(data=dat, aes(x=RPM)) + geom_histogram(bins=80) + theme_few(12) +
labs(x='Telomeric Reads per 1M reads', y='N') +
facet_grid(minRepeats~facetLabel) +
theme(text=element_text(family="Atkinson Hyperlegible", size = 20))

# dat[, bin := cut(RPM, breaks=seq(0,10000,10))]
# dat[, bin := as.numeric(as.factor(bin))]
# dat.ag <- dat[, .N, by=list(SRA, minRepeats, bin, facetLabel)]

# g <- ggplot(data=dat.ag, aes(x=bin, y=N)) + geom_smooth() + theme_few(12) +
# labs(x='Telomeric Reads per 1M reads', y='N') +
# facet_grid(minRepeats~facetLabel) +
# theme(text=element_text(family="Atkinson Hyperlegible", size = 20))