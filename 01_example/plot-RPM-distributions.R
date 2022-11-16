#!/usr/bin/env R

library(data.table)
library(ggplot2)
library(foreach)
library(ggthemes)

dat <- foreach(file=list.files(pattern='*.telo-counts.txt'), .combine='rbind') %do% {
    sra_id <- strsplit(file, '_')[[1]][1]
    tmp <- fread(file)
    setnames(tmp, c('id', 'nReads', 'nTeloReads'))
    tmp[, SRA := sra_id][]
}
dat[, RPM := 1e6 * (nTeloReads / nReads) ]


dat[SRA=='SRR9330808', facetLabel := 'YJM981 x CBS2888']
dat[SRA=='SRR9330831', facetLabel := 'CLIB219 x CBS2888']
dat[SRA=='SRR9330809', facetLabel := 'YJM981 x 273614']
dat[SRA=='SRR9330830', facetLabel := 'PW5 x 273614']

ggplot(data=dat, aes(x=RPM)) + geom_histogram(bins=80) + theme_few() +
labs(x='Telomeric Reads per 1M reads', y='N') +
facet_grid(facetLabel~.)

