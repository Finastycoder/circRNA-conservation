.script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NA_character_)
if (is.null(.script_file) || !length(.script_file) || is.na(.script_file)) {
  .file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  .script_file <- if (length(.file_arg)) sub("^--file=", "", .file_arg[1]) else NA_character_
}
if (!is.na(.script_file) && nzchar(.script_file)) {
  setwd(dirname(normalizePath(.script_file, winslash = "/", mustWork = TRUE)))
}
wkdir <- getwd()

# 读取本地数据 ------------------------------------------------------------------

library(dplyr)
#读取上游分析的原始数据
#读取find_circ结果
#物种：hsa
fc_hsa1 <- read.table('re_analysis/hsa/find_circ/hsa1.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_hsa1, 'fc_hsa1.csv')
fc_hsa1 <- fc_hsa1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))  #提取染色体位置作为circ_id产生新的列
fc_hsa1 <- fc_hsa1 %>% select(circ_id, hsa1 = n_reads)  #提取circ_id和n_reads列作为表达矩阵
fc_hsa2 <- read.table('re_analysis/hsa/find_circ/hsa2.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_hsa2, 'fc_hsa2.csv')
fc_hsa2 <- fc_hsa2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_hsa2 <- fc_hsa2 %>% select(circ_id, hsa2 = n_reads)
fc_hsa3 <- read.table('re_analysis/hsa/find_circ/hsa3.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_hsa3, 'fc_hsa3.csv')
fc_hsa3 <- fc_hsa3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_hsa3 <- fc_hsa3 %>% select(circ_id, hsa3 = n_reads)
fc_hsa <- merge(fc_hsa1, fc_hsa2, by = 'circ_id', all = T)  #将3个生物学重复样本合并
fc_hsa <- merge(fc_hsa, fc_hsa3, by = 'circ_id', all = T)  #将3个生物学重复样本合并
fc_hsa <- fc_hsa[!duplicated(fc_hsa$circ_id),]  #以circ_id为标准去重
write.csv(fc_hsa, 'fc_hsa.csv')  #输出表达矩阵
#物种：mfu
fc_mfu1 <- read.table('re_analysis/mfu/find_circ/mfu1.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mfu1, 'fc_mfu1.csv')
fc_mfu1 <- fc_mfu1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mfu1 <- fc_mfu1 %>% select(circ_id, mfu1 = n_reads)
fc_mfu2 <- read.table('re_analysis/mfu/find_circ/mfu2.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mfu2, 'fc_mfu2.csv')
fc_mfu2 <- fc_mfu2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mfu2 <- fc_mfu2 %>% select(circ_id, mfu2 = n_reads)
fc_mfu3 <- read.table('re_analysis/mfu/find_circ/mfu3.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mfu3, 'fc_mfu3.csv')
fc_mfu3 <- fc_mfu3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mfu3 <- fc_mfu3 %>% select(circ_id, mfu3 = n_reads)
fc_mfu <- merge(fc_mfu1, fc_mfu2, by = 'circ_id', all = T)  
fc_mfu <- merge(fc_mfu, fc_mfu3, by = 'circ_id', all = T)  
fc_mfu <- fc_mfu[!duplicated(fc_mfu$circ_id),]  
write.csv(fc_mfu, 'fc_mfu.csv')
#物种：mma
fc_mma1 <- read.table('re_analysis/mma/find_circ/mma1.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mma1, 'fc_mma1.csv')
fc_mma1 <- fc_mma1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mma1 <- fc_mma1 %>% select(circ_id, mma1 = n_reads)
fc_mma2 <- read.table('re_analysis/mma/find_circ/mma2.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mma2, 'fc_mma2.csv')
fc_mma2 <- fc_mma2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mma2 <- fc_mma2 %>% select(circ_id, mma2 = n_reads)
fc_mma3 <- read.table('re_analysis/mma/find_circ/mma3.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mma3, 'fc_mma3.csv')
fc_mma3 <- fc_mma3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mma3 <- fc_mma3 %>% select(circ_id, mma3 = n_reads)
fc_mma <- merge(fc_mma1, fc_mma2, by = 'circ_id', all = T)  
fc_mma <- merge(fc_mma, fc_mma3, by = 'circ_id', all = T)  
fc_mma <- fc_mma[!duplicated(fc_mma$circ_id),]  
write.csv(fc_mma, 'fc_mma.csv')
#物种：mmu
fc_mmu1 <- read.table('re_analysis/mmu/find_circ/mmu1.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mmu1, 'fc_mmu1.csv')
fc_mmu1 <- fc_mmu1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mmu1 <- fc_mmu1 %>% select(circ_id, mmu1 = n_reads)
fc_mmu2 <- read.table('re_analysis/mmu/find_circ/mmu2.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mmu2, 'fc_mmu2.csv')
fc_mmu2 <- fc_mmu2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mmu2 <- fc_mmu2 %>% select(circ_id, mmu2 = n_reads)
fc_mmu3 <- read.table('re_analysis/mmu/find_circ/mmu3.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mmu3, 'fc_mmu3.csv')
fc_mmu3 <- fc_mmu3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mmu3 <- fc_mmu3 %>% select(circ_id, mmu3 = n_reads)
fc_mmu <- merge(fc_mmu1, fc_mmu2, by = 'circ_id', all = T)  
fc_mmu <- merge(fc_mmu, fc_mmu3, by = 'circ_id', all = T)  
fc_mmu <- fc_mmu[!duplicated(fc_mmu$circ_id),]  
write.csv(fc_mmu, 'fc_mmu.csv')
#物种：mpi
fc_mpi1 <- read.table('re_analysis/mpi/find_circ/mpi1.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mpi1, 'fc_mpi1.csv')
fc_mpi1 <- fc_mpi1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mpi1 <- fc_mpi1 %>% select(circ_id, mpi1 = n_reads)
fc_mpi2 <- read.table('re_analysis/mpi/find_circ/mpi2.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mpi2, 'fc_mpi2.csv')
fc_mpi2 <- fc_mpi2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mpi2 <- fc_mpi2 %>% select(circ_id, mpi2 = n_reads)
fc_mpi3 <- read.table('re_analysis/mpi/find_circ/mpi3.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_mpi3, 'fc_mpi3.csv')
fc_mpi3 <- fc_mpi3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mpi3 <- fc_mpi3 %>% select(circ_id, mpi3 = n_reads)
fc_mpi <- merge(fc_mpi1, fc_mpi2, by = 'circ_id', all = F)  
fc_mpi <- merge(fc_mpi, fc_mpi3, by = 'circ_id', all = F)  
fc_mpi <- fc_mpi[!duplicated(fc_mpi$circ_id),]  
write.csv(fc_mpi, 'fc_mpi.csv')
#物种：rsi
fc_rsi1 <- read.table('re_analysis/rsi/find_circ/rsi1.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_rsi1, 'fc_rsi1.csv')
fc_rsi1 <- fc_rsi1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_rsi1 <- fc_rsi1 %>% select(circ_id, rsi1 = n_reads)
fc_rsi2 <- read.table('re_analysis/rsi/find_circ/rsi2.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_rsi2, 'fc_rsi2.csv')
fc_rsi2 <- fc_rsi2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_rsi2 <- fc_rsi2 %>% select(circ_id, rsi2 = n_reads)
fc_rsi3 <- read.table('re_analysis/rsi/find_circ/rsi3.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_rsi3, 'fc_rsi3.csv')
fc_rsi3 <- fc_rsi3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_rsi3 <- fc_rsi3 %>% select(circ_id, rsi3 = n_reads)
fc_rsi4 <- read.table('re_analysis/rsi/find_circ/rsi4.candidates.bed',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','n_reads','strand','n_uniq','uniq_bridges',
                                    'best_qual_left','best_qual_right','tissues','tiss_counts','edits','anchor_overlap',
                                    'breakpoints','signal','strandmatch','category'))
write.csv(fc_rsi4, 'fc_rsi4.csv')
fc_rsi4 <- fc_rsi4 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_rsi4 <- fc_rsi4 %>% select(circ_id, rsi4 = n_reads)
fc_rsi <- merge(fc_rsi1, fc_rsi2, by = 'circ_id', all = T)  
fc_rsi <- merge(fc_rsi, fc_rsi3, by = 'circ_id', all = T)
fc_rsi <- merge(fc_rsi, fc_rsi4, by = 'circ_id', all = T) 
fc_rsi <- fc_rsi[!duplicated(fc_rsi$circ_id),]  
write.csv(fc_rsi, 'fc_rsi.csv')
#读取CIRIquant结果
#在CIRIquant分析的结果中，染色体location的Start比其他两种分析方法的得到多1，因此需要用start=start - 1，把多的1减掉，才能与其他两种分析方法联合分析
#物种：hsa
library(rtracklayer)
cq_hsa1 <- import('re_analysis/hsa/CIRIquant/hsa1.gtf')
cq_hsa1 <- as.data.frame(cq_hsa1)
write.csv(cq_hsa1, 'cq_hsa1.csv')
cq_hsa1 <- cq_hsa1 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_hsa1 <- cq_hsa1 %>% select(circ_id = circ_ID, hsa1 = bsj)
cq_hsa2 <- import('re_analysis/hsa/CIRIquant/hsa2.gtf')
cq_hsa2 <- as.data.frame(cq_hsa2)
write.csv(cq_hsa2, 'cq_hsa2.csv')
cq_hsa2 <- cq_hsa2 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_hsa2 <- cq_hsa2 %>% select(circ_id = circ_ID, hsa2 = bsj)
cq_hsa3 <- import('re_analysis/hsa/CIRIquant/hsa3.gtf')
cq_hsa3 <- as.data.frame(cq_hsa3)
write.csv(cq_hsa3, 'cq_hsa3.csv')
cq_hsa3 <- cq_hsa3 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_hsa3 <- cq_hsa3 %>% select(circ_id = circ_ID, hsa3 = bsj)
cq_hsa <- merge(cq_hsa1, cq_hsa2, by = 'circ_id', all = T)
cq_hsa <- merge(cq_hsa, cq_hsa3, by = 'circ_id', all = T)
cq_hsa <- cq_hsa[!duplicated(cq_hsa$circ_id),]
write.csv(cq_hsa, 'cq_hsa.csv')
#物种：mfu
library(rtracklayer)
cq_mfu1 <- import('re_analysis/mfu/CIRIquant/mfu1.gtf')
cq_mfu1 <- as.data.frame(cq_mfu1)
write.csv(cq_mfu1, 'cq_mfu1.csv')
cq_mfu1 <- cq_mfu1 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mfu1 <- cq_mfu1 %>% select(circ_id = circ_ID, mfu1 = bsj)
cq_mfu2 <- import('re_analysis/mfu/CIRIquant/mfu2.gtf')
cq_mfu2 <- as.data.frame(cq_mfu2)
write.csv(cq_mfu2, 'cq_mfu2.csv')
cq_mfu2 <- cq_mfu2 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mfu2 <- cq_mfu2 %>% select(circ_id = circ_ID, mfu2 = bsj)
cq_mfu3 <- import('re_analysis/mfu/CIRIquant/mfu3.gtf')
cq_mfu3 <- as.data.frame(cq_mfu3)
write.csv(cq_mfu3, 'cq_mfu3.csv')
cq_mfu3 <- cq_mfu3 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mfu3 <- cq_mfu3 %>% select(circ_id = circ_ID, mfu3 = bsj)
cq_mfu <- merge(cq_mfu1, cq_mfu2, by = 'circ_id', all = T)
cq_mfu <- merge(cq_mfu, cq_mfu3, by = 'circ_id', all = T)
cq_mfu <- cq_mfu[!duplicated(cq_mfu$circ_id),]
write.csv(cq_mfu, 'cq_mfu.csv')
#物种：mma
library(rtracklayer)
cq_mma1 <- import('re_analysis/mma/CIRIquant/mma1.gtf')
cq_mma1 <- as.data.frame(cq_mma1)
write.csv(cq_mma1, 'cq_mma1.csv')
cq_mma1 <- cq_mma1 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mma1 <- cq_mma1 %>% select(circ_id = circ_ID, mma1 = bsj)
cq_mma2 <- import('re_analysis/mma/CIRIquant/mma2.gtf')
cq_mma2 <- as.data.frame(cq_mma2)
write.csv(cq_mma2, 'cq_mma2.csv')
cq_mma2 <- cq_mma2 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mma2 <- cq_mma2 %>% select(circ_id = circ_ID, mma2 = bsj)
cq_mma3 <- import('re_analysis/mma/CIRIquant/mma3.gtf')
cq_mma3 <- as.data.frame(cq_mma3)
write.csv(cq_mma3, 'cq_mma3.csv')
cq_mma3 <- cq_mma3 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mma3 <- cq_mma3 %>% select(circ_id = circ_ID, mma3 = bsj)
cq_mma <- merge(cq_mma1, cq_mma2, by = 'circ_id', all = T)
cq_mma <- merge(cq_mma, cq_mma3, by = 'circ_id', all = T)
cq_mma <- cq_mma[!duplicated(cq_mma$circ_id),]
write.csv(cq_mma, 'cq_mma.csv')
#物种：mmu
library(rtracklayer)
cq_mmu1 <- import('re_analysis/mmu/CIRIquant/mmu1.gtf')
cq_mmu1 <- as.data.frame(cq_mmu1)
write.csv(cq_mmu1, 'cq_mmu1.csv')
cq_mmu1 <- cq_mmu1 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mmu1 <- cq_mmu1 %>% select(circ_id = circ_ID, mmu1 = bsj)
cq_mmu2 <- import('re_analysis/mmu/CIRIquant/mmu2.gtf')
cq_mmu2 <- as.data.frame(cq_mmu2)
write.csv(cq_mmu2, 'cq_mmu2.csv')
cq_mmu2 <- cq_mmu2 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mmu2 <- cq_mmu2 %>% select(circ_id = circ_ID, mmu2 = bsj)
cq_mmu3 <- import('re_analysis/mmu/CIRIquant/mmu3.gtf')
cq_mmu3 <- as.data.frame(cq_mmu3)
write.csv(cq_mmu3, 'cq_mmu3.csv')
cq_mmu3 <- cq_mmu3 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mmu3 <- cq_mmu3 %>% select(circ_id = circ_ID, mmu3 = bsj)
cq_mmu <- merge(cq_mmu1, cq_mmu2, by = 'circ_id', all = T)
cq_mmu <- merge(cq_mmu, cq_mmu3, by = 'circ_id', all = T)
cq_mmu <- cq_mmu[!duplicated(cq_mmu$circ_id),]
write.csv(cq_mmu, 'cq_mmu.csv')
#物种：mpi
library(rtracklayer)
cq_mpi1 <- import('re_analysis/mpi/CIRIquant/mpi1.gtf')
cq_mpi1 <- as.data.frame(cq_mpi1)
write.csv(cq_mpi1, 'cq_mpi1.csv')
cq_mpi1 <- cq_mpi1 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mpi1 <- cq_mpi1 %>% select(circ_id = circ_ID, mpi1 = bsj)
cq_mpi2 <- import('re_analysis/mpi/CIRIquant/mpi2.gtf')
cq_mpi2 <- as.data.frame(cq_mpi2)
write.csv(cq_mpi2, 'cq_mpi2.csv')
cq_mpi2 <- cq_mpi2 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mpi2 <- cq_mpi2 %>% select(circ_id = circ_ID, mpi2 = bsj)
cq_mpi3 <- import('re_analysis/mpi/CIRIquant/mpi3.gtf')
cq_mpi3 <- as.data.frame(cq_mpi3)
write.csv(cq_mpi3, 'cq_mpi3.csv')
cq_mpi3 <- cq_mpi3 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mpi3 <- cq_mpi3 %>% select(circ_id = circ_ID, mpi3 = bsj)
cq_mpi <- merge(cq_mpi1, cq_mpi2, by = 'circ_id', all = F)
cq_mpi <- merge(cq_mpi, cq_mpi3, by = 'circ_id', all = F)
cq_mpi <- cq_mpi[!duplicated(cq_mpi$circ_id),]
write.csv(cq_mpi, 'cq_mpi.csv')
#物种：rsi
library(rtracklayer)
cq_rsi1 <- import('re_analysis/rsi/CIRIquant/rsi1.gtf')
cq_rsi1 <- as.data.frame(cq_rsi1)
write.csv(cq_rsi1, 'cq_rsi1.csv')
cq_rsi1 <- cq_rsi1 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_rsi1 <- cq_rsi1 %>% select(circ_id = circ_ID, rsi1 = bsj)
cq_rsi2 <- import('re_analysis/rsi/CIRIquant/rsi2.gtf')
cq_rsi2 <- as.data.frame(cq_rsi2)
write.csv(cq_rsi2, 'cq_rsi2.csv')
cq_rsi2 <- cq_rsi2 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_rsi2 <- cq_rsi2 %>% select(circ_id = circ_ID, rsi2 = bsj)
cq_rsi3 <- import('re_analysis/rsi/CIRIquant/rsi3.gtf')
cq_rsi3 <- as.data.frame(cq_rsi3)
write.csv(cq_rsi3, 'cq_rsi3.csv')
cq_rsi3 <- cq_rsi3 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_rsi3 <- cq_rsi3 %>% select(circ_id = circ_ID, rsi3 = bsj)
cq_rsi4 <- import('re_analysis/rsi/CIRIquant/rsi4.gtf')
cq_rsi4 <- as.data.frame(cq_rsi4)
write.csv(cq_rsi4, 'cq_rsi4.csv')
cq_rsi4 <- cq_rsi4 %>% mutate(start = start - 1, circ_ID=paste(seqnames, ':', start, '|', end, sep = ''))
cq_rsi4 <- cq_rsi4 %>% select(circ_id = circ_ID, rsi4 = bsj)
cq_rsi <- merge(cq_rsi1, cq_rsi2, by = 'circ_id', all = T)
cq_rsi <- merge(cq_rsi, cq_rsi3, by = 'circ_id', all = T)
cq_rsi <- merge(cq_rsi, cq_rsi4, by = 'circ_id', all = T)
cq_rsi <- cq_rsi[!duplicated(cq_rsi$circ_id),]
write.csv(cq_rsi, 'cq_rsi.csv')
#读取CIRCexplorer3结果
#物种：hsa
ce_hsa1 <- read.table('re_analysis/hsa/CIRCexplorer3/hsa1.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_hsa1, 'ce_hsa1.csv')
ce_hsa1 <- ce_hsa1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_hsa1 <- ce_hsa1 %>% select(circ_id, hsa1 = readNumber)
ce_hsa2 <- read.table('re_analysis/hsa/CIRCexplorer3/hsa2.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_hsa2, 'ce_hsa2.csv')
ce_hsa2 <- ce_hsa2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_hsa2 <- ce_hsa2 %>% select(circ_id, hsa2 = readNumber)
ce_hsa3 <- read.table('re_analysis/hsa/CIRCexplorer3/hsa3.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_hsa3, 'ce_hsa3.csv')
ce_hsa3 <- ce_hsa3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_hsa3 <- ce_hsa3 %>% select(circ_id, hsa3 = readNumber)
ce_hsa <- merge(ce_hsa1, ce_hsa2, by = 'circ_id', all = T)
ce_hsa <- merge(ce_hsa, ce_hsa3, by = 'circ_id', all = T)
ce_hsa <- ce_hsa[!duplicated(ce_hsa$circ_id),]
write.csv(ce_hsa, 'ce_hsa.csv')
#物种：mfu
ce_mfu1 <- read.table('re_analysis/mfu/CIRCexplorer3/mfu1.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mfu1, 'ce_mfu1.csv')
ce_mfu1 <- ce_mfu1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mfu1 <- ce_mfu1 %>% select(circ_id, mfu1 = readNumber)
ce_mfu2 <- read.table('re_analysis/mfu/CIRCexplorer3/mfu2.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mfu2, 'ce_mfu2.csv')
ce_mfu2 <- ce_mfu2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mfu2 <- ce_mfu2 %>% select(circ_id, mfu2 = readNumber)
ce_mfu3 <- read.table('re_analysis/mfu/CIRCexplorer3/mfu3.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mfu3, 'ce_mfu3.csv')
ce_mfu3 <- ce_mfu3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mfu3 <- ce_mfu3 %>% select(circ_id, mfu3 = readNumber)
ce_mfu <- merge(ce_mfu1, ce_mfu2, by = 'circ_id', all = T)
ce_mfu <- merge(ce_mfu, ce_mfu3, by = 'circ_id', all = T)
ce_mfu <- ce_mfu[!duplicated(ce_mfu$circ_id),]
write.csv(ce_mfu, 'ce_mfu.csv')
#物种：mma
ce_mma1 <- read.table('re_analysis/mma/CIRCexplorer3/mma1.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mma1, 'ce_mma1.csv')
ce_mma1 <- ce_mma1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mma1 <- ce_mma1 %>% select(circ_id, mma1 = readNumber)
ce_mma2 <- read.table('re_analysis/mma/CIRCexplorer3/mma2.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mma2, 'ce_mma2.csv')
ce_mma2 <- ce_mma2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mma2 <- ce_mma2 %>% select(circ_id, mma2 = readNumber)
ce_mma3 <- read.table('re_analysis/mma/CIRCexplorer3/mma3.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mma3, 'ce_mma3.csv')
ce_mma3 <- ce_mma3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mma3 <- ce_mma3 %>% select(circ_id, mma3 = readNumber)
ce_mma <- merge(ce_mma1, ce_mma2, by = 'circ_id', all = T)
ce_mma <- merge(ce_mma, ce_mma3, by = 'circ_id', all = T)
ce_mma <- ce_mma[!duplicated(ce_mma$circ_id),]
write.csv(ce_mma, 'ce_mma.csv')
#物种：mmu
ce_mmu1 <- read.table('re_analysis/mmu/CIRCexplorer3/mmu1.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mmu1, 'ce_mmu1.csv')
ce_mmu1 <- ce_mmu1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mmu1 <- ce_mmu1 %>% select(circ_id, mmu1 = readNumber)
ce_mmu2 <- read.table('re_analysis/mmu/CIRCexplorer3/mmu2.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mmu2, 'ce_mmu2.csv')
ce_mmu2 <- ce_mmu2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mmu2 <- ce_mmu2 %>% select(circ_id, mmu2 = readNumber)
ce_mmu3 <- read.table('re_analysis/mmu/CIRCexplorer3/mmu3.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mmu3, 'ce_mmu3.csv')
ce_mmu3 <- ce_mmu3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mmu3 <- ce_mmu3 %>% select(circ_id, mmu3 = readNumber)
ce_mmu <- merge(ce_mmu1, ce_mmu2, by = 'circ_id', all = T)
ce_mmu <- merge(ce_mmu, ce_mmu3, by = 'circ_id', all = T)
ce_mmu <- ce_mmu[!duplicated(ce_mmu$circ_id),]
write.csv(ce_mmu, 'ce_mmu.csv')
#物种：mpi
ce_mpi1 <- read.table('re_analysis/mpi/CIRCexplorer3/mpi1.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mpi1, 'ce_mpi1.csv')
ce_mpi1 <- ce_mpi1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mpi1 <- ce_mpi1 %>% select(circ_id, mpi1 = readNumber)
ce_mpi2 <- read.table('re_analysis/mpi/CIRCexplorer3/mpi2.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mpi2, 'ce_mpi2.csv')
ce_mpi2 <- ce_mpi2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mpi2 <- ce_mpi2 %>% select(circ_id, mpi2 = readNumber)
ce_mpi3 <- read.table('re_analysis/mpi/CIRCexplorer3/mpi3.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_mpi3, 'ce_mpi3.csv')
ce_mpi3 <- ce_mpi3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mpi3 <- ce_mpi3 %>% select(circ_id, mpi3 = readNumber)
ce_mpi <- merge(ce_mpi1, ce_mpi2, by = 'circ_id', all = F)
ce_mpi <- merge(ce_mpi, ce_mpi3, by = 'circ_id', all = F)
ce_mpi <- ce_mpi[!duplicated(ce_mpi$circ_id),]
write.csv(ce_mpi, 'ce_mpi.csv')
#物种：rsi
ce_rsi1 <- read.table('re_analysis/rsi/CIRCexplorer3/rsi1.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_rsi1, 'ce_rsi1.csv')
ce_rsi1 <- ce_rsi1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_rsi1 <- ce_rsi1 %>% select(circ_id, rsi1 = readNumber)
ce_rsi2 <- read.table('re_analysis/rsi/CIRCexplorer3/rsi2.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_rsi2, 'ce_rsi2.csv')
ce_rsi2 <- ce_rsi2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_rsi2 <- ce_rsi2 %>% select(circ_id, rsi2 = readNumber)
ce_rsi3 <- read.table('re_analysis/rsi/CIRCexplorer3/rsi3.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_rsi3, 'ce_rsi3.csv')
ce_rsi3 <- ce_rsi3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_rsi3 <- ce_rsi3 %>% select(circ_id, rsi3 = readNumber)
ce_rsi4 <- read.table('re_analysis/rsi/CIRCexplorer3/rsi4.txt',
                      header = F,
                      sep = '\t',
                      col.names = c('chrom','start','end','name','score','strand','thickStart','thickEnd','itemRgb',
                                    'exonCount','exonSize','exonOffsets','readNumber','circType','geneName',
                                    'isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore'))
write.csv(ce_rsi4, 'ce_rsi4.csv')
ce_rsi4 <- ce_rsi4 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_rsi4 <- ce_rsi4 %>% select(circ_id, rsi4 = readNumber)
ce_rsi <- merge(ce_rsi1, ce_rsi2, by = 'circ_id', all = T)
ce_rsi <- merge(ce_rsi, ce_rsi3, by = 'circ_id', all = T)
ce_rsi <- merge(ce_rsi, ce_rsi4, by = 'circ_id', all = T)
ce_rsi <- ce_rsi[!duplicated(ce_rsi$circ_id),]
write.csv(ce_rsi, 'ce_rsi.csv')

# 按照物种将三种分析方法分的得到的circRNA进行Venn图绘制 ----------------------------------------

#按照每个物种，将三种方法分析得到的环状RNA进行Venn图绘制
#物种：hsa
library(VennDiagram)
Venn_hsa <- venn.diagram(x = list(ce_hsa$circ_id, cq_hsa$circ_id, fc_hsa$circ_id), category.names = c('CIRCexplorer3', 'CIRIquant', 'find_circ'),
                             disable.logging = T, lwd = 1, col = 'black', fill=c('red','blue','yellow'), alpha = 0.50, cex = 1, fontfamily = 'serif',
                             scaled = F, margin = 0.1, cat.pos = c(330,30,180), cat.dist = c(0.07,0.06,0.04), cat.cex = 1, cat.col = rep('black',3),
                             filename = NULL, height = 1600, width = 1600, resolution = 400, compression = 'lzw')
pdf('Venn_hsa.pdf')
grid.draw(Venn_hsa)
dev.off()
#物种：mfu
library(VennDiagram)
Venn_mfu <- venn.diagram(x = list(ce_mfu$circ_id, cq_mfu$circ_id, fc_mfu$circ_id), category.names = c('CIRCexplorer3', 'CIRIquant', 'find_circ'),
                         disable.logging = T, lwd = 1, col = 'black', fill=c('red','blue','yellow'), alpha = 0.50, cex = 1, fontfamily = 'serif',
                         scaled = F, margin = 0.1, cat.pos = c(330,30,180), cat.dist = c(0.07,0.06,0.04), cat.cex = 1, cat.col = rep('black',3),
                         filename = NULL, height = 1600, width = 1600, resolution = 400, compression = 'lzw')
pdf('Venn_mfu.pdf')
grid.draw(Venn_mfu)
dev.off()
#物种：mma
library(VennDiagram)
Venn_mma <- venn.diagram(x = list(ce_mma$circ_id, cq_mma$circ_id, fc_mma$circ_id), category.names = c('CIRCexplorer3', 'CIRIquant', 'find_circ'),
                         disable.logging = T, lwd = 1, col = 'black', fill=c('red','blue','yellow'), alpha = 0.50, cex = 1, fontfamily = 'serif',
                         scaled = F, margin = 0.1, cat.pos = c(330,30,180), cat.dist = c(0.07,0.06,0.04), cat.cex = 1, cat.col = rep('black',3),
                         filename = NULL, height = 1600, width = 1600, resolution = 400, compression = 'lzw')
pdf('Venn_mma.pdf')
grid.draw(Venn_mma)
dev.off()
#物种：mmu
library(VennDiagram)
Venn_mmu <- venn.diagram(x = list(ce_mmu$circ_id, cq_mmu$circ_id, fc_mmu$circ_id), category.names = c('CIRCexplorer3', 'CIRIquant', 'find_circ'),
                         disable.logging = T, lwd = 1, col = 'black', fill=c('red','blue','yellow'), alpha = 0.50, cex = 1, fontfamily = 'serif',
                         scaled = F, margin = 0.1, cat.pos = c(330,30,180), cat.dist = c(0.07,0.06,0.04), cat.cex = 1, cat.col = rep('black',3),
                         filename = NULL, height = 1600, width = 1600, resolution = 400, compression = 'lzw')
pdf('Venn_mmu.pdf')
grid.draw(Venn_mmu)
dev.off()
#物种：mpi
library(VennDiagram)
Venn_mpi <- venn.diagram(x = list(ce_mpi$circ_id, cq_mpi$circ_id, fc_mpi$circ_id), category.names = c('CIRCexplorer3', 'CIRIquant', 'find_circ'),
                         disable.logging = T, lwd = 1, col = 'black', fill=c('red','blue','yellow'), alpha = 0.50, cex = 1, fontfamily = 'serif',
                         scaled = F, margin = 0.1, cat.pos = c(330,30,180), cat.dist = c(0.07,0.06,0.04), cat.cex = 1, cat.col = rep('black',3),
                         filename = NULL, height = 1600, width = 1600, resolution = 400, compression = 'lzw')
pdf('Venn_mpi.pdf')
grid.draw(Venn_mpi)
dev.off()
#物种：rsi
library(VennDiagram)
Venn_rsi <- venn.diagram(x = list(ce_rsi$circ_id, cq_rsi$circ_id, fc_rsi$circ_id), category.names = c('CIRCexplorer3', 'CIRIquant', 'find_circ'),
                         disable.logging = T, lwd = 1, col = 'black', fill=c('red','blue','yellow'), alpha = 0.50, cex = 1, fontfamily = 'serif',
                         scaled = F, margin = 0.1, cat.pos = c(330,30,180), cat.dist = c(0.07,0.06,0.04), cat.cex = 1, cat.col = rep('black',3),
                         filename = NULL, height = 1600, width = 1600, resolution = 400, compression = 'lzw')
pdf('Venn_rsi.pdf')
grid.draw(Venn_rsi)
dev.off()


# 以gene_name为桥梁跨物种分析 ------------------------------------------------------

#建立以gene_name为桥接的物种间分析连接，即将每个物种中表达的circRNA以host gene的gene_name表示
#以CIRIquant的分析结果为后续分析基础，原因：(1)CIRIquant是最新发表的分析软件；(2)CIRIquant提供了分析得到的circRNA全长序列
#在hsa，mma，mmu，rsi的注释文件中包含gene_name，因此CIRIquant分析时能够产生gene_id列和gene_name列，可直接进行转化
#读取上述4个物种中，三种分析方法merge后的circRNA的circ_id
hsa_id <- merge(ce_hsa, cq_hsa, by = 'circ_id', all = F)
hsa_id <- merge(hsa_id, fc_hsa, by = 'circ_id', all = F)
hsa_id <- data.frame(circ_id = hsa_id$circ_id)
mma_id <- merge(ce_mma, cq_mma, by = 'circ_id', all = F)
mma_id <- merge(mma_id, fc_mma, by = 'circ_id', all = F)
mma_id <- data.frame(circ_id = mma_id$circ_id)
mmu_id <- merge(ce_mmu, cq_mmu, by = 'circ_id', all = F)
mmu_id <- merge(mmu_id, fc_mmu, by = 'circ_id', all = F)
mmu_id <- data.frame(circ_id = mmu_id$circ_id)
rsi_id <- merge(ce_rsi, cq_rsi, by = 'circ_id', all = F)
rsi_id <- merge(rsi_id, fc_rsi, by = 'circ_id', all = F)
rsi_id <- data.frame(circ_id = rsi_id$circ_id)
#读取CIRIquant分析得到的上述4个物种的每个样本的结果
cq_hsa1 <- read.csv('cq_hsa1.csv', header = T, row.names = 1)
cq_hsa2 <- read.csv('cq_hsa2.csv', header = T, row.names = 1)
cq_hsa3 <- read.csv('cq_hsa3.csv', header = T, row.names = 1)
cq_mma1 <- read.csv('cq_mma1.csv', header = T, row.names = 1)
cq_mma2 <- read.csv('cq_mma2.csv', header = T, row.names = 1)
cq_mma3 <- read.csv('cq_mma3.csv', header = T, row.names = 1)
cq_mmu1 <- read.csv('cq_mmu1.csv', header = T, row.names = 1)
cq_mmu2 <- read.csv('cq_mmu2.csv', header = T, row.names = 1)
cq_mmu3 <- read.csv('cq_mmu3.csv', header = T, row.names = 1)
cq_rsi1 <- read.csv('cq_rsi1.csv', header = T, row.names = 1)
cq_rsi2 <- read.csv('cq_rsi2.csv', header = T, row.names = 1)
cq_rsi3 <- read.csv('cq_rsi3.csv', header = T, row.names = 1)
cq_rsi4 <- read.csv('cq_rsi4.csv', header = T, row.names = 1)
#因为在分析时，除mpi外每个物种中的3个重复样本是取并集，以充分囊括所有可能的在该物种中表达的circRNA
#因此将每个物种的3或4个重复样本进行行合并，通过circ_id进行去重
#去掉错误的circ_id列（start多1），将start-1后产生新的circ_id，以使其能够与三种方法交集得到的circ_id进行merge
#在CIRIquant分析的结果中，一个circ_id会对应多个gene_id或gene_name，这是因为circ_id是染色体位置，这个位置上可能包含protein gene，pseudogene，nc gene
#因此使用separate_rows将逗号分隔的这些gene_name单独成行，再通过gene_name去重，得到三种分析方法交集circRNA的host gene_name
#物种：hsa
combined_hsa <- do.call(rbind, list(cq_hsa1, cq_hsa2, cq_hsa3))
combined_hsa <- combined_hsa[!duplicated(combined_hsa$circ_id),]
combined_hsa <- subset(combined_hsa, select = -circ_id)
library(dplyr)
combined_hsa <- combined_hsa %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
hsa.cq <- merge(hsa_id, combined_hsa, by = 'circ_id', all = F)
write.csv(hsa.cq, 'hsa.cq.csv', row.names = F)
library(tidyr)
hsa.cq2 <- separate_rows(hsa.cq, gene_name, sep = ',')
hsa.cq2 <- unique(hsa.cq2)
hsa.cq2 <- hsa.cq2[!duplicated(hsa.cq2$gene_name),]
#物种：mma
combined_mma <- do.call(rbind, list(cq_mma1, cq_mma2, cq_mma3))
combined_mma <- combined_mma[!duplicated(combined_mma$circ_id),]
combined_mma <- subset(combined_mma, select = -circ_id)
library(dplyr)
combined_mma <- combined_mma %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
mma.cq <- merge(mma_id, combined_mma, by = 'circ_id', all = F)
write.csv(mma.cq, 'mma.cq.csv', row.names = F)
library(tidyr)
mma.cq2 <- separate_rows(mma.cq, gene_name, sep = ',')
mma.cq2 <- unique(mma.cq2)
mma.cq2 <- mma.cq2[!duplicated(mma.cq2$gene_name),]
#物种：mmu
combined_mmu <- do.call(rbind, list(cq_mmu1, cq_mmu2, cq_mmu3))
combined_mmu <- combined_mmu[!duplicated(combined_mmu$circ_id),]
combined_mmu <- subset(combined_mmu, select = -circ_id)
library(dplyr)
combined_mmu <- combined_mmu %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
mmu.cq <- merge(mmu_id, combined_mmu, by = 'circ_id', all = F)
write.csv(mmu.cq, 'mmu.cq.csv', row.names = F)
library(tidyr)
mmu.cq2 <- separate_rows(mmu.cq, gene_name, sep = ',')
mmu.cq2 <- unique(mmu.cq2)
mmu.cq2 <- mmu.cq2[!duplicated(mmu.cq2$gene_name),]
#物种：rsi
combined_rsi <- do.call(rbind, list(cq_rsi1, cq_rsi2, cq_rsi3, cq_rsi4))
combined_rsi <- combined_rsi[!duplicated(combined_rsi$circ_id),]
combined_rsi <- subset(combined_rsi, select = -circ_id)
library(dplyr)
combined_rsi <- combined_rsi %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
rsi.cq <- merge(rsi_id, combined_rsi, by = 'circ_id', all = F)
write.csv(rsi.cq, 'rsi.cq.csv', row.names = F)
library(tidyr)
rsi.cq2 <- separate_rows(rsi.cq, gene_name, sep = ',')
rsi.cq2 <- unique(rsi.cq2)
rsi.cq2 <- rsi.cq2[!duplicated(rsi.cq2$gene_name),]
#对于mfu和mpi两个物种来说，上游分析中使用的gtf文件不包含gene_name，只包含gene_id
#但gff3文件包含gene_id和gene_name，因此需要手动依据gff3文件构建gene_id=gene_name的索引，在gff3文件中分别为ID和Name
#读取mfu和mpi的gff3文件转换为数据框格式
library(rtracklayer)
gffmfu <- import('mfu.gff3')
gffmfu <- as.data.frame(gffmfu)
gffmpi <- import('mpi.gff3')
gffmpi <- as.data.frame(gffmpi)
#读取mfu和mpi中，三种分析方法merge后的circRNA的circ_id
mfu_id <- merge(ce_mfu, cq_mfu, by = 'circ_id', all = F)
mfu_id <- merge(mfu_id, fc_mfu, by = 'circ_id', all = F)
mfu_id <- data.frame(circ_id = mfu_id$circ_id)
mpi_id <- merge(ce_mpi, cq_mpi, by = 'circ_id', all = F)
mpi_id <- merge(mpi_id, fc_mpi, by = 'circ_id', all = F)
mpi_id <- data.frame(circ_id = mpi_id$circ_id)
#读取CIRIquant分析得到的mfu和mpi的每个样本的结果
cq_mfu1 <- read.csv('cq_mfu1.csv', header = T, row.names = 1)
cq_mfu2 <- read.csv('cq_mfu2.csv', header = T, row.names = 1)
cq_mfu3 <- read.csv('cq_mfu3.csv', header = T, row.names = 1)
cq_mpi1 <- read.csv('cq_mpi1.csv', header = T, row.names = 1)
cq_mpi2 <- read.csv('cq_mpi2.csv', header = T, row.names = 1)
cq_mpi3 <- read.csv('cq_mpi3.csv', header = T, row.names = 1)
#对mfu和mpi每个物种中的重复样本进行行合并，去掉错误的circ_id列，产生新的circ_id列（start-1），与三种分析方法交集的circ_id进行merge
#同样1个circ_id会对应多个gene_id，使用separate_rows进行拆分并产生新的行，再进行去重
#只取gff3文件的ID和Name列，将去重后的gene_id与gff3的ID列进行merge，得到gene_id对应的Name（gene_name）列，修改好列名为gene_name
#物种：mfu
combined_mfu <- do.call(rbind, list(cq_mfu1, cq_mfu2, cq_mfu3))
combined_mfu <- combined_mfu[!duplicated(combined_mfu$circ_id),]
combined_mfu <- subset(combined_mfu, select = -circ_id)
library(dplyr)
combined_mfu <- combined_mfu %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
mfu.cq <- merge(mfu_id, combined_mfu, by = 'circ_id', all = F)
write.csv(mfu.cq, 'mfu.cq.csv', row.names = F)
library(tidyr)
mfu.cq2 <- separate_rows(mfu.cq, gene_id, sep = ',')
mfu.cq2 <- unique(mfu.cq2)
mfu.cq2 <- mfu.cq2[!duplicated(mfu.cq2$gene_id),]
mfu.gff <- gffmfu[,c('ID', 'Name')]
mfu.cq2 <- merge(mfu.cq2, mfu.gff, by.x = 'gene_id', by.y = 'ID', all = F)
colnames(mfu.cq2)[16] <- 'gene_name'
#物种：mpi
combined_mpi <- do.call(rbind, list(cq_mpi1, cq_mpi2, cq_mpi3))
combined_mpi <- combined_mpi[!duplicated(combined_mpi$circ_id),]
combined_mpi <- subset(combined_mpi, select = -circ_id)
library(dplyr)
combined_mpi <- combined_mpi %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
mpi.cq <- merge(mpi_id, combined_mpi, by = 'circ_id', all = F)
write.csv(mpi.cq, 'mpi.cq.csv', row.names = F)
library(tidyr)
mpi.cq2 <- separate_rows(mpi.cq, gene_id, sep = ',')
mpi.cq2 <- unique(mpi.cq2)
mpi.cq2 <- mpi.cq2[!duplicated(mpi.cq2$gene_id),]
mpi.gff <- gffmpi[,c('ID', 'Name')]
mpi.cq2 <- merge(mpi.cq2, mpi.gff, by.x = 'gene_id', by.y = 'ID', all = F)
colnames(mpi.cq2)[16] <- 'gene_name'
#以上分析过程已经获得每个物种中三种方法交集表达的circRNA对应的host gene name
#接下来以gene_name为桥接，绘制6个物种的Venn图，得到455个在mpi中特异性表达的host gene_name
library(venn)
library(ggplot2)
venn_list <- list(hsa.cq2$gene_name,
                  mfu.cq2$gene_name,
                  mma.cq2$gene_name,
                  mmu.cq2$gene_name,
                  mpi.cq2$gene_name,
                  rsi.cq2$gene_name)
names(venn_list) <- c('hsa', 'mfu', 'mma', 'mmu', 'mpi', 'rsi')
venn_list = purrr::map(venn_list, na.omit)
pdf('venn_multispecies.pdf')
venn_multispecies <- venn(venn_list,
                          zcolor = 'style',
                          ilabels = F,
                          ellipse = F,
                          opacity = 0.5,
                          box = F,
                          borders = F,
                          ilcs = 0.8,
                          sncs = 1.5)
dev.off()
#输出在mpi中特异存在的host gene name
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ',')
inter <- subset(inter, select = -..values.. )
inter <- subset(inter, select = -..set.. )
write.table(inter, "venn_multispecies.csv", row.names = FALSE, sep = ',', quote = FALSE)



# 联合转录组数据分析主要以环状形式表达的host gene --------------------------------------------

#联合mRNA转录组分析，找到这455个基因中哪些主要以环状形式表达，哪些主要以线性mRNA形式表达，筛选那些主要以环状形式表达的基因及其产生的circrna
#手动从venn_multispecies.csv中摘出mpi中特异性存在的455个host gene name，以mpi_specific_name读入
#通过mpi.gff文件将gene_name转换为gene_id，因为1个gene_name包含多个转录本，每个转录本有一个gene_id，因此会产生很多gene_id
#按照gene_name去重，删除gene_name列，将ID列名改为gene_id列名，此时只有一列gene_id
#将gene_id输出到本地，用excel的分列选项将gene_id单独分出来（排除t1,t2等转录本），重新读入
#将gene_id与mpi.cq2进行merge，得到包含gene_id，gene_name的455个特异性表达基因的转换列表
mpi_specific_name <- read.csv('mpi_specific_name.csv', header = T)
mpi_specific_geneid <- merge(mpi_specific_name, mpi.gff, by.x = "gene_name", by.y = 'Name', all = F)
mpi_specific_geneid <- mpi_specific_geneid[!duplicated(mpi_specific_geneid$gene_name),]
mpi_specific_geneid <- subset(mpi_specific_geneid, select = -gene_name)
colnames(mpi_specific_geneid)[1] <- 'gene_id'
write.csv(mpi_specific_geneid, 'mpi_specific_geneid.csv', row.names = F) #输出后手动分列
mpi_specific_geneid <- read.csv('mpi_specific_geneid.csv', header = T)
mpi_specific_convert <- merge(mpi_specific_geneid, mpi.cq2, by = 'gene_id', all = F)
mpi_specific_convert <- mpi_specific_convert[,c(1,16)]
mpi_specific_convert <- mpi_specific_convert[!duplicated(mpi_specific_convert$gene_name),]
#455个基因有455个gene_name，以及对应的455个gene_id，但对应的应该会有>455个circ_id
#产生包含所有选择性剪接信息的mpi中circrna信息表——mpi.cq.all
#用包含选择性剪接信息的circnra信息表与455个基因的gene_id做merge，这样保留了circrna的选择性剪接信息，455个gene对应着651个circrna（包含重复circrna）
#但在这455个gene中，存在着某些gene对应的是同一条circrna，用circ_id去重后，455个gene对应着549个circrna
mpi.cq.all <- separate_rows(mpi.cq, gene_id, sep = ',')
mpi_specific_circid <- merge(mpi.cq.all, mpi_specific_geneid, by = 'gene_id', all = F)
mpi_specific_circid <- mpi_specific_circid[,c(1,2)]
mpi_specific_allcircid <- mpi_specific_circid
mpi_specific_circid <- mpi_specific_circid[!duplicated(mpi_specific_circid$circ_id),]  
write.csv(mpi_specific_circid, 'mpi_specific_circid.csv', row.names = F)  
#对于这455个基因中，存在1个gene_id对应多个circ_id的情况，对于这样的gene，它产生circrna的bsj reads应该是它产生的所有circrna的bsj reads的和
#在CIRIquant的分析结果中，还会存在1个circ_id对应多个gene_id的情况，事实上我们不清楚这个circrna究竟是哪个gene产生的
#因此将这个circrna的表达值算在每个产生它的gene头上，这样让每个gene对应的所有circrna的表达值都按样本相加，这样的gene共有651个
#基于CIRIquant检测到的bsj reads，mpi物种CIRIquant分析得到的表达谱是cq_mpi，将circ_id与表达谱进行merge，得到651个gene（包含重复）对应的circ的表达值
#将bsj reads表达值转化为数值型
#对于circ_id不同，gene_id相同的行，将mpi1,2,3中对应的bsj reads相加并对gene_id去重，得到455个gene所有对应circrna的bsj reads总和（按样本）
#比如，对于两个circrna，它们都对应1个gene_id A，则将gene_id A的两个mpi1的bsj reads相加，作为其在mpi1样本中产生的总的circrna bsj reads
mpi_specific_circ_matrix <- merge(mpi_specific_allcircid, cq_mpi, by = 'circ_id', all = F)
library(dplyr)
mpi_specific_circ_matrix <- mpi_specific_circ_matrix %>% mutate_at(vars(mpi1,mpi2,mpi3), as.numeric)
mpi_specific_circ_matrix <- subset(mpi_specific_circ_matrix, select = -circ_id)
mpi_specific_circ_matrix <- mpi_specific_circ_matrix %>% group_by(gene_id) %>% 
  summarise(circ_mpi1 = sum(mpi1), circ_mpi2 = sum(mpi2), circ_mpi3 = sum(mpi3))
#从上游分析的featureCounts结果中提取所有样本的基因表达矩阵
#将455个gene_id的mpi_specific_circ_matrix与mRNA表达矩阵进行merge，得到包含gene_id的circ和mrna表达矩阵
#把gene_id转换为gene_name并作为行名，修改mrna表达值的列名，
mpi_all <- read.table('mpi_raw.txt', header = T,quote = '\t',skip = 1)
names(mpi_all)[7:9] <- c('mpi1','mpi2','mpi3')  
all_counts <- as.data.frame(mpi_all[7:9])  
all_counts$gene_id <- mpi_all$Geneid
all_counts <- all_counts[!duplicated(all_counts$gene_id),]
mpi_specific_rna_matrix <- merge(mpi_specific_circ_matrix, all_counts, by = 'gene_id', all = F)
write.csv(mpi_specific_rna_matrix, 'mpi_specific_rna_matrix.csv', row.names = F)
mpi_specific_rna_matrix <- merge(mpi_specific_rna_matrix, mpi_specific_convert, by = 'gene_id', all = F)
rownames(mpi_specific_rna_matrix) <- mpi_specific_rna_matrix$gene_name
mpi_specific_rna_matrix <- mpi_specific_rna_matrix[,-c(1,8)]
colnames(mpi_specific_rna_matrix) <- c('circ_mpi1','circ_mpi2','circ_mpi3',
                                       'linear_mpi1','linear_mpi2','linear_mpi3')
library(dplyr)
mpi_specific_rna_matrix <- mpi_specific_rna_matrix %>% mutate_at(vars(circ_mpi1,circ_mpi2,circ_mpi3,linear_mpi1,linear_mpi2,linear_mpi3), as.numeric)
#编写用于热图绘制的注释文件
annotation_col <- data.frame(RNA = c(rep('circRNA', 3), rep('mRNA', 3)))
row.names(annotation_col) <- colnames(mpi_specific_rna_matrix)
#绘制400个基因的circrna和mrna表达热图
library(pheatmap)
save_pheatmap_pdf <- function(x, filename, width=4, height=6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()}
heatmap_mpi_rna <- pheatmap(mpi_specific_rna_matrix,
                       scale="row",
                       color = colorRampPalette(colors = c("blue","white","red"))(100),
                       annotation_col = annotation_col,
                       annotation_colors = list(RNA=c(circRNA='#fb0e22',mRNA='#00cfff')),
                       #border_color = NA,
                       cluster_rows = T,
                       #treeheight_row = 40,
                       #cluster_cols = col_cluster,
                       #treeheight_col = 20,
                       cellwidth = 40,
                       cellheight = 0.5,
                       #labels_row = T,
                       show_rownames=F,
                       angle_col=45)
save_pheatmap_pdf(heatmap_mpi_rna, "heatmap_mpi_rna.pdf",width=8, height=6)
#将热图拉长，展示每个row的gene_name，以从中手动选择出主要以circrna形式表达的基因列表
#手动建立96个主要以circrna形式表达的基因的gene_name列表输入
#将96个gene_name与包含选择性剪接信息的mpi.cq.all进行merge，获得192个circrna的具体信息，只从中提取gene_id，gene_name和circ_id
#其中gene_id，gene_name有重复，circ_id是唯一的，包含gene的选择性剪接信息，即一条gene对应的多条circrna
heatmap_mpi_rna2 <- pheatmap(mpi_specific_rna_matrix,
                            scale="row",
                            color = colorRampPalette(colors = c("blue","white","red"))(100),
                            annotation_col = annotation_col,
                            annotation_colors = list(RNA=c(circRNA='#fb0e22',mRNA='#00cfff')),
                            #border_color = NA,
                            cluster_rows = T,
                            #treeheight_row = 40,
                            #cluster_cols = col_cluster,
                            #treeheight_col = 20,
                            cellwidth = 10,
                            cellheight = 3,
                            #labels_row = T,
                            show_rownames=T,
                            fontsize_row = 3,
                            angle_col=45)
save_pheatmap_pdf(heatmap_mpi_rna2, "heatmap_mpi_rna2.pdf",width=8, height=30)
mpi_specific_highcircrna_name <- read.csv('mpi_specific_highcircrna_name.csv', header = T) #手动创建该表
mpi_specific_highcircrna_name <- unique(mpi_specific_highcircrna_name)
mpi_specific_highcircrna <- merge(mpi_specific_highcircrna_name, mpi_specific_convert, by = 'gene_name', all = F)
mpi_specific_highcircrna <- merge(mpi.cq.all, mpi_specific_highcircrna, by = 'gene_id', all = F)
mpi_specific_highcircrna <- mpi_specific_highcircrna[!duplicated(mpi_specific_highcircrna$circ_id),]
mpi_specific_highcircrna_convert <- mpi_specific_highcircrna[,c(1,2,16)]
write.csv(mpi_specific_highcircrna_convert, 'mpi_specific_highcircrna_convert.csv', row.names = F)



# 筛选mpi中高表达的circRNA -------------------------------------------------------

#基于CIRIquant分析结果，筛选mpi物种每个样本中bsj reads大于100的circrna作为在mpi中高表达的circrna
#根据高表达circrna的circ_id与CIRIquant表达谱cq_mpi进行merge，得到高表达circrna的表达谱矩阵
#取bsj reads > 100的前50个circrna绘制高表达circrna热图
cq_mpi <- cq_mpi %>% mutate_at(vars(mpi1, mpi2, mpi3), as.numeric)
mpi1_readsover100 <- subset(cq_mpi, mpi1 > 100, select = circ_id)
mpi2_readsover100 <- subset(cq_mpi, mpi2 > 100, select = circ_id)
mpi3_readsover100 <- subset(cq_mpi, mpi3 > 100, select = circ_id)
mpi_readsover100_id <- merge(mpi1_readsover100, mpi2_readsover100, by = 'circ_id', all = F)
mpi_readsover100_id <- merge(mpi_readsover100_id, mpi3_readsover100, by = 'circ_id', all = F)
mpi_readsover100_matrix <- merge(mpi_readsover100_id, cq_mpi, by = 'circ_id', all = F)
rownames(mpi_readsover100_matrix) <- mpi_readsover100_matrix[,1]
mpi_readsover100_matrix <- mpi_readsover100_matrix[,-1]
library(dplyr)
mpi_readsover100_matrix <- arrange(mpi_readsover100_matrix, desc(mpi1))
mpi_readsover100_matrix <- mpi_readsover100_matrix[1:20,]
library(pheatmap)
heatmap_mpi_readsover100 <- pheatmap(mpi_readsover100_matrix,
                            scale="none",
                            color = colorRampPalette(colors = c("blue","white","red"))(100),
                            #annotation_col = annotation_col,
                            #annotation_colors = list(RNA=c(circRNA='#fb0e22',mRNA='#00cfff')),
                            #border_color = NA,
                            cluster_cols = F,
                            cluster_rows = T,
                            #treeheight_row = 40,
                            #cluster_cols = col_cluster,
                            #treeheight_col = 20,
                            cellwidth = 15,
                            cellheight = 5,
                            #labels_row = T,
                            show_rownames=T,
                            fontsize_row = 5,
                            angle_col=45)
save_pheatmap_pdf(heatmap_mpi_readsover100, "heatmap_mpi_readsover100_TOP20.pdf",width=8, height=8)



# 合并高表达、特异性、主要以环状形式表达的三种筛选标准 ----------------------------------------------

#合并三个思路分析到的的circrna(circ_id)
#H：Highly-expressed，mpi的3个样本中均高表达的circrna，即bsj reads > 100的circrna
#S：Specifically-expressed，455个mpi特异性gene对应的549个circrna
#P：Primarily-expressed，455个mpi特异性gene中的149个主要以环状形式表达的gene对应的194个circrna
#通过circ_id的merge，获得29个HS-circrna和其中11个HSP-circrna
library(VennDiagram)
Venn_mpicandidate <- venn.diagram(x = list(mpi_readsover100_id$circ_id, 
                                           mpi_specific_circid$circ_id,
                                           mpi_specific_highcircrna_convert$circ_id),
                                   category.names = c('Highly-\nexpressed', 'Specifically-\nexpressed', 'Primarily-\nexpressed'),
                                   disable.logging = T,
                                   #rotation = 1,
                                   #圆圈
                                   lwd = 1,
                                   col = 'black',  #线条色
                                   fill=c('red','blue','yellow'),  #填充色
                                   alpha = 0.50,  #透明度
                                   cex = 0.8,  #圈内字号
                                   fontfamily = 'serif',  #圈内字体
                                   scaled = F, #圆圈一样大
                                   margin = 0.1,
                                   #标签
                                   cat.pos = c(330,30,150),  #标签位置
                                   cat.dist = c(0.2,0.2,0.2),
                                   cat.cex = 0.8,  #标签字号
                                   cat.col = rep('black',3),  #标签字体颜色
                                   #输出
                                   filename = NULL,
                                   height = 1000,
                                   width = 1000,
                                   resolution = 400,
                                   compression = 'lzw')
pdf('Venn_mpicandidate.pdf')
grid.draw(Venn_mpicandidate)
dev.off()
mpi_HS_id <- merge(mpi_readsover100_id, mpi_specific_circid, by = 'circ_id', all = F)
mpi_HSP_id <- merge(mpi_readsover100_id, mpi_specific_highcircrna_convert, by = 'circ_id', all = F)
#通过重新构建CIRIquant分析的mpi的每个样本的circrna信息表，与HS和HSP的circrna进行merge，获得这些candidate circrna在每个样本的信息
#构建CIRIquant分析的mpi样本信息总表，只保留每个样本的width，strand，circ_type信息
#将mpi_HS_id和mpi_HSP_id的start位置恢复为+1，与总表进行merge
cq_mpi1_merge <- cq_mpi1[,c(4,5,10,11)]
colnames(cq_mpi1_merge) <- c('mpi1_width', 'mpi1_strand', 'circ_id', 'mpi1_circtype')
cq_mpi2_merge <- cq_mpi2[,c(4,5,10,11)]
colnames(cq_mpi2_merge) <- c('mpi2_width', 'mpi2_strand', 'circ_id', 'mpi2_circtype')
cq_mpi3_merge <- cq_mpi3[,c(4,5,10,11)]
colnames(cq_mpi3_merge) <- c('mpi3_width', 'mpi3_strand', 'circ_id', 'mpi3_circtype')
cq_mpi_merge <- merge(cq_mpi1_merge, cq_mpi2_merge, by = 'circ_id', all = F)
cq_mpi_merge <- merge(cq_mpi_merge, cq_mpi3_merge, by = 'circ_id', all = F)
library(dplyr)
mpi_HS_id <- separate(mpi_HS_id, circ_id, into = c('seqnames', 'start', 'end'), sep = ':|\\|')
mpi_HS_id <- mpi_HS_id %>% mutate_at(vars(start,end), as.numeric)
mpi_HS_id <- mpi_HS_id %>% mutate(start = start + 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
mpi_HS_id <- mpi_HS_id[,-c(1:3)]
mpi_candidate_HS <- merge(mpi_HS_id, cq_mpi_merge, by = 'circ_id', all = F)
mpi_candidate_HS <- merge(mpi_specific_convert, mpi_candidate_HS, by = 'gene_id', all = F)
write.csv(mpi_candidate_HS, 'mpi_candidate_HS.csv', row.names = F)
mpi_HSP_id <- separate(mpi_HSP_id, circ_id, into = c('seqnames', 'start', 'end'), sep = ':|\\|')
mpi_HSP_id <- mpi_HSP_id %>% mutate_at(vars(start,end), as.numeric)
mpi_HSP_id <- mpi_HSP_id %>% mutate(start = start + 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
mpi_HSP_id <- mpi_HSP_id[,-c(1:3)]
mpi_candidate_HSP <- merge(mpi_HSP_id, cq_mpi_merge, by = 'circ_id', all = F)
write.csv(mpi_candidate_HSP, 'mpi_candidate_HSP.csv', row.names = F)



# 使用CIRI2对mpi样本circRNA进行分析以获得预测全长序列 ---------------------------------------

#基于使用CIRI2的上游分析结果，找到上述candidate circrna序列
#CIRI2分析得到sample.list文件，需要转换为csv后，按照染色体位置id查找candidate circrna
#将mpi_candidate_HS_id与ciri2_mpi1，ciri2_mpi2，ciri2_mpi3的id列分别merge，得到id对应在各个样本中的Image_ID
#确定后在与样本中对应的.fasta文件中查找序列
#29个candidate circrna对应在CIRI2分析结果中，有125个circrna，有的circ_id相同，但在不同样本或同一样本中存在不同长度的转录本
#手动筛选，每个长度的转录本只保留1条，最终有26个circrna的40个转录本，有3个circrna在CIRI2结果中不存在
ciri2_mpi1 <- read.table(file = 'ciri2_mpi1.list', row.names = NULL, header = T)
write.csv(ciri2_mpi1, 'ciri2_mpi1.csv', row.names = F)
ciri2_mpi1_id <- ciri2_mpi1[,c(1,2,9,10)]
colnames(ciri2_mpi1_id) <- c('mpi1_ID', 'circ_id', 'mpi1_length', 'mpi1_state')
ciri2_mpi2 <- read.table(file = 'ciri2_mpi2.list', row.names = NULL, header = T)
write.csv(ciri2_mpi2, 'ciri2_mpi2.csv', row.names = F)
ciri2_mpi2_id <- ciri2_mpi2[,c(1,2,9,10)]
colnames(ciri2_mpi2_id) <- c('mpi2_ID', 'circ_id', 'mpi2_length', 'mpi2_state')
ciri2_mpi3 <- read.table(file = 'ciri2_mpi3.list', row.names = NULL, header = T)
write.csv(ciri2_mpi3, 'ciri2_mpi3.csv', row.names = F)
ciri2_mpi3_id <- ciri2_mpi3[,c(1,2,9,10)]
colnames(ciri2_mpi3_id) <- c('mpi3_ID', 'circ_id', 'mpi3_length', 'mpi3_state')
mpi_candidate_HS_id <- mpi_candidate_HS[,c(1:6)]
colnames(mpi_candidate_HS_id) <- c('gene_id','gene_name','circ_id','width','strand','type')
candidate_list_ciri2 <- merge(mpi_candidate_HS_id, ciri2_mpi1_id, by = 'circ_id', all = F)
candidate_list_ciri2 <- merge(candidate_list_ciri2, ciri2_mpi2_id, by = 'circ_id', all = F)
candidate_list_ciri2 <- merge(candidate_list_ciri2, ciri2_mpi3_id, by = 'circ_id', all = F)
write.csv(candidate_list_ciri2, 'candidate_list_ciri2.csv', row.names = F)



#绘制出CIRI2方法分析得到的3个样本中共表达的circrna
library(VennDiagram)
Venn_ciri2<- venn.diagram(x = list(ciri2_mpi1_id$circ_id,
                                   mpi_HS_id$circ_id,
                                   ciri2_mpi2_id$circ_id,
                                   ciri2_mpi3_id$circ_id),
                                  category.names = c('mpi1', 'HS_circ', 'mpi2','mpi3'),
                                  disable.logging = T,
                                  #rotation = 1,
                                  #圆圈
                                  lwd = 1,
                                  col = 'black',  #线条色
                                  fill=c('red','blue','yellow','green'),  #填充色
                                  alpha = 0.50,  #透明度
                                  cex = 0.65,  #圈内字号
                                  fontfamily = 'serif',  #圈内字体
                                  scaled = F, #圆圈一样大
                                  margin = 0.1,
                                  #标签
                                  cat.pos = c(340,20,0,0),  #标签位置
                                  cat.dist = c(0.23,0.23,0.12,0.12),
                                  cat.cex = 0.8,  #标签字号
                                  cat.col = rep('black',4),  #标签字体颜色
                                  #输出
                                  filename = NULL,
                                  height = 1000,
                                  width = 1000,
                                  resolution = 400,
                                  compression = 'lzw')
pdf('Venn_ciri2.pdf')
grid.draw(Venn_ciri2)
dev.off()

# ciri2分析得到的mpi三个样本交集 -----------------------------------------------------

# 取Circle_ID和total_exp列生成表达谱
ciri2_mpi1_exp <- ciri2_mpi1[,c('Circle_ID', 'total_exp')]
colnames(ciri2_mpi1_exp) <- c('circ_id', 'mpi1')
ciri2_mpi2_exp <- ciri2_mpi2[,c('Circle_ID', 'total_exp')]
colnames(ciri2_mpi2_exp) <- c('circ_id', 'mpi2')
ciri2_mpi3_exp <- ciri2_mpi3[,c('Circle_ID', 'total_exp')]
colnames(ciri2_mpi3_exp) <- c('circ_id', 'mpi3')
ciri2_mpi_exp <- merge(ciri2_mpi1_exp, ciri2_mpi2_exp, by = 'circ_id', all = F)
ciri2_mpi_exp <- merge(ciri2_mpi_exp, ciri2_mpi3_exp, by = 'circ_id',  all = F)
ciri2_mpi_exp <- unique(ciri2_mpi_exp)
write.csv(ciri2_mpi_exp, 'ciri2_mpi_exp.csv', row.names = F)

# 重新分析：只用MPI样本筛选高表达的候选circRNA ---------------------------------------------
# 第1步：读取4种分析方法获得的MPI的表达谱
ce_mpi1 <- read.csv('ce_mpi1.csv', header = T, row.names = 1)
ce_mpi2 <- read.csv('ce_mpi2.csv', header = T, row.names = 1)
ce_mpi3 <- read.csv('ce_mpi3.csv', header = T, row.names = 1)
ciri2_mpi1 <- read.csv('ciri2_mpi1.csv', header = T)
ciri2_mpi2 <- read.csv('ciri2_mpi2.csv', header = T)
ciri2_mpi3 <- read.csv('ciri2_mpi3.csv', header = T)
cq_mpi1 <- read.csv('cq_mpi1.csv', header = T, row.names = 1)
cq_mpi2 <- read.csv('cq_mpi2.csv', header = T, row.names = 1)
cq_mpi3 <- read.csv('cq_mpi3.csv', header = T, row.names = 1)
fc_mpi1 <- read.csv('fc_mpi1.csv', header = T, row.names = 1)
fc_mpi2 <- read.csv('fc_mpi2.csv', header = T, row.names = 1)
fc_mpi3 <- read.csv('fc_mpi3.csv', header = T, row.names = 1)

# 第2步：绘制Venn图交集
# 对齐不同分析方法的染色质location
library(dplyr)
ce_mpi1 <- ce_mpi1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mpi2 <- ce_mpi2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mpi3 <- ce_mpi3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ciri2_mpi1 <- ciri2_mpi1 %>% mutate(start = start - 1, circ_id=paste(Chr, ':', start, '|', end, sep = ''))
ciri2_mpi2 <- ciri2_mpi2 %>% mutate(start = start - 1, circ_id=paste(Chr, ':', start, '|', end, sep = ''))
ciri2_mpi3 <- ciri2_mpi3 %>% mutate(start = start - 1, circ_id=paste(Chr, ':', start, '|', end, sep = ''))
cq_mpi1 <- cq_mpi1 %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mpi2 <- cq_mpi2 %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mpi3 <- cq_mpi3 %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
fc_mpi1 <- fc_mpi1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mpi2 <- fc_mpi2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mpi3 <- fc_mpi3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))

# 生物学重复取交集
ce_mpi <- merge(ce_mpi1, ce_mpi2, by = 'circ_id', all = F)
ce_mpi <- merge(ce_mpi, ce_mpi3, by = 'circ_id', all = F)
ciri2_mpi <- merge(ciri2_mpi1, ciri2_mpi2, by = 'circ_id', all = F)
ciri2_mpi <- merge(ciri2_mpi, ciri2_mpi3, by = 'circ_id', all = F)
cq_mpi <- merge(cq_mpi1, cq_mpi2, by = 'circ_id', all = F)
cq_mpi <- merge(cq_mpi, cq_mpi3, by = 'circ_id', all = F)
fc_mpi <- merge(fc_mpi1, fc_mpi2, by = 'circ_id', all = F)
fc_mpi <- merge(fc_mpi, fc_mpi3, by = 'circ_id', all = F)

# 绘制Venn图
library(VennDiagram)
Venn_mpi <- venn.diagram(x = list(ce_mpi$circ_id, ciri2_mpi$circ_id, cq_mpi$circ_id, fc_mpi$circ_id),
                         category.names = c('CIRCexplorer3', 'CIRI2','CIRIquant', 'find_circ'),
                         disable.logging = T,
                         lwd = 1,
                         col = 'black',
                         fill=c('red','blue','yellow','green'),
                         alpha = 0.50,
                         cex = 1,
                         fontfamily = 'serif',
                         scaled = F,
                         margin = 0.1,
                         cat.pos = c(330,30,0,0),
                         cat.dist = c(0.22,0.22,0.1,0.1),
                         cat.cex = 1,
                         cat.col = rep('black',4),
                         filename = NULL,
                         height = 1600,
                         width = 1600,
                         resolution = 400,
                         compression = 'lzw')
pdf('Venn_mpi.pdf')
grid.draw(Venn_mpi)
dev.off()

# 获得4种方法交集的circ_id
mpi_exp_id <- merge(ce_mpi, ciri2_mpi, by = 'circ_id', all = F)
mpi_exp_id <- data.frame(circ_id = mpi_exp_id$circ_id)
mpi_exp_id <- merge(mpi_exp_id, cq_mpi, by = 'circ_id', all = F)
mpi_exp_id <- data.frame(circ_id = mpi_exp_id$circ_id)
mpi_exp_id <- merge(mpi_exp_id, fc_mpi, by = 'circ_id', all = F)
mpi_exp_id <- data.frame(circ_id = mpi_exp_id$circ_id)
mpi_exp_id <- unique(mpi_exp_id)

# 第3步：

ce_mpi1_exp <- merge(mpi_exp_id, ce_mpi1, by = 'circ_id', all = F)
ce_mpi1_exp <- data.frame(circ_id = ce_mpi1_exp$circ_id,
                          mpi1 = ce_mpi1_exp$FPBcirc)

# 重新分析：在高表达中找共表达的（Gemini3） ------------------------------------------------
# 第1步：读取4种分析方法获得的MPI的表达谱
ce_mpi1 <- read.csv('ce_mpi1.csv', header = T, row.names = 1)
ce_mpi2 <- read.csv('ce_mpi2.csv', header = T, row.names = 1)
ce_mpi3 <- read.csv('ce_mpi3.csv', header = T, row.names = 1)
ciri2_mpi1 <- read.csv('ciri2_mpi1.csv', header = T)
ciri2_mpi2 <- read.csv('ciri2_mpi2.csv', header = T)
ciri2_mpi3 <- read.csv('ciri2_mpi3.csv', header = T)
cq_mpi1 <- read.csv('cq_mpi1.csv', header = T, row.names = 1)
cq_mpi2 <- read.csv('cq_mpi2.csv', header = T, row.names = 1)
cq_mpi3 <- read.csv('cq_mpi3.csv', header = T, row.names = 1)
fc_mpi1 <- read.csv('fc_mpi1.csv', header = T, row.names = 1)
fc_mpi2 <- read.csv('fc_mpi2.csv', header = T, row.names = 1)
fc_mpi3 <- read.csv('fc_mpi3.csv', header = T, row.names = 1)

# 对齐不同分析方法的染色质location
library(dplyr)
ce_mpi1 <- ce_mpi1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mpi2 <- ce_mpi2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mpi3 <- ce_mpi3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ciri2_mpi1 <- ciri2_mpi1 %>% mutate(start = start - 1, circ_id=paste(Chr, ':', start, '|', end, sep = ''))
ciri2_mpi2 <- ciri2_mpi2 %>% mutate(start = start - 1, circ_id=paste(Chr, ':', start, '|', end, sep = ''))
ciri2_mpi3 <- ciri2_mpi3 %>% mutate(start = start - 1, circ_id=paste(Chr, ':', start, '|', end, sep = ''))
cq_mpi1 <- cq_mpi1 %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mpi2 <- cq_mpi2 %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mpi3 <- cq_mpi3 %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
fc_mpi1 <- fc_mpi1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mpi2 <- fc_mpi2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mpi3 <- fc_mpi3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))

# 数据整理
# 定义统一的清洗函数：提取 circ_id 和指定的表达量列，并按表达量降序排序
format_circ_data <- function(df, expr_col) {
  data.frame(circ_id = df$circ_id, expr = as.numeric(df[[expr_col]])) %>%
    filter(!is.na(expr)) %>%
    arrange(desc(expr))
}

# 将 12 个数据框整理成结构化的列表，并提取对应的定量列
data_list <- list(
  # mpi1 = list(
  #   ce = format_circ_data(ce_mpi1, "FPBcirc"),
  #   ciri2 = format_circ_data(ciri2_mpi1, "total_exp"),
  #   cq = format_circ_data(cq_mpi1, "score"),
  #   fc = format_circ_data(fc_mpi1, "n_reads")
  # ),
  mpi2 = list(
    ce = format_circ_data(ce_mpi2, "FPBcirc"),
  #  ciri2 = format_circ_data(ciri2_mpi2, "total_exp"),
    cq = format_circ_data(cq_mpi2, "score"),
    fc = format_circ_data(fc_mpi2, "n_reads")
  ),
  mpi3 = list(
    ce = format_circ_data(ce_mpi3, "FPBcirc"),
  #  ciri2 = format_circ_data(ciri2_mpi3, "total_exp"),
    cq = format_circ_data(cq_mpi3, "score"),
    fc = format_circ_data(fc_mpi3, "n_reads")
  )
)

# 第2步：在高表达中找共表达的
# 设定每种方法提取的 TopN 数量（根据实际结果灵活调整，如果最后结果为0，请调大至 300 或 500）
TopN <- 300 

sample_intersections <- list()

# A. 样本内取 4 种方法的交集
for (s_name in names(data_list)) {
  # 提取该样本下 4 种方法的 TopN circ_id
  method_top_list <- lapply(data_list[[s_name]], function(df) head(df$circ_id, TopN))
  
  # 取交集
  sample_intersections[[s_name]] <- Reduce(intersect, method_top_list)
  
  cat(sprintf("样本 %s 中 4 种方法共有的 Top%d circRNA 数量为: %d\n", 
              s_name, TopN, length(sample_intersections[[s_name]])))
}

# B. 在 3 个生物学重复之间取交集，找到最稳健的高表达 circRNA
final_candidates <- Reduce(intersect, sample_intersections)

cat("\n最终在 3 个样本、4 种方法中均高表达的 circRNA 数量为：", length(final_candidates), "\n")

if (length(final_candidates) > 0) {
  # 提取 3 个样本的 CIRIquant score 作为最终定量参考
  final_table <- data.frame(circ_id = final_candidates)
  for (s_name in names(data_list)) {
    tmp <- data_list[[s_name]]$cq %>% filter(circ_id %in% final_candidates)
    final_table[[paste0(s_name, "_CQ_score")]] <- tmp$expr[match(final_candidates, tmp$circ_id)]
  }
  write.csv(final_table, "Final_Consensus_High_Expression.csv", row.names = FALSE)
  cat("已保存最终稳健结果表。\n")
}




library(tidyverse)
library(UpSetR)
library(ComplexHeatmap)
library(RColorBrewer)
# A. 保存 UpSet Plot (以 mpi2 为代表)
pdf("Upset_Plot_mpi2.pdf", width = 7, height = 7, onefile = F)
upset_input <- list(
  CIRCexplorer3 = head(data_list$mpi2$ce$circ_id, TopN),
 # CIRI2 = head(data_list$mpi2$ciri2$circ_id, TopN),
  CIRIquant = head(data_list$mpi2$cq$circ_id, TopN),
  find_circ = head(data_list$mpi2$fc$circ_id, TopN)
)
upset(fromList(upset_input), order.by = "freq", main.bar.color = "#2c3e50",
      # --- 字号设置 (比例 = 目标字号 / 12) ---
      # text.scale 向量的 6 个值依次代表：
      # c(柱状图标题, 柱状图刻度, 条形图标题, 条形图刻度, 方法名称, 柱子上方数字)
       text.scale = c(2.75, 2.8, 2.75, 2.8, 2.8, 2.8), 
      
      # --- 点图（Matrix）修饰 ---
      point.size = 4,        # 点的大小 (默认是 3，这里加大以协调)
      line.size = 1.2,       # 连线的粗细 (默认是 0.7)
      
      # --- 柱子数字设置 ---
      show.numbers = "yes",  # 确保显示上方数字
      
      # --- 布局协调 ---
      # mb.ratio 设置上方柱状图和下方点图的比例 (0.6 代表上方占 60%)
      mb.ratio = c(0.6, 0.4))
dev.off()


pdf("Upset_Plot_mpi3.pdf", width = 7, height = 7, onefile = F)
upset_input <- list(
  CIRCexplorer3 = head(data_list$mpi3$ce$circ_id, TopN),
  # CIRI2 = head(data_list$mpi2$ciri2$circ_id, TopN),
  CIRIquant = head(data_list$mpi3$cq$circ_id, TopN),
  find_circ = head(data_list$mpi3$fc$circ_id, TopN)
)
upset(fromList(upset_input), order.by = "freq", main.bar.color = "#2c3e50",
      # --- 字号设置 (比例 = 目标字号 / 12) ---
      # text.scale 向量的 6 个值依次代表：
      # c(柱状图标题, 柱状图刻度, 条形图标题, 条形图刻度, 方法名称, 柱子上方数字)
      text.scale = c(2.75, 2.8, 2.75, 2.8, 2.8, 2.8), 
      
      # --- 点图（Matrix）修饰 ---
      point.size = 4,        # 点的大小 (默认是 3，这里加大以协调)
      line.size = 1.2,       # 连线的粗细 (默认是 0.7)
      
      # --- 柱子数字设置 ---
      show.numbers = "yes",  # 确保显示上方数字
      
      # --- 布局协调 ---
      # mb.ratio 设置上方柱状图和下方点图的比例 (0.6 代表上方占 60%)
      mb.ratio = c(0.6, 0.4))
dev.off()



# 绘制Venn图
library(ggVennDiagram)
# 提取 mpi2 中 4 种方法共有的 TopN
mpi2_top_list <- lapply(data_list$mpi2, function(df) head(df$circ_id, TopN))
mpi2_consensus <- Reduce(intersect, mpi2_top_list)

# 提取 mpi3 中 4 种方法共有的 TopN (顺便做了，备用)
mpi3_top_list <- lapply(data_list$mpi3, function(df) head(df$circ_id, TopN))
mpi3_consensus <- Reduce(intersect, mpi3_top_list)

# 这里我们对比 mpi1 和 mpi2
venn_data <- list(
  #Sample_mpi1 = mpi1_consensus,
  mpi2 = mpi2_consensus,
  mpi3 = mpi3_consensus
)

p_venn <- ggVennDiagram(venn_data, 
                        label_alpha = 0, 
                        edge_size = 0.1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#4292C6") +
  labs(title = paste0("Consistency of Top ", TopN, " circRNAs (3-Tool Intersection)"),
       subtitle = "Comparison between mpi2 and mpi3",
       caption = "Calculated by CIRCexplorer3, CIRIquant, find_circ") +
  theme(legend.position = "none")

print(p_venn)
# 保存为 PDF
ggsave("Venn_mpi2_vs_mpi3.pdf", plot = p_venn, width = 9, height = 6)









# B. 保存 Heatmap
if (length(final_candidates) > 1) {
  mat_log <- log2(as.matrix(final_table[,-1]) + 0.001)
  rownames(mat_log) <- final_table$circ_id
  
  pdf("Heatmap_Final_Candidates.pdf", width = 8, height = 10)
  p <- Heatmap(mat_log, 
               name = "log2(CQ_Score)", 
               column_title = "Final Consensus circRNAs",
               col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
               cluster_columns = FALSE)
  draw(p)
  dev.off()
}






if (length(final_candidates) > 1) {
  
  # 准备矩阵
  mat_log <- log2(as.matrix(final_table[,-1]) + 0.001)
  rownames(mat_log) <- final_table$circ_id
  
  # 定义你需要标注的特定 circRNA ID
  labels_to_show <- c("scaffold0007:91656590|91656785", "scaffold0031:2255403|2256630")
  
  # 找到这些 ID 在矩阵中的行索引
  row_indices <- which(rownames(mat_log) %in% labels_to_show)
  
  # 检查是否找到了 ID
  if(length(row_indices) == 0) {
    warning("在最终列表中未找到指定的 ID，请检查坐标格式是否完全一致。")
  }
  
  # 创建行名标注对象 (anno_mark)
  row_anno <- rowAnnotation(
    foo = anno_mark(
      at = row_indices, 
      labels = rownames(mat_log)[row_indices],
      labels_gp = gpar(fontsize = 10, fontface = "bold"),
      padding = unit(1, "mm")
    )
  )
  
  # 保存为 PDF
  pdf("Heatmap_Marked_Candidates.pdf", width = 7, height = 10)
  
  p <- Heatmap(mat_log, 
               name = "log2(CQ_Score)", 
               column_title = "Bat circRNA Expression (Marked)",
               col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
               cluster_columns = FALSE, 
               show_row_names = FALSE,      # 隐藏所有原始行名
               right_annotation = row_anno,  # 添加指引线标注
               border = TRUE)
  
  draw(p)
  dev.off()
  
  cat("热图已保存，已标注指定 ID 的指引线。\n")
}
# 再次重新分析：合并样本后找3种方法各自的高表达然后merge ------------------------------------------
ce_mpi1 <- read.csv('ce_mpi1.csv', header = T, row.names = 1)
ce_mpi2 <- read.csv('ce_mpi2.csv', header = T, row.names = 1)
ce_mpi3 <- read.csv('ce_mpi3.csv', header = T, row.names = 1)
ciri2_mpi1 <- read.csv('ciri2_mpi1.csv', header = T)
ciri2_mpi2 <- read.csv('ciri2_mpi2.csv', header = T)
ciri2_mpi3 <- read.csv('ciri2_mpi3.csv', header = T)
cq_mpi1 <- read.csv('cq_mpi1.csv', header = T, row.names = 1)
cq_mpi2 <- read.csv('cq_mpi2.csv', header = T, row.names = 1)
cq_mpi3 <- read.csv('cq_mpi3.csv', header = T, row.names = 1)
fc_mpi1 <- read.csv('fc_mpi1.csv', header = T, row.names = 1)
fc_mpi2 <- read.csv('fc_mpi2.csv', header = T, row.names = 1)
fc_mpi3 <- read.csv('fc_mpi3.csv', header = T, row.names = 1)

# 对齐不同分析方法的染色质location
library(dplyr)
ce_mpi1 <- ce_mpi1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mpi2 <- ce_mpi2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ce_mpi3 <- ce_mpi3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
ciri2_mpi1 <- ciri2_mpi1 %>% mutate(start = start - 1, circ_id=paste(Chr, ':', start, '|', end, sep = ''))
ciri2_mpi2 <- ciri2_mpi2 %>% mutate(start = start - 1, circ_id=paste(Chr, ':', start, '|', end, sep = ''))
ciri2_mpi3 <- ciri2_mpi3 %>% mutate(start = start - 1, circ_id=paste(Chr, ':', start, '|', end, sep = ''))
cq_mpi1 <- cq_mpi1 %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mpi2 <- cq_mpi2 %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
cq_mpi3 <- cq_mpi3 %>% mutate(start = start - 1, circ_id=paste(seqnames, ':', start, '|', end, sep = ''))
fc_mpi1 <- fc_mpi1 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mpi2 <- fc_mpi2 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))
fc_mpi3 <- fc_mpi3 %>% mutate(circ_id=paste(chrom, ':', start, '|', end, sep = ''))

fc_mpi_norm <- merge(fc_mpi1, fc_mpi2, by = 'circ_id', all = F)
fc_mpi_norm <- merge(fc_mpi_norm, fc_mpi3, by = 'circ_id', all = F)
fc_mpi_norm <- data.frame(circ_id = fc_mpi_norm$circ_id,
                          mpi1_norm = fc_mpi_norm$n_reads.x,
                          mpi2_norm = fc_mpi_norm$n_reads.y,
                          mpi3_norm = fc_mpi_norm$n_reads)
rownames(fc_mpi_norm) <- fc_mpi_norm[,1]
fc_mpi_norm <- fc_mpi_norm[,-1]
library(edgeR)
fc_mpi_norm <- cpm(fc_mpi_norm)
write.csv(fc_mpi_norm, 'fc_mpi_norm.csv', row.names = T)

# 选取TOP300/样本/方法，ce用FPBcirc，cq用score，fc用n_reads
# TopN <- 300
# 
# top300_ce_mpi2 <- ce_mpi2 %>%
#   arrange(desc(FPBcirc)) %>%
#   filter(FPBcirc > 19) %>%
#   select(circ_id, FPBcirc)
# 
# top300_ce_mpi3 <- ce_mpi3 %>%
#   arrange(desc(FPBcirc)) %>%
#   slice_head(n = TopN) %>%
#   select(circ_id, FPBcirc)
# 
# top300_cq_mpi2 <- cq_mpi2 %>%
#   arrange(desc(score)) %>%
#   slice_head(n = TopN) %>%
#   select(circ_id, score)
# 
# top300_cq_mpi3 <- cq_mpi3 %>%
#   arrange(desc(score)) %>%
#   slice_head(n = TopN) %>%
#   select(circ_id, score)
# 
# top300_fc_mpi2 <- fc_mpi2 %>%
#   arrange(desc(n_reads)) %>%
#   slice_head(n = TopN) %>%
#   select(circ_id, n_reads)
# 
# top300_fc_mpi3 <- fc_mpi3 %>%
#   arrange(desc(n_reads)) %>%
#   slice_head(n = TopN) %>%
#   select(circ_id, n_reads)

# --- 处理 CIRCexplorer3 (ce) ---
# 定义样本列表和阈值
ce_list <- list(ce_mpi1, ce_mpi2, ce_mpi3)
ce_ids_list <- lapply(ce_list, function(df) {
  df %>% filter(FPBcirc > 19) %>% pull(circ_id)
})
# 取 3 个样本的交集
ce_mpi19 <- Reduce(intersect, ce_ids_list)

# --- 处理 CIRIquant (cq) ---
cq_list <- list(cq_mpi1, cq_mpi2, cq_mpi3)
cq_ids_list <- lapply(cq_list, function(df) {
  df %>% filter(score > 1.3) %>% pull(circ_id)
})
# 取 3 个样本的交集
cq_mpi1.3 <- Reduce(intersect, cq_ids_list)

# --- 处理 find_circ (fc) ---
fc_list <- list(fc_mpi1, fc_mpi2, fc_mpi3)
fc_ids_list <- lapply(fc_list, function(df) {
  df %>% filter(n_reads > 49) %>% pull(circ_id)
})
# 取 3 个样本的交集
fc_mpi49 <- Reduce(intersect, fc_ids_list)

library(VennDiagram)
Venn_mpi_filter <- venn.diagram(x = list(ce_mpi19, cq_mpi1.3, fc_mpi49),
                         category.names = c('CIRCexplorer3\nFPBcirc>19',
                                            'CIRIquant\nscore>1.3',
                                            'find_circ\nn_reads>49'),
                         disable.logging = T,
                         lwd = 1, col = 'black',
                         fill=c('red','blue','yellow'),
                         alpha = 0.50,
                         cex = 1,
                         fontfamily = 'serif',
                         scaled = F,
                         margin = 0.1,
                         cat.pos = c(330,30,180),
                         cat.dist = c(0.07,0.07,0.06),
                         cat.cex = 1,
                         cat.col = rep('black',3),
                         filename = NULL,
                         height = 1600,
                         width = 1600,
                         resolution = 400,
                         compression = 'lzw')
pdf('Venn_mpi_filter.pdf')
grid.draw(Venn_mpi_filter)
dev.off()

venn_input <- list(
  CIRCexplorer3 = ce_mpi19,
  CIRIquant     = cq_mpi1.3,
  find_circ     = fc_mpi49
)
final_consensus_ids <- Reduce(intersect, venn_input)
write.csv(data.frame(circ_id = final_consensus_ids), 
          "Final_Consensus_IDs_from_Venn.csv", row.names = FALSE)

# 绘制热图
mat_data <- data.frame(row.names = final_consensus_ids)
mat_data$mpi1 <- cq_mpi1$score[match(final_consensus_ids, cq_mpi1$circ_id)]
mat_data$mpi2 <- cq_mpi2$score[match(final_consensus_ids, cq_mpi2$circ_id)]
mat_data$mpi3 <- cq_mpi3$score[match(final_consensus_ids, cq_mpi3$circ_id)]

# 转换为矩阵并进行 log2 转换 (加 0.01 防止 log0)
# 转置矩阵：使行成为样本 (mpi1, 2, 3)，列成为 circ_id
mat_final <- t(as.matrix(mat_data))
mat_log <- log2(mat_final + 0.01)

# 定义需要标注的 4 个 ID
special_ids <- c(
  "scaffold0007:91656590|91656785", 
  "scaffold0029:6092716|6092982", 
  "scaffold0002:7550372|7550683", 
  "scaffold0031:2255403|2256630"
)

# 找到这些 ID 在矩阵列中的位置索引
# 注意：此时 circ_id 是列名 (colnames)
col_indices <- which(colnames(mat_log) %in% special_ids)

# 创建列标注对象
col_anno <- columnAnnotation(
  mark = anno_mark(
    at = col_indices, 
    labels = colnames(mat_log)[col_indices],
    labels_gp = gpar(fontsize = 7, fontface = "italic"), # 蝙蝠物种名通常斜体，这里ID也可用斜体区分
    labels_rot = 0,
    side = "bottom" # 标注放在下方
  )
)

pdf("Transposed_Heatmap_Consensus.pdf", width = 12, height = 4) # 横向图建议宽度设大

p <- Heatmap(mat_log, 
             name = "log2(Score)", 
             column_title = "Robust High-Expression circRNAs (Transposed)",
             
             # 颜色设置：蓝-黄-红
             col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
             
             # 聚类设置
             cluster_rows = FALSE,    # 样本顺序固定为 mpi1, 2, 3
             cluster_columns = TRUE,  # 对 circRNA 进行聚类
             
             # 名字显示设置
             show_row_names = TRUE,   # 显示 mpi1, 2, 3
             show_column_names = FALSE, # 隐藏下方密集的几百个 ID
             
             # 关键：添加顶部的指示线标注
             top_annotation = col_anno,
             
             # 视觉美化
             border = F,
             #rect_gp = gpar(col = "white", lwd = 0.5) # 给格子加极细的白边，增加质感
)

draw(p)
dev.off()

cat("横向热图已生成，已标注 4 个特定 circRNA。\n")
