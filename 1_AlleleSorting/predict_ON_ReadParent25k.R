# This will load results of counts of specific Mother and Father kmers for a batch and compute an indicative pvalue
# It is used to predict if a read is of mother or father origin (or if uncertain should be used for both assemblies)
# since the various types of reads (normal, microsatellite, etc...) will have different distributions of unique kmers,
# we adjust the pvalue independently for each kind of read.
#
# (c) N.Guex 2024-2025
# Bioinformatics Competence Center (BICC)
# University of Lausanne
# Switzerland
# https://bix.unil.ch
#
# R --vanilla --quiet --args /home/nguex/GOUDET/allele_sort/jf49o/ON_25k/  /data2/chris/tyto/ON_25k_table.txt  ./predict_ON_ReadParent25k_rslts < predict_ON_ReadParent25k.R
#

options(width=280)
args = commandArgs(TRUE)
indir       <- args[1] # indir <- '/home/nguex/GOUDET/allele_sort/jf49o/ON_25k/'
infosfile   <- args[2] # infosfile <- '/data2/chris/tyto/ON_25k_table.txt'
outfn       <- args[3] # outfn <- './predict_ON_ReadParent25k_rslts'
dir <- './'
cutoff <- 0.05

infos <- read.delim(file=infosfile,sep="\t",header=FALSE)
colnames(infos) <- c('id','id_ori','kind','len')

#######

alldata <- NULL
fl <- list.files(indir,pattern="x..$")
for (fn in fl)
{
    # read results block and add header
    data <- read.delim(file=paste0(indir,fn,'.txt'),sep=",",header=FALSE)
    hdr <- read.delim(file=paste0(dir,'predictReadParentHeaders','.txt'),sep="\t",header=FALSE)
    hdr <- hdr[1:18,]  # WARNING: redundant, consider only unique
    colnames(hdr) <- c('test','FatherSet','MotherSet')
    colnames(data) <- c('id',paste(rep(hdr[,1],each=4),c('Fcnt','Mcnt','Fscore','Mscore'),sep='_'))
    rm(hdr)
    
    data$id <- substr(data$id,2,999) # remove leading #
    
    # build result table and predict parent
    FcntIdx <- grep('Fcnt',colnames(data))
    McntIdx <- grep('Mcnt',colnames(data))
    FcntMean <- apply(data[,FcntIdx],1,mean)
    McntMean <- apply(data[,McntIdx],1,mean)
    data <- cbind(data,FcntMean=NA,McntMean=NA,highest_parent_mean=NA,parent=NA,pval=NA,kind_padj=NA)
    data[,'FcntMean'] <- FcntMean
    data[,'McntMean'] <- McntMean
    idx <- which(FcntMean > McntMean)
    data[ idx,'highest_parent_mean'] <- 'Father'
    data[-idx,'highest_parent_mean'] <- 'Mother'
    idx <- which(FcntMean == McntMean)
    data[idx,'highest_parent_mean'] <- NA
    
    # estimate strength of prediction
    data[,'pval'] <- unlist(lapply(1:nrow(data),function(i) {  pval <- try(t.test(as.integer(data[i,FcntIdx]),as.integer(data[i,McntIdx]),paired=TRUE)$p.val ,silent=TRUE); if (is.finite(pval)) {return(pval)} else return(NA)     }))
    alldata <- rbind(alldata,data)
    rm(data)
}

alldata <- merge(alldata,infos,by='id',sort=FALSE)

kind_list <- unique(alldata$kind) # c('other','musat','centro','kineto','telo')
# adjust pvalues for each kind
for (kind in kind_list)
{
    idx <- which(!is.na(alldata[,'pval']) & alldata[,'kind'] == kind)
    alldata[idx,'kind_padj']  <- p.adjust(alldata[idx,'pval'], method = 'BH')
}

# make prediction at padj=0.05
alldata[,'parent'] <- alldata[,'highest_parent_mean']
idx <- which(alldata$kind_padj > cutoff | is.na(alldata$kind_padj))
alldata$parent[idx] <- 'unknown'


# order columns and save results
rslt <- alldata[,c(1,82,81,77,79,80,74,75,76,78,2:73)]
save(rslt,file=paste0(outfn,'.Rdata'))

write.table(rslt[,1:10],file=paste0(outfn,'.txt'),row.names=FALSE,quote=FALSE,sep="\t")


tbl <- table(rslt$kind,rslt$parent)
rs <- rowSums(tbl)
tbl.norm <- 100.0 * sweep(tbl,1,rs,FUN="/")
print(paste("predictions raw counts with padj=",cutoff))
tbl
print(paste("predictions percentages with padj=",cutoff))
tbl.norm

#################


