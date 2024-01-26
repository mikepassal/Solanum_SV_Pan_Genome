parse_counts <- function(args){    
    # check stranded or unstranded
    # count total read counts on genes in the 3 and 4 columns 
    # (these column represent different library strandedness). 
    # For stranded data one of the columns should be much larger than the other.
    library(stringr)
    SRA = args
    path1 = paste('/data/passala/Collaborator_Data/Marchantia_data/Full_marchan_network/', SRA, '/', sep = '')
    files1 = list.files(path1, pattern = "\\ReadsPerGene.out.tab$")  # all sub-folders
    filepath1 = paste(path1, files1[1], sep = '')
    
    strand_mat = read.delim(filepath1, sep = '\t', header = T)
    strand_mat = strand_mat[-c(1:3),]
    a1 = sum(strand_mat[,3])
    a2 = sum(strand_mat[,4])
    stranded = 0
    if ((a1>10*a2) | (a2>10*a1)){
      stranded = 1
    } 
    print(paste(SRA, ':  ', a1, '  ', a2, '  ', 'stranded: ', stranded, sep = ''))
    
    
    genecounts = "ReadsPerGene.out.tab"
    splicejunctions = "SJ.out.tab"
    logname = "Log.final.out"
    
    files = unlist(read.table(paste0(SRA, "_accList_1.txt"),stringsAsFactors=F))
    ind_rm = NULL
    Ns = list()
    i = 1
    for (fil in files){
      filedir = paste0(SRA, "/", fil)
      countfile = paste0(filedir, genecounts)
      logfile = paste0(filedir, logname)
      N = list()
      if( file.exists(countfile) ) {
        #             print(countfile)
        
        counts  =read.table(countfile)
        log1    =read.table(logfile, sep="\t", nrows=6)
        log2    =read.table(logfile, sep="\t", skip=8, nrows=14)
        log3    =read.table(logfile, sep="\t", skip=23, nrows=4)
        log4    =read.table(logfile, sep="\t", skip=28, nrows=3)
        
        N$mapinfo      = rbind(log1,log2,log3,log4)
        N$unmapped     = counts[1,]
        N$multimapping = counts[2,]
        N$noFeature    = counts[3,]
        N$ambiguous    = counts[4,]
        N$length       = dim(counts)[1]-4
        N$genes        = counts[(1:N$length)+4,1]
        N$counts1      = counts[(1:N$length)+4,2]
        N$counts2      = counts[(1:N$length)+4,3]
        N$counts3      = counts[(1:N$length)+4,4]
      } else {        
        print("uh oh..." )
        print(fil)
        ind_rm = c(ind_rm, fil)
      }
      
      # Stranded or unstranded? Spot check this before running the code. 
      # Should use the last column if stranded (or max of 2nd/3rd). 
      # Otherwise first column (unstranded). 
      # N$counts3 = rep(0, length(attr$ensemblID ) )
      if (stranded == 1){
        if (sum(N$counts2) > sum(N$counts3)){
          counts_this = N$counts2
        } else {
          counts_this = N$counts3
        }
      } else {
        counts_this = N$counts1
      }
      
      
      # remove 0 counts files
      if (sum(counts_this)>0){
        if( i > 1  ){
          counts_exp = cbind(counts_exp, counts_this)
        } else {
          counts_exp = counts_this
        }
      } else {
        print(paste('0 counts:  ', fil, ' ', sep = ''))
        ind_rm = c(ind_rm, fil)
      }
      
      Ns[[i]] = N
      print(paste('sample ', i, '  total counts: ', sum(counts_this), sep = ''))
      i = i + 1
      
    }
    
    rownames(counts_exp) = N$genes
    colnames(counts_exp) = setdiff(files, ind_rm)
    write.table(counts_exp, file = paste0('metaExpData',SRA,'.csv'), sep = ',',col.names = TRUE, row.names = TRUE, quote = FALSE)
    print(paste('Saved! ', dim(counts_exp), sep = ''))
}

args = commandArgs(trailingOnly=TRUE)
parse_counts(args)