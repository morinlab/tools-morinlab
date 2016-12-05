require(maftools);
library(argparse);
require(data.table);

###

parser <- ArgumentParser(description="Create a Gene Lollipop using Maftools");

parser$add_argument(
    "--input_maf", "-maf",
    required="True",
    help="Input Variants in MAF format"
    );

parser$add_argument(
    "--gene_blacklist", "-gl",
    help="Input gene list with separated by newline"
    );

parser$add_argument(
   "--min_mut", "-mm",
      default=5,
         help="Minimum number of mutations seen in the gene for it to be included in the calculation");

parser$add_argument(
   "--fdr", "-f",
         default=0.1,
	          help="FDR threshold to use in plots and returned gene list");

parser$add_argument(
   "--aacol", "-ac",
   help="Optionally provide the name of the column that contains the amino acid annotation in your MAF file");

parser$add_argument(
   "--output_detail", "-o",
   required="True",
   help="Output text file for oncodriveclust detail"
   )

parser$add_argument(
   "--output_plot", "-p",
      required="True",
         help="Output pdf file for oncodriveclust detail"
	    )

args <- parser$parse_args();

###


aacol = 'HGVSp_Short'
if(!is.null(args$aacol)){
aacol = args$aacol
}

min_mut = as.integer(args$min_mut)


#--------------------- based on binaomial distribution, estimate threshhold.
get_threshold = function(gene_muts, gene_length){
  th = which(unlist(lapply(X = 2:gene_muts, FUN = function(x) dbinom(x = x, size = gene_muts, prob = 1/gene_length) )) < 0.01)[1]
    return(th+1)
    }
    #-------------------- end of function.

parse_prot_fix = function(dat, AACol, gl, m, calBg = FALSE, nBg){

  if(is.null(AACol)){
      if(! 'AAChange' %in% colnames(dat)){
            message('Available fields:')
	          print(colnames(dat))
		        stop('AAChange field not found in MAF. Use argument AACol to manually specifiy field name containing protein changes.')
			    }
			      }else{
			          colnames(dat)[which(colnames(dat) == AACol)] = 'AAChange'
				    }

  all.prot.dat = dat[,.(Hugo_Symbol, Variant_Classification, AAChange)]
    all.prot.dat = all.prot.dat[Variant_Classification != 'Splice_Site']
      #parse AAchanges to get postion
        prot.spl = strsplit(x = as.character(all.prot.dat$AAChange), split = '.', fixed = TRUE)
	  prot.conv = sapply(prot.spl, function(x) x[length(x)])

  all.prot.dat[,conv := prot.conv]
    all.prot.dat = all.prot.dat[!conv == 'NULL']

  #If conversions are in HGVSp_long (default HGVSp) format, we will remove strings Ter followed by anything (e.g; p.Asn1986GlnfsTer13)
    pos = gsub(pattern = 'Ter.*', replacement = '',x = all.prot.dat$conv)

  #Following parsing takes care of most of HGVSp_short and HGVSp_long format
    pos = gsub(pattern = '[[:alpha:]]', replacement = '', x = pos)
      pos = gsub(pattern = '\\*$', replacement = '', x = pos) #Remove * if nonsense mutation ends with *
        pos = gsub(pattern = '^\\*', replacement = '', x = pos) #Remove * if nonsense mutation starts with *
	  pos = gsub(pattern = '\\*.*', replacement = '', x = pos) #Remove * followed by position e.g, p.C229Lfs*18

  pos = suppressWarnings( as.numeric(sapply(strsplit(x = pos, split = '_', fixed = TRUE), '[', 1)) )
    all.prot.dat[,pos := pos]

  if(nrow( all.prot.dat[is.na(all.prot.dat$pos),]) > 0){
      #message(paste('Removed', nrow( all.prot.dat[is.na(all.prot.dat$pos),]), 'mutations for which AA position was not available', sep = ' '))
          #print(prot.dat[is.na(prot.dat$pos),])
	      all.prot.dat = all.prot.dat[!is.na(all.prot.dat$pos),]
	        }

  gene.sum = summarizeMaf_fix(maf = dat)$gene.summary
    #gene.sum = merge.data.frame(x = gene.sum, y = gl, by = 'Hugo_Symbol', all.x = TRUE)
      gene.sum = merge(x = gene.sum, y = gl, by = 'Hugo_Symbol', all.x = TRUE)
        #gene.sum = gene.sum[!is.na(gene.sum$aa.length),]
	  gene.sum = gene.sum[!is.na(gene.sum$aa.length)]

  num_mut_colIndex = which(colnames(gene.sum) == 'total')
    aalen_colIndex = which(colnames(gene.sum) == 'aa.length')

  #Get background threshold
    gene.sum$th = apply(gene.sum, 1, function(x) get_threshold(gene_muts = as.numeric(x[num_mut_colIndex]), gene_length = as.numeric(x[aalen_colIndex])))
      #use only genes with atleast 2 (or m ) mutations.
        gene.sum = gene.sum[total >= m]

  if(calBg){
      if(nrow(gene.sum) < nBg){
            #message("Not enough genes to build background. Using predefined values. (Mean = 0.279; SD = 0.13)")
	          return(NULL)
		      } else{
		            syn.res = c()
			          pb <- txtProgressBar(min = 0, max = nrow(gene.sum), style = 3) #progress bar

      for(i in 1:nrow(gene.sum)){
              prot.dat = all.prot.dat[Hugo_Symbol == gene.sum[i, "Hugo_Symbol"]]
	              syn.res = rbind(syn.res, cluster_prot_fix(prot.dat = prot.dat, gene = gene.sum[i, "Hugo_Symbol"], th = gene.sum[i,"th"], protLen = gene.sum[i,"aa.length"]))
		              setTxtProgressBar(pb, i)
			            }
				          return(syn.res)
					      }
					        } else{
						    nonsyn.res = c()
						        pb <- txtProgressBar(min = 0, max = nrow(gene.sum), style = 3) #progress bar

    for(i in 1:nrow(gene.sum)){
          hs = gene.sum[i, Hugo_Symbol]
	        #print(hs)
		      prot.dat = all.prot.dat[Hugo_Symbol %in% hs]
		            nonsyn.res = rbind(nonsyn.res, cluster_prot_fix(prot.dat = prot.dat, gene = hs, th = gene.sum[Hugo_Symbol %in% hs, th], protLen = gene.sum[Hugo_Symbol %in% hs, aa.length]))
			          setTxtProgressBar(pb, i)
				      }
				          return(nonsyn.res)
					    }
					    }

cluster_prot_fix = function(prot.dat, gene, th, protLen){

  mergeDist = 5 #hard coded inter event distance.
    #prot.dat = all.prot.dat[Hugo_Symbol == gene]

  #Summarise counts per position
    pos.counts = prot.dat[,.N,pos]
      pos.counts = pos.counts[order(pos)]

  #classify position as meaningful if its greater than background threshhold.
    pos.counts$cluster = ifelse(test = pos.counts$N >= th, yes = 'meaningful', no = 'nonMeaningful')

  #Just choose meaningful positions
    clust.tbl = pos.counts[cluster %in% 'meaningful']
      nonclust.tbl = pos.counts[cluster %in% 'nonMeaningful']

  if(nrow(clust.tbl) == 0){
      #message(paste('No meaningful positions found for', gene, sep=' '))
          return(NULL)
	    }

  clust.tbl$distance = c(0,diff(clust.tbl$pos)) #calculate inter event distance.

  #If more than one meaningful positions are found within a 5 aa distance, join them to form a cluster.
    if(nrow(clust.tbl) > 1){

    #initialize variables.
        cstart = end = clust.tbl[1,pos]
	    n = clust.tbl[1,N]
	        cdf = c()
		    cluster = 1

    #Go through entire table and update variables.
        for(i in 2:nrow(clust.tbl)){
	      pos = clust.tbl[i,pos]

      d = clust.tbl[i,distance]

      if(d < mergeDist){
              end = pos
	              n = n + clust.tbl[i,N]
		            }else{
			            tempdf = data.frame(cluster = paste('cluster', cluster, sep='_'), start = cstart, end = end ,N = n)
				            cdf = rbind(cdf, tempdf)
					            cstart = end = pos
						            n = clust.tbl[i,N]
							            cluster = cluster + 1
								          }
									      }
									          cdf = rbind(cdf, data.frame(cluster = paste('cluster', cluster, sep='_'), start = cstart, end = end ,N = n))
										    } else {
										        cdf = data.frame(cluster = 'cluster_1', start = clust.tbl$pos, end = clust.tbl$pos ,N = clust.tbl$N)
											  }

  #merge adjacent variants to clusters.
    for(i in 1:nrow(cdf)){
        tempcdf = cdf[i,]
	    nonclust.tbl$startDist = nonclust.tbl$pos - tempcdf$start
	        nonclust.tbl$endDist = nonclust.tbl$pos - tempcdf$end

    merge.adj.to.start = nonclust.tbl[startDist >= -5 & startDist <= 0]
        if(nrow(merge.adj.to.start) > 0){
	      tempcdf$start = merge.adj.to.start[which(merge.adj.to.start$startDist == min(merge.adj.to.start$startDist)),pos]
	            tempcdf$N = tempcdf$N + sum(merge.adj.to.start$N)
		        }

    merge.adj.to.end = nonclust.tbl[endDist <= 5 & endDist >= 0]
        if(nrow(merge.adj.to.end) > 0){
	      tempcdf$end = merge.adj.to.end[which(merge.adj.to.end$endDist == max(merge.adj.to.end$endDist)),pos]
	            tempcdf$N = tempcdf$N + sum(merge.adj.to.end$N)
		        }
			    cdf[i,] = tempcdf
			      }
			        cdf$Hugo_Symbol = gene

  #Calcluate cluster score.

  total.muts = nrow(prot.dat) #total variants for this gene.
    clusterScores = c()

  for(i in 1:nrow(cdf)){
      temp.prot.dat = prot.dat[pos >= as.numeric(cdf$start[i]) & pos <= as.numeric(cdf$end[i])]
          temp.prot.dat.summary = temp.prot.dat[,.N, pos]
	      temp.prot.dat.summary[,fraction:= N/total.muts]

    peak = temp.prot.dat.summary[N == max(N), pos]

    posVector = as.numeric(temp.prot.dat.summary[,pos])
        fractionMutVector = unlist(lapply(posVector, FUN = function(x) temp.prot.dat.summary[pos == x, fraction]))
	    distanceVector = suppressWarnings(abs(posVector - peak))

    clusterScores = c(clusterScores,  sum( fractionMutVector / (sqrt(2)^ distanceVector)))

  }

  cdf$clusterScore = clusterScores

  gene.clust.res = data.frame(Hugo_Symbol = gene, clusters = nrow(cdf), muts_in_clusters = sum(cdf$N), clusterScores = sum(cdf$clusterScore), protLen = protLen)
    return(gene.clust.res)
    }
    




createOncoMatrix<- function(maf){

    message('Creating oncomatrix (this might take a while)..')

     oncomat = data.table::dcast(data = maf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                                      fun.aggregate = function(x) {ifelse(test = length(as.character(x))>1 ,
				                                      no = as.character(x), yes = vcr(x, gis = FALSE))
								                                       }, value.var = 'Variant_Classification', fill = '')

    #If maf contains only one sample converting to matrix is not trivial.
        if(ncol(oncomat) == 2){
	      genes = oncomat[,Hugo_Symbol]
	            sampleId = colnames(oncomat)[2]
		          oncomat = as.matrix(data.frame(row.names = genes, sample = oncomat[,2, with =FALSE]))
			      }else if(nrow(oncomat) == 1){
			            #If MAF has only one gene
				          gene = oncomat[,Hugo_Symbol]
					        oncomat[,Hugo_Symbol:= NULL]
						      oncomat = as.matrix(oncomat)
						            rownames(oncomat) = gene
							          sampleID = colnames(oncomat)
								        }else{
									      oncomat = as.matrix(oncomat)
									            rownames(oncomat) = oncomat[,1]
										          oncomat = oncomat[,-1]
											        }

     variant.classes = as.character(unique(maf[,Variant_Classification]))
          variant.classes = c('',variant.classes, 'Multi_Hit')
	       names(variant.classes) = 0:(length(variant.classes)-1)

     #Complex variant classes will be assigned a single integer.
          vc.onc = unique(unlist(apply(oncomat, 2, unique)))
	       vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
	            names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
		         variant.classes2 = c(variant.classes, vc.onc)

     oncomat.copy <- oncomat
         #Make a numeric coded matrix
	     for(i in 1:length(variant.classes2)){
	           oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
		       }

    #If maf has only one gene
        if(nrow(oncomat) == 1){
	      mdf  = t(matrix(as.numeric(oncomat)))
	            rownames(mdf) = gene
		          colnames(mdf) = sampleID
			        return(list(oncomat = oncomat.copy, nummat = mdf, vc = variant.classes))
				    }

    #convert from character to numeric
        mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
	    rownames(mdf) = rownames(oncomat.copy)

    message('Sorting..')

    #If MAF file contains a single sample, simple sorting is enuf.
        if(ncol(mdf) == 1){
	      mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
	            colnames(mdf) = sampleId

      oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
            colnames(oncomat.copy) = sampleId

      return(list(oncomat = oncomat.copy, nummat = mdf, vc = variant.classes))
          } else{
	        #Sort by rows as well columns if >1 samples present in MAF
		      #Add total variants per gene
		            mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
			            length(x[x != "0"])
				          }))
					        #Sort by total variants
						      mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
						            colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
							          nMut = mdf[, ncol(mdf)]

      mdf = mdf[, -ncol(mdf)]

      mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

      mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
            tmdf = t(mdf) #transposematrix
	          mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

      mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
            mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
	          mdf = mdf.temp.copy

      #organise original character matrix into sorted matrix
            oncomat.copy <- oncomat.copy[,colnames(mdf)]
	          oncomat.copy <- oncomat.copy[rownames(mdf),]

      return(list(oncomat = oncomat.copy, nummat = mdf, vc = variant.classes))
          }
	  }

validateMaf<-function(maf, rdup = TRUE, isTCGA = isTCGA){

  #necessary fields.
    required.fields = c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
                          'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode')

  #Change column names to standard names; i.e, camel case
    for(i in 1:length(required.fields)){
        colId = suppressWarnings(grep(pattern = required.fields[i], x = colnames(maf), ignore.case = TRUE))
	    if(length(colId) > 0){
	          colnames(maf)[colId] = required.fields[i]
		      }
		        }

  missing.fileds = required.fields[!required.fields %in% colnames(maf)] #check if any of them are missing

  if(length(missing.fileds) > 0){
      missing.fileds = paste(missing.fileds[1], sep = ',', collapse = ', ')
          stop(paste('missing required fields from MAF:', missing.fileds)) #stop if any of required.fields are missing
	    }

  #convert "-" to "." in "Tumor_Sample_Barcode" to avoid complexity in naming
    maf$Tumor_Sample_Barcode = gsub(pattern = '-', replacement = '.', x = as.character(maf$Tumor_Sample_Barcode))

  if(rdup){
      maf = maf[, variantId := paste(Chromosome, Start_Position, Tumor_Sample_Barcode, sep = ':')]
          if(nrow(maf[duplicated(variantId)]) > 0){
	        message("NOTE: Removed ",  nrow(maf[duplicated(variantId)]) ," duplicated variants")
		      maf = maf[!duplicated(variantId)]
		          }
			      maf[,variantId := NULL]
			        }

  if(nrow(maf[Hugo_Symbol %in% ""]) > 0){
      message('NOTE: Found ', nrow(maf[Hugo_Symbol %in% ""]), ' variants with no Gene Symbols.')
          print(maf[Hugo_Symbol %in% "", required.fields, with = FALSE])
	      message("Annotating them as 'UnknownGene' for convenience")
	          maf$Hugo_Symbol = ifelse(test = maf$Hugo_Symbol == "", yes = 'UnknownGene', no = maf$Hugo_Symbol)
		    }

  if(nrow(maf[is.na(Hugo_Symbol)]) > 0){
      message('NOTE: Found ', nrow(maf[is.na(Hugo_Symbol) > 0]), ' variants with no Gene Symbols.')
          print(maf[is.na(Hugo_Symbol), required.fields, with =FALSE])
	      message("Annotating them as 'UnknownGene' for convenience")
	          maf$Hugo_Symbol = ifelse(test = is.na(maf$Hugo_Symbol), yes = 'UnknownGene', no = maf$Hugo_Symbol)
		    }

  if(isTCGA){
      maf$Tumor_Sample_Barcode = substr(x = maf$Tumor_Sample_Barcode, start = 1, stop = 12)
        }

  return(maf)
  }

read.maf_fix = function(maf, removeSilent = TRUE, useAll = TRUE, gisticAllLesionsFile = NULL, gisticAmpGenesFile = NULL,
                    gisticDelGenesFile = NULL, cnTable = NULL, removeDuplicatedVariants = TRUE, isTCGA = FALSE){
  
  message('reading maf..')
  
  if(as.logical(length(grep(pattern = 'gz$', x = maf, fixed = FALSE)))){
    #If system is Linux use fread, else use gz connection to read gz file.
    if(Sys.info()[['sysname']] == 'Windows'){
      maf.gz = gzfile(description = maf, open = 'r')
      suppressWarnings(maf <- data.table(read.csv(file = maf.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE)))
      close(maf.gz)
    } else{
      maf = suppressWarnings(data.table::fread(input = paste('zcat <', maf), sep = '\t', stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE, showProgress = TRUE, header = TRUE))
    }
  } else{
    suppressWarnings(maf <- data.table::fread(input = maf, sep = "\t", stringsAsFactors = FALSE, verbose = FALSE, data.table = TRUE, showProgress = TRUE, header = TRUE))
  }
  
  #validate MAF file
  maf = validateMaf(maf = maf, isTCGA = isTCGA, rdup = removeDuplicatedVariants)
  
  #validation check for variants classified as Somatic in Mutation_Status field.
  if(length(colnames(maf)[colnames(x = maf) %in% 'Mutation_Status']) > 0){
    if(!useAll){
      message('Using only Somatic variants from Mutation_Status. Switch on useAll to include everything.')
      maf = maf[Mutation_Status %in% "Somatic"]
      
      if(nrow(maf) == 0){
        stop('No more Somatic mutations left after filtering for Mutation_Status! Maybe set useAll to TRUE ?')
      }
      
      #maf = subset(maf, Mutation_Status == 'Somatic')
    }else {
      message('Using all variants.')
    }
  }else{
    message('Mutation_Status not found. Assuming all variants are Somatic and validated.')
  }
  #Variant Classification with Low/Modifier variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  silent = c("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron",
             "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA")
  #Variant Classification with High/Moderate variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                   "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                   "In_Frame_Ins", "Missense_Mutation")
  
  maf.silent = maf[Variant_Classification %in% silent]
  
  if(removeSilent){
    
    if(nrow(maf.silent) > 0){
      maf.silent.vc = maf.silent[,.N, .(Tumor_Sample_Barcode, Variant_Classification)]
      maf.silent.vc.cast = data.table::dcast(data = maf.silent.vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N') #why dcast is not returning it as data.table ?
      summary.silent = data.table(ID = c('Samples',colnames(maf.silent.vc.cast)[2:ncol(maf.silent.vc.cast)]),
                                  N = c(nrow(maf.silent.vc.cast), colSums(maf.silent.vc.cast[,2:ncol(maf.silent.vc.cast), with = FALSE])))
      
      maf = maf[!Variant_Classification %in% silent] #Remove silent variants from main table
      message(paste('Excluding',nrow(maf.silent), 'silent variants.'))
      print(summary.silent)
    } else{
      message(message(paste('Excluding',nrow(maf.silent), 'silent variants.')))
    }
  }else{
    message('Silent variants are being kept!')
  }
  
  if(!is.null(gisticAllLesionsFile)){
    gisticIp = readGistic(gisticAllLesionsFile = gisticAllLesionsFile, gisticAmpGenesFile = gisticAmpGenesFile,
                          gisticDelGenesFile = gisticDelGenesFile, isTCGA = isTCGA)
    gisticIp = gisticIp@data
    
    gisticIp[, id := paste(Hugo_Symbol, Tumor_Sample_Barcode, sep=':')]
    gisticIp = gisticIp[!duplicated(id)]
    gisticIp[,id := NULL]
    
    maf = rbind(maf, gisticIp, fill =TRUE)
    oncomat = createOncoMatrix(maf)
  }else if(!is.null(cnTable)){
    message('Processing copy number data..')
    cnDat = data.table::fread(input = cnTable, sep = '\t', stringsAsFactors = FALSE, header = TRUE, colClasses = 'character')
    colnames(cnDat) = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'Variant_Classification')
    cnDat$Variant_Type = 'CNV'
    suppressWarnings(cnDat[, id := paste(Hugo_Symbol, Tumor_Sample_Barcode, sep=':')])
    cnDat = cnDat[!duplicated(id)]
    cnDat[,id := NULL]
    maf = rbind(maf, cnDat, fill =TRUE)
    oncomat = createOncoMatrix(maf)
  }else{
    oncomat = createOncoMatrix(maf)
  }
  
  #convert to factors
  maf$Variant_Type = as.factor(as.character(maf$Variant_Type))
  maf$Variant_Classification = as.factor(as.character(maf$Variant_Classification))
  maf$Tumor_Sample_Barcode = as.factor(as.character(maf$Tumor_Sample_Barcode))
  
  message('Summarizing..')
  mafSummary = summarizeMaf_fix(maf = maf)
  
  #Create MAF object
  m = MAF(data = maf, variants.per.sample = mafSummary$variants.per.sample, variant.type.summary = mafSummary$variant.type.summary,
          variant.classification.summary = mafSummary$variant.classification.summary,gene.summary = mafSummary$gene.summary,
          oncoMatrix = oncomat$oncomat, numericMatrix = oncomat$nummat, summary = mafSummary$summary,
          classCode = oncomat$vc, maf.silent = maf.silent)
  
  
  message('Done !')
  return(m)
}


#' Class MAF
#' @description S4 class for storing summarized MAF.
#' @slot data data.table of original MAF file.
#' @slot variants.per.sample table containing variants per sample
#' @slot variant.type.summary table containing variant types per sample
#' @slot variant.classification.summary table containing variant classification per sample
#' @slot gene.summary table containing variant classification per gene
#' @slot oncoMatrix character matrix of dimension n*m where n is number of genes and m is number of variants
#' @slot numericMatrix numeric matrix of dimension n*m where n is number of genes and m is number of variants
#' @slot summary table with basic MAF summary stats
#' @slot classCode mapping between numeric values in numericMatrix and Variant Classification
#' @slot maf.silent subset of main MAF containing only silent variants
#' @exportClass MAF
#' @import methods
#' @seealso \code{\link{getGeneSummary}} \code{\link{getSampleSummary}} \code{\link{getFields}}

## MAF object
MAF <- setClass(Class = 'MAF', slots =  c(data = 'data.table', variants.per.sample = 'data.table', variant.type.summary = 'data.table',
                                         variant.classification.summary = 'data.table', gene.summary = 'data.table', oncoMatrix = 'matrix',
					                                          numericMatrix = 'matrix', summary = 'data.table', classCode = 'character',
										                                           maf.silent = 'data.table'))

setMethod(f = 'show', signature = 'MAF', definition = function(object){
  cat(paste('An object of class ', class(object), "\n"))
    print(object@summary)
    })


summarizeMaf_fix = function(maf){
  
  if('NCBI_Build' %in% colnames(maf)){
    NCBI_Build = unique(maf[!Variant_Type %in% 'CNV', NCBI_Build])
    NCBI_Build = NCBI_Build[!is.na(NCBI_Build)]
    
    if(length(NCBI_Build) > 1){
      message('NOTE: Mutiple reference builds found!')
      NCBI_Build = do.call(paste, c(as.list(NCBI_Build), sep=";"))
      message(NCBI_Build)
    }
  }else{
    NCBI_Build = NA
  }
  
  if('Center' %in% colnames(maf)){
    Center = unique(maf[!Variant_Type %in% 'CNV', Center])
    #Center = Center[is.na(Center)]
    if(length(Center) > 1){
      message('Mutiple centers found.')
      Center = do.call(paste, c(as.list(Center), sep=";"))
      print(Center)
    }
  }else{
    Center = NA
  }
  
  #nGenes
  nGenes = length(unique(maf[,Hugo_Symbol]))
  
  
  
  #Top 20 FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/
  flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
            "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
            "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17")
  
  #Variants per TSB
  tsb = maf[,.N, Tumor_Sample_Barcode]
  colnames(tsb)[2] = 'Variants'
  tsb = tsb[order(tsb$Variants, decreasing = TRUE),]
  
  #summarise and casting by 'Variant_Classification'
  vc = maf[,.N, .(Tumor_Sample_Barcode, Variant_Classification )]
  vc.cast = data.table::dcast(data = vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N')
  
  if(any(colnames(vc.cast) %in% c('Amp', 'Del'))){
    vc.cast.cnv = vc.cast[,colnames(vc.cast)[colnames(vc.cast) %in% c('Amp', 'Del')], with =FALSE]
    vc.cast.cnv$CNV_total = rowSums(x = vc.cast.cnv)
    
    vc.cast = vc.cast[,!colnames(vc.cast)[colnames(vc.cast) %in% c('Amp', 'Del')], with =FALSE]
    vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = FALSE])]
    
    vc.cast = cbind(vc.cast, vc.cast.cnv)
    vc.cast = vc.cast[order(total, CNV_total, decreasing = TRUE)]
    
    vc.mean = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, mean))))
    vc.median = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, median))))
  }else{
    vc.cast[,total:=rowSums(vc.cast[,2:ncol(vc.cast), with = FALSE])]
    vc.cast = vc.cast[order(total, decreasing = TRUE)]
    
    vc.mean = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, mean))))
    vc.median = as.numeric(as.character(c(NA, NA, NA, NA, apply(vc.cast[,2:ncol(vc.cast), with = FALSE], 2, median))))
  }
  
  #summarise and casting by 'Variant_Type'
  vt = maf[,.N, .(Tumor_Sample_Barcode, Variant_Type )]
  vt.cast = data.table::dcast(data = vt, formula = Tumor_Sample_Barcode ~ Variant_Type, value.var = 'N', fill = 0)
  if(any(colnames(vt.cast) %in% c('CNV'))){
    vt.cast.cnv = vt.cast[,colnames(vt.cast)[colnames(vt.cast) %in% c('CNV')], with =FALSE]
    
    vt.cast = vt.cast[,!colnames(vt.cast)[colnames(vt.cast) %in% c('CNV')], with =FALSE]
    vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = FALSE])]
    vt.cast = vt.cast[order(total, decreasing = TRUE)]
    
    vt.cast = cbind(vt.cast, vt.cast.cnv)
    vt.cast[order(total, CNV, decreasing = TRUE)]
  }else{
    vt.cast[,total:=rowSums(vt.cast[,2:ncol(vt.cast), with = FALSE])]
    vt.cast = vt.cast[order(total, decreasing = TRUE)]
  }
  
  #summarise and casting by 'Hugo_Symbol'
  hs = maf[,.N, .(Hugo_Symbol, Variant_Classification)]
  hs.cast = data.table::dcast(data = hs, formula = Hugo_Symbol ~Variant_Classification, fill = 0, value.var = 'N')
  #----
  if(any(colnames(hs.cast) %in% c('Amp', 'Del'))){
    hs.cast.cnv = hs.cast[,colnames(hs.cast)[colnames(hs.cast) %in% c('Amp', 'Del')], with =FALSE]
    hs.cast.cnv$CNV_total = rowSums(x = hs.cast.cnv)
    
    hs.cast = hs.cast[,!colnames(hs.cast)[colnames(hs.cast) %in% c('Amp', 'Del')], with =FALSE]
    hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE])]
    
    hs.cast = cbind(hs.cast, hs.cast.cnv)
    hs.cast = hs.cast[order(total, CNV_total, decreasing = TRUE)]
    
  }else{
    hs.cast[,total:=rowSums(hs.cast[,2:ncol(hs.cast), with = FALSE])]
    hs.cast = hs.cast[order(total, decreasing = TRUE)]
    
  }
  #Get in how many samples a gene ismutated
  numMutatedSamples = maf[!Variant_Type %in% 'CNV', .(MutatedSamples = length(unique(Tumor_Sample_Barcode))), by = Hugo_Symbol]
  #Merge and sort
  hs.cast = merge(hs.cast, numMutatedSamples, by = 'Hugo_Symbol', all = TRUE)
  hs.cast = hs.cast[order(MutatedSamples, total, decreasing = TRUE)]
  #Make a summarized table
  summary = data.table::data.table(ID = c('NCBI_Build', 'Center','Samples', 'nGenes',colnames(vc.cast)[2:ncol(vc.cast)]),
                                   summary = c(NCBI_Build, Center, nrow(vc.cast), nGenes, colSums(vc.cast[,2:ncol(vc.cast), with =FALSE])))
  summary[,Mean := vc.mean]
  summary[,Median := vc.median]
  
  print(summary)
  
  message("Frequently mutated genes..")
  print(hs.cast)
  
  #Check for flags.
  if(nrow(hs.cast) > 10){
    topten = hs.cast[1:10, Hugo_Symbol]
    topten = topten[topten %in% flags]
    if(length(topten) > 0){
      message('NOTE: Possible FLAGS among top ten genes:')
      print(topten)
    }
  }
  
  return(list(variants.per.sample = tsb, variant.type.summary = vt.cast, variant.classification.summary = vc.cast,
              gene.summary = hs.cast, summary = summary))
}

oncodrive_fix = function(maf, AACol = NULL, minMut = 5, pvalMethod = 'zscore', nBgGenes = 100, bgEstimate = TRUE, ignoreGenes = NULL){

  #Proetin Length source
    gl = system.file('extdata', 'prot_len.txt.gz', package = 'maftools')

  if(Sys.info()[['sysname']] == 'Windows'){
      gl.gz = gzfile(description = gl, open = 'r')
          gl <- suppressWarnings( data.table(read.csv( file = gl.gz, header = TRUE, sep = '\t', stringsAsFactors = FALSE)) )
	      close(gl.gz)
	        } else{
		    gl = data.table::fread(input = paste('zcat <', gl), sep = '\t', stringsAsFactors = FALSE)
		      }

  pval.options = c('zscore', 'poisson', 'combined')

  if(!pvalMethod %in% pval.options){
      stop('pvalMethod can only be either zscore, poisson or combined')
        }

  if(length(pvalMethod) > 1){
      stop('pvalMethod can only be either zscore, poisson or combined')
        }


  #syn variants for background
    syn.maf = maf@maf.silent
      #number of samples in maf
        numSamples = as.numeric(maf@summary[3,summary])
	  #Perform clustering and calculate background scores.
	    if(bgEstimate){
	        if(nrow(syn.maf) == 0){
		      message('No syn mutations found! Skipping background estimation. Using predefined values. (Mean = 0.279; SD = 0.13)')
		            bg.mean = 0.279
			          bg.sd = 0.13
				      }else{
				            message('Estimating background scores from synonymous variants..')
					          syn.bg.scores = parse_prot_fix(dat = syn.maf, AACol = AACol, gl, m = minMut, calBg = TRUE, nBg = nBgGenes)

      #If number of genes to calculate background scores is not enough, use predefined scores.
            if(is.null(syn.bg.scores)){
	            message("Not enough genes to build background. Using predefined values. (Mean = 0.279; SD = 0.13)")
		            bg.mean = 0.279
			            bg.sd = 0.13
				          }else {
					          if(nrow(syn.bg.scores) < nBgGenes){
						            message("Not enough genes to build background. Using predefined values. (Mean = 0.279; SD = 0.13)")
							              bg.mean = 0.279
								                bg.sd = 0.13
										        }else{
											          bg.mean = mean(syn.bg.scores$clusterScores)
												            bg.sd = sd(syn.bg.scores$clusterScores)
													              message(paste('Estimated background mean: ', bg.mean))
														                message(paste('Estimated background SD: ', bg.sd))
																        }
																	      }
																	          }
																		    }else{
																		        message("Using predefined values for background. (Mean = 0.279; SD = 0.13)")
																			    bg.mean = 0.279
																			        bg.sd = 0.13
																				  }



  #non-syn variants
    non.syn.maf = maf@data
      #Variant Classification with Low/Modifier variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
        silent = c("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron",
	             "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA", "Amp", "Del")
		       non.syn.maf = non.syn.maf[!Variant_Classification %in% silent] #Remove silent variants from main table

  #Remove genes to ignore
    if(!is.null(ignoreGenes)){
        ignoreGenes.count = nrow(non.syn.maf[Hugo_Symbol %in% ignoreGenes])
	    message(paste('Removed', ignoreGenes.count, 'variants belonging to', paste(ignoreGenes, collapse = ', ', sep=',')))
	        non.syn.maf = non.syn.maf[!Hugo_Symbol %in% ignoreGenes]
		  }

  #Perform clustering and calculate cluster scores for nonsyn variants.
    message('Estimating cluster scores from non-syn variants..')
      nonsyn.scores = parse_prot_fix(dat = non.syn.maf, AACol = AACol, gl = gl, m = minMut, calBg = FALSE, nBg = nBgGenes)

  if(pvalMethod == 'combined'){
      message('Comapring with background model and estimating p-values..')
          nonsyn.scores$zscore = (nonsyn.scores$clusterScores - bg.mean) / bg.sd
	      nonsyn.scores$tPval = 1- pnorm(nonsyn.scores$zscore)
	          nonsyn.scores$tFdr = p.adjust(nonsyn.scores$tPval, method = 'fdr')

    nonsyn.scores = merge(getGeneSummary(maf), nonsyn.scores, by = 'Hugo_Symbol')
        nonsyn.scores[,fract_muts_in_clusters := muts_in_clusters/total]

    counts.glm = glm(formula = total ~ protLen+clusters, family = poisson(link = identity), data = nonsyn.scores) #Poisson model
        nonsyn.scores$Expected = counts.glm$fitted.values #Get expected number of events (mutations) from the model

    observed_mut_colIndex = which(colnames(nonsyn.scores) == 'total')
        expected_mut_colIndex = which(colnames(nonsyn.scores) == 'Expected')

    #Poisson test to caluclate difference (p-value)
        nonsyn.scores$poissonPval = apply(nonsyn.scores, 1, function(x) {
	      poisson.test(as.numeric(x[observed_mut_colIndex]), as.numeric(x[expected_mut_colIndex]))$p.value
	          })

    nonsyn.scores$poissonFdr = p.adjust(nonsyn.scores$poissonPval)
        nonsyn.scores = nonsyn.scores[order(poissonFdr)]

    nonsyn.scores$fdr = apply(nonsyn.scores[,.(tFdr, poissonFdr)], MARGIN = 1, FUN = min)

  } else if(pvalMethod == 'zscore'){
      #Oncodrive clust way of caluclating pvalues
          #Calculate z scores; compare it to bg scores and estimate z-score, pvalues, corrected pvalues (fdr) (assumes normal distribution)
	      message('Comapring with background model and estimating p-values..')
	          nonsyn.scores$zscore = (nonsyn.scores$clusterScores - bg.mean) / bg.sd
		      nonsyn.scores$pval = 1- pnorm(nonsyn.scores$zscore)
		          nonsyn.scores$fdr = p.adjust(nonsyn.scores$pval, method = 'fdr')

    nonsyn.scores = merge(getGeneSummary(maf), nonsyn.scores, by = 'Hugo_Symbol')
        nonsyn.scores[,fract_muts_in_clusters := muts_in_clusters/total]
	    #nonsyn.scores[,fract_MutatedSamples := MutatedSamples/numSamples]
	        nonsyn.scores = nonsyn.scores[order(fdr)]
		  }else{
		      #Assuming poisson distribution of mutation counts
		          #Now model observed number of mutations as a function of number of clusters and protein length. Calculate expected number of events based on poisson distribution.
			      nonsyn.scores = merge(getGeneSummary(maf), nonsyn.scores, by = 'Hugo_Symbol')
			          nonsyn.scores[,fract_muts_in_clusters := muts_in_clusters/total]

    counts.glm = glm(formula = total ~ protLen+clusters, family = poisson(link = identity), data = nonsyn.scores) #Poisson model
        nonsyn.scores$Expected = counts.glm$fitted.values #Get expected number of events (mutations) from the model

    observed_mut_colIndex = which(colnames(nonsyn.scores) == 'total')
        expected_mut_colIndex = which(colnames(nonsyn.scores) == 'Expected')

    #Poisson test to caluclate difference (p-value)
        nonsyn.scores$pval = apply(nonsyn.scores, 1, function(x) {
	      poisson.test(as.numeric(x[observed_mut_colIndex]), as.numeric(x[expected_mut_colIndex]))$p.value
	          })

    nonsyn.scores$fdr = p.adjust(nonsyn.scores$pval)
        nonsyn.scores = nonsyn.scores[order(fdr)]
	  }
	    message('Done !')
	      return(nonsyn.scores)
	      }


laml = read.maf(maf = args$input_maf, removeSilent = F, useAll = T)

if(is.null(args$gene_blacklist)){
        laml.sig = oncodrive(maf =laml, AACol = aacol, pvalMethod = 'zscore',minMut = min_mut)
	write.table(laml.sig,file=args$output_detail, quote=FALSE,row.names=FALSE,sep="\t")
	pdf(args$output_plot)
	plotOncodrive(res=laml.sig,fdrCutOff=as.numeric(args$fdr),useFraction=TRUE)
	dev.off()
	}else{
        all_genes <- read.table(args$gene_blacklist, stringsAsFactors=FALSE)[,1]
        laml.sig = oncodrive(maf =laml, AACol = aacol, pvalMethod = 'zscore',minMut = min_mut,ignoreGenes=all_genes)
	write.table(laml.sig,file=args$output_detail, quote=FALSE,row.names=FALSE,sep="\t")
	pdf(args$output_plot)
	plotOncodrive(res=laml.sig,fdrCutOff=as.numeric(args$fdr),useFraction=TRUE)
	dev.off()
}