
library(ExomeDepth)
library(GenomeInfoDb)
library(Rsamtools)

#setwd("/mapped/")

getwd()
setwd("./final/")

if (!dir.exists("../results")) {
  dir.create("../results")
}

exons.hg19 <- read.table("../100bp_exon.bed", sep="\t", as.is=TRUE)
segments <- cbind(rownames(exons.hg19), exons.hg19)
rownames(exons.hg19) <- NULL
colnames(exons.hg19) <- c("chromosome","start","end","name")


exons.hg19.GRanges <- GRanges(seqnames = exons.hg19$chromosome,
								IRanges(start=exons.hg19$start,end=exons.hg19$end),
								names = exons.hg19$name)

reference.file <- "../index/hg19.fa"

analysisConfig <- read.csv('../config_single.csv',
							              header = TRUE,
							              fill = TRUE)

list_of_bam_files <- as.vector(analysisConfig$list_of_bam_files)

#list_of_bam_files

my.counts <- getBamCounts(bed.frame = exons.hg19,
                          bam.files = list_of_bam_files,
                          include.chr = FALSE,
						  referenceFasta = reference.file
						  )

#print(head(my.counts))

# Create dataframe

ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')
ExomeCount.dafr$chr <- gsub(as.character(ExomeCount.dafr$chromosome),
pattern = 'chr',
replacement = '')

head(ExomeCount.dafr)

# Create matrix of the bam counts

ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr),
                                                        pattern = '*.final.bam')])

head(ExomeCount.mat)

nsamples <- ncol(ExomeCount.mat)

for (i in 1:nsamples) {

		

		# Combine multiple samples to optimize the reference set in order to maximize the power to detect CNV
		### start looping over each sample

		my.test.data <- as.matrix(ExomeCount.mat[, i])

		my.reference.set <- as.matrix(ExomeCount.mat[, -i])

		head(my.test.data)
		head(my.reference.set)	

		my.choice <- select.reference.set(test.counts = my.test.data,
											reference.counts = my.reference.set,
											bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
											n.bins.reduced = 10000)



		head(my.choice)

		my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE])

		my.reference.set <- apply(X = my.matrix, MAR = 1, FUN = sum)
		
		## CNV calling

		all.exons <- new('ExomeDepth',
                  test = ExomeCount.mat[,i],
                  reference = my.reference.set,
                  formula = 'cbind(test, reference) ~ 1')


		all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = ExomeCount.dafr$chr,
                      start = ExomeCount.dafr$start,
                      end = ExomeCount.dafr$end,
                      name = ExomeCount.dafr$exon)
					
		#check output
		head(all.exons@CNV.calls)

		# ranking by BF
		all.exons@CNV.calls <- all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),]

		#Now save it in an easily readable format

		# Now save it in BED format
		output.file <- paste('../results/', list_of_bam_files[i], '.bed', sep = '')
		write.table(file = output.file, 
					x = data.frame(chrom = all.exons@CNV.calls$chromosome,
                    Start = all.exons@CNV.calls$start,
                    End = all.exons@CNV.calls$end,
                    Variant_type = ifelse(all.exons@CNV.calls$type == "duplication", "DUP", "DEL"),
                    all.exons@CNV.calls), 
                    sep = "\t", 
                    quote = FALSE, 
                    row.names = FALSE,
                    col.names = TRUE)
	    

}


q()
