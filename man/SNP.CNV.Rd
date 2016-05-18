\name{SNP.CNV}
\alias{SNP.CNV}
\title{
CNV Analysis Pipeline for SNP array Data
}
\description{
All analysis steps are integrate into a pipeline. The results, 
including visualization plots are placed in a directory as specified by user.
}
\usage{
SNP.CNV(snp, output.dir, sample.id, 
    do.GC.adjust = FALSE,
    gc.file = system.file("extdata","GC_1kb_hg19.txt.gz",package="saasCNV"),
    min.chr.probe = 100, min.snps = 10, 
    joint.segmentation.pvalue.cutoff = 1e-04, max.chpts = 30, 
    do.merge = TRUE, use.null.data = TRUE, 
    num.perm = 1000, maxL = NULL, 
    merge.pvalue.cutoff = 0.05, 
    do.cnvcall.on.merge = TRUE, 
    cnvcall.pvalue.cutoff = 0.05, 
    do.boundary.refine = FALSE, 
    do.plot = TRUE, cex = 0.3, 
    ref.num.probe = NULL, 
    do.gene.anno = FALSE, 
    gene.anno.file = NULL, 
    seed = NULL, verbose = TRUE)
}
\arguments{
  \item{snp}{a data frame constructed from a text file with LRR and BAF information. 
  	%See \code{\link{vcf2txt}
}
  \item{output.dir}{the directory to which all the results will be located.
}
  \item{sample.id}{sample ID to be displayed in the data frame of the results and 
  	the title of some diagnosis plots.
}
  \item{do.GC.adjust}{logical. If GC content adjustment on \code{log2ratio} to be carried out. 
  	Default is \code{FALSE}. See \code{\link{GC.adjust}} for details.
}
  \item{gc.file}{the location of tab-delimit file with GC content (averaged per 1kb window) 
  	information. See \code{\link{GC.adjust}} for details.
}
  \item{min.chr.probe}{the minimum number of probes tagging a chromosome for it 
to be passed to the subsequent analysis.
}
  \item{min.snps}{the minimum number of probes a segment needs to span.
}
  \item{joint.segmentation.pvalue.cutoff}{the p-value cut-off one (or a pair) of change points to
  	be determined as significant in each cycle of joint segmentation.
}
  \item{max.chpts}{the maximum number of change points to be detected for each 
  	chromosome.
}
  \item{do.merge}{logical. If segments merging step to be carried out. 
  	Default is \code{TRUE}.
}
  \item{use.null.data}{logical. If only data for probes located in normal copy 
  	segments to be used for bootstrapping. Default is \code{TRUE}. If a more 
  	aggressive merging is needed, it can be switched to \code{FALSE}.
}
  \item{num.perm}{the number of replicates drawn by bootstrap.
}
  \item{maxL}{integer. The maximum length in terms of number of probes a bootstrapped segment
  	may span. Default is \code{NULL}. If \code{NULL}, It will be automatically specified 
  	as 1/100 of the number of data points.
}
  \item{merge.pvalue.cutoff}{a p-value cut-off for merging. If the empirical p-value is 
  	greater than the cut-off value, the two adjacent segments under consideration will 
  	be merged.
}
  \item{do.cnvcall.on.merge}{logical. If CNV call to be done for the segments after 
  	merging step. Default is \code{TRUE}. If \code{TRUE}, CNV call will be done on the 
  	segments resulting directly from joint segmentation without merging step.
}
  \item{cnvcall.pvalue.cutoff}{a p-value cut-off for CNV calling.
}
  \item{do.boundary.refine}{logical. If the segment boundaries based on the grid of 
  	heterozygous probes to be refined by all probes with LRR data. Default is \code{FALSE}. 
  	We do not recommend to perform this step except in the case that the segment 
  	boundaries need to be aligned well on the same grid of probes for downstream analysis.
}
  \item{do.plot}{logical. If diagnosis plots to be output. Default is \code{TRUE}.
}
  \item{cex}{a numerical value giving the amount by which plotting text and 
  	symbols should be magnified relative to the default. It can be adjusted 
  	in order to make the plot legible.
}
  \item{ref.num.probe}{integer. The reference number of probes against which 
  	a segment is compared in order to determine the cex of the segment 
  	to be displayed. Default is \code{NULL}. If \code{NULL}, It will be automatically 
  	specified as 1/100 of the number of data points.
}
  \item{do.gene.anno}{logical. If gene annotation step to be performed. Default is \code{FALSE}.
}
  \item{gene.anno.file}{a tab-delimited file containing gene annotation information. 
  	For example, RefSeq annotation file which can be found at UCSC genome browser.
}
  \item{seed}{integer. Random seed can be set for reproducibility of results.
}
  \item{verbose}{logical. If more details to be output. Default is \code{TRUE}.
}
}
\details{
See the vignettes of the package for more details.
}
\value{
The results, including visualization plots are placed in 
subdirectories of the output directory \code{output.dir}
as specified by user.
}
\references{
Zhongyang Zhang and Ke Hao. (2015) SAAS-CNV: A Joint Segmentation Approach on Aggregated and Allele Specific Signals 
for the Identification of Somatic Copy Number Alterations with Next-Generation Sequencing Data. 
PLoS Computational Biology, 11(11):e1004618.
}
\author{
Zhongyang Zhang <zhongyang.zhang@mssm.edu>
}
%% \note{
%%  ~~further notes~~
%% }

\seealso{
\code{\link{NGS.CNV}}, \code{\link{snp.cnv.data}},
\code{\link{joint.segmentation}}, \code{\link{merging.segments}}
\code{\link{cnv.call}}, \code{\link{diagnosis.seg.plot.chr}},
\code{\link{genome.wide.plot}}, \code{\link{diagnosis.cluster.plot}}, 
\code{\link{snp.refine.boundary}}
}
\examples{
\dontrun{
## the pipeline for SNP array analysis
## download snp_table.txt.gz
url <- "https://zhangz05.u.hpc.mssm.edu/saasCNV/data/snp_table.txt.gz"
tryCatch({download.file(url=url, destfile="snp_table.txt.gz")
         }, error = function(e) {
          download.file(url=url, destfile="snp_table.txt.gz", method="curl")
         })
## If download.file fails to download the data, please manually download it from the url.

snp_table <- read.delim(file="snp_table.txt.gz", as.is=TRUE)

## download refGene_hg19.txt.gz
url <- "https://zhangz05.u.hpc.mssm.edu/saasCNV/data/refGene_hg19.txt.gz"
tryCatch({download.file(url=url, destfile="refGene_hg19.txt.gz")
         }, error = function(e) {
          download.file(url=url, destfile="refGene_hg19.txt.gz", method="curl")
         })
## If download.file fails to download the data, please manually download it from the url.

sample.id <- "SNP_0116"
output.dir <- file.path(getwd(), "test_saasCNV")
	
SNP.CNV(snp=snp_table, output.dir=output.dir, sample.id=sample.id,
        min.chr.probe=100,
        min.snps=10,
        joint.segmentation.pvalue.cutoff=1e-4,
        max.chpts=30,
        do.merge=TRUE, use.null.data=TRUE, num.perm=1000, maxL=5000, 
        merge.pvalue.cutoff=0.05,
        do.cnvcall.on.merge=TRUE, 
        cnvcall.pvalue.cutoff=0.05,
        do.boundary.refine=TRUE,
        do.plot=TRUE, cex=0.3, ref.num.probe=5000,
        do.gene.anno=TRUE,
        gene.anno.file="refGene_hg19.txt.gz",
        seed=123456789,
        verbose=TRUE)
}

}
\keyword{ SNP array }
\keyword{ CNV }
\keyword{ pipeline }
