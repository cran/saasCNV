\name{snp.cnv.data}
\alias{snp.cnv.data}
\title{
Construct Data Frame for CNV Inference with SNP Array Data
}
\description{
Transform LRR and BAF information into log2ratio and log2mBAF
that we use for joint segmentation and CNV calling.
}
\usage{
snp.cnv.data(snp, min.chr.probe = 100, verbose = FALSE)
}
\arguments{
  \item{snp}{a data frame with LRR and BAF information from SNP array. See the 
  	example below for details.
}
  \item{min.chr.probe}{the minimum number of probes tagging a chromosome for it 
to be passed to the subsequent analysis.
}
  \item{verbose}{logical. If more details to be output. Default is \code{FALSE}.
}
}
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
A data frame containing the log2raio and log2mBAF values for each probe site.
}
\references{
Staaf, J., Vallon-Christersson, J., Lindgren, D., Juliusson, G., Rosenquist, R., 
Hoglund, M., Borg, A., Ringner, M. (2008) Normalization of Illumina Infinium 
whole-genome SNP data improves copy number estimates and allelic intensity ratios. 
\emph{BMC bioinformatics}, \bold{9}:409.
}
\author{
Zhongyang Zhang <zhongyang.zhang@mssm.edu>
}
%% \note{
%%  ~~further notes~~
%% }
\seealso{
\code{\link{cnv.data}}
}
\examples{
\dontrun{
## an example data with LRR and BAF information
url <- "https://zhangz05.u.hpc.mssm.edu/saasCNV/data/snp_table.txt.gz"
tryCatch({download.file(url=url, destfile="snp_table.txt.gz")
         }, error = function(e) {
          download.file(url=url, destfile="snp_table.txt.gz", method="curl")
         })
## If download.file fails to download the data, please manually download it from the url.

snp_table <- read.delim(file="snp_table.txt.gz", as.is=TRUE)
snp.data <- snp.cnv.data(snp=snp_table, min.chr.probe=100, verbose=TRUE) 

## see how seq.data looks like
url <- "https://zhangz05.u.hpc.mssm.edu/saasCNV/data/snp.data.RData"
tryCatch({download.file(url=url, destfile="snp.data.RData")
         }, error = function(e) {
          download.file(url=url, destfile="snp.data.RData", method="curl")
         })
## If download.file fails to download the data, please manually download it from the url.

load("snp.data.RData")
head(snp.data)
}

}
\keyword{ CNV }
