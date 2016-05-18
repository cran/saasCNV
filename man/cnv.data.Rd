\name{cnv.data}
\alias{cnv.data}
\title{
Construct Data Frame for CNV Inference with NGS Data
}
\description{
Transform read depth information into log2ratio and log2mBAF
that we use for joint segmentation and CNV calling.
}
\usage{
cnv.data(vcf, min.chr.probe = 100, verbose = FALSE)
}
\arguments{
  \item{vcf}{a data frame constructed from a vcf file. See \code{\link{vcf2txt}}.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{vcf2txt}}
}
\examples{
## load a data frame constructed from a vcf file with vcf2txt

\dontrun{
## download vcf_table.txt.gz
url <- "https://zhangz05.u.hpc.mssm.edu/saasCNV/data/vcf_table.txt.gz"
tryCatch({download.file(url=url, destfile="vcf_table.txt.gz")
         }, error = function(e) {
          download.file(url=url, destfile="vcf_table.txt.gz", method="curl")
         })
## If download.file fails to download the data, please manually download it from the url.

vcf_table <- read.delim(file="vcf_table.txt.gz", as.is=TRUE)
seq.data <- cnv.data(vcf=vcf_table, min.chr.probe=100, verbose=TRUE)
}

## see how seq.data looks like
data(seq.data)
head(seq.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ CNV }
