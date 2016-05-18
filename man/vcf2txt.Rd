\name{vcf2txt}
\alias{vcf2txt}
\title{
Covert VCF File to A Data Frame
}
\description{
It parses a VCF file and extract necessary information for CNV analysis.
}
\usage{
vcf2txt(vcf.file, normal.col = 10, tumor.col = 11, MQ.cutoff = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vcf.file}{vcf file name.
}
  \item{normal.col}{the number of the column in which the genotype
and read depth information of normal tissue are located in the vcf file.
}
  \item{tumor.col}{the number of the column in which the genotype
and read depth information of tumor tissue are located in the vcf file.
}
  \item{MQ.cutoff}{the minimum criterion of mapping quality.
}
}
\details{
Note that the first 9 columns in vcf file are mandatory, followed by the
information for called variant starting from the 10th column.
}
\value{
A data frame of detailed information about each variant, including 
chrosome position, reference and alternative alleles, genotype and read depth 
carrying reference and alternative alleles for normal and tumor respectively.
}
\references{
Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., 
Handsaker, R. E., Lunter, G., Marth, G. T., Sherry, S. T., et al. (2011) 
The variant call format and VCFtools. \emph{Bioinformatics}, \bold{27}:2156--2158.

http://www.1000genomes.org/node/101
}
\author{
Zhongyang Zhang <zhongyang.zhang@mssm.edu>
}
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
\dontrun{
## an example VCF file from WES
## download WES_example.vcf.gz
url <- "https://zhangz05.u.hpc.mssm.edu/saasCNV/data/WES_example.vcf.gz"
tryCatch({download.file(url=url, destfile="WES_example.vcf.gz")
         }, error = function(e) {
          download.file(url=url, destfile="WES_example.vcf.gz", method="curl")
         })
## If download.file fails to download the data, please manually download it from the url.

## convert VCf file to a data frame
vcf_table <- vcf2txt(vcf.file="WES_example.vcf.gz", normal.col=9+1, tumor.col=9+2)

## see how vcf_table looks like
## download vcf_table.txt.gz
url <- "https://zhangz05.u.hpc.mssm.edu/saasCNV/data/vcf_table.txt.gz"
tryCatch({download.file(url=url, destfile="vcf_table.txt.gz")
         }, error = function(e) {
          download.file(url=url, destfile="vcf_table.txt.gz", method="curl")
         })
## If download.file fails to download the data, please manually download it from the url.

vcf_table <- read.delim(file="vcf_table.txt.gz", as.is=TRUE)
head(vcf_table)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ vcf }
\keyword{ VCF }
