\name{reannotate.CNV.res}
\alias{reannotate.CNV.res}
\title{
Gene Annotation
}
\description{
An optional function to add gene annotation to each CNV segment.
}
\usage{
reannotate.CNV.res(res, gene, only.CNV = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{a data frame resultingfrom \code{\link{cnv.call}}.
}
  \item{gene}{a data frame containing gene annotation information.
}
  \item{only.CNV}{logical. If only segment assigned to gain/loss/LOH 
  	to be annotated and output. Default is \code{FALSE}.
}
}
\details{
The RefSeq gene annotation file can be downloaded from UCSC Genome Browser.
}
\value{
A gene annotation column have been add to the data frame resulting from 
\code{\link{cnv.call}}.
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
\author{
Zhongyang Zhang <zhongyang.zhang@mssm.edu>
}
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{joint.segmentation}}, \code{\link{cnv.call}}
}
\examples{

\dontrun{
## An example of RefSeq gene annotation file,
## the original version of which can be downloaded from UCSC Genome Browser
url <- "https://zhangz05.u.hpc.mssm.edu/saasCNV/data/refGene_hg19.txt.gz"
tryCatch({download.file(url=url, destfile="refGene_hg19.txt.gz")
         }, error = function(e) {
          download.file(url=url, destfile="refGene_hg19.txt.gz", method="curl")
         })
## If download.file fails to download the data, please manually download it from the url.

gene.anno <- read.delim(file="refGene_hg19.txt.gz", as.is=TRUE, comment.char="")
data(seq.cnv)
seq.cnv.anno <- reannotate.CNV.res(res=seq.cnv, gene=gene.anno, only.CNV=TRUE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ CNV }
\keyword{ gene }
\keyword{ annotation }
