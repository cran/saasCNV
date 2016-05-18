\name{joint.segmentation}
\alias{joint.segmentation}
\title{
Joint Segmentation on log2ratio and log2mBAF Dimensions
}
\description{
We employ the algorithm developed by (Zhang et al., 2010) to perform joint segmentation on log2ratio 
and log2mBAF dimensions. The function outputs the starting and ending points of each 
CNV segment as well as some summary statistics.
}
\usage{
joint.segmentation(data, min.snps = 10, global.pval.cutoff = 1e-04, 
    max.chpts = 30, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame containing log2ratio and log2mBAF data generated 
  	by \code{\link{cnv.data}}.
}
  \item{min.snps}{the minimum number of probes a segment needs to span.
}
  \item{global.pval.cutoff}{the p-value cut-off a (or a pair) of change points to
  	be determined as significant in each cycle of joint segmentation.
}
  \item{max.chpts}{the maximum number of change points to be detected for each 
  	chromosome.
}
  \item{verbose}{logical. If more details to be output. Default is \code{TRUE}.
}
}
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
\value{
A data frame containing the starting and ending points of each 
CNV segment as well as some summary statistics.
}
\references{
Zhang, N. R., Siegmund, D. O., Ji, H., Li, J. Z. (2010) 
Detecting simultaneous changepoints in multiple sequences. 
\emph{Biometrika}, \bold{97}:631--645.
}
\author{
Zhongyang Zhang <zhongyang.zhang@mssm.edu>
}
%% \note{
%%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cnv.data}}
}
\examples{
data(seq.data)

\dontrun{
seq.segs <- joint.segmentation(data=seq.data, min.snps=10, 
                               global.pval.cutoff=1e-4, max.chpts=30, 
                               verbose=TRUE)
}

## how the joint segmentation results look like
data(seq.segs)
head(seq.segs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ joint segmentation }
\keyword{ CBS }
\keyword{ segmentation }
