\name{diagnosis.cluster.plot}
\alias{diagnosis.cluster.plot}
\title{
Visualize Genome-Wide SCNA Profile in 2D Cluster Plot
}
\description{
An optional function to visualize genome-wide SCNA Profile in 
log2mBAF-log2ratio 2D cluster plot.
}
\usage{
diagnosis.cluster.plot(segs, chrs, min.snps, max.cex = 3, ref.num.probe = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segs}{a data frame containing segment location, summary statistics 
  	and SCNA status resulting from \code{\link{cnv.call}}.}
  \item{chrs}{the chromosomes to be visualized. For example, 1:22.
}
  \item{min.snps}{the minimum number of probes a segment span.
}
  \item{max.cex}{the maximum of cex a circle is associated with. See details.
}
  \item{ref.num.probe}{integer. The reference number of probes against which 
  	a segment is compared in order to determine the cex of the segment 
  	to be displayed. Default is \code{NULL}. If \code{NULL}, It will be automatically 
  	specified as 1/100 of the number of data points.
}
}
\details{
on the main log2mBAF-log2ratio panel, each circle corresponds to a segment, 
with the size reflecting the length of the segment; the color code is specified 
in legend; the dashed gray lines indicate the adjusted baselines. 
The side panels, corresponding to log2ratio and log2mBAF dimension respectively, 
show the distribution of median values of each segment weighted by its length.
}
\value{
An R plot will be generated.
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
\code{\link{joint.segmentation}}, \code{\link{cnv.call}}, 
\code{\link{diagnosis.seg.plot.chr}}, \code{\link{genome.wide.plot}}
}
\examples{
data(seq.data)
data(seq.cnv)

diagnosis.cluster.plot(segs=seq.cnv, 
                       chrs=sub("^chr","",unique(seq.cnv$chr)), 
                       min.snps=10, max.cex=3, ref.num.probe=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SCNA }
\keyword{ cluster }
\keyword{ visualization }
\keyword{ diagnosis }

