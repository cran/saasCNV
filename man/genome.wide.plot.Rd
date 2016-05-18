\name{genome.wide.plot}
\alias{genome.wide.plot}
\title{
Visualize Genome-Wide SCNA Profile
}
\description{
An optional function to visualize genome-wide SCNA Profile.
}
\usage{
genome.wide.plot(data, segs, sample.id, chrs, cex = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a data frame containing log2ratio and log2mBAF data generated 
  	by \code{\link{cnv.data}}.
}
  \item{segs}{a data frame containing segment location, summary statistics 
  	and SCNA status resulting from \code{\link{cnv.call}}.
}
  \item{sample.id}{sample ID to be displayed in the title of the plot.
}
  \item{chrs}{the chromosomes to be visualized. For example, 1:22.
}
  \item{cex}{a numerical value giving the amount by which plotting text and 
  	symbols should be magnified relative to the default. It can be adjusted 
  	in order to make the plot legible.
}
}
\details{On the top panel, the log2ratio signal is plotted against 
	chromosomal position and on the panels blow, the log2mBAF, 
	tumor mBAF, and normal mBAF signals. The dots, each representing 
	a probe data point, are colored alternately to distinguish chromosomes. 
	The segments, each representing a DNA segment resulting from the 
	joint segmentation, are colored based on inferred copy number status. 
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
\code{\link{diagnosis.seg.plot.chr}}, \code{\link{diagnosis.cluster.plot}}
}
\examples{
data(seq.data)
data(seq.cnv)

genome.wide.plot(data=seq.data, segs=seq.cnv, 
                 sample.id="PT116", 
                 chrs=sub("^chr","",unique(seq.cnv$chr)), 
                 cex=0.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SCNA }
\keyword{ visualization }
\keyword{ diagnosis }


