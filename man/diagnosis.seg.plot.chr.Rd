\name{diagnosis.seg.plot.chr}
\alias{diagnosis.seg.plot.chr}
\title{
Visualize Segmentation Results for Diagnosis
}
\description{
The results from joint segmentation and segments merging are visualized 
for the specified choromosome.
}
\usage{
diagnosis.seg.plot.chr(data, segs, sample.id = "Sample", chr = 1, cex = 0.3)
}
\arguments{
  \item{data}{a data frame containing log2ratio and log2mBAF data generated 
  	by \code{\link{cnv.data}}.
}
  \item{segs}{a data frame containing segment locations and summary statistics 
  	resulting from \code{\link{joint.segmentation}} or \code{\link{merging.segments}}.
}
  \item{sample.id}{sample ID to be displayed in the title of the plot.
}
  \item{chr}{the chromosome number (e.g. 1) to be visualized.
}
  \item{cex}{a numerical value giving the amount by which plotting text and 
  	symbols should be magnified relative to the default. It can be adjusted 
  	in order to make the plot legible.
}
}
%% \details{
%%  ~~ If necessary, more details than the description above ~~
%% }
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
\code{\link{joint.segmentation}}, \code{\link{merging.segments}}, \code{\link{cnv.data}}
}
\examples{
## visual diagnosis of joint segmentation results
data(seq.data)
data(seq.segs)
diagnosis.seg.plot.chr(data=seq.data, segs=seq.segs, 
                       sample.id="Joint Segmentation", 
                       chr=1, cex=0.3)

## visual diagnosis of results from merging step
data(seq.segs.merge)
diagnosis.seg.plot.chr(data=seq.data, segs=seq.segs.merge, 
                       sample.id="After Segments Merging Step", 
                       chr=1, cex=0.3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ segmentation }
\keyword{ diagnosis }
