% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdb2wave.R
\docType{data}
\name{chop}
\alias{chop}
\title{Tentative parameters to reproduce Yu, Chi-Hua, et al.ACS nano (2019).
sustain, decay, offset, and volume for secondary structure elements
shorten prepares a wave based on secondary structure}
\format{An object of class \code{numeric} of length 1.}
\usage{
chop
}
\arguments{
\item{wave}{the wave in}

\item{ss}{the secondary structure label (s for sheet, h for helix, l for loop)}

\item{pad}{shift for L/R delay}

\item{volfactor}{scale volume}

\item{decayfactor}{scale decay (e.g. for last residue)}
}
\value{
the enveloped wave
}
\description{
Tentative parameters to reproduce Yu, Chi-Hua, et al.ACS nano (2019).
sustain, decay, offset, and volume for secondary structure elements
shorten prepares a wave based on secondary structure
}
\keyword{datasets}
