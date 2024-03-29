% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdb2wave.R
\name{gen_tune}
\alias{gen_tune}
\title{Generates a tune from sequence, secondary structure
if distances to two points are provided, L/R waves are attenuated and delayed
to try and give stereo impression.}
\usage{
gen_tune(seq, ses, distL = NULL, distR = NULL)
}
\arguments{
\item{seq}{an array of 1-aa codes}

\item{ses}{an array of 1-letter lowercase DSSP codes}

\item{distL}{distances to left ref}

\item{distR}{distances to right ref}
}
\value{
stereo Wave object of the whole tune
}
\description{
Generates a tune from sequence, secondary structure
if distances to two points are provided, L/R waves are attenuated and delayed
to try and give stereo impression.
}
