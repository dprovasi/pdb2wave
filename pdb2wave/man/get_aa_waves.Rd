% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pdb2wave.R
\name{get_aa_waves}
\alias{get_aa_waves}
\title{Reads WAV files for residues. 
data from SI of Yu, Chi-Hua, et al.ACS nano (2019).
Assumes alphabetical order matches order of dir
Uniforms to shortest lenght, normalises to 16 bit}
\usage{
get_aa_waves(path)
}
\arguments{
\item{path}{The directory with the wav files.}
}
\value{
A function of the aa1 code that returns the numerical base-16 wave
}
\description{
Reads WAV files for residues. 
data from SI of Yu, Chi-Hua, et al.ACS nano (2019).
Assumes alphabetical order matches order of dir
Uniforms to shortest lenght, normalises to 16 bit
}
\examples{
aas = get_aa_waves()
ala_wave = aas("A") 

}
