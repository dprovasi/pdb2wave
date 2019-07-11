## ------------------------------------------------------------------------
require(pdb2wave)
require(tidyverse)

## ------------------------------------------------------------------------
setwd("~/Desktop/protein_sounf")
wavpath = "FINAL SI/DFT BASED AA SOUNDS/"
aa_waves = get_aa_waves(wavpath)

