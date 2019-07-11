require(tidyverse)

require(bio3d)
require(seewave)
require(tuneR)

 
#' Reads WAV files for residues. 
#' data from SI of Yu, Chi-Hua, et al.ACS nano (2019).
#' Assumes alphabetical order matches order of dir
#' Uniforms to shortest lenght, normalises to 16 bit
#'    
#' @param path The directory with the wav files.
#' @return A function of the aa1 code that returns the numerical base-16 wave
#' @export
#' @examples
#' aas = get_aa_waves()
#' ala_wave = aas("A") 
#' 
get_aa_waves = function(path){
  files = dir(path)
  getwav = function(filen,chop=0){
    wave=tuneR::readWave(glue::glue("{wavpath}/{filen}"))
    wave_norm=tuneR::normalize(wave,"16",level=1)
    #print(wave_norm)
    if(chop>0){
      unlist(wave_norm@left)[1:chop]
    } else {
      unlist(wave_norm@left)
    }
  }
  # read once without chopping to get sample lengths, then read again chopping
  ws = lapply(files, getwav)
  chop = min(sapply(X = ws, length))
  ws = lapply(files, getwav, chop=chop)
  
  aa3 = c("ALA", "ARG", "ASN", "ASP", 
          "CYS", "GLU", "GLN", "GLY", 
          "HIS", "ILE", "LEU", "LYS", 
          "MET", "PHE", "PRO", "SER", 
          "THR", "TRP", "TYR", "VAL")
  
  getaa = function(label){
    i=which(aa3==bio3d::aa123(label))
    #print(bio3d::aa123(label))
    ws[[i]]
  }
  
  getaa
}



#' Add two waves, offsetting the starting point of the second from the end of the first
#' @param a The first wave (vector of numbers)
#' @param b The second wave (vector of numbers)
#' @param offset how many indices to move the second wave with respect to the end of the first. 
#'  Negative offset results in overlap
#' @return The vector of numbers of the joint waves
#' @export
#' @examples
#' aas = get_aa_waves()
#' AH_wave = over(aas("A"),aas("H"))
#'  
over = function(a,b, offset=0){
  la = length(a)
  lb = length(b)
  c1 = 0*(1:(la+lb+offset))
  c2 = 0*(1:(la+lb+offset))
  c1[1:la] = a
  c2[(la+offset+1):(la+offset+lb)]=b
  z=c1+c2
  z[z>(2^15-1)]=2^15-1
  z[z<(-2^15)]=-2^15
  z
}

#' Envelopes a wave with a simple envelope.
#' @param a The wave (vector of numbers)
#' @param sus How many sampling points to sustain (envelope = 1) from the first
#' @param decay sampling points scale for exponential decay after sustain
#' @param predecay second exponential summed (results in envelope > 1), for attack
#' @param w overall weight
#' @param pad shift the whole wave ahead of this many points (attempt to delay for location sensitivity)
#' @export
#' @return The vector of numbers of the enveloped wave
#' @examples
#' aas = get_aa_waves()
#' A_wave = taper(aas("A"),sus=10000, decay=100)
taper = function(a, sus, decay, predecay=0, w=1, pad=0){
  #print(c(sus,decay,w))
  la = length(a)
  enve = 0*(1:la)
  enve[pad:(sus+pad)]=1
  enve[(sus+pad):la] = exp(-((sus:(la-pad))-sus)/decay)
  if(predecay>0){
    enve[(sus+pad):la] = enve[(sus+pad):la] + exp(-((sus:(la-pad))-sus)/predecay)
  }
  a*enve*w
}



#' Tentative parameters to reproduce Yu, Chi-Hua, et al.ACS nano (2019).
#' sustain, decay, offset, and volume for secondary structure elements
#' shorten prepares a wave based on secondary structure
#' @param wave the wave in
#' @param ss the secondary structure label (s for sheet, h for helix, l for loop)
#' @param pad shift for L/R delay
#' @param volfactor scale volume
#' @param decayfactor scale decay (e.g. for last residue)
#' @return the enveloped wave
#' @export
#' @example 
#' shorten(aawav("A"),"h")
chop=99975
sustain  = list(`s`=100, `h`=100, `l`=100)
decay = list( `s`=8000, `h`=8000, `l`=8000)
offset = list(`s`= -round(chop-4*3500),  
              `h`=-round(chop-4*1750), 
              `l`=-round(chop-4*6000))
vol = list(`s`=1,`h`=.5,`l`=.25)

shorten = function(wave, ss, pad=0, volfactor=1, decayfactor=1){
  list(w=taper(wave,sustain[[ss]],
               decayfactor*decay[[ss]], predecay = 200, 
               w=vol[[ss]]*volfactor,pad=pad),
       off = offset[[ss]])
}



#' Read pdb extract sequence, call dssp for secondary structure, 
#' calculate distances to two points for sound source localisation
#' @param filename the pdb filename, only calphas are used
#' @param refL_resno residue number for left reference point
#' @param refR_resno residue number for rigth reference point
#' @export
#' @return a list with sequence, secondary structure and distances to R/L refs
get_structure_pdb = function(filename, 
                             refL_resno=NULL, 
                             refR_resno=NULL){
  
  coords = read.pdb(filename)
  atom.select(coords,"calpha")->indexca
  xxx=matrix(coords$xyz[indexca$xyz], ncol=3, byrow=TRUE)
  
  if(!is.null(refL_resno) & !is.null(refR_resno)){
    refL=atom.select(coords,"calpha", resno=refL_resno)
    refR=atom.select(coords,"calpha", resno=refR_resno)
    xxxL = matrix(coords$xyz[refL$xyz], ncol=3, byrow=TRUE)
    xxxR = matrix(coords$xyz[refR$xyz], ncol=3, byrow=TRUE)
    print(xxxL)
    print(xxxR)
    uuuL = t(t(xxx) - c(xxxL))
    uuuR = t(t(xxx) - c(xxxR))
    distL = sqrt(rowSums(uuuL^2))
    distR = sqrt(rowSums(uuuR^2))
  } else {
    distL = NULL
    distR = NULL
  }
  # plot(dist)
  seq2 = bio3d::aa321(coords$atom[indexca$atom,"resid"])
  ss =   pdb.dssp(coords)
  ses2 = ss$sse 
  ses3 = ses2
  ses3[ses2=="S"] = "l"
  ses3[ses2=="T"] = "l"
  ses3[ses2=="B"] = "l"
  ses3[ses2==" "] = "l"
  ses3[ses2=="E"] = "s"
  ses3[ses2=="H"] = "h"
  ses3[ses2=="G"] = "h"
  ses3[ses2=="I"] = "h"
  
  list(ses=ses3, seq=seq2, distR=distR, distL=distL)
}

#' Read sequence and seconday structure from csv file
#' @param filename the csv filename, two columns with no header for three-letter aa code and ss
#' @export
#' @return a list with sequence, secondary structure
get_structure_nopdb = function(filename){
  protein = read_csv(filename, col_names = c("aa3","ss")) %>% 
    mutate(aa1 = bio3d::aa321(aa3), ss = tolower(ss)) %>%
    replace_na(list(ss = "l"))
  seq = protein$aa1
  ses = protein$ss
  list(ses=ses,seq=seq)
}


#' Generates a tune from sequence, secondary structure
#' if distances to two points are provided, L/R waves are attenuated and delayed
#' to try and give stereo impression.
#' @param seq an array of 1-aa codes
#' @param ses an array of 1-letter lowercase DSSP codes
#' @param distL distances to left ref
#' @param distR distances to right ref
#' @export
#' @return stereo Wave object of the whole tune
#' 
gen_tune = function(seq, ses, distL=NULL, distR=NULL){
  if(!is.null(distL) & !is.null(distR)){
    Lp = seewave::attenuation(0, dref = 40, n=100)
    attL = (10^(Lp/10))[round(distL)+1]
    attR = (10^(Lp/10))[round(distR)+1]
    plot(attL,type='b', col='purple',pch=20)
    lines(attR,type='b', col='orange',pch=20)
    offL = round(distL)
    offR = round(distR)
    
  } else {
    attL = 0*(1:length(seq))+1
    attR = 0*(1:length(seq))+1
    offL = 0*(1:length(seq))
    offR = 0*(1:length(seq))
  }
  for (i in 1:length(seq)){ 
    #fac = facRL[[ses[i]]]
    fac = c(1,1)
    if(i==1){
      l1L = shorten(aa_waves(seq[i]),ses[i],pad=offL[i])
      l1R = shorten(aa_waves(seq[i]),ses[i],pad=offR[i])
      
      tuneL = attL[i]*fac[1]*l1L$w; 
      tuneR = attR[i]*fac[2]*l1R$w;
    } else {
      if(i == length(seq)) decayfactor=1.5 else decayfactor=1
      off = l1R$off
      l1L = shorten(aa_waves(seq[i]),ses[i],pad=offL[i],decayfactor=decayfactor)
      l1R = shorten(aa_waves(seq[i]),ses[i],pad=offR[i],decayfactor=decayfactor)
      tuneL = tuneL %>% over(attL[i]*fac[1]*l1L$w,off)
      tuneR = tuneR %>% over(attR[i]*fac[2]*l1R$w,off)
    }
  }
  
  allaaWstereo = tuneR::Wave(
    left=tuneL,
    right=tuneR,
    bit=16,samp.rate = 44100)
  allaaWstereo
  
}




