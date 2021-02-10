# PolyHaplotyper

# A haploblock (or block) is a set of bi-allelic markers (SNPs) that are assumed
# to have no recombinations within the population studied. A good example is
# (the set of SNPs within) a contig.

# A haplotype is an allele of a haploblock: a unique combination of marker
# alleles.

# PolyHaplotyper aims to infer the genetic composition of individuals in
# terms of haplotype dosages for a haploblock, based on marker dosages, for any
# ploidy level.

# PolyHaplotyper should not be confused with PediHaplotyper, which is aimed
# at diploid individuals that are all part of a pedigree, and where the
# (bi- or multi-allelic) marker alleles have already been phased.
# In contrast, PolyHaplotyper can handle any ploidy level and uses unphased
# bi-allelic marker genotypes (i.e. marker allele dosages). It can work with
# unrelated individuals and/or independent or linked FS families but
# (so far) does not use or need an entire pedigree, except (optionally) to
# check for haplotype assignments that violate Mendelian inheritance, after the
# haplotyping has been completed.

# A note on the representation of genotypes:
# 1. genotypes in terms of biallelic marker dosages: these can be represented as
#    a. a vector with one element for each marker in the haploblock,
#       with the dosages of one of the marker alleles. Each marker dosage
#       is in 0:ploidy, or NA. This type of data often has mrkdos in the name
#    b. a single integer number: the index of the marker dosage combination
#       in an ordered list of all possible marker dosage combinations. This
#       is NA when at least one of the marker dosages in the haploblock is NA.
#       The index is 1-based and ranges from 1 to (ploidy+1)^nmrk, where nmrk
#       is the number of markers in the haploblock.
#       This index number is often called did (dosage ID) or mrkdid.
# 2. genotypes in terms of haplotype dosages: these can be represented as
#    a. a vector with one element for each of the 2^nmrk possible haplotypes
#       (where nmrk is the number of markers in the haploblock), with the
#       dosages of the haplotypes. The dosages sum to ploidy or all dosages
#       are NA. This type of data often has hapdos or had in the name.
#    b. an ordered vector of length ploidy with the haplotype numbers present
#       in the genotype (e.g. c(1, 5, 8, 8) : 1 copy of haplotypes 1 and 5,
#       and 2 copies of haplotype 8). This type of data often has hapcomb or
#       hac in the name (but comb is also used for other types of
#       combinations, e.g. different combinations of parental genotypes)
#    c. a single integer number: the index of the haplotype combination
#       in an ordered list of all possible haplotype combinations. The index
#       is 1-based and ranges from 1 to choose(nhap+ploidy-1, ploidy), where
#       nhap is the number of possible haplotypes: nhap = 2^nmrk.
#       This index number is often called hapindex or hapix.
#       This indexing fails with more than 8 markers (tetraploid) or 6 markers
#       (hexaploid) because the number of possible haplotype combinations
#       exceeds the capacity of a 32-bit integer.
#    d. an index to the hapcomb in ahc(complete)list: the column number in
#       the matrix at ahccompletelist[[nmrk]][[mrkdid]] (with mrkdid a numeric
#       or integer) or ahclist[[nmrk]][[mrkdid]] (with mrkdid a string).
#       This requires the ploidy, nmrk (the number of markers in the
#       haploblock) and mrkdid (the marker dosage ID, see above).
#       The hapcomb is obtained by
#       getAllHapcomb(did=mrkdid, allhap=allhap, ploidy=ploidy)[, colnum]

#'@importFrom utils combn flush.console read.table write.table
#'@importFrom stats chisq.test pchisq pbinom qbinom

# documentation of the objects in PolyHaplotyper_demo.RData:

#'@title dosages of SNP alleles
#'@description A data.frame with the allele dosages for 30 SNPs in 661
#'samples; used in vignette, stored in PolyHaplotyper_demo.RData
#'@docType data
#'@keywords datasets
#'@name demo_snpdos
#'@format a data.frame with 30 rows (SNPs) and 663 columns. The first
#'two columns contain the SNP name and contig; the remaining columns
#'contain the allele dosages (integers in 0..6 or NA) of the SNPs in
#'661 samples
NULL

#'@title pedigree
#'@description A data.frame with a pedigree; used in vignette, stored in
#'PolyHaplotyper_demo.RData
#'@docType data
#'@keywords datasets
#'@name demo_ped
#'@format a data.frame with 661 rows (samples) and
#'4 columns: genotype, mother, father, sample_nr
NULL

# documentation of the objects in PolyHaplotyper_small.RData:

#'@title dosages of SNP alleles
#'@description A matrix with the allele dosages for 8 SNPs in 661
#'samples; used in manuals, stored in PolyHaplotyper_small.RData
#'@docType data
#'@keywords datasets
#'@name phdos
#'@format a data.frame with 30 rows (SNPs) and 663 columns. The first
#'two columns contain the SNP name and contig; the remaining columns
#'contain the allele dosages (integers in 0..6 or NA) of the SNPs in
#'661 samples
NULL

#'@title pedigree
#'@description A data.frame with a pedigree; used in manuals, stored in
#'PolyHaplotyper_small.RData
#'@docType data
#'@keywords datasets
#'@name phped
#'@format a data.frame with 661 rows (samples) and
#'4 columns: genotype, mother, father, sample_nr
NULL

#'@title List of markers per haploblock
#'@description A list with for each haploblock the names of the markers it
#'contains; used in manuals, stored in PolyHaplotyper_small.RData
#'@docType data
#'@keywords datasets
#'@name phblocks
#'@format a list with 2 character vectors
NULL

#'@title parents of FS families
#'@description A matrix with the names of the parents of each FS family;
#'used in manuals, stored in PolyHaplotyper_small.RData
#'@docType data
#'@keywords datasets
#'@name phpar
#'@format a matrix with 4 rows, 1 per FS family, and 2 columns, 1 for each
#'parent
NULL

#'@title members of FS families
#'@description A list with for each FS family the names of the individuals it
#'contains; used in manuals, stored in PolyHaplotyper_small.RData
#'@docType data
#'@keywords datasets
#'@name phFS
#'@format a list with 4 vectors, each with the (here numeric) names of
#'the FS members
NULL

#'@title haplotyping results
#'@description A list with the results of function inferHaplotypes; used in
#'manuals, stored in PolyHaplotyper_small.RData
#'@docType data
#'@keywords datasets
#'@name phresults
#'@format a list with 2 elements, one per haploblock; each element itself
#'a list with multiple components
NULL

# Start of function definitions

#'@title get all haplotypes for the given markers
#'@description Given a set of bi-allelic (SNP) marker names, generate all
#'possible haplotypes
#'@usage allHaplotypes(mrknames)
#'@param mrknames the names of the (bi-allelic) markers in the
#'haploblock (contig)
#'@return a matrix with markers in columns and all possible (2 ^ nmrk)
#'haplotypes in rows. 0: haplotype contains the non-dosage-counted marker
#'allele (the reference allele); 1: haplotype contains the dosage-counted
#'(alternative) marker allele. The colnames are the marker names.
#'@examples
#'# show the 8 possible haplotypes with 3 bi-allelic markers:
#'allHaplotypes(mrknames=c("mrkA", "mrkB", "mrkC"))
#'@export
allHaplotypes <- function(mrknames) {
  lst <- list()
  for (s in 1:length(mrknames)) {
    lst[[s]] <- 0:1
  }
  eg <- as.matrix(expand.grid(lst, stringsAsFactors=FALSE))
  eg <- eg[, ncol(eg):1, drop=FALSE]
  colnames(eg) <- mrknames
  eg
} #allHaplotypes

getahcdir <- function(ahcdir=NULL) {
  # for internal use only; get ahcdir from (1) param ahcdir or (2) working dir
  # result ends always with "/" so a filename can be pasted onto it
  if (is.null(ahcdir)) ahcdir <- getwd() else {
    if (!is.character(ahcdir) || length(ahcdir) != 1)
      stop(paste("invalid ahcdir:", ahcdir))
    if (ahcdir=="") ahcdir <- getwd()
  }
  if (!dir.exists(ahcdir))
    stop(paste("non-existing ahcdir:", ahcdir))
  if (substr(ahcdir, nchar(ahcdir), nchar(ahcdir)) != "/")
    ahcdir <- paste0(ahcdir, "/")
  ahcdir
} #getahcdir

loadAllHapcombLists <- function(ploidy, nmrk, ahcdir) {
  # This function is called only at start of inferHaplotypes
  # ploidy: 1 integer
  # nmrk: 1 integer, the number of markers in the haploblock
  # result: a list with elements ahccompletelist, ahclist, complete_nmrk,
  #         ploidy, ahcdir

  checklist <- function(lst, ploidy, namespresent) {
    if (!is.list(lst) || length(lst) < 1) return(FALSE)
    if (!all(sapply(lst, class) == "list")) return(FALSE)
    if (length(lst[[length(lst)]]) == 0) return(FALSE) # in highest nmrk at least one mrkdid should be present
    # cl <- lapply(lst, FUN=function(x) sapply(x, class))
    # if (!all(unlist(cl) == "matrix")) return(FALSE)
    # due to change in class of matrix in R version 4:
    cl <- lapply(lst, FUN=function(x) sapply(x, inherits, "matrix"))
    if (!all(unlist(cl))) return(FALSE)
    cl <- lapply(lst, FUN=function(x) sapply(x, nrow))
    if (!all(unlist(cl) == ploidy)) return(FALSE)
    lnames <- names(lst[[length(lst)]])
    if (xor(!namespresent, (is.null(lnames) || any(lnames == "")))) return(FALSE)
    return(TRUE)
  } #checklist

  loadlist <- function(listname, ploidy, ahcdir, namespresent, message) {
    # try to load:
    if (message) cat(paste("loading ", listname, "..."))
    x <- tryCatch({
      fname <- paste0(ahcdir, listname, "_", ploidy, "x.RData")
      suppressWarnings(n <- load(fname))
      if (!(listname %in% n))
        return(paste("file", fname, "does not contain", listname))
    }, error=function(e) "error")
    if (identical(x, "error")) {
      if (message) cat (" not found!\n")
      return("") # not an error but no list file
    } else {
      #checks:
      if (!checklist(get(listname), ploidy, namespresent))
        return(paste(listname, " has incorrect format!\n"))
      if (message) cat (" done!\n")
      return(get(listname))
    }
  } #loadlist


  ahcinfo <- list(ahccompletelist=NULL, ahclist=NULL, complete_nmrk=0,
                  ploidy=ploidy, ahcdir=getahcdir(ahcdir))
  if (exists("ahccompletelist") && # in .GlobalEnv, loaded by user
      checklist(lst=ahccompletelist, ploidy=ploidy, namespresent=FALSE)) {
    ahcinfo$ahccompletelist <- ahccompletelist
  } else {
    ahccompletelist <-
      loadlist(listname="ahccompletelist", ploidy=ploidy, ahcdir=ahcinfo$ahcdir,
               namespresent=FALSE, message=TRUE)
    if (is.list(ahccompletelist)) {
      ahcinfo$ahccompletelist <- ahccompletelist
    } else if (ahccompletelist != "") stop(ahccompletelist)
  }
  if (!is.null(ahcinfo$ahccompletelist))
    ahcinfo$complete_nmrk <- length(ahcinfo$ahccompletelist)
  if (ahcinfo$complete_nmrk < nmrk) {
    # we need an additional ahclist for the mrkdids with larger nmrk
    if (exists("ahclist") && # in .GlobalEnv, loaded by user
        checklist(ahclist, ploidy, namespresent=TRUE)) {
      ahcinfo$ahclist <- ahclist
    } else {
      ahclist <- loadlist("ahclist", ploidy, ahcinfo$ahcdir,
                          namespresent=TRUE, message=FALSE)
      if (is.list(ahclist)) {
        ahcinfo$ahclist <- ahclist
      } else {
        if (ahclist != "") stop(ahclist)
        ahcinfo$ahclist <- list() # new empty ahclist
      }
    }
  }
  ahcinfo
} # loadAllHapcombLists

# not exported function (relies on ahcinfo)
#_#@title get all haplotype combinations that result in the given marker
#_#'dosages
#_#'@description get all haplotype combinations that result in the given
#_#'marker dosages
#_#'@usage getAllHapcomb(mrkDosage=NULL, mrkdid=NULL, nmrk, ahcinfo)
#_#'@param mrkDosage an integer vector with the dosages of each marker in one
#_#'individual; all dosages must be in 0:ahcinfo$ploidy or NA
#_#'Either mrkDosage or mrkdid must be specified, not both
#_#'@param mrkdid marker dosage ID: in combination with nmrk and ploidy specifying
#_#'the mrkDosage. Can be NA or an integer, numeric or character which must
#_#'specify an integer > 0
#_#'Either mrkDosage or mrkdid must be specified, not both
#_#'@param nmrk the number of markers (should match mrkDosage if that is specified)
#_#'@param ahcinfo a list as returned by loadAllHapcombLists
#_#'@details Each column of the return value represents one combination of
#_#'haplotypes that yields the observed marker dosages. If any of the marker
#_#'dosages, or mrkdid, are NA a matrix with 0 columns is returned.\cr
#_#'This function makes use of precalculated lists ahcinfo$ahccompletelist
#_#'and ahcinfo$ahclist; the mrkdid must be present in one of them.
#_#'@return an integer matrix with <ploidy> rows and one column for each haplotype
#_#'combination that matches the marker dosages; the haplotype numbers
#_#'within a column are in increasing order
getAllHapcomb <- function(mrkDosage=NULL, mrkdid=NULL, nmrk, ahcinfo) {
  # mrkDosage has already been checked at start of inferHaplotypes,
  # mrkdid are computed
  # so we don't check either, only find out of did or mrkDosage specified:
  if (is.null(mrkdid)) {
    if (anyNA(mrkDosage)) return(matrix(integer(0), nrow=ahcinfo$ploidy))
    mrkdid <- mrkdidfun(mrkDosage, ploidy=ahcinfo$ploidy)
  }
  if (is.na(mrkdid)) return(matrix(integer(0), nrow=ahcinfo$ploidy))
  if (nmrk <= ahcinfo$complete_nmrk)
    return(ahcinfo$ahccompletelist[[nmrk]][[as.integer(mrkdid)]])
  #debug check:
  # if (nmrk>length(ahcinfo$ahclist) ||
  #     !(mrkdid %in% names(ahcinfo$ahclist[[nmrk]])))
  #   stop("mrkdid not in ahc(complete)list in getAllHaplocomb")
  return(ahcinfo$ahclist[[nmrk]][[as.character(mrkdid)]])
} # getAllHapcomb

calcAllhapcomb4mrkdid <- function(mrkDosage=NULL, mrkdid=NULL, ploidy, allhap) {
  # mrkDosage: an integer vector with the dosages (0:ploidy) of nmrk markers
  # mrkdid: 1 integer (marker dosage id); if NULL, mrkDosage is used
  # ploidy: 1 integer
  # allhap: a matrix as returned by allHaplotypes, for nmrk markers
  # in the same order as mrkDosage
  # result: a matrix with <ploidy> rows and one column for each haplotype
  #         combination that matches mrkdid
  #
  # check if mrkdid or mrkDosage is valid (we assume ploidy and nmrk are valid)
  nmrk <- ncol(allhap)
  if (is.null(mrkdid)) {
    # check mrkDosage:
    if (length(mrkDosage) != nmrk || !all(mrkDosage %in% 0:ploidy))
      stop("invalid mrkDosage")
  } else {
    if (length(mrkdid)!=1 || is.na(mrkdid) || as.integer(mrkdid) != mrkdid ||
        mrkdid<1 || mrkdid > (ploidy+1)^nmrk)
      stop("invalid mrkdid")
    mrkDosage <- mrkdid2mrkdos(mrkdid, nmrk, ploidy)
  }
  #TODO: the next algorithm is less efficient if the SNP dosages are high than if they are
  #low (because with low dosages a haplotype combination will sooner exceed the dosage of
  #at least one of the SNPs). Therefore in these cases (with average SNP dosage > ploidy/2)
  #it may be more efficient to find a solution for the did with all SNP dosages set to
  #<ploidy> - actual dosage, and afterwards replacing the haplotypes by their inverses.
  allmat1 <- matrix(NA_integer_, nrow=ploidy, ncol=0) #matrix with haplotype sets in columns
  allset <- rep(nrow(allhap), ploidy) #current haplotype set: a full set has length ploidy
  while(length(allset) > 0) {
    # if the last allele is 0, change the next:
    if (length(allset) > 0 && allset[length(allset)] <= 0) {
      if (length(allset) == 1) {
        allset <- integer(0)
      } else allset <- allset[1:(length(allset) - 1)]
      if (length(allset) > 0)
        allset[length(allset)] <- allset[length(allset)] - 1
    } else if (all(colSums(allhap[allset,, drop=FALSE]) == mrkDosage)) {
      #solution found; fill up with nr 1 (dosage 0 for all markers) and add to allmat1
      allset <- c(allset, rep(1, ploidy-length(allset)))
      allmat1 <- cbind(allmat1, allset)
      allset[ploidy] <- 0 #make sure to change the next level
    } else {
      if (length(allset) < ploidy  &&
          !any(colSums(allhap[allset,, drop=FALSE]) > mrkDosage) &&
          allset[length(allset)] > 1) {
        # only add the next allele if none of the dosages too high;
        # adding allele 1 (all zeroes) won't lead to a solution
        allset <- c(allset, allset[length(allset)]) #allele nrs never increase
      } else {
        # change last allele
        allset[length(allset)] <- allset[length(allset)] - 1
      }
    }
  }
  if (ncol(allmat1) == 0) stop("getAllHapcomb: 0 solutions found")
  colnames(allmat1) <- NULL
  #reverse the row order:
  allmat1[nrow(allmat1):1,, drop=FALSE]
} # calcAllhapcomb4mrkdid

getHapcombCount_1mrk <- function(mrkdosage, ploidy, nhap) {
  # mrkdosage: the dosage of the marker (integer of length 1)
  # ploidy
  # nhap = total nr of different haplotypes for the haploblock: 2 ^ nmrk
  #returns the number of haplotype combinations for the given dosage at one
  #marker
  #currently not used
  d <- mrkdosage
  e <- ploidy - d
  h <- nhap %/% 2
  if (d == ploidy) C0 <- 1 else
    C0 <- cumprod(seq(h, h+e-1))[e] / cumprod(seq_len(e))[e]
  if (d == 0) C1 <- 1 else
    C1 <- cumprod(seq(h, h+d-1))[d] / cumprod(seq_len(d))[d]
  C0 * C1
} #getHapcombCount_1mrk

#'@title calculate the total nr of possible haplotype combinations
#'@description calculate the total nr of possible haplotype combinations
#'@usage totHapcombCount(ploidy, nmrk)
#'@param ploidy a vector of 1 or more ploidy levels
#'@param nmrk a vector of 1 or more numbers of markers per haploblock
#'@details nmrk is used to calculate the total number of possible haplotypes
#'= 2^nmrk. The shorter vector of ploidy and nmrk is recycled.
#'@return a vector with the number of possible haplotype combinations for
#'each pair of ploidy and nmrk values
#'@examples
#'totHapcombCount(ploidy=4, nmrk=c(1:8))
#'@export
totHapcombCount <- function(ploidy, nmrk) {
  nhap <- 2^nmrk
  choose(nhap + ploidy - 1, ploidy)
} #totHapcombCount

# not exported function - users should use build_ahccompletelist()
#_#'@title generate all haplotype combinations
#_#'@description generate a list which contains for each marker dosage combination
#_#'at a given ploidy all matching haplotype combinations
#_#'@usage completeAllHaploComb(ploidy, nmrk, savesec=1800, printsec=60,
#_#'fname=NULL)
#_#'@param ploidy ploidy (may be even or odd)
#_#'@param nmrk number of markers in the haploblock
#_#'@param savesec default 1800: number of seconds between successive saves of
#_#'the intermediate results
#_#'@param printsec default 60: number of seconds between printout of
#_#'the current set of haplotypes (the last two are always equal at the time of
#_#'printing). NA or 0 suppresses printing.
#_#'@param fname filename (to which the extension .RData will be added) to store
#_#'the results. Intermediate saves will go to fname + the current set of
#_#'haplotypes; these intermediate files are temporary and will be deleted.\cr
#_#'If NULL (default), fname will be set to e.g. ahc_4x_nmrk6, where 4 is the
#_#'ploidy and nmrk = 6
#_#'@details This function is used by build_ahccompletelist and would not
#_#'normally be called by the user
#_#'@examples
#_#'# Generate an ahccompletelist for a given ploidy, for 1 .. nmrk markers:
#_#'\dontrun{
#_#'ahccompletelist <- list(1:nmrk)
#_#'for (n in seq_len(nmrk)) {
#_#'  ahccompletelist[[n]] <-
#_#'    completeAllHaploComb(ploidy=ploidy, nmrk=n)
#_#'  save(ahccompletelist, file=paste0("ahccompletelist_", ploidy, "x.RData"))
#_#'}
#_#'Note that this takes lots of time and memory already for
#_#'ploidy 4, nmrk = 7 and ploidy 6, nmrk = 5
#_#'}
#_#'@return a list with for each mrkdid (each combination of marker dosages) a
#_#'matrix with one row per haplotype and one column per haplotype combination,
#_#'containing the dosages of the haplotypes in each haplotype combination.
completeAllHaploComb <- function(ploidy, nmrk, savesec=1800, printsec=60,
                                 fname=NULL) {
  if (is.null(fname)) fname <- paste0("ahc_", ploidy, "x_nmrk", nmrk)
  if (is.null(printsec) || length(printsec) !=1 || is.na(as.numeric(printsec)) ||
      printsec == 0) printsec <- NA

  ndids <- as.integer((ploidy + 1) ^ nmrk + 0.001)
  ahc <- list()
  for (did in ndids:1) ahc[[did]] <- matrix(NA_integer_, nrow=ploidy, ncol=1)
  lasthap <- rep(0, ndids) # last stored haplotype column for each did
  allhap <- allHaplotypes(paste0("m", 1:nmrk))
  lastfile <- ""
  lastsave <- proc.time()
  lastprint <- lastsave - printsec - 1
  # generate all hapcombs in turn and calculate their mrkdids:
  starttime <- proc.time()
  nhap <- nrow(allhap)
  allset <- rep(nhap, ploidy) #current haplotype set: a full set has length ploidy
  while(TRUE) {
    mrkdos <- colSums(allhap[allset,, drop=FALSE])
    did <- mrkdidfun(mrkdos, ploidy=ploidy)
    lasthap[did] <- lasthap[did] + 1
    if ((nc<-ncol(ahc[[did]])) < lasthap[did]) {
      #we don't want to reallocate for each new column but increase size
      #when needed by 50%
      extmat <- matrix(NA_integer_, nrow=ploidy, ncol= nc %/% 2 + 2)
      ahc[[did]] <- cbind(ahc[[did]], extmat)
    }
    ahc[[did]][,lasthap[did]] <- allset
    # get the next allset:
    len <- ploidy
    allset[len] <- allset[len] - 1
    while (len > 0 && allset[len] <= 0) {
      len <- len - 1
      if (len > 0) allset[len] <- allset[len] - 1
    }
    if (len <= 0) break
    if (len < ploidy) {
      allset[(len+1):ploidy] <- allset[len] #nhap
      #progress and save: only check when len < ploidy
      if (!is.na(printsec) && (proc.time()-lastprint)[3] > printsec) {
        cat(paste0("ploidy: ", ploidy, " nmrk: ", nmrk,
                   " allset: ", paste(allset, collapse=" "), " sec: ",
                   paste(round((proc.time()-starttime)[1:3], 3), collapse=" "), "\n"))
        lastprint <- proc.time()
        if ((lastsave-lastprint)[3] > savesec) {
          save(ahc, file=paste0(fname, "_", paste(allset, collapse="-"),
                                ".RData"),
               version=2)
          if (file.exists(lastfile)) file.remove(lastfile)
          lastfile <- fname
          lastsave <- proc.time()
        }
      }
    } # len < ploidy
  }
  # all elements of ahc are now matrices with one solution per column,
  # in the form of sets like 1-1-2-4 (i.e. all occurring haplotypes listed).
  # We still only need to cut off the unused but allocated columns
  # of the matrices if any; and also we reverse the order of the rows, so that
  # the haplotype numbers are increasing rather than decreasing:
  for (d in seq_along(ahc)) {
    ahc[[d]] <- ahc[[d]][ploidy:1, 1:lasthap[d], drop=FALSE]
  }
  save(ahc, file=paste0(fname, ".RData"), version=2)
  if (file.exists(lastfile)) file.remove(lastfile)
  invisible(ahc)
} #completeAllHaploComb

#'@title generate a list with all haplotype combinations
#'@description generate a list which contains for each marker dosage combination
#'at a given ploidy all matching haplotype combinations
#'@usage build_ahccompletelist(ploidy, maxmrk, savesec=1800, printsec=300,
#'overwrite, shorten=FALSE)
#'@param ploidy ploidy; may be even or odd, but inferHaplotypes only works
#'with even ploidy if FS families are present, so building an ahccompletelist
#'for odd ploidy is probably not useful
#'@param maxmrk the list will countain all haplotype combinations for
#'haploblocks of 1...maxmrk markers
#'@param savesec default 1800: number of seconds between successive saves of
#'the intermediate results
#'@param printsec default 300: number of seconds between printout of
#'the current set of haplotypes. NA or 0 suppresses printing.
#'@param overwrite the new file ahccompletelist_<ploidy>x.RData is written in
#'the working directory. If a file with that name already exists and would be
#'changed, this is aborted and the old file maintained if overwrite is FALSE
#'@param shorten if file ahccompletelist_<ploidy>x.RData already exists and
#'covers haploblocks of lengths longer than maxmrk,
#'this file is left unchanged if shorten is FALSE (default), but is shortened
#'to maxmrk if shorten is TRUE
#'@details
#'An ahccompletelist reduces the processing time of inferHaplotypes enormously
#'but takes a long time to build. This should therefore be done when
#'PolyHaplotyper will be used multiple times.\cr
#'If an ahccompletelist file already exists in the working directory, this file is
#'used as starting point; if maxmrk is larger than the length of the existing
#'list only the additional items are calculated.\cr
#'If this function crashes (due to exceeding the memory limits of the computer
#'or of R) the file ahccompletelist_<ploidy>x.RData in the working directory
#'contains an intact version of the ahccompletelist for the last finished
#'marker count. Also, the temporary file
#'"buildAHCcompletelist_<datetime>"_"<nmrk>.RData" contains a part of the
#'list item for the number of markers being calculated at that moment.\cr
#'Note that this function takes lots of time and memory already for
#'ploidy 4 / maxmrk 7 and ploidy 6 / maxmrk 5, and is not practicable above
#'ploidy 4 / maxmrk 8, ploidy 6 / maxmrk 6, ploidy 8 / maxmrk 5. Also, the
#'files are large and take up considerable memory, and most of the size is
#'used for the largest haploblock size. Therefore maxmrk should
#'not be taken larger than necessary.
#'@return NULL (invisible); the actual result is a file in the working directory
#'named ahccompletelist_<ploidy>x.RData
#'@examples
#'\donttest{
#'# this example will create a file ahccompletelist_4x.RData
#'# in the working directory, for ploidy=4x and 1 - 4 markers:
#'build_ahccompletelist(ploidy=4, maxmrk=4, overwrite=FALSE)
#'}
#'@export
build_ahccompletelist <- function(ploidy, maxmrk, savesec=1800, printsec=300,
                                 overwrite, shorten=FALSE) {
  fname <- paste0("ahccompletelist_", ploidy, "x.RData")
  tmpfname <- paste0("buildAHCcompletelist_", gsub("[^0-9]", "", Sys.time()), "_")
  ahccompletelist <- NULL
  if (missing(overwrite) || length(overwrite) != 1 || is.na(overwrite) ||
      !is.logical(overwrite))
    stop("overwrite must be FALSE or TRUE")
  if (file.exists(fname)) { # in working directory
    x <- load(fname) # contains ahccompletelist
    if (length(x) != 1 || x !="ahccompletelist")
      stop(paste("the existing file", fname, "does not contain (only) ahccompletelist"))
    if (length(ahccompletelist) == maxmrk ||
        (length(ahccompletelist) > maxmrk && !shorten)) {
      # file in working dir already ok
      return(invisible(NULL))
    }
    if (!overwrite)
      stop(paste("file", fname, "already exists in working directory; aborted"))
  } else {
    ahccompletelist <- vector(0, mode="list")
  }
  newmaxmrk <- ifelse(shorten, maxmrk, length(ahccompletelist))
  if (length(ahccompletelist) >= maxmrk) {
    ahccompletelist <- ahccompletelist[1:newmaxmrk]
    if (file.exists(fname)) file.rename(from=fname, to=paste0(fname,".backup"))
    save(ahccompletelist, file=fname, version=2)
    if (file.exists(paste0(fname,".backup"))) file.remove(paste0(fname,".backup"))
  } else {
    if (totHapcombCount(ploidy=ploidy, nmrk=maxmrk) > 2e9) {
      ans <- readline("this will result in a very large list and may crash the computer; proceed (y/n)? ")
      if (!(ans %in% c("Y", "y")))
        stop("aborted")
    }
    for (nmrk in (length(ahccompletelist)+1):maxmrk) {
      cat(paste0("starting block size ", nmrk, " markers\n"))
      if (file.exists(fname)) file.rename(from=fname, to=paste0(fname,".backup"))
      ahccompletelist[[nmrk]] <-
        completeAllHaploComb(ploidy=ploidy, nmrk=nmrk, savesec=savesec,
                             printsec=printsec,
                             fname=paste0(tmpfname, nmrk))
      save(ahccompletelist, file=fname, version=2)
      cat(paste("saved", fname, "for block sizes up to", nmrk, "markers\n"))
      fn <- paste0(tmpfname, nmrk, ".RData")
      if (file.exists(fn)) file.remove(fn) #left by completeAllHaploComb
      fn <- paste0(fname,".backup")
      if (file.exists(fn)) file.remove(fn)
    }
  }
  invisible(NULL)
} #buildAHCcompletelist

#'@title calculate the marker dosages resulting from haplotype dosage
#'combinations
#'@description calculate the marker dosages resulting from haplotype dosage
#'combinations
#'@usage hapdos2mrkdos(hapdos, allhap)
#'@param hapdos a matrix with one column per combination of haplotypes
#'and one row for each possible haplotype (corresponding to the rows of allhap)
#'with dosage of the haplotypes in each combination. A vector is interpreted
#'as a one-column matrix; all columns must sum to ploidy. The rownames of
#'the matrix (or names of the vector) must contain the haplotype numbers
#'@param allhap a matrix as returned by allHaplotypes
#'@return a matrix with columns corresponding to the columns of hapdos
#'and one row for each marker, with the dosages of each marker in each combination;
#'colnames are the mrkdids (marker dosage IDs), rownames are the marker names taken
#'from allhap
#'@details if hapdos contains NA values, all values in the corresponding
#'column of the result will also be NA
#'@examples
#'# get a matrix of all haplotypes with the 3 specified markers:
#'ah <- allHaplotypes(mrknames=c("mrkA", "mrkB", "mrkC"))
#'# specify haplotype dosages of 4 tetraploid individuals,
#'# only the 3 occurring haplotypes (1, 5 and 6) are given:
#'haplodosg <-
#'  matrix(c(1,2,1, 4,0,0, 0,4,0, 0,0,4), nrow=3,
#'         dimnames=list(paste0("demohap_", c(1,5,6)), paste0("indiv", 1:4)))
#'# calculate the corresponding marker (SNP) dosages:
#'hapdos2mrkdos(hapdos=haplodosg, allhap=ah)
#'@export
hapdos2mrkdos <- function(hapdos, allhap) {
  if (is.vector(hapdos)) hapdos <- as.matrix(hapdos)
  hapdos <- expandHapdos(hapdos, nhap=nrow(allhap))
  res <- t(t(hapdos) %*% allhap)
  colnames(res) <- mrkdos2mrkdid(res, ploidy=sum(hapdos[, 1]), check=FALSE)
  res
} #hapdos2mrkdos

hapcomb2mrkdos <- function(hapcomb, allhap) {
  # similar to hapdos2mrkdos, but with haplotype combinations specified as
  # <ploidy> haplotype id's instead of the dosages of each haplotype
  if (is.vector(hapcomb)) hapcomb <- as.matrix(hapcomb)
  if (!is.matrix(allhap)) stop("allhap must be a matrix")
  if (!all(hapcomb %in% seq_len(nrow(allhap))))
    stop("hapcomb2mrkdos: input error")
  res <- matrix(NA_integer_, nrow=ncol(allhap), ncol=ncol(hapcomb),
                dimnames=list(colnames(allhap), NULL))
  for (i in seq_len(ncol(hapcomb))) {
    tmp <- allhap[hapcomb[, i],, drop=FALSE]
    res[, i] <- colSums(tmp)
  }
  colnames(res) <- mrkdos2mrkdid(res, ploidy=nrow(hapcomb), check=FALSE)
  res
} #hapcomb2mrkdos

split_hapnames <- function(hapnames) {
  # hapnames: a character vector of haplotype names, all starting with the same
  #           contig name, with the haplotype number appended after an
  #           underscore. Or a numeric vector with just the haplotype numbers.
  # return value: a list with 2 elements:
  # $ prefix: a character of length 1 with the entire (common) haplotype name
  #           up to and including the underscore separating the haplotype
  #           number; "x_" in case the hapnames consist of only the haplotype
  #           numbers
  # $ hapnrs: an integer vector with the haplotype numbers
  if (is.null(hapnames)) hapnames <- integer(0)
  if (length(hapnames) == 0) return(list(prefix="x_", hapnrs=integer(0)))
  #are the names just integers?
  suppressWarnings( hapnrs <- as.integer(hapnames) )
  if (!anyNA(hapnrs)) {
    prefix <- ""
  } else {
    #assume that all haplotype names are composed as
    # <haploblock name>_<hapnr>
    pos <- gregexpr("_", hapnames)
    num <- sapply(pos, length)
    if (any(num==0) || max(num) != min(num)) stop("haplotype names invalid")
    # same nr of underscores (1 or more) in all haplotype names
    lastpos <- sapply(pos, FUN="[", num[1])
    if (!min(lastpos) == max(lastpos)) stop("haplotype names invalid")
    # last underscore at same position in all haplotype names
    prefix <- substring(hapnames[1], 1, lastpos[1])
    hapnrs <- as.integer(substring(hapnames, lastpos[1]+1, nchar(hapnames)))
  }
  if (anyNA(hapnrs) || !all(hapnrs > 0) || !all(diff(hapnrs) > 0))
    stop("haplotype names invalid")
  return(list(prefix=prefix, hapnrs=hapnrs))
} #split_hapnames

#'@title Add dropped rows back to haplotype dosage matrix
#'@description Haplotype dosage matrices generated with the dropUnused=TRUE
#'lack some haplotypes that are required by pedigreeHapCheck; this
#'function adds them back.
#'@usage expandHapdos(hapdos, nhap, hbname)
#'@param hapdos a haplotype dosage matrix as returned by a.o. inferHaplotypes.
#'The rownames are assumed to have all the same length:
#'they are composed of the haploblock name to which the 0-padded haplotype number
#'is appended, separated by "_"
#'@param nhap the total number of haplotypes for the haploblock:
#'nhap = 2^nmrk with nmrk the number of markers in the haploblock
#'@param hbname the haploblock name (only used when hapdos has 0 rows)
#'@return a matrix similar to hapdos with rows for the dropped haplotypes
#'re-inserted, with the correct rownames and containing only 0 and NA values
#'@examples
#'# specify haplotype dosages of 4 tetraploid individuals,
#'# only the 3 occurring haplotypes (1, 5 and 6) are given:
#'haplodosg <- matrix(c(1,2,1, 4,0,0, 0,4,0, 0,0,4), nrow=3,
#'                    dimnames=list(paste0("demohap_", c(1,5,6)), paste0("indiv", 1:4)))
#'# add the rows for the absent haplotypes, assuming the haploblock consists
#'# of 3 markers, so 8 haplotypes are possible:
#'expandHapdos(hapdos=haplodosg, nhap=8)
#'@export
expandHapdos <- function(hapdos, nhap, hbname="x") {
  if (!is.matrix(hapdos)) stop("invalid hapdos, must be matrix")
  if (nrow(hapdos) == 0)
    return(matrix(NA_integer_, nrow=nhap, ncol=ncol(hapdos),
                  dimnames=list(
                    getHaplotypeNames(haploblock=hbname, hapcount=nhap),
                    colnames(hapdos))))
  hapnrs <- split_hapnames(rownames(hapdos)) # error if no rownames and other checks
  if (nhap == nrow(hapdos)) return(hapdos)
  hbname <- hapnrs$prefix
  hapnrs <- hapnrs$hapnrs
  if (hapnrs[length(hapnrs)] > nhap) stop("invalid hapdos")
  had <- matrix(NA_integer_, nrow=nhap, ncol=ncol(hapdos),
                dimnames=list(
                  getHaplotypeNames(haploblock=hbname, hapcount=nhap),
                  colnames(hapdos)))
  NAcols <- is.na(colSums(hapdos))
  had[hapnrs,] <- hapdos
  newrows <- setdiff(1:nhap, hapnrs)
  had[newrows, !NAcols] <- 0
  had
} #expandHapdos


#'@title convert haplotype dosages to haplotype combinations
#'@description convert a vector or matrix of haplotype dosages to a vector
#'or matrix of haplotype combinations
#'@usage hapdos2hapcomb(hapdos, ploidy)
#'@param hapdos a vector with one item for some or all possible haplotypes,
#'with the dosage of the haplotypes, summing to ploidy; or a matrix where
#'each column is such a vector, each representing an individual. Missing data
#'are allowed (but any NA makes the entire vector or matrix column unknown);
#'if a matrix, 0 rows indicate missing values for all individuals.
#'When hapdos is a matrix, the haplotype numbers must be given as rownames;
#'when a vector, as names. Not all possible haplotype numbers need to be present
#'but the ones that are present must be in ascending order, no duplicates
#'or missing values allowed.
#'@param ploidy the ploidy (one number)
#'@return a sorted vector of length ploidy with all haplotype numbers present,
#'or a matrix where each column is such a vector, with the same
#'colnames as hapdos
#'@examples
#'# specify haplotype dosages of 4 tetraploid individuals,
#'# only the 3 occurring haplotypes (1, 5 and 6) are given:
#'haplodosg <- matrix(c(1,2,1, 4,0,0, 0,4,0, 0,0,4), nrow=3,
#'                    dimnames=list(paste0("demohap_", c(1,5,6)), paste0("indiv", 1:4)))
#'# usage with hapdos as matrix:
#'hapdos2hapcomb(hapdos=haplodosg, ploidy=4)
#'# usage with hapdos as vector:
#'hapdos2hapcomb(hapdos=haplodosg[, 1], ploidy=4)
#'@export
hapdos2hapcomb <- function(hapdos, ploidy) {
  fn <- function(x, hapnrs) unlist(mapply(rep, hapnrs, x))

  if (missing(ploidy)) stop("ploidy not specified")
  vec <- !is.matrix(hapdos)
  if (vec) hapdos <- as.matrix(hapdos) #names become rownames
  result <- matrix(NA_integer_, nrow=ploidy, ncol=ncol(hapdos),
                   dimnames=list(NULL, colnames(hapdos)))
  cs <- colSums(hapdos) #all sum to ploidy or NA
  if (nrow(hapdos) > 0 && !all(is.na(cs))) {
    if (!all(cs %in% c(NA, ploidy)))
      stop("not all columns of hapdos sum to ploidy")
    if (is.null(rownames(hapdos))) {
      # we include this to allow more general use of hapdos2hapcomb
      hapnrs <- seq_len(nrow(hapdos))
    } else hapnrs <- split_hapnames(rownames(hapdos))$hapnrs
    NAcols <- is.na(cs)
    tmp <- apply(hapdos[, !NAcols, drop=FALSE], MARGIN=2, FUN=fn, hapnrs)
    result[, !NAcols] <- tmp
  }
  if (vec) result <- result[, 1, drop=TRUE] #rownames become names
  result
} #hapdos2hapcomb

#'@title convert haplotype combinations to haplotype dosages
#'@description converts matrices that contain haplotype combinations in columns
#'(as the matrices in ahccompletelist or ahclist) to matrices with the haplotype
#'dosages
#'@usage hapcomb2hapdos(hapcomb, nhap)
#'@param hapcomb matrix with <ploidy> rows and any number of columns. Each column
#'contains a set of haplotypes (numbers in 1 ... nhap). NAs allowed but no
#'values outside this range (not checked). E.g. a column with numbers
#'1-1-3-12 means two copies of haplotype 1 and one of haplotypes 3 and 12 each.
#'@param nhap the total number of possible haplotypes: 2^nmrk where nmrk is the
#'number of bi-allelic markers in the haploblock
#'@return matrix with nhap rows (one row for each possible haplotype) and as many
#'columns as in hapmat, giving the dosages of each haplotype in each column
#'of hapmat.
#'@examples
#'# specify haplotype combinations of 2 individuals:
#'haplocomb <- matrix(c(1,5,5,6, 5,6,6,6), ncol=2,
#'                    dimnames=list(NULL, c("indiv1", "indiv2")))
#'# convert to dosage matrix,
#'# assuming haploblock has 3 markers, so 8 possible haplotypes:
#'hapcomb2hapdos(hapcomb=haplocomb, nhap=8)
#'@export
hapcomb2hapdos <- function(hapcomb, nhap) {
  if (!is.matrix(hapcomb)) hapcomb <- matrix(hapcomb, ncol=1)
  cs <- colSums(hapcomb)
  hapdos <- matrix(NA_integer_, nrow=nhap, ncol=ncol(hapcomb),
                   dimnames=list(1:nhap, colnames(hapcomb)))
  for (col in which(!is.na(cs)))
    hapdos[, col] <- tabulate(hapcomb[, col], nbins=nhap)
  # version with apply is about 5% slower:
  # hapcombs <- as.integer(hapcomb) #prevent tabulate error if all are NA (logical)
  # hapdos <- apply(hapcombs, MARGIN=2, FUN=tabulate, nbins=nhap)
  # also this version has 0's instead of NA's in columns with missing data
  return(hapdos)
} #hapcomb2hapdos

hapindex2hapcomb <- function(hapindex, nhap, ploidy) {
  # hapindex: the index number of the desired haplotype combination (hapcomb)
  #           in a sorted sequence of all possible hapcombs at this ploidy
  #           and nhap
  # nhap: the number of possible haplotypes: 2^nmrk, with nmrk the number of
  #       markers in the haploblock
  # ploidy: the ploidy
  # return: an ordered vector of length ploidy with all haplotype numbers
  #         present, e.g. c(1, 5, 8, 8) : 1 copy of 1 and 5, 2 of 8
  # See function findcomb_wReplacement for details
  findcomb_wReplacement(n=nhap, k=ploidy, p=hapindex, check=FALSE)
}

hapcomb2hapindex <- function(hapcomb, nhap) {
  # hapcomb: an ordered vector of length ploidy with all haplotype numbers
  #          present, e.g. c(1, 5, 8, 8) : 1 copy of 1 and 5, 2 of 8
  # nhap: the total number of possible haplotypes: 2^nmrk
  # return: the index number of this hapcomb in a sorted sequence of all
  #         possible hapcombs at this ploidy and nmrk. See the
  #         findpos_wReplacement function
  findpos_wReplacement(n=nhap, comb=hapcomb, check=FALSE)
}

##!title Find a specific combination-with-replacement
##!description Find a specific combination-with-replacement assuming
##!the combinations are ordered
##!usage findcomb_wReplacement(n, k, p, check=TRUE)
##!param n a single integer > 0: the number of different object types
##!(which are numbered 1 to n) available
##!param k a single integer >= 0: the number of objects in each combination
##!param p a single integer: the 1-based index of the combination (see details)
##!param check if TRUE (default) the input data are checked against
##!incorrect length or values. If FALSE and n<=0 or k<0 or p outside the range
##!of possible combinations an incorrect value (even of incorrect length) may be
##!returned without warning.
##!details Assume that all possible combinations-with-replacement are ordered like
##!(for n=3, k=3) 111, 112, 113, 122, 123, 133, 222, 223, 233, 333 (there are
##!choose(n+k-1, k) such combinations). p indices one of these combinations
##!and this is returned as an integer vector with k elements.\cr
##!If k == 0 there is one combination: integer(0).
##!return an integer vector with the selected combination; if p does not index
##!a combination (p <= 0 or p > total number of combinations) NA is returned
##!(in that case the length of the returned vector is 1, not k)
findcomb_wReplacement <- function(n, k, p, check=TRUE) {

  findcomb_as_Delimiterpos <- function(n, k, p) {
    # n = nr of object types (colors, letters etc)
    # k = number of objects (balls) to select
    # p = 1-based index of target combination
    # return = positions of delimiters at index p
    p <- p - 1
    if (n <= 1 || k < 1) return(integer(0))
    if (n == 2) return(k+1-p)
    nd <- n-1 #nr of delimiters: 1 - nr of colors
    pos <- seq(n+k-nd, n+k-1) #original positions of delimiters, all at right
    for (j in 1:(nd-1)) {
      s <- 0 #cumulative nr of accounted-for combinations with this delimiter
      while (TRUE) {
        N <- choose(nd+k-pos[j], nd-j)
        if (s + N <= p) {
          pos[j] <- pos[j] - 1
          s <- s + N
        } else break
      }
      p <- p - s
    }
    #last delimiter:
    pos[nd] <- pos[nd] - p
    pos
  } #findcomb_as_Delimiterpos

  delimiterpos2comb <- function(pos, k) {
    # pos = vector of delimiter positions
    # k = number of objects (balls) to select
    # return = selected balls (ordered by increasing number, 1-based)
    if (length(pos) == 0) return(rep(1, k))
    pos <- pos - seq_along(pos)
    result <- integer(k)
    last <- 0
    for (j in 1:length(pos)) if (pos[j] > last) {
      result[(last+1):pos[j]] <- j
      last <- pos[j]
    }
    if (k > last) result[(last+1):k] <- length(pos) + 1
    result
  } # delimiterpos2comb

  if (check) {
    n <- as.integer(n); k <- as.integer(k); p <- as.integer(p)
    if (length(n) != 1 || n <= 0) stop ("n must be a single integer > 0")
    if (length(k) != 1 || n < 0) stop ("k must be a single integer >= 0")
    if (length(p) != 1) stop("p must be a single integer")
    combcount <- choose(n+k-1, k)
    if (p <= 0 || p > combcount) return(NA)
  }
  if (k == 0) {
    if (p == 1) return(integer(0)) else return(NA)
  }
  delimpos <- findcomb_as_Delimiterpos(n, k, p)
  delimiterpos2comb(delimpos, k)
} # findcomb_wReplacement

# don't export this function
#_#'@title Find the index number of a combination-with-replacement
#_#'@description Find the index number of a combination-with-replacement assuming
#_#'the combinations are ordered
#_#'@usage findpos_wReplacement(n, comb, check=TRUE)
#_#'@param n a single integer > 0: the number of different object types
#_#'(which are numbered 1 to n) available
#_#'@param comb an integer vector representing the combination of objects;
#_#'sorted from small to large; all items must be in 1:n
#_#'@param check if TRUE (default) the input data are checked against
#_#'incorrect length or values and an eror is generated if needed; also
#_#'comb is sorted first. If FALSE and n<=0 or comb not sorted correctly
#_#'or any elements of comb not in 1:n an incorrect value
#_#'may be returned without warning.
#_#'@details Assume that all possible combinations-with-replacement are ordered like
#_#'(for n=3, k=3) 111, 112, 113, 122, 123, 133, 222, 223, 233, 333 (there are
#_#'choose(n+k-1, k) such combinations). This function returns the index number
#_#'of a given combination in this sequence.\cr
#_#'If comb == integer(0), i.e. an empty combination, the index number is 1
#_#'as there is just one way of selecting 0 elements from any collection.\cr
#_#'This function relies on integer accuracy. Above 2^31-1 (2.1e9) R uses doubles
#_#'instead of integers and at some point exact
#_#'accuracy is lost. This can be overcome by using integer64 from package bit64
#_#'for the intermediate (ix) and returned results (and an integer64 - adapted
#_#'version of choose) - the maximum value is then 2^63-1 or 9.2e18.
#_#'@return an integer: the position (index number) of the given combination,
#_#'1-based.
findpos_wReplacement <- function(n, comb, check=TRUE) {
  #this function is used to assign an ID (index nr) to haplotype combinations
  #for easy comparing and countig in solveOneBlock
  if (check) {
    if (length(n) != 1 || n <= 0) stop("n must be a single integer > 0")
    if (anyNA(comb) || any(!(comb %in% 1:n)))
      stop("all elements in comb must be in 1:n")
    comb <- sort(comb)
  }
  k <- length(comb)
  if (k == 0) return(1)
  if (k == 1) return(comb)
  ix <- 1
  for (p in seq_len(length(comb)-1)) {
    #p is the current position within the (ordered) combination
    #find the index ix of the first combination with comb[p] at position p:
    #preceding colors at pos p to calculate nr of combs for:
    if (p == 1) start <- 1 else start <- comb[p-1]
    end <- comb[p] - 1
    x <- start # x loops through the preceding colors at position p
    while (x <= end) {
      tmpn <- n - x + 1
      tmpk <- k - p
      ix <- ix + choose(tmpn + tmpk - 1, tmpk)
      x <- x + 1
    }
  }
  # last position:
  ix <- ix + comb[length(comb)] - comb[length(comb)-1]
  ix
} #findpos_wReplacement

mrkdidfun <- function(dos, ploidy) {
  # for use within mrkdos2mrkdid and getAllHapcomb
  # dos: vector with marker dosages (length is nr of markers in the haploblock)
  # Fun fact: by using ploidy==1 and only dosages 0 and 1, this converts
  # a haplotype to its hapnr
  # return: the mrkdid (marker dosage combination ID)
  id <- 0
  fac <- 1
  for (i in length(dos):1) {
    id <- id + fac * dos[i]
    fac <- fac * (ploidy + 1)
  }
  id +1
} #mrkdidfun

#'@title get marker dosage IDs from marker dosages
#'@description get marker dosage IDs (mrkdid) from marker dosages
#'@usage mrkdos2mrkdid(mrkDosage, indiv=NULL, ploidy, check=TRUE)
#'@param mrkDosage matrix or data.frame. Markers are in rows, individuals in
#'columns, each cell has a marker dosage. Names of individuals are the column
#'names, marker names are the row names or (if a data.frame) in a column named
#'MarkerNames. All marker dosages must be in 0:ploidy or NA. If a data.frame,
#'additional columns may be present.
#'@param indiv NULL (default) or a character vector with names of individuals
#'to be selected. If NULL, all columns are selected;
#'if mrkDosage is a data.frame, that is probably not what is intended.
#'@param ploidy all marker dosages are checked to be in 0:ploidy or NA
#'@param check if TRUE (default) checkmrkDosage is called. If FALSE it is
#'assumed that mrkDosage is a matrix (not a data.frame) and it is not checked.
#'@details with ploidy==1 and (of course) all dosages 0 or 1 this function
#'returns the haplotype numbers for the haplotype specified by each column
#'@return a vector of marker dosage IDs, one for each column of mrkDosage:
#'each a number in 1:((ploidy+1)^nrow(mrkDosage)), NA for each column in
#'dosages where any of the dosages are NA
#'@examples
#'# dosages of 3 markers in 3 tetraploid individuals:
#'mrkdosg <-
#'  matrix(c(1,2,2, 4,0,0, 3,0,2), nrow=3,
#'         dimnames=list(c("mrkA", "mrkB", "mrkC"), c("indiv1", "indiv2", "indiv3")))
#'# get the "marker dosage IDs":
#'dids <- mrkdos2mrkdid(mrkDosage=mrkdosg, ploidy=4)
#'# convert dids back to marker dosages:
#'mrkdid2mrkdos(dosageIDs=dids, nmrk=3, ploidy=4, mrknames=c("mrkA", "mrkB", "mrkC"))
#'@export
mrkdos2mrkdid <- function(mrkDosage, indiv=NULL, ploidy, check=TRUE) {
  if (check) {
    mrkDosage <- checkmrkDosage(mrkDosage=mrkDosage, ploidy=ploidy,
                                indiv=indiv)
    if (is.character(mrkDosage)) stop(mrkDosage)
  }
  apply(mrkDosage, MARGIN=2, FUN=mrkdidfun, ploidy=ploidy)
} #mrkdos2mrkdid

#'@title get the marker dosages from mrkdids (marker dosage IDs)
#'@description get the marker dosages from mrkdids (marker dosage IDs)
#'@usage mrkdid2mrkdos(dosageIDs, nmrk, ploidy, mrknames=NULL)
#'@param dosageIDs vector of marker-dosage-combination IDs (mrkdid)
#'@param nmrk character vector of marker names in the haploblock
#'@param ploidy the ploidy level, a single positive integer
#'@param mrknames a vector of nmrk marker names (default NULL): if not NULL
#'these are used as rownames of the returned matrix
#'@return a matrix with in columns the marker dosages corresponding to the
#'marker dosageIDs, with these mrkdids as colnames, and one row per marker,
#'with marker names as rownames if mrknames are specified
#'
#'@examples
#'# dosages of 3 markers in 3 tetraploid individuals:
#'mrkdosg <-
#'  matrix(c(1,2,2, 4,0,0, 3,0,2), nrow=3,
#'         dimnames=list(c("mrkA", "mrkB", "mrkC"), c("indiv1", "indiv2", "indiv3")))
#'# get the "marker dosage IDs":
#'dids <- mrkdos2mrkdid(mrkDosage=mrkdosg, ploidy=4)
#'# convert dids back to marker dosages:
#'mrkdid2mrkdos(dosageIDs=dids, nmrk=3, ploidy=4, mrknames=c("mrkA", "mrkB", "mrkC"))
#'@export
mrkdid2mrkdos <- function(dosageIDs, nmrk, ploidy, mrknames=NULL) {
  mrkdos <- matrix(NA_integer_, nrow=nmrk, ncol=length(dosageIDs))
  if (!is.null(mrknames)) {
    if (length(mrknames) != nmrk) stop ("nmrk and mrknames don't match")
    rownames(mrkdos) <- mrknames
  }
  colnames(mrkdos) <- dosageIDs # no problem with NA
  if (is.character(dosageIDs)) dosageIDs <- as.integer(dosageIDs)
  uID <- unique(dosageIDs)
  pp1 <- ploidy + 1
  for (i in seq_along(uID)) {
    if (!is.na(uID[i]) && uID[i] >= 1 && uID[i] <= pp1^nmrk) {
      mdos <- integer(nmrk)
      di <- uID[i] - 1
      fac <- as.integer(round(pp1 ^ (nmrk-1), 0))
      m <- 1
      while (m <= nmrk) {
        mdos[m] <- di %/% fac
        di <- di %% fac
        fac <- fac %/% pp1
        m <- m + 1
      }
      mrkdos[,colnames(mrkdos)==uID[i]] <- mdos
    }
  }
  mrkdos
} #mrkdid2mrkdos

##!title get the min and max dosage of each haplotype over all
##!haplotype combinations
##!description get the min and max dosage of each haplotype over all
##!haplotype combinations
##!usage haplofrqMinMax(hapdos)
##!param hapdos a matrix with one row for each haplotype and one column per
##!combination of haplotypes (each column summing to ploidy);
##!a vector is interpreted as a one-column matrix
##!return a matrix with haplotypes in rows, and columns min and max with the
##!minimum and maximum times each haplotype occurs over all different
##!haplotype combinations
haplofrqMinMax <- function(hapdos) {
  if (is.vector(hapdos)) hapdos <- as.matrix(hapdos)
  minmax <- matrix(NA_integer_, nrow=nrow(hapdos), ncol=2)
  colnames(minmax) <- c("min", "max")
  minmax[,1] <- apply(hapdos, 1, min)
  minmax[,2] <- apply(hapdos, 1, max)
  minmax
} #haplofrqMinMax

#'@title pad an integer (prefix with zeroes to a fixed length)
#'@description pad an integer (prefix with zeroes to a fixed length)
#'@usage padded(x, maxx=0)
#'@param x vector of non-negative integers
#'@param maxx a single integer to whose nchar all x elements will be padded;
#'if 0 (default) the largest value in x will be used
#'@return a character vector representing the values of x left-padded with 0's
#'to the length of integer maxx or of max(x)
#'@examples
#'padded(c(21, 1, 121, NA, 0))
#'padded(c(21, 1, 121, NA, 0), maxx=1000)
#'@export
padded <- function(x, maxx=0) {
  formatC(x, width=nchar(max(x, maxx, na.rm=TRUE)), flag="0")
} #padded

getHaplotypeNames <- function(haploblock, hapcount, sep="_") {
  # haploblock: the name of the haploblock
  # hapcount: the number of possible haplotypes for this haploblock (= the
  #           number of the last haplotype)
  # result: a character vector with the haplotype names: all of the same
  #         length, with the appended haplotype number padded with zeroes
  paste0(haploblock, sep, padded(1:hapcount))
}

##!title Assign haplotype combinations based on output of inferHaps_noFS
##!description Assign haplotype combinations based on output of inferHaps_noFS
##!usage hapcomb_from_IHlist(ihl, mrkdids, ploidy, nocombs=FALSE)
##!param haploblockname Prefix to which the (zero-padded) haplotype numbers
##!are appended with separator '_'.
##!param ihl A list as returned by inferHaps_noFS
##!param mrkdids the mrkdids of all individuals
##!param ahcinfo a list as returned by loadAllHapcombLists
##!param nocombs if FALSE (default) the returned list does not have the matrix
##!hapcombs, only the vector ahccols
##!details The results of inferHaps_noFS, as passed in parameter ihl, determine
##!the haplotype combinations that match each mrkdid. This function translates
##!these results in a set of haplotype combinations for each individual.
##!return a list with
##!2 element: $hapcombs is the matrix as before, $ahccols is an integer vector
##!with for each indiv the column number in its mrkdid element of
##!ahc(complete)list (or NA if more than 1 ahccol possible).
hapcomb_from_IHlist <- function(ihl, mrkdids, ahcinfo, nocombs=FALSE) {
  #haplotypeNames <- getHaplotypeNames(haploblockname, nrow(ihl$allhap))
  nind <- sum(ihl$mrkdidsTable)
  hapcombs <- matrix(NA_integer_, nrow=ahcinfo$ploidy, #length(haplotypeNames),
                   ncol=length(mrkdids),
                   dimnames=list(NULL, names(mrkdids)))
  # we also return for each indiv the column nr in ahc(complete)list
  # for its mrkdid
  ahccols <- rep(NA_integer_, ncol(hapcombs))
  names(ahccols) <- colnames(hapcombs)
  for (i in seq_along(ihl$mrkdidsTable)) {
    did <- as.integer(names(ihl$mrkdidsTable)[i])
    mrkdidind <- which(mrkdids == did)
    #if (ncol(ihl$hclist[[i]]) == 1) {
    if (length(ihl$hclist[[i]]) == 1) {
      #we don't assign a haplotype combination if there is more than
      #one combination possible
      ahccols[mrkdidind] <- ihl$hclist[[i]]
      if (!nocombs) hapcombs[, mrkdidind] <-
         getAllHapcomb(mrkdid=did, nmrk=ncol(ihl$allhap),
                       ahcinfo=ahcinfo)[,ihl$hclist[[i]]]
    }
  }
  res <- list(ahccols=ahccols)
  if (!nocombs) res$hapcombs <- hapcombs
  res
} #hapcomb_from_IHlist

#'@title check a marker dosages matrix or data.frame
#'@description check a marker dosages matrix or data.frame, select columns and rows,
#'convert to matrix
#'@usage checkmrkDosage(mrkDosage, ploidy, indiv=NULL, markers=NULL,
#'generateMarkernames=TRUE)
#'@param mrkDosage matrix or data.frame. Markers are in rows, individuals in
#'columns, each cell has a marker dosage. If mrkDosage is a matrix,
#'the colnames are the individual names and the rownames are the marker names.
#'If mrkDosage is a data.frame and the name of the first column starts with
#'"marker" (upper/lowercase not relevant) the contents of this column are used
#'as marker names, else the rownames are used as marker names. The (other)
#'column names are the individual names.
#'All marker dosages must be in 0:ploidy or NA.
#'@param ploidy an integer value; all dosages are checked to be in 0:ploidy or NA
#'@param indiv NULL (default) or a character vector with names of individuals
#'to be selected. If NULL, all individuals are selected
#'@param markers NULL (default) or a character vector with names of markers
#'to be selected. If NULL, all markers are selected.
#'@param generateMarkernames if TRUE (default) and mrkDosage has no markernames
#'specified, markernames are generated automatically. Markernames
#'must be either present or be generated, and none may be missing.
#'@details This function is called by inferHaplotypes and by mergeReplicates, so
#'normally there is no need for the user to call this function directly.
#'@return a matrix with the selected columns in the order of indiv and the
#'selected rows in order of markers (or all columns or rows in the original
#'order if indiv or markers are not specified), with
#'names of individuals as column names, marker names as row names
#'@examples
#'mrkdos <- data.frame(
#'  marker=paste0("mrk", 1:3),
#'  indiv1=c(1, 2, 2),
#'  indiv2=c(4, 0, 0),
#'  indiv3=c(3, 0, 2))
#'# use all rows and columns:
#'checkmrkDosage(mrkDosage=mrkdos, ploidy=4)
#'# use only first and last row and column and change the order:
#'checkmrkDosage(mrkDosage=mrkdos, ploidy=4, indiv=c("indiv3", "indiv1"),
#'               markers=c("mrk3", "mrk1"))
#'@export
checkmrkDosage <- function(mrkDosage, ploidy, indiv=NULL, markers=NULL,
                           generateMarkernames=TRUE) {

  checkvector <- function(v, vname) {
    # v is not NULL
    if (anyNA(v)) stop(paste(vname, "may not contain missing values"))
    if (is.factor(v)) v <- as.character(v)
    if (is.character(v)) {
      v <- gsub("[[:blank:]]", "", v)
      if (any(v == "")) stop(paste(vname, "may not contain blank names"))
    }
    v
  }

  if (!is.null(indiv)) indiv <- checkvector(indiv, "indiv")
  if (!is.null(markers)) markers <- checkvector(markers, "markers")
  if (is.vector(mrkDosage)) {
    mrkDosage <- as.matrix(mrkDosage)
    #             1 col; names -> rownames
    colnames(mrkDosage) <- "1" #name of one individual
  }
  if (nrow(mrkDosage) == 0 || ncol(mrkDosage) == 0)
    stop("No data in mrkDosage")
  if (is.data.frame(mrkDosage)) {
    mrkcol <- which(tolower(substr(names(mrkDosage), 1, 6)) == "marker")
    if (length(mrkcol) == 1 && mrkcol==1) {
      rownames(mrkDosage) <- checkvector(mrkDosage[, 1], "markernames column")
      mrkDosage <- mrkDosage[, -1]
    }
    if (is.null(names(mrkDosage))) {
      names(mrkDosage) <- seq_len(length(mrkDosage))
    } else {
      names(mrkDosage) <- checkvector(names(mrkDosage), "column names of mrkDosage")
      if (anyDuplicated(names(mrkDosage)))
      stop("mrkDosage contains duplicated individual names")
    }
    if (is.null(indiv)) indiv <- names(mrkDosage)
    indcol <- match(indiv, names(mrkDosage))
    if (anyNA(indcol)) stop("Not all indiv occur in mrkDosage")
    isnum <- sapply(mrkDosage[,indiv], is.numeric)
    if (!all(isnum))
      stop("Not all selected columns of mrkDosage are numeric")
    mrkDosage <- as.matrix(mrkDosage[,indcol])
    #             rownames (markernames) and colnames (indiv) are retained
  }
  if (!is.matrix(mrkDosage))
    stop("mrkDosage is not a data.frame or matrix")
  # mrkDosage is (now) a numeric matrix
  # we will check it, possibly duplicating some earlier checks if mrkDosage
  # originally was a data.frame
  if (is.null(colnames(mrkDosage))) {
    colnames(mrkDosage) <- seq_len(ncol(mrkDosage))
  } else {
    colnames(mrkDosage) <-
      checkvector(colnames(mrkDosage), "column names of mrkDosage")
    if (anyDuplicated(colnames(mrkDosage)))
      stop("mrkDosage contains duplicated individual names")
  }
  if (is.null(indiv)) indiv <- colnames(mrkDosage)
  indcol <- match(indiv, colnames(mrkDosage))
  if (anyNA(indcol)) stop("Not all indiv occur in mrkDosage")
  mrkDosage <- mrkDosage[, indcol, drop=FALSE]
  if (is.null(rownames(mrkDosage))) {
    if (generateMarkernames) {
      rownames(mrkDosage) <- paste0("mrk", padded(seq_len(nrow(mrkDosage))))
    } else stop("no marker names specified")
  }
  rownames(mrkDosage) <- checkvector(rownames(mrkDosage), "marker names")
  if (is.null(markers)) markers <- rownames(mrkDosage)
  mrkrow <- match(markers, rownames(mrkDosage))
  if (anyNA(mrkrow)) stop("Not all markers occur in mrkDosage")
  mrkDosage <- mrkDosage[mrkrow,, drop=FALSE]
  if (!missing(ploidy) && !all(mrkDosage %in% c(NA, 0:ploidy)))
    stop("mrkDosage contains values not in 0:ploidy")
  mrkDosage
} #checkmrkDosage


single_cycle_infer <- function(allhap, mrkdidsTable, ahcinfo, knownHap,
                               useKnownDosage, progress) {
  #function used only by inferHaps_noFS
  #For each mrkdid, select the hapcombs that require the least new haplotypes
  #(not present in knownHap); if useKnownDosages is TRUE, among these select the
  #hapcombs with the maximum dosage of known haplotypes.
  #return: a list with 3 elements:
  # - nPresent: a vector with for each possible haplotype the minimum number of
  #             individuals in which it is inferred to be present (based on the
  #             selected hapcombs for each mrkdid)
  # - nAbsent: a vector with for each possible haplotype the minimum number of
  #            individuals in which it is inferred to be absent (based on the
  #            selected hapcombs for each mrkdid)
  # - hclist: a list with for each mrkdid the selected hapcombs as column
  #           numbers in ahc(complete)list
  abpr <- array(NA_integer_, dim=c(length(mrkdidsTable), 2, nrow(allhap)),
                dimnames=list(names(mrkdidsTable), c("nabs", "npres"),
                              1:nrow(allhap)))
  #       abpr stands for absent-present; contains for each mrkdid the min and
  #       max dosages of each haplotype, taken over all (remaining) haplotype
  #       dosage combinations
  hclist <- vector(length=length(mrkdidsTable), mode="list")
  for (i in seq_along(mrkdidsTable)) {
    #   for the current mrkdid, get the matrix with all possible
    #   haplotype combinations (different combinations in columns):
    ahc <- getAllHapcomb(mrkdid=names(mrkdidsTable)[i],
                        nmrk=ncol(allhap), ahcinfo=ahcinfo)
    ahccols <- 1:ncol(ahc)
    if (length(knownHap) > 0) {
      # Based on comparisons of 20171023, testscheme 1 is clearly worse, so:
      testscheme <- 2 # for now
      if (testscheme==1) {
        # versions before 20171006:
        # select the haplotype combination(s) with the highest total dosage
        # of known haplotypes:
        knowndos <- colSums(ahc[knownHap,, drop=FALSE])
        ahccols <- which(knowndos == max(knowndos))
        #ahc <- ahc[, knowndos==max(knowndos), drop=FALSE]
        ahc <- ahc[, ahccols, drop=FALSE]
      } else {
        #testscheme==2: version 20171006 (minimize the nr of haplotypes,
        #instead of maximize the dosage of known haplotypes, and if
        #useKnownDosage is TRUE, maximize the dosage of known haplotypes
        #among the selected hapcombs)
        ahc_new <- ahc
        ahc_new[ahc_new %in% knownHap] <- NA
        fn <- function(x) {length(unique(x[!is.na(x)]))}
        newcount <- apply(ahc_new, MARGIN=2, FUN=fn)
        ahccols <- which(newcount == min(newcount)) #least nr of new haplotypes
        if (useKnownDosage & length(ahccols) > 1) {
          subahc <- ahc_new[, ahccols, drop=FALSE]
          knowndos <- colSums(apply(subahc, MARGIN=2, FUN=is.na))
          ahccols <- ahccols[knowndos == max(knowndos)] #max dosage of known haplotypes
        }
        ahc <- ahc[, ahccols, drop=FALSE]
      }
      # (a different approach could be: for each possible set of haplotypes in
      # addition to the known ones, from small to large, see how many
      # individuals have one or more solutions, and stop adding haplotypes at
      # some point (when extra haplotypes would probably only be used to
      # fit errors in the SNP dosages);
      # but then, no selection is made between solutions for a specific did)
    }
    hclist[[i]] <- ahccols #was: ahc; now a vector of column nrs rather than the entire matrix
    mima <- haplofrqMinMax(hapdos=hapcomb2hapdos(ahc, nhap=nrow(allhap)))
    abpr[i,1,] <- mrkdidsTable[i] * (mima[,2] == 0) #nr of ind with this mrkdid
    #             that cannot have the haplotypes
    abpr[i,2,] <- mrkdidsTable[i] * (mima[,1] > 0) # nr of ind with this mrkdid
    #             that must have the haplotypes
  } # for i in mrkdidsTable
  nPresent <- apply(abpr[,2,,drop=FALSE], MARGIN=3, FUN=sum)
  #           vector with the nr of indiv in which each haplotype must be present
  nAbsent <- apply(abpr[,1,,drop=FALSE], MARGIN=3, FUN=sum)
  #           vector with the nr of indiv in which each haplotype must be absent
  names(hclist) <- names(mrkdidsTable)
  list(nPresent=nPresent,
       nAbsent=nAbsent,
       hclist=hclist)
} #single_cycle_infer


##!title infer haplotypes from marker dosages
##!@description infer haplotypes from marker dosages for a single haploblock,
##!treating all individuals as unrelated
##!usage inferHaps_noFS(mrkdids, ahcinfo, markers,
##minfrac=c(0.1, 0.01), knownHap=integer(0), useKnownDosage=TRUE, progress=TRUE)
##!param mrkdids a named integer vector with for each indiv its mrkdid; the
##names are the individual names.
##!param ahcinfo a list produced by loadAllHapcombFiles
##!param markers a character vector with the names of the markers in the
##haploblock, in the same order as used for calculating the mrkdids
##!param minfrac a vector of 1 or 2 fractions, the second smaller than the
##!first. A haplotype is considered to be certainly present if it must occur
##!in at least a fraction minfrac[1] of all individuals; default 0.1. For the
##!meaning of the optional second value in minfrac see Details.
##!param knownHap selected haplotypes (haplotypes that must be present according
##!to prior inference or knowledge, numbers refer to rows of matrix produced by
##!allHaplotypes); default integer(0), i.e. no known haplotypes.
##!param useKnownDosage a single logical, intended for internal use only.
##!; if TRUE (default), haplotype combinations for a mrkdid are selected not
##!only based on the smallest number of new haplotypes but also on the highest
##!total dosage of known haplotypes
##!param progress if TRUE (default), and new haplotype combinations need to be
##!calculated, and the number of markers and the ploidy are both >= 6, progress
##!is indicated by printed messages
##!details This function is for internal use; it is called only from
##!inherHaplotypes. It performs haplotype inference ignoring the presence of
##!FS families.
##The returned list contains some general calculations and statistics
##!based on mrkdids and ploidy, and some parts that are the result of the
##!inference. The primary of these is hclist: this contains for each mrkdid in
##!the population the most likely combination of haplotypes (sometimes
##!more than one, if several are equally likely = all have the maximum dosage
##!of haplotypes inferred to be present).\cr
##!Principle: first the haplotypes are derived that must be present in at least
##!a fraction minfrac of the individuals (because those haplotypes occur in all
##!possible haplotype configurations that result in the observed marker dosages).
##!In subsequent iterations, for individuals not completely explained by these
##!known haplotypes we select the hapcombs that require the smallest number of
##!extra haplotypes and we determine which haplotypes are required in all of
##!these solutions. For each haplotype we count in how many individuals it is
##!required, and all haplotypes that are required in a minimum fraction of
##individuals (minfrac), together with any haplotypes known from other
##!information (knownHap) become the new set of known haplotypes. We repeat this
##!step until the set of known haplotypes does not change any more.\cr
##!During this process there may be additional haplotypes needed for some
##!mrkdids, that are not frequent enough to reach the minfrac[1] criterion. If
##!there is a second (smaller) fraction in minfrac, in a final round we also consider
##!those extra haplotypes as must-be-present, provided they occur in at least
##!a fraction minfrac[2] of all individuals (minimum 2 individuals). If that
##!leads to more mrkdids having a unique haplotype combination AND not having
##!other mrkdids getting more combination the new assignments are used.\cr
##!If knownHap in the function call already contains haplotypes, these are
##!considered to be certainly present, even if they don't meet the minfrac
##!criterion. This may help to reduce the number of haplotype configurations
##!for some individuals, but may also increase them.
##!return a list with elements:
##!\itemize{
##!\item{nPresent: a vector with for each haplotype the number of
##!individuals in which it must occur}
##!\item{nAbsent: a vector with for each haplotype the number of
##!individuals in which it cannot occur}
##!\item{minfrac: same as parameter minfrac}
##!\item{allhap: a matrix as returned by allHaplotypes}
##!\item{mrkDosage: a  matrix with marker dosages for all individuals in indiv
##!(produced by calling checkmrkDosage on mrkDosage)}
##!\item{mrkdids: a vector of the mrkdid (marker dosage ID) for each individual in
##!mrkDosage (each combination of marker dosages has its own ID; if any of the
##!markers has an NA dosage the corresponding mrkdid is also NA)}
##!\item{mrkdidsTable: a table of the counts of all non-NA mrkdids in the
##!population}
##!\item{hclist: a list in which each element has the name of one of the mrkdids
##!in mrkdidsTable, and contains a matrix with <ploidy> rows and one column
##!for each of the remaining haplotype combinations: (a subset of the columns
##!of) the matrix returned by getAllHapcomb for that mrkdid}
##!}
inferHaps_noFS <- function(mrkdids, ahcinfo, markers,
                            minfrac=c(0.1, 0.01),
                            knownHap=integer(0), useKnownDosage=TRUE,
                            progress=TRUE) {
  if (anyNA(knownHap) || !all(knownHap %in% 1:(2^length(markers))) ||
      anyDuplicated(knownHap)) stop("invalid knownHap")
  if (useKnownDosage) useKnownDosage <- c(FALSE, TRUE) #else just FALSE
  origHap <- sort(knownHap)
  mrkdidsTable <- table(mrkdids)
  allhap <- allHaplotypes(markers)
  # haplocombs are already calculated in inferHaplotypes
  nind <- sum(mrkdidsTable) #includes only individuals with no missing marker data
  prev_knownHap <- list() #stores the knownHap after each cycle
  knownHap_ukd <- list(NULL)
  cycle <- 0
  for (ukd in seq_along(useKnownDosage)) {
    while (TRUE) {
      cycle <- cycle + 1
      res <- single_cycle_infer(allhap=allhap, mrkdidsTable=mrkdidsTable,
                                ahcinfo=ahcinfo, knownHap=knownHap,
                                useKnownDosage=useKnownDosage[ukd],
                                progress=progress)
      knownNew <- sort(union(origHap, which(res$nPresent >= minfrac[1] * nind)))
      haploLost <- setdiff(knownHap, knownNew)
      if (length(haploLost) > 0) {
        #The initial expectation was that within these cycles known haplotypes can
        #only be gained, not lost (except the origHap which need not be supported
        #by the current data, but these have been added back to knownNew).
        #Therefore the following warning was generated if haplotypes were lost:
        #warning(paste0("inferHaps_noFS: lost in cycle ", cycle, ": ",
        #              paste(haploLost, collapse=" ")))
        #However it seems that it is possible to lose previous must-have
        #haplotypes, although this does not happen often.
        #We must prevent infinite looping, so we store the knownHap and
        #if the situation occurs again with the same knownHap we break with
        #the original knownHap and the result of one single_cycle_infer with that
        #knownHap
        if (length(prev_knownHap) == 0 ||
            !any(sapply(prev_knownHap, identical, knownHap))) {
          prev_knownHap[[length(prev_knownHap) + 1]] <- knownHap
        } else {
          #recalculate the result of the first cycle - very fast and avoids
          #having to store it:
          res <- single_cycle_infer(allhap=allhap, mrkdidsTable=mrkdidsTable,
                                    ahcinfo=ahcinfo, knownHap=origHap,
                                    useKnownDosage=useKnownDosage[ukd],
                                    progress=FALSE)

          knownHap <- origHap
          break #break out of while loop
        }
      }
      if (identical(knownNew, knownHap)) {
        knownHap_ukd[[ukd]] <- knownHap
        break
      }
      knownHap <- knownNew
    } #while TRUE
  } #for ukd
  # if useKnownDosage has length 2 (FALSE, TRUE) and we had convergence after
  # the first ukd but not the second we should fall back on the first:
  if (length(useKnownDosage) == 2 && length(knownHap_ukd) == 1 &&
      !is.null(knownHap_ukd[[1]])) {
    res <- single_cycle_infer(allhap=allhap, mrkdidsTable=mrkdidsTable,
                              ahcinfo=ahcinfo, knownHap=knownHap_ukd[[1]],
                              useKnownDosage=FALSE,
                              progress=FALSE)
    knownHap <- sort(union(origHap, which(res$nPresent >= minfrac[1] * nind)))
  }
  #Now we have haplotype combinations selected for all mrkdids,
  #stored in res$hclist, based on the final knownHap set of
  #haplotypes that must be present.
  #The selected hapcombs for some mrkdids may include additional haplotypes
  #that did not make it into knownHap because they did not reach the minfrac
  #criterion.
  #We may now consider also those extra haplotypes as must-be-present
  #(or some of them, above a certain threshold: occurring in at least 2 and
  #at least 1% of the individuals) and see if that leads to more mrkdids
  #having a unique hapcomb (and not having other mrkdids getting additional
  #hapcombs)
  if (length(minfrac) > 1 && (minfrac[2] < minfrac[1])) {
    mrkdidHaccount <- vapply(res$hclist[names(mrkdidsTable)],
                             FUN=length, FUN.VALUE=0)
    onecomb <- mrkdidHaccount == 1
    if (any(!onecomb) && any(onecomb)) {
      # any(!onecomb): else all mrkdid have already one hapcomb, no need
      # for further selection
      # any(onecomb): else there are no mrkdid to generate additional knownHap
      hac1comb <- unlist(res$hclist[names(mrkdidsTable)[onecomb]])
      names(hac1comb) <- names(mrkdidsTable)[onecomb]
      #           a vector with the mrkdids as names and the ahc column nr as values
      #           for the mrkdids with just 1 hapcomb
      hapcomb <- matrix(NA_integer_, nrow=ahcinfo$ploidy, ncol=length(hac1comb))
      for (d in seq_along(hac1comb)) {
        hapcomb[, d] <-
          getAllHapcomb(mrkdid=names(hac1comb)[d], nmrk=ncol(allhap),
                        ahcinfo=ahcinfo)[, hac1comb[d]]
      }
      #calculate the nr of individuals in which each haplotype occurs:
      #select the mrkdids for which only 1 haplotype combination is possible,
      #and count the nr of indiv in which they occur:
      hapfrq <- colSums(as.vector(mrkdidsTable[onecomb]) *
                          t(hapcomb2hapdos(hapcomb, nhap=nrow(allhap)) > 0))
      final_hap <-
        sort(union(origHap,
                   which (hapfrq >= 2 & hapfrq > minfrac[2]*sum(mrkdidsTable))))
      resfinal <-
        single_cycle_infer(allhap=allhap, mrkdidsTable=mrkdidsTable,
                           ahcinfo=ahcinfo, knownHap=final_hap,
                           useKnownDosage=useKnownDosage[length(useKnownDosage)],
                           progress=progress)
      mrkdidHaccountfinal <-
        vapply(resfinal$hclist[names(mrkdidsTable)], FUN=length, FUN.VALUE=0)
      #keep resfinal instead of res only if all mrkdids that had only one
      #haplotype combi still have only one (and no dids have more columns
      #than before)
      if (all(mrkdidHaccountfinal <= mrkdidHaccount)) {
        res <- resfinal
        if (FALSE) { # set to TRUE for checking this
          #print comparison:
          lesscol <- mrkdidHaccountfinal < mrkdidHaccount
          m <- matrix(NA_integer_, nrow=2, ncol=sum(lesscol))
          colnames(m) <- names(mrkdidsTable)[lesscol]
          rownames(m) <- c("before", "after")
          m[1,] <- mrkdidHaccount[lesscol]
          m[2,] <- mrkdidHaccountfinal[lesscol]
          print(m)
        }
        # perhaps now there are more mrkdids with a unique haplotype combination
        # (which is the reason for this final test). If so, these may again have
        # extra haplotypes so we could do the final step again. We will not do
        # that because these extra haplotypes are selected on a weak basis
        # (weaker than minfrac[1])
      }
    }
  }

  #finally we add some extra data to the output list:
  res$allhap <- allhap
  res$mrkdidsTable <- mrkdidsTable
  res
} #inferHaps_noFS


chisq.probtest <- function(x, p=rep(1/length(x), length(x)),
                           correct=TRUE, rescale.p = FALSE) {
  #version of the standard chisq.test function, only for given probabilities,
  #that applies the Yates continuity correction in this case if correct=TRUE.
  #The original function applies Yates only in the case of a contingency test,
  #not for given probabilities.
  DNAME <- deparse(substitute(x))
  if (!is.null(dim(x))) {
    if (length(dim(x)) == 1) x <- as.vector(x) else
      stop("'x' must be a vector")
  }
  if (any(x < 0) || anyNA(x))
    stop("all entries of 'x' must be nonnegative and finite")
  if ((n <- sum(x)) == 0)
    stop("at least one entry of 'x' must be positive")
  if (length(x) == 1L)
    stop("'x' must at least have 2 elements")
  if (length(x) != length(p))
    stop("'x' and 'p' must have the same number of elements")
  if (any(p < 0))
    stop("probabilities must be non-negative.")
  if (abs(sum(p) - 1) > sqrt(.Machine$double.eps)) {
    if (rescale.p)
      p <- p/sum(p)
    else stop("probabilities must sum to 1.")
  }
  METHOD <- "Chi-squared test for given probabilities"
  E <- n * p
  if (correct) {
    YATES <- min(0.5, abs(x - E))
    if (YATES > 0)
      METHOD <- paste(METHOD, "with Yates' continuity correction")
  } else YATES <- 0
  V <- n * p * (1 - p)
  STATISTIC <- sum((abs(x - E) - YATES)^2/E) #old: sum((x - E)^2/E)
  names(E) <- names(x)
  PARAMETER <- length(x) - 1
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)

  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  if (any(E < 5) && is.finite(PARAMETER))
    warning("Chi-squared approximation may be incorrect")
  structure(list(statistic = STATISTIC, parameter = PARAMETER,
                 p.value = PVAL, method = METHOD, data.name = DNAME,
                 observed = x, expected = E,
                 residuals = (x - E)/sqrt(E), stdres = (x - E)/sqrt(V)),
            class = "htest")
} #chisq.probtest

testFSseg <- function(parhac, DRrate, errfrac,
                      FSmrkdidsTable, allhap) {
  # This function performs a test of the observed segregation of mrkdids in a
  # FS family (supplied in FSmrkdidsTable) against the segregation expected
  # with parents having the haplotype combinations (given in parhac) under
  # assumed DRrate and errfrac.
  # This function does most of the work in solveOneFS::calcParcombStats().
  # parameters:
  # parhac: a matrix of 2 columns and ploidy rows with the haplotype combinations
  #         of the parents. Does not contain NAs!
  # DRrate, errfrac: the assumed frequency of DR and of SNP genotyping errors
  # FSmrkdidsTable: a table of the *observed* counts of each mrkdid in the FS
  # allhap: matrix of nmrk columns and 2^nmrk rows with all possible haplotypes
  # return value: a list with 2 components:
  # $stats: list with statistics for the given parental combination (parhac):
  #   $nhap: the number of haplotypes
  #   $P: P-value of the (chi-squared ...) test of the observed segregation
  #   $toosmall: logical vector with one value per pco, TRUE if the FS is
  #              too small and/or too fragmented to calculate a P-value
  # $didNhapcombLst: a list with the following components:
  #   $mrkdid: vector with names of expected mrkdids in the FS, in order of
  #            hapcomb (multiple hacs may have the same mrkdid!)
  #   $hapcomb: matrix of hacs expected in the FS (1 col per hac); with
  #             order of columns sorted so that all hacs with the same mrkdid
  #             are adjacent
  #   $expHCfreq: vector of the expected hac frequencies in order of hapcomb
  #   $dupdid: logical vector, TRUE where mrkdid[i] == mrkdid[i-1], else FALSE

  calcErrprob <- function(origP, ninv, invObs, FSsize, maxErrprob) {
    # Function within testFSseg
    # origP: vector, the error-free probabilities of all the non-invalid
    #        genotypes; the last value is the error-free total probability of
    #        the Invalid group (for the segregation check we may add some
    #        valid but low-probability genotypes to the Invalid group,
    #        hence it can have a non-zero probability). origP must sum to 1
    # ninv: scalar, the nr of genotypes in the Invalid group
    # invObs: the observed nr of individuals in the Invalid group
    # FSsize: the total nr of (fully genotyped) FS indiv
    # maxErrprob: the specified (maximum) probability of marker genotype
    #             errors (the probability that at least one of the markers in
    #             the haploblock is mis-genotyped)
    # Return: the inferred error probability: the errorprob that best
    #         explains the observed invObs, with the following restrictions:
    #         errprob<=maxErrprob, errprob >=
    #         the value leading to an expected invObs=1, except if that
    #         is higher than maxErrprob
    # The calculation is based on the same principles as in
    # solveOneFS::redistErr()
    if (invObs == 0) invObs <- 1
    invFrac <- invObs / FSsize
    ntot <- ninv + length(origP) - 1
    invP <- origP[length(origP)]
    errprob <- (invFrac - invP) * (ntot-1) / (ninv - ntot*invP)
    errprob <- max(errprob, 0)
    errprob <- min(errprob, maxErrprob)
    errprob
  } #calcErrprob in testFSseg

  redistErr <- function(origP, ninv, errP) {
    # Function within testFSseg
    # calculate the resulting probabilities from the original probabilities,
    # taking the error rate and the number of genotypes in the Invalid group
    # into account.
    # origP: vector, the error-free probabilities of all the non-invalid
    #        genotypes; the last value is the error-free total probability of
    #        the Invalid group (for the segregation check we may add some
    #        valid but low-probability genotypes to the Invalid group,
    #        hence it can have a non-zero probability). origP must sum to 1
    # ninv: scalar, the nr of different genotypes (*not individuals*) in
    #       the Invalid group
    # errP: the probability of genotyping errors
    # Principle: the errors are uniformly distributed over all genotypes except
    # the correct one. Depending on the number of genotypes in the Invalid group
    # more or less of the errors will go to that group.
    if (length(origP) < 2) return(origP)
    origPinv <- origP[length(origP)]
    origP <- origP[-length(origP)]
    nval <- length(origP) # nr of valid genotypes
    ntot <- ninv + nval
    Plost <- origP * errP
    PlostInv <- origPinv * errP * (nval / (ntot-1))
    Pgained <- matrix(0, nrow=nval, ncol=nval)
    # rows: valid genotypes to contribute those errors
    # cols: valid genotypes to gain from errors;
    PgainedInv <- rep(0, nval)
    for (rw in 1:nval) {
      Pgained[rw, ] <- Plost[rw] / (ntot-1)
      Pgained[rw, rw] <- 0
      PgainedInv[rw] <- Plost[rw] * ninv / (ntot-1)
    }
    c(origP - Plost + colSums(Pgained) + PlostInv/nval,
      origPinv - PlostInv + sum(PgainedInv))
  } #redistErr within testFSseg

  nmrk <- ncol(allhap)
  ploidy <- nrow(parhac)
  FScount <- sum(FSmrkdidsTable) # count of fully genotyped FS indiv
  stats <- list()
  stats$nhap <- length(unique(as.vector(parhac))) #length(table(parhac)) nr of different haplotypes
  FSseg <- getFSfreqs(parhac=parhac, DRrate=DRrate)
  #we will now get the expected freq of unique mrkDosage in the FS (note that
  #(FSseg$FShac has unique *haplotype* combinations, but that different
  #haplotype combinations may result in the same *marker dosage*
  #combinations)
  FSmrkdos <- hapcomb2mrkdos(FSseg$FShac, allhap=allhap)
  #find the order to sort the FSmrkdos columns (on all rows) = in order of mrkdid
  #(in this way columns with the same *marker dosages* are next to each other):
  o <- do.call(order, lapply(1:nrow(FSmrkdos), function(i) FSmrkdos[i,]))
  FSmrkdos <- FSmrkdos[, o, drop=FALSE]
  expFSfreq <- FSseg$freq [o]
  du <- duplicated(FSmrkdos, MARGIN=2)
  didNhapcombLst <- list(mrkdid=colnames(FSmrkdos), #expected mrkdids (may contain duplicates)
                         hapcomb=FSseg$FShac[, o, drop=FALSE], #expected hacs in same order
                         expHCfreq=expFSfreq, #expected hac freqs in same order
                         dupdid=du) #whether each hapcomb col has same did as previous
  uniq <- which(!du)
  #sum the freqs per set of duplicates:
  uniqfreq <- numeric(length(uniq)); names(uniqfreq) <- colnames(FSmrkdos)[uniq]
  #uhx <- c(uniq, length(du)+1)
  #for (u in seq_along(uniqfreq))
  #  uniqfreq[u] <- sum(expFSfreq[uhx[u]:(uhx[u+1]-1)])
  uniqfreq <- tapply(expFSfreq, as.integer(didNhapcombLst$mrkdid), FUN=sum) # expected frequencies of mrkdids in FS
  FSmrkdos <- FSmrkdos[, uniq, drop=FALSE]
  expFSmrkdids <- colnames(FSmrkdos)
  # the expected expFSmrkdids are now unique,
  # their (error-free) probabilities are in uniqfreq

  # The chisquare test of the segregation is very sensitive to mrkdids with
  # a low probability where still a plant is observed. Therefore we group
  # those mrkdids with the invalid mrkdids and test the total nr of FS indiv
  # in that group against the total expectation. So:
  # - group all the hapcombs with an expected nr of individuals < 1 together
  #   with all invalid hapcombs into one group and get its total probability
  # - adjust the probabilities of the valid mrkdids and the invalid group
  #   with the error rate.
  # - Test the observed segregation (incl the invalid group) against the
  #   expected segregation
  lowfrqMrkdids <- uniqfreq * FScount < 1
  # and those that are observed but unexpected (invalid):
  invalidMrkdids <- which(!(names(FSmrkdidsTable) %in% expFSmrkdids))
  # calculate the probabilities:
  invalidProb <- sum(uniqfreq[lowfrqMrkdids]) #the really invalid have prob=0
  validProbs <- uniqfreq[!lowfrqMrkdids]
  # get the observed frequencies in the same order:
  invalidObs <- sum(FSmrkdidsTable[invalidMrkdids]) +
    sum(FSmrkdidsTable[names(FSmrkdidsTable) %in% expFSmrkdids[lowfrqMrkdids]])
  validObs <- FSmrkdidsTable[match(expFSmrkdids[!lowfrqMrkdids],
                                   names(FSmrkdidsTable))]
  validObs[is.na(validObs)] <- 0
  # redistribute the error probability over the mrkdids:
  ninv <- (ploidy+1)^nmrk - length(validProbs) #the number of
  #       mrkdids in the Invalid group
  errP <- calcErrprob(origP=c(validProbs, invalidProb), ninv=ninv,
                      invObs=invalidObs, FSsize=FScount, maxErrprob=errfrac)
  probs <- redistErr(origP=c(validProbs, invalidProb),
                     ninv =ninv, errP = errP)
  # We always have at least one expected mrkdid (in case of
  # no segregation) and a category Invalid
  # Four situations:
  # (1) no valid mrkdid remain (small FS family, highly segregating, large
  #     errfrac, so that no mrkdid has an expectation of >= 1 FS indiv):
  #     toosmall <- TRUE
  # (2) at least one invalid FS expected (including low-freq mrkdids and the
  #     errfrac): P <- chisquare test on valid mrkdids + invalid group
  # (3) less than 1 invalid expected:
  #     calculate Pinvalid <- binomial P for invalidObs
  # (3a) only one valid mrkdid: Pseg <- 1.0
  # (3b) more than one valid mrkdid: Pseg <- chisquare test on valid mrkdids
  # (3a and 3b) P <- Pinvalid * Pseg
  if (length(probs) == 1) {
    # (1): only the invalid group, no valid mrkdids with >= 1 expected FS indiv
    stats$toosmall <- TRUE
    stats$P <- NA
  } else {
    stats$toosmall <- FALSE
    if (probs[length(probs)] * FScount >= 1) {
      # (2): invalid group large enough for inclusion in chiquared test
      stats$P <- suppressWarnings({
        #suppress warnings about approximation for small values
          chisq.probtest(c(validObs, invalidObs), p=probs)$p.value
      })
    } else {
      # (3): less than 1 "invalid" FS expected, don't include invalid group
      #      in chisquared test:
      Pinvalid <- 1 - pbinom(invalidObs, FScount, prob=probs[length(probs)])
      if (length(validObs) <= 1) {
        Pseg <- 1.0 # (3a): only one valid mrkdid
      } else {
        # (3b): more than one valid mrkdid
        suppressWarnings({
          #suppress warnings about approximation for small values
          p <- probs[-length(probs)]; p <- p / sum(p)
          Pseg <- chisq.probtest(validObs, p=p)$p.value
        })
      }
      stats$P <- Pinvalid * Pseg
    }
  }
  list(stats=stats, didNhapcombLst=didNhapcombLst)
} #testFSseg

##!title find the optimal combination of parental haplotype combinations
##!description For one FS family and one haploblock, find the optimal
##!combination of parental haplotype combinations
##!usage solveOneFS(mrkDosage, ihl, FS, P, priorPhac, hacXmat, maxparcombs,
##!ploidy, minfrac, errfrac, DRrate, minPseg)
##!param mrkDosage dosage matrix with only the markers in the current haploblock
##!param ihl a list returned by inferHaps_noFS
##!param FS character vector with sample names of the current FS
##!param P vector with 2 sample names, one for parent1 and one for parent2.
##!There should be corresponding columns in mrkDosage , and orig_ihl should be
##!based on that mrkDosage
##!param knownPhac: a vector with two numbers (or NA) indexing the column in
##ahc(complete)list for the Pdids that has already be inferred to be the solution
##for that parent
##!triedPhacs: a list of two vectors of numbers (possibly empty) with the Phacs
##!already tried with no success for firstcombs in an earlier call for this FS,
##!before the FS was skipped because of too many extracombs. Numbers are colnrs
##!in ahclist
##!param maxparcombs Parent 1 and 2 both may have multiple possible haplotype
##!combinations. For each pair of haplotype combinations (one from P1 and one
##!from P2) the expected FS segregation must be checked against the observed.
##!This may take a long time if many such combinations need to be checked.
##!This parameter sets a limit to the number of allowed combinations.
##!param ploidy all marker dosages should be in 0:ploidy or NA
##!param minfrac vector of two fractions, default 0.1 and 0.01. A haplotype is
##!considered to be certainly present if it occurs in at least a fraction
##!minfrac[1] of all individuals; in the final stage for the "other"
##!individuals (those that do not belong to the FS or its parents) this fraction
##!is lowered to minfrac[2]; see also inferHaps_noFS
##!param errfrac the assumed fraction marker genotypes with an error (over all
##!markers in the haploblock). The errors are assumed to be uniformly distributed
##!over all except the original marker dosage combinations (mrkdids)
##!param DRrate the rate of double reduction per meiosis (NOT per allele!); e.g.
##!with a DRrate of 0.04, a tetraploid parent with genotype ABCD will produce
##!a fraction of 0.04 of DR gametes AA, BB, CC and DD (each with a frequency of
##!0.01), and a fraction of 0.96 of the non-DR gametes AB, AC, AD, BC, BD, CD
##!(each with a frequency of 0.16)
##!param minPseg candidate solutions with a chisquared P-value for segregation
##!smaller than this threshold are rejected; if none of the candidates remain
##!no solution is found (and in the return list parcombok has 0 rows and
##!NoSolution is TRUE)
##!return a list with two items\cr
##! $FSdat: itself a list with items
##!   $message\cr
##!   $parcombok: a matrix with 0 or more rows and two columns: each row represents
##!             one solution, with the two numbers on that row the hapixes
##!             of the parental hacs (Phacs). parcombok has 0 rows if
##!             no solution found. NEW: stores the ahclist column numbers, not hapixes
##!   $didNhapcomb: a list with one item per selected solution; each item is itself
##!               a list with (a.o.) expected FS frequencies of mrkdids and the
##!               corresponding haplotype combinations\cr
##!   $NoSolution: TRUE or FALSE: TRUE if one or more solutions selected,
##!              FALSE if no solution attempted (one or both parent(s) or all
##!              FS indiv with missing mrk data, or no fitting
##!              solutions found)\cr
##!              If NoSolution is TRUE didNhapcomb and Pseg may not be
##!              present
##!   $Pseg: a vector of P-values, one for each row of parcombok, or one if
##!             no solutions (P of the best rejected solution)
##! and a second item:
##! $hacXmat: matrix as parameter hacXmat, possibly extended with additional hacs
##! NEW: haxXmat not used or output anymore
# solveOneFS is not exported: internal function
solveOneFS <- function(mrkdids, ihl, FS, P, knownPhac, triedPhacs, maxparcombs,
                       ahcinfo, minfrac, errfrac, DRrate, minPseg) {

  calcParcombok <- function(potPhacs, oldPhacs, allhap, Pdids, FSfrac,
                            P.reject=c(0.05, 0.001), ahcinfo) {
    # Function within solveOneFS
    # Calculates the combinations of parental haplotype combinations
    # that can explain (almost) all FS individuals
    # potPhacs: list of 2 vectors with the columns of potential haplotype
    #           combinations for parents P1 and P2 (indexing columns in the
    #           ahclist items for their mrkdids)
    # oldPhacs: list of 2 vectors with columns of parental haplotype combinations
    #           that were already calculated and whose P1/P2 combinations were
    #           no solutions (may have length 0 if none)
    # allhap:   a matrix with all possible haplotypes with the nr of markers
    #           in the currect haploblock
    # FSfrac:   the fraction of FS individuals for which it is accepted
    #           that they are not explained by a parental haplotype combination
    # Pdids:    vector of 2 integers: the parental mrkdids
    # P.reject: two P values, decending: used in a binomial test for the
    #           nr of unexplained FS indiv
    # ahcinfo: a list returned by loadAllHapcombLists
    #TODO: parallellize; see also getAllHapcomb
    # Phacs are column numbers in ahc(complete)list
    if (all(potPhacs[[1]] %in% oldPhacs[[1]]) &&
        all(potPhacs[[2]] %in% oldPhacs[[2]])) {
       #no new parental combinations, so no solutions
      parcombok <- matrix(integer(0), ncol=2)
      setaside_pcs <- parcombok
    } else {
      missedFSs <- matrix(NA_integer_,
                          nrow=length(potPhacs[[1]]), ncol=length(potPhacs[[2]]))
      # missedFSs is the number of FS indivs with a mrkdid that cannot be
      # produced with each parental combination (assuming no DR)
      for (pc1 in seq_along(potPhacs[[1]])) {
        for (pc2 in seq_along(potPhacs[[2]])) {
          isnewpc <- !(potPhacs[[1]][pc1] %in% oldPhacs[[1]] &&
                         potPhacs[[2]][pc2] %in% oldPhacs[[2]])
          if (isnewpc) {
            #this is a new combination, check:
            FScombs <- getFScombs(cbind(
              getAllHapcomb(mrkdid=Pdids[1], nmrk=ncol(allhap),
                            ahcinfo=ahcinfo)[, potPhacs[[1]][pc1]],
              getAllHapcomb(mrkdid=Pdids[2], nmrk=ncol(allhap),
                            ahcinfo=ahcinfo)[, potPhacs[[2]][pc2]]))
            FSmrkdos <- hapcomb2mrkdos(FScombs, allhap=allhap)
            FScombmrkdids <- colnames(FSmrkdos)
            missedmrkdids <- setdiff(names(FSmrkdidsTable), FScombmrkdids)
            missedFSs[pc1, pc2] <-
              sum(FSmrkdidsTable[names(FSmrkdidsTable) %in% missedmrkdids])
          } else {
            #this is an old parcomb, either already rejected or set aside.
            #in both cases we don't need to check it again, so have it discarded:
            missedFSs[pc1, pc2] <- length(FS)
          }
        }
      }
      # explanation of the two levels of P.reject: see large comment at
      # the start of the main code of solveOneFS
      FSsize <- sum(FSmrkdidsTable) #nr of non-missing FS indiv
      T.reject <- qbinom(1-P.reject[1], FSsize, FSfrac) #reject  if missedFSs >= T.reject
      if (T.reject <= 1) {
        #always allow 1 missedFS:
        T.reject <- 2
      }
      parcombok <- missedFSs < T.reject
      parcombok <- which(parcombok, arr.ind=TRUE) #now a matrix with one row for
      # each valid combination of potPhacs[[1]] and potPhacs[[2]]
      # (indexing the PotPhacs). Change to indexing in ahclist:
      for (p in 1:2) parcombok[, p] <- potPhacs[[p]][parcombok[, p]]
      # now the indices to the ahc(complete)list columns for the Pdids
      # instead of indices to the potPhacs
      colnames(parcombok) <- c("P1col", "P2col")
      # also produce data for the parcombs passing P.reject[2]
      # but not P.reject[1]:
      T.rej2 <- qbinom(1-P.reject[2], FSsize, FSfrac) #reject  if missedFSs >= T.reject
      if (T.rej2 <= 1) {
        #always allow 1 missedFS:
        T.rej2 <- 2
      }
      # get the setaside_pcs analogous to parcombok
      setaside_pcs <- missedFSs >= T.reject & missedFSs < T.rej2
      setaside_pcs <- which(setaside_pcs, arr.ind=TRUE) #indexing PotPhacs
      for (p in 1:2) setaside_pcs[, p] <- potPhacs[[p]][setaside_pcs[, p]] #indexing ahclist
      colnames(setaside_pcs) <- c("P1col", "P2col")
    }
    list(parcombok=parcombok, setaside_pcs=setaside_pcs)
  } #calcParcombok within solveOneFS

  calcParcombstats <- function(pco, parcombok, Pdids, parcombstats,
                               didNhapcombLst, allhap, ahcinfo) {
    # pco: number, parental combination (index to parcombok rows)
    # Pdids: vector of 2 integers: the parental mrkdids
    # parcombstats: list with statistics for each parental combination, see
    #               return value
    # didNhapcombLst: list as in return value
    # return value: a list with 2 components:
    # $parcombstats: list with statistics for each parental combination:
    #   $nhap: a vector with the number of haplotypes for each pco
    #   $P: vector with the P-value  for each pco
    #   $toosmall: logical vector with one value per pco, TRUE if the FS is
    #              too small and/or too fragmented to calculate a P-value
    # $didNhapcombLst: a list with one component for each pco; each of these
    #   has the following components:
    #   $mrkdid: vector with names of mrkdids in the FS, in order of hapcomb
    #            (multiple hacs may have the same mrkdid!)
    #   $hapcomb: matrix of hacs expected in the FS (1 col per hac); with
    #             order of columns sorted so that all hacs with the same mrkdid
    #              are adjacent
    #   $expHCfreq: vector of the expected hac frequencies in order of hapcomb
    #   $dupdid: logical vector, TRUE where mrkdid[i] == mrkdid[i-1], else FALSE
    parhac <- cbind(
      getAllHapcomb(mrkdid=Pdids[1], nmrk=ncol(allhap),
                    ahcinfo=ahcinfo)[, parcombok[pco, 1]],
      getAllHapcomb(mrkdid=Pdids[2], nmrk=ncol(allhap),
                    ahcinfo=ahcinfo)[, parcombok[pco, 2]])
    # parhac does not contain NAs!
    segdat <- testFSseg(parhac=parhac, DRrate=DRrate, errfrac=errfrac,
                        FSmrkdidsTable=FSmrkdidsTable, allhap=allhap)
    parcombstats$nhap[pco] <- segdat$stats$nhap
    parcombstats$toosmall[pco] <- segdat$stats$toosmall
    parcombstats$P[pco] <- segdat$stats$P
    didNhapcombLst[[pco]] <- segdat$didNhapcombLst
    list(parcombstats=parcombstats, didNhapcombLst=didNhapcombLst)
  } #calcParcombstats in solveOneFS

  selFromParcombok <- function(parcombok, allhap, minPseg) {
    # parcombok: matrix with 2 columns (1 for each parent) and one row
    #      per parental combination; numbers index the columns of the
    #      parental mrkdid entry in ahc(complete)list
    # allhap: matrix of all possible haplotypes
    # minPseg: minimum chi-squared P-value to test the observed segregation
    # return: list of
    # $parcombstats for all parcombs in parcombok
    # $didNhapcombLst for all parcombs in parcombok
    # $sel: vector of indices to selected rows of parcombok (and to
    #       $parcombstats and $didNhapcombLst)
    if (nrow(parcombok) == 0) {
      return(list(sel=integer(0)))
    }
    parcombstats <-
      data.frame(P = rep(NA_real_, nrow(parcombok)),
                 nhap = rep(NA_integer_, nrow(parcombok)),
                 toosmall = rep(FALSE, nrow(parcombok)))
    didNhapcombLst <- list() #save calculated data that can be used later
    #                        to convert FS mrkdids into FS hapcombs
    for (pco in seq_len(nrow(parcombok))) {
      #calc parcombstats[pco,] and didNhapcombLst[[pco]] :
      tmp <- calcParcombstats(pco=pco, parcombok=parcombok, Pdids=Pdids,
                              parcombstats=parcombstats,
                              didNhapcombLst=didNhapcombLst,
                              allhap=allhap, ahcinfo=ahcinfo)
      parcombstats <- tmp$parcombstats
      didNhapcombLst <- tmp$didNhapcombLst
    } # for pco

    # now we have the parcombstats calculated for all pco.
    # we reject all that have P <= minPseg and
    # from the remaining we select the one with the "best" combination
    # of nhap and P
    firstP <- max(minPseg, 1e-4)
    if (any(!is.na(parcombstats$P) & parcombstats$P >= firstP)) {
      #if possible we limit our selection to those with P >= 1e-4:
      sel <- which(!is.na(parcombstats$P) & parcombstats$P >= firstP)
    } else {
      #else we select all with P >= minPseg:
      sel <- which(!is.na(parcombstats$P) & parcombstats$P >= minPseg)
    }
    list(sel=sel, parcombstats=parcombstats, didNhapcombLst=didNhapcombLst)
  } # selFromParcombok within solveOneFS

  # START of solveOneFS
  FSfrac <- errfrac + DRrate # the naive fraction "unexpected" FS individuals
  NoSolution <- FALSE
  maxP <- -1
  PahccolsTried <- list(integer(0), integer(0)) # to which oldPhacs must be added at the end
  FSmrkdidsTable <- table(mrkdids[names(mrkdids) %in% FS], useNA="no")
  FScount <- sum(FSmrkdidsTable)
  Pdids <- mrkdids[match(P, names(mrkdids))] #Pdids: Parental marker dosage IDs
  Phacs <- vector("list", 2) #Phacs: Parental haplotype combinations (column nrs in ahc(complete)list)
  for (p in 1:2) {
    if (is.na(Pdids[p])) {
      Phacs[[p]] <- integer(0)
    } else {
      #get the parental haplotype combinations based on ihl:
      if (is.na(knownPhac[p])) {
        i <- which(names(ihl$hclist) == Pdids[p])
        Phacs[[p]] <- ihl$hclist[[i]] # a subset of all possible hacs
      } else {
        Phacs[[p]] <- knownPhac[p]
      }
    }
  }
  sel <- integer(0)
  if (anyNA(Pdids) || FScount == 0) {
    message <- ("missing marker data in parents or no FS data; skipped")
    NoSolution <- TRUE
  } else {
    message <- ""
    # We consider two sets of parental hapcombs: first the
    # ones based on the ihl obtained ignoring any population structure:
    # this leads to a subset of possible hapcombs for each parent. And
    # second, if no suitable parental hapcombs are found in the first set,
    # we consider all possible parental hapcombs (from the ahc(complete)list);
    # of course we skip the ones that were already considered in the first set.
    # In the first and the second set, we first only check which fraction
    # of the FS cannot be explained by each parcomb (ignoring DR) - this
    # is a relatively fast check. We do a binomial test and only keep the
    # parcombs where the probability of that amount of "missingFS" is less
    # than P.reject[1] (=0.05). The parcombs that pass this first test
    # are next checked further; the segregation of mrkdids in the FS is
    # tested against the expected segregation using a chi-squared test, and
    # this is combined with a P-value for unexpected mrkdids. Here DR is
    # taken into account.
    # We use two P.reject levels: a strict (0.05 as before) and
    # a more lenient (0.001) in the first filtering. Now we consider
    # successively the segregation of the parcombs: (1) that result from
    # the restricted ihl set and pass P.reject 0.05, (2) that result from
    # the complete ahclist set and pass P.reject 0.05, (3) that result from the
    # restricted ihl set and pass P.reject 0.001, (4) that result from the
    # complete ahclist set and pass P.reject 0.001. As soon as one of these
    # 4 sets yields acceptable parcombs the next are all skipped.
    # This procedure is quite complex, but it will allow more often a solution
    # while not increasing the computations too much, and not accepting a bad
    # (incorrect) solution from the ihl set if a better one exists.
    firstcombs <- length(Phacs[[1]]) * length(Phacs[[2]]) # from ihl
    if (firstcombs == 0) stop("firstcombs==0, should not happen!")
    if (firstcombs > maxparcombs) {
      # too many potential combinations (it takes about 45 min to check
      # 150000 combinations)
      message <- paste("more than", maxparcombs, "parental combinations; skipped")
      parcombok <- matrix(integer(0), ncol=2)
    } else {
      parcombok <- calcParcombok(potPhacs=Phacs,
                                 oldPhacs=PahccolsTried,
                                 allhap=ihl$allhap, Pdids=Pdids,
                                 FSfrac=FSfrac, ahcinfo=ahcinfo)
      PahccolsTried <- Phacs
      step1pc <- parcombok$setaside_pcs
      parcombok <- parcombok$parcombok
      # note that in this first try the oldPhacs have 0 columns as
      # PahccolsTried is created empty at start of SolveOneFS
      tmp <- selFromParcombok(parcombok=parcombok, allhap=ihl$allhap,
                              minPseg=minPseg)
      sel <- tmp$sel; parcombstats <- tmp$parcombstats;
      didNhapcombLst <- tmp$didNhapcombLst
      maxP <- max(maxP, parcombstats$P, na.rm=TRUE)
      #print(paste("firstcombs, firstsel: P=", parcombstats$P, "length(sel)=", length(sel)))
      if (length(sel) == 0) {
        #No parental combinations from inferHaps_noFS, preselected with
        #stringent P.reject, fit
        #1. See how many parental combinations are possible without first
        #   inferHaps_noFS (i.e. with getAllHapcomb).
        #   If p1*p2 < maxparcombs we try them all.
        nwPhacs <- vector("list", 2)
        for (p in 1:2) {
          if (is.na(knownPhac[p])) {
            nwPhacs[[p]] <-
              1:ncol(getAllHapcomb(mrkdid=Pdids[p], nmrk=ncol(ihl$allhap),
                                   ahcinfo=ahcinfo))
          } else {
            nwPhacs[[p]] <- knownPhac[p]
          }
        }
        extracombs <-
          length(nwPhacs[[1]]) * length(nwPhacs[[2]]) - firstcombs
        if (extracombs > maxparcombs) {
          message <- paste("more than", maxparcombs,
                           "extra parental combinations; skipped")
          parcombok <- step2pc <- matrix(nrow=0, ncol=2,
                                         dimnames=list(NULL,c("col1", "col2")))
          sel <- integer(0)
        } else {
          # we try all possible parcombs
          # calculate for extra parcombs:
          parcombok <- calcParcombok(potPhacs=nwPhacs,
                                     oldPhacs=PahccolsTried,
                                     allhap=ihl$allhap, Pdids=Pdids,
                                     FSfrac=FSfrac, ahcinfo=ahcinfo)
          for (p in 1:2)
            PahccolsTried[[p]] <- sort(union(PahccolsTried[[p]], nwPhacs[[p]]))
          step2pc <- parcombok$setaside_pcs
          parcombok <- parcombok$parcombok
          tmp <- selFromParcombok(parcombok=parcombok, allhap=ihl$allhap,
                                  minPseg=minPseg)
          sel <- tmp$sel; parcombstats <- tmp$parcombstats;
          didNhapcombLst <- tmp$didNhapcombLst
          maxP <- max(maxP, parcombstats$P, na.rm=TRUE)
        }
        if (length(sel) == 0) {
          # now we add two more sets of candidate parcombs: those that were
          # set aside in the firstcombs and in the extracombs check.
          # (i.e. those that failed the preselection with P.reject[1] but
          # passed with P.reject[2])
          parcombok <- step1pc # the setaside_pcs from the firstcombs set
          tmp <- selFromParcombok(parcombok=parcombok, allhap=ihl$allhap,
                                  minPseg=minPseg)
          sel <- tmp$sel; parcombstats <- tmp$parcombstats;
          didNhapcombLst <- tmp$didNhapcombLst
          maxP <- max(maxP, parcombstats$P, na.rm=TRUE)
          if (length(sel) > 0) {
            # remove a possible message about too many extracombs,
            # else the FS will be considered as not fitted:
            message <- ""
          } else {
            # length(sel) == 0
            parcombok <- step2pc
            tmp <- selFromParcombok(parcombok=parcombok, allhap=ihl$allhap,
                                    minPseg=minPseg)
            sel <- tmp$sel; parcombstats <- tmp$parcombstats;
            #print(paste("extracombs, secondsel: P=", parcombstats$P, "length(sel)=", length(sel)))
            didNhapcombLst <- tmp$didNhapcombLst
            maxP <- max(maxP, parcombstats$P, na.rm=TRUE)
          } # no firstcombs remain if preselected with P.reject[2]
        } # no extracombs remain if preselected with P.reject[1]
      } # no firstcombs remain if preselected with P.reject[1]
    } # not too many firstcombs
  } # parental and FS data available

  # now we have a parcombok with matching sel, parcombstats and didNhapcombLst
  if (maxP == -1) maxP <- NA
  if (length(sel) == 0) {
    if (message=="") message <- "no valid FS solutions found"
    NoSolution <- TRUE #after analyzing all possible parental combinations
  }
  if (message != "") {
    return(list(message=message,
                Pdids=Pdids,
                parcombok=matrix(integer(0), ncol=2),
                NoSolution=NoSolution,
                SolFound=FALSE,
                Pseg=NA,
                maxPseg=maxP,
                PahccolsTried=PahccolsTried))
  }
  # now we have acceptable solutions: sel refers to the acceptable rows of
  # parcombok

  # SECOND STEP: find the best fitting of the acceptable parental haplotype
  # combinations.
  # parcombok indicates which combinations of these match (almost) all FS's.
  # Now, any of these combinations could be the correct one. We check all
  # and select the best, based on match with expected segregation (including
  # the expected error and DR rates) and the number of haplotypes.
  if (length(sel) == 1) {
    sel3 <- sel
  } else {
    #among these we need to select the best solution(s), based on P and nhap.
    #If we take one extra haplotype to be as bad as a 5-fold lower P,
    #we can use the selection criterion:
    selcrit <- parcombstats$nhap[sel] - log(parcombstats$P[sel], base=5)
    if (sum(!is.na(selcrit)) == 0) stop("error in selcrit")
    toler <- log(0.5, base=5) #tolerance: select all that have a P-value >= 0.5 * highest P
    # old: sel2 <- which.min(selcrit) # sel2 indexes the best among the sel
    sel2 <- which(selcrit >= min(selcrit) + toler) #multiple equivalent solutions possible
    sel3 <- sel[sel2] #sel3 indexes the best among the original parcombs
  }
  return(list(message="",
              Pdids=Pdids,
              parcombok=parcombok[sel3,, drop=FALSE],
              didNhapcomb=didNhapcombLst[sel3],
              NoSolution=FALSE,
              SolFound=TRUE,
              Pseg=parcombstats$P[sel3],
              maxPseg=max(parcombstats$P[sel3]),
              PahccolsTried=PahccolsTried))
} #solveOneFS

checkpops <- function(parents, FS, mrkDosage) {
  message <- ""
  if (missing(parents) || is.null(parents) ||
      (is.matrix(parents) && nrow(parents)==0))
    parents <- matrix(character(0), ncol=2)
  if (missing(FS) || is.null(FS)) FS <- list()
  if (!is.list(FS))
    message <- "FS must be a list"
  if (!is.matrix(parents) || ncol(parents) != 2)
    message <- "parents must be a matrix with 2 columns"
  if (length(FS) != nrow(parents))
    message <- "FS does not match parents"
  FS <- lapply(FS, as.character)
  allFSind <- unlist(FS)
  if (anyNA(allFSind) || !all(allFSind %in% colnames(mrkDosage)))
    message <- "All FS indiv must occur in mrkDosage"
  parents <- matrix(as.character(parents), ncol=2)
  if (anyNA(parents) || !all(parents %in% colnames(mrkDosage)))
    message <- "All parents must occur in mrkDosage"
  list(parents=parents, FS=FS, message=message)
} #checkpops

checkHaploblock <- function(haploblock, mrkDosage) {
  # haploblock is a list of haploblocks, or a vector of marker names in
  # a single haploblock.
  # Returns the haploblock list with empty haploblocks deleted and with
  # names added to nameless haploblocks.
  # Produces an error if not all markers occur in matrix mrkDosage (the
  # reverse is allowed: mrkDosage may contain markers not in haploblock)
  if (is.list(haploblock)) {
    haploblock <- haploblock[sapply(haploblock, length) > 0]
    mrknames <- unlist(haploblock)
    if (anyNA(mrknames) || !all(mrknames %in% rownames(mrkDosage)))
      stop("Not all markers in haploblock occur in mrkDosage")
  } else {
    #haploblock not a list: then must be the name of the only
    #haploblock present, which contains all markers in mrkDosage
    if (!is.vector(haploblock) || length(haploblock) != 1)
      stop("Haploblock invalid")
    hbname <- as.character(haploblock)
    haploblock <- list(rownames(mrkDosage))
    names(haploblock) <- hbname
  }
  if (is.null(names(haploblock)))
    names(haploblock) <- paste0("hb", padded(seq_len(length(haploblock))))
  nonames <- which(is.na(names(haploblock)) | names(haploblock)=="")
  if (length(nonames) > 0)
    names(haploblock)[nonames] <- paste0("hb",padded(nonames))
  haploblock
} #checkHaploblock

# hapOneBlock approach
# A - For all FSs calculate the possible solutions, but dont apply any of them
#     From large to small, without grouping
# B - use all possible parental haplotypes as "known" to re-do
#     the skipped FSs
# C - Per group of related FSs with one or more solutions:
#     loop: find and resolve conflicts between FSs:
#     where both FSs have one or more solutions but the parent(s) cannot have
#     the same hac
#     For each parental solution, see for how many total FS progeny (over
#     multiple FSs) it is acceptable. Discard the one with the smallest total
#     FS progeny;
#     check for which FSs this was the only solution for that parent and
#     remove these FSs from further consideration (i.e. treat as unrelated).
#     Continue until no further conflicts: all remaining parental solutions
#     fit all remaining FSs where they are a parent
#     We do this over all parents at the same time, i.e. each time for
#     all parents calculate the solution applicable to the smallest total
#     FS progeny and reject that
# D - While there are parents with multiple solutions left (these fit all
#     remaining FSs of which they are parents, so these FSs also have multiple
#     solutions):
#     For each combination of remaining parental solutions over all remaining
#     parents calculate an overall P-value (by multiplying the P-values of all
#     remaining FSs for that combination).
#     If any single one of these combinations of two parental hacs has the
#     largest P-value, select that and apply it to all these remaining FSs;
#     else treat all remaining FSs as unrelated.
# E - solve all unrelated material (incl FSs for which no solutions found,
#     skipped FSs and FSs set aside during C, including their parents if not
#     solved yet;
#     for this use the original knownhap and the haplotypes of the ACCEPTED
#     FS solutions (i.e. not all the ones present at the end of step C,
#     which might be conflicting between FSs in the same group)

hapOneBlock <- function(mrkDosage, hbname, ahcinfo, parents, FS,
                        minfrac, errfrac, DRrate,
                        dropUnused, maxparcombs, minPseg,
                        knownHap, progress) {
  # called by inferHaplotypes
  # mrkDosage is now checked, a matrix with only the rows for the markers
  #           in the current haploblock
  # hbname: the name of the haploblock
  # all other parameters as in inferHaplotypes, plus ahcinfo: a list
  # returned by loadAllHapcombLists at start of inferHaplotypes
  # parents and FS are checked, and FS has length 0 if no FS families defined
  # Return value is a list with elements:
  # $hapdos: a matrix with one column per individual and one row per
  #          haplotype, with the dosages och each haplotypes in each indiv
  #          (if dropUnused, only the rows for the haplotypes that occur
  #           at least once in the population are retained). The columns are
  #           named as the individuals, the rows have the haplotype names.
  # $mrkdids: the observed mrkdids for all individuals, in same order as the
  #           columns of $hapdos; names are individual names
  # $allhap: a matrix with the composition of all haplotypes
  # $FSfit: logical vector with one element per FS family; TRUE if a (or
  #          more than one) acceptable solution for the FS is found (although
  #          if more than one solution they might not be used if unclear
  #          which is the best solution)
  # $FSmessages: character vector with one message per FS. If not FSfit then
  #          an error message; if FSfit either "" or a report of mrkdids that
  #          were not expected or that were not assigned a haplotype combination
  #          because multiple possible (and likely) combinations would give the
  #          same mrkdid
  # $FSpval: numeric vector with a P value for each FS. If one solution was
  #         chosen, the chi-squared P-value of that solution; else the maximum
  #         P-value over all considered solutions

  getFSgroupStructure <- function(FSgrouping, grp, FSdata, pcos, parents) {
    # FSgrouping: a list as produced by
    # grp: index to one group in FSgrouping
    # FSdata: list in hapOneBlock with one element per FS, we use here
    #         the $parcombok
    # pcos: a list with one element per FS, each of which is a vector indexing
    #       the still valid rows of the parcombok for that FS
    # parents: a matrix with the nams of the two parents of each FS
    # for one group of related FSs produce a list with item:
    # $availFSs: integer vector with FS numbers:
    #         which FSs belong to the group and still have at least
    #         1 solution (at least 1 row in FSsof[[fs]]$parcombok)
    # $availParents: character vector with the names of all parents of the
    #         availFSs
    # $FSmatlist: list of matrices, one for each availFS, with one row for each
    #         parent occurring in any of the availFSs, and one column for each
    #         parcomb (combination of 2 parental hacs).
    #         The matrices are filled with the ahccols of the hacs.
    #         (i.e. in each matrix only the 2 rows of the parents of that FS
    #         are not NAs)
    # $parhapcombs: list of vectors, one for each availParent, with all
    #         hacs (haplotype combinations) occurring in any of the
    #         solutions for any availFS for that parent.
    availFSs <- integer(0) #all FSs in the grp that still have some solutions
    # for (fs in FSgrouping$FSgroups[[grp]]) {
    #   #if (nrow(FSdata[[fs]]$parcombok) > 0) availFSs <- c(availFSs, fs)
    #   if (length(pcos[[fs]]) > 0) availFSs <- c(availFSs, fs)
    # }
    pcocounts <- sapply(pcos, length)
    availFSs <- intersect(FSgrouping$FSgroups[[grp]], which(pcocounts > 0))
    #           all FSs in group with some solutions left
    # the rest also works if no FSs left:
    availparents <- sort(unique(as.vector(parents[availFSs,]))) # names of
    #         parents of available (non-conflicting) FSs
    FSmatlist <- vector("list", length(availFSs)) #for each FS still
    #         available in the group, we make a matrix with one row for each
    #         parent occurring in any of them, and one column for each
    #         parcomb (combination of 2 parental hacs).
    #         The matrices are filled with the ahccols of the hacs,
    for (avfs in seq_along(availFSs)) {
      fs <- availFSs[avfs]
      #pccount <- nrow(FSdata[[fs]]$parcombok) # count of parental
      #       combinations for the current fs
      FSmatlist[[avfs]] <- matrix(integer(0), nrow=length(availparents),
                                  ncol=pcocounts[fs])
      rownames(FSmatlist[[avfs]]) <- availparents
      for (pc in seq_len(pcocounts[fs])) {
        prows <- match(parents[fs,], availparents)
        FSmatlist[[avfs]][prows, pc] <-
          FSdata[[fs]]$parcombok[pcos[[fs]][pc],] #already ahccols
      }
    } #for avfs

    # now we have the FSmatlist filled and we can directly
    # compare the ahccols for each parent between FSs
    # the FSmatlist matrices have NAs for the parents that don't occur
    # in that FS (i.e. only 2 rows non-NA)
    # Now we create a list with for each availparent the unique hapcombs
    # over all availFSs:
    parhapcombs <- vector("list", length(availparents))
    for (prow in seq_along(availparents)) {
      parhapcombs[[prow]] <- integer(0)
      for (avfs in seq_along(availFSs)) {
        parhapcombs[[prow]] <-
          c(parhapcombs[[prow]], FSmatlist[[avfs]][prow,])
      }
      parhapcombs[[prow]] <-
        sort(unique(parhapcombs[[prow]][!is.na(parhapcombs[[prow]])]))
    } # for prow
    names(parhapcombs) <- availparents
    list(availFSs=availFSs, availparents=availparents,
         FSmatlist=FSmatlist, parhapcombs=parhapcombs)
  } #getFSgroupStruct within hapOneBlock

  getAllParCombs <- function(parents, FSdata, FSgrpstruct) {
    # produce a list with
    # $allparcombs: a matrix with 1 parhapcomb per row with ahccols of all
    #               available parents, one column per available parent
    # $combPval: ("combined P-value") the product of the P-values over all FSs
    #            for the corresponding row of allparcomb
    allparcombs <- as.matrix(expand.grid(FSgrpstruct$parhapcomb)) #matrix with one row
    #              for each possible combination of ahccols of all parents;
    #              parents in columns, in order of parhapcomb
    Pval <- rep(1.0, nrow(allparcombs))
    for (fs in FSgrpstruct$availFSs) {
      ok <- rep(FALSE, nrow(allparcombs))
      parcols <- match(parents[fs,], colnames(allparcombs))
      for (cmb in seq_len(nrow(FSdata[[fs]]$parcombok))) {
        parahccols <- FSdata[[fs]]$parcombok[cmb,]
        okrow <- which(allparcombs[, parcols[1]] == parahccols[1] &
                       allparcombs[, parcols[2]] == parahccols[2])
        ok[okrow] <- TRUE
        Pval[okrow] <- Pval[okrow] * FSdata[[fs]]$Pseg[cmb]
      }
      allparcombs <- allparcombs[ok,, drop=FALSE]
      Pval <- Pval[ok]
    }
    list(allparcombs=allparcombs, Pval=Pval)
  } #getAllParCombs within hapOneBlock

  findParcombokRow <- function(allparcombs, selcomb, parents, FSsof, fs) {
    # allparcombs: matrix with all parents in group in columns and
    #             group solutions in rows, containing the parental hapixes
    # selcomb: one number: a row in allparcombs
    # parents: the 2-column matrix defining the parents of each FS
    # FSsof: the list with all solveOneFS results
    # fs: the fs we want
    # return value: the row number of FSsof$parcombok corresponding to
    #               the parental solution selcomb
    apccols <- match(parents[fs,], colnames(allparcombs))
    # for debugging:
    if (length(apccols) != 2 || anyNA(apccols)) stop("apccols error")
    pcorow <- which(FSsof[[fs]]$parcombok[, 1] ==
                      allparcombs[selcomb, apccols[1]]
                       &
                      FSsof[[fs]]$parcombok[, 2] ==
                      allparcombs[selcomb, apccols[2]])
    # for debugging:
    if (length(pcorow) != 1) stop ("pcorow error")
    pcorow
  } #findParcombokRow within hapOneBlock

  selOneMrkdidHac <- function(FSsof, pco, expdid) {
    # FSsof: list returned by solveOneFS
    # pco: the selected solution (row in FSsof$parcombok)
    # expdid: vector of all unique expected mrkdids in this FS
    # For the indicated solution FSsof$didNhapcomb[[pco]],
    # find out which hacs (haplotype combinations) correspond to the
    # given mrkdid (names of marker dosage combinations). If more than one,
    # select the one with the highest probability if much more probable than
    # the others, else none
    # return value: a matrix with ploidy rows and one column for each did in
    # expdid, will contain the unique solution (hac) for that did, or NA if
    # not one most likely hac. The names are the expdid.
    # TODO debug check:
    if (anyNA(expdid) || length(unique(expdid)) != length(expdid))
      stop("invalid expdid")
    hacmatrix <- matrix(integer(0),
                         nrow=nrow(FSsof$didNhapcomb[[1]]$hapcomb),
                         ncol=length(expdid),
                         dimnames=list(NULL, expdid))
    for (md in seq_along(expdid)) {
      haccols <- which(FSsof$didNhapcomb[[pco]]$mrkdid == expdid[md])
      if (length(haccols) == 0) {
        #results$didproblems[1, md] <- TRUE
      } else if (length(haccols) == 1) {
        #only one hapcomb corresponds to mrkdid, assign:
        hacmatrix[, md] <- FSsof$didNhapcomb[[pco]]$hapcomb[, haccols]
      } else { # length(haccols) > 1
        #multiple hapcomb correspond to same mrkdid. Assign only if one
        #hapcomb has a very high probability compared with the rest (to avoid
        #e.g. that DR solutions get in the way of assigning a good solution)
        didprobs <- FSsof$didNhapcomb[[pco]]$expHCfreq[haccols]
        didprobs <- didprobs / sum(didprobs)
        srtprobs <- sort(didprobs, decreasing=TRUE)
        if ((srtprobs[1] >= 0.8 && srtprobs[2] <= 0.05) ||
            (srtprobs[1] >= 0.5 && srtprobs[2] <= 0.01)) {
          # one most probable hac for this mrkdid, select:
          sel <- which(didprobs >= 0.5)
          hacmatrix[, md] <- FSsof$didNhapcomb[[pco]]$hapcomb[, haccols[sel]]
        } else {
          # not one single most likely hapcomb for this did,
          # leave result$hacmatrix column as NAs
        }
      } # multiple hacs for mrkdid
    } # for md
    hacmatrix
  } # selOneMrkdidHac within hapOneBlock

  newSolstats <- function(knownHap, FS) {
    stats <- list( # the initial "results" before cycle 1
      knownHap=knownHap, #the haplotypes present in the assigned parental
      # hapcombs and/or in the specified knownhaps
      Pahccols=matrix(NA_integer_, nrow=nrow(parents), ncol=2), # the assigned
      # parental hapcombs that are identical in all solutions, as ahccols -
      # these are considered certain and used as basis for next cycle
      selPahccols=matrix(NA_integer_, nrow=nrow(parents), ncol=2), # as
      # Pahccols, but now including the selected solutions for those that
      # were not the same over all solutions - these are not necessarily
      # certain and only used if this cycle is chosen as bestcycle
      pcos=vector("list", length(FS)), # for each FS a vector with the remaining
      # pcos (parental haplotype combinations), indexing FSdata[[fs]]$parcombok
      # and $Pseg
      selpco=rep(NA_integer_, length(FS)), # the selected pco, if any
      solvedFSindiv=0 # the number of FS indiv with mrkdids (i.e. with complete
      # marker data) in the FSs for which a solution is assigned
      # TODO: consider if we can use an extra statistic to select between
      # cycles with the same maximum solvedFSindiv, if the cycles don't converge
      # to a single solution
    )
    for (fs in seq_along(FS)) stats$pcos[[fs]] <- integer(0)
    stats
  } # function newSolstats

  # hapOneBlock START ####
  allmrkdids <- mrkdos2mrkdid(mrkDosage, ploidy=ahcinfo$ploidy, check=FALSE)
  if (length(FS) == 0) {
    ihl <- inferHaps_noFS(mrkdids=allmrkdids, ahcinfo=ahcinfo,
                           markers=rownames(mrkDosage),
                           minfrac=minfrac,
                           knownHap=knownHap, useKnownDosage=TRUE,
                           progress=progress)
    hapdos <- hapcomb2hapdos(
      hapcomb=hapcomb_from_IHlist(ihl=ihl, mrkdids=allmrkdids,
                                  ahcinfo=ahcinfo, nocombs=FALSE)$hapcomb,
      nhap=nrow(ihl$allhap))
    result <- list(
      hapdos=hapdos, mrkdids=allmrkdids, markers=colnames(ihl$allhap))
  } else {
    # one or more FS
    #first calculate ihl without using FS info:
    # - for elements of the returned list
    # - for a pre-selection of likely parental hapcomb
    ihl <- inferHaps_noFS(mrkdids=allmrkdids, ahcinfo=ahcinfo,
                           markers=rownames(mrkDosage),
                           minfrac=minfrac[1], #without final inference
                           knownHap=knownHap, useKnownDosage=FALSE,
                           progress=progress)
    ahccols <- hapcomb_from_IHlist(ihl=ihl, mrkdids=allmrkdids,
                                   ahcinfo=ahcinfo, nocombs=TRUE)$ahccols
    ihl <- ihl[names(ihl) != "hclist"] #sometimes big, not used here, delete
    nhap <- nrow(ihl$allhap)
    # we calculate the FS sizes per haploblock because of differences
    # in missing data:
    indiv <- colnames(mrkDosage)
    mrkcomplete <- !is.na(colSums(mrkDosage)) #for each indiv: complete marker dosage info?
    FSgrouping <- getLinkedFSgroups(parents=parents)

    # We perform a looping where in each cycle we find the optimal solution
    # per group of linked FSs, which yields
    # - the haplotypes present in the parents of the known FSs
    # - the solutions for the parents of the solved FSs
    # - the number of FS indiv with complete marker genotypes over all solved
    #   FSs
    FSdata <- vector(mode="list", length(FS)) #some accumulated data for each FS
    # over all cycles, as follows:
    for (fs in seq_along(FSdata)) {
      FSdata[[fs]]$mrksize <-
        sum(!is.na(allmrkdids[names(allmrkdids) %in% FS[[fs]]]))
      #  the number of indiv in the FS with complete marker info
      FSdata[[fs]]$message <- "" # will be about skipped or about FS-mrkdids with
      #  multiple possible hapcombs
      FSdata[[fs]]$NoSolution <- FALSE # will only be TRUE
      #  when all possible parcombs were checked and failed, i.e. when
      #  firstcombs and extracombs not skipped but no acceptable solution
      FSdata[[fs]]$SolFound <- FALSE # will remain FALSE until at least one
      #  acceptable solution found for this FS. Will then always remain TRUE
      #  even if all solutions in parcombok are discarded later due to conflicts
      #  with other FSs in the group
      FSdata[[fs]]$parcombok <- matrix(NA_integer_, nrow=0, ncol=2) # will remain
      #  empty (0 rows) until at least one acceptable solution found; and then
      #  not recomputed or changed anymore in subsequent cycles
      FSdata[[fs]]$Pseg <- numeric(0) # the chi-squared P-values for the mrkdid
      #  segregation for each row of parcombok; will be filled in at the same
      #  time as parcombok.
      FSdata[[fs]]$maxPseg <- NA # the maximum of $Pseg or if no solutions found
      #  the maximum of the discarded solutions (insufficient for acceptance).
      #  NA as long as no solutions considered
      FSdata[[fs]]$PahccolsTried <- list(P1=integer(0), P2=integer(0)) # the
      #  Pahccols of which the parcombs were checked; empty except when no
      #  acceptable solution found and some but but not all possible parhapcomb
      #  checked (i.e. firstcombs done but extracombs skipped; can be extended
      #  over multiple cycles but emptied when NoSolution TRUE or parcombok not
      #  empty anymore
    }

    solstats <- list() # for each cycle this will store the solutions and their
    # statistics over all FSs. "cycle 1" stores the initial data.
    solstats[[1]] <- newSolstats(knownHap, FS) # the initial "results" are
    # stored as the results of cycle 1; the first true cycle with computations
    # is cycle 2
    cycle <- 1 # "cycle 1" is the start situation
    while (TRUE) {
      cycle <- cycle + 1
      # did we converge to a stable solution or to a repeating set of different
      # solutions?
      if (cycle > 2) {
        bestcycle <- NA
        if (identical(solstats[[cycle-1]], solstats[[cycle-2]])) {
          bestcycle <- cycle - 1
        } else {
          oldcycle <- cycle - 3
          while (oldcycle > 0 &&
                 !identical(solstats[[cycle - 1]], solstats[[oldcycle]]))
            oldcycle <- oldcycle - 1
          if (oldcycle > 0) {
            # oldcycle had same solution as cycle-1: find the best cycle
            # in [oldcycle ...cycle-2]:
            maxSolvedFSindiv <- solstats[[oldcycle]]$solvedFSindiv
            minhapcount <- length(solstats[[oldcycle]]$knownhap)
            bestcycle <- oldcycle
            while (oldcycle < cycle - 2) {
              oldcycle <- oldcycle + 1
              if (solstats[[oldcycle]]$solvedFSindiv > maxSolvedFSindiv ||
                  (solstats[[oldcycle]]$solvedFSindiv == maxSolvedFSindiv &&
                   length(solstats[[oldcycle]]$knownhaps) < minhapcount)) {
                # we find a cycle with a better solution in terms of
                # solvedFSindiv and haplotype count
                # if two of the repeating cycles have the same solvedFSindiv
                # and hapcount we use the first (arbitrarily).
                # TODO: we could also discard the entire FS-based solution of
                # this haploblock but that seems a waste. Or we may think of
                # an additional statistic. How often will this occur anyway?
                maxSolvedFSindiv <- solstats[[oldcycle]]$solvedFSindiv
                minhapcount <- length(solstats[[oldcycle]]$knownhap)
                bestcycle <- oldcycle
              }
            }
          }
        }
        if (!is.na(bestcycle)) break
      }
      # No break, go on with new cycle
      solstats[[cycle]] <- newSolstats(knownHap, FS) #only the specified haplotypes
      # 1: for each FS find a set of acceptable solutions
      ih <- inferHaps_noFS(mrkdids=allmrkdids, ahcinfo=ahcinfo,
                            markers=rownames(mrkDosage),
                            minfrac=minfrac[1],
                            knownHap=solstats[[cycle-1]]$knownHap,
                            useKnownDosage=FALSE)
      # find acceptable solutions per FS in this cycle:
      for (fs in seq_along(FS)) {
        if (!FSdata[[fs]]$NoSolution && !FSdata[[fs]]$SolFound) {
          sof <- solveOneFS(mrkdids=allmrkdids, ihl=ih, FS=FS[[fs]],
                            P=parents[fs,],
                            knownPhac=solstats[[cycle-1]]$Pahccols[fs,],
                            triedPhacs=FSdata[[fs]]$PahccolsTried,
                            maxparcombs=maxparcombs,
                            ahcinfo=ahcinfo, minfrac=minfrac, errfrac=errfrac,
                            DRrate=DRrate, minPseg=minPseg)
          # we update FSdata[[fs]]:
          FSdata[[fs]]$message <- sof$message
          FSdata[[fs]]$NoSolution <- sof$NoSolution
          FSdata[[fs]]$SolFound <- sof$SolFound
          FSdata[[fs]]$parcombok <- sof$parcombok
          FSdata[[fs]]$Pseg <- sof$Pseg
          FSdata[[fs]]$maxPseg <- sof$maxPseg #different from max(Pseg) if no solution
          #FSdata[[fs]]$didNhapcombLst <- sof$didNhapcombLst #20151115 incorrect? sof has didNhapcomb, not List
          # a list with for each pco of this fs a list with for each FS mrkdid
          # a vector of the possible ahccols
          solstats[[cycle]]$pcos[[fs]] <- seq_len(nrow(FSdata[[fs]]$parcombok))
          for (p in 1:2) FSdata[[fs]]$PahccolsTried[[p]] <-
            sort(union(FSdata[[fs]]$PahccolsTried[[p]],
                       sof$PahccolsTried[[p]]))
        } else {
          solstats[[cycle]]$pcos[[fs]] <- solstats[[cycle-1]]$pcos[[fs]]
        }
      } # for fs
      # find the best solution over all FSs in this cycle:

      # Per group of linked FSs with one or more solutions:
      # loop: find and resolve conflicts between FSs:
      # where both FS have one or more solutions but the parent(s) cannot
      # have the same haplotype combination
      # For each parental solution, see for how many total FS progeny
      # (over multiple FSs) it is acceptable. Discard the one with the
      # smallest total FS progeny; check for which FSs this was the only
      # solution for that parent and remove these FSs from further
      # consideration (i.e. treat as unrelated).
      # Continue until no further conflicts: all remaining parental
      # solutions fit all remaining FSs where they are a parent.
      # We do this over all parents at the same time, i.e. each time
      # for all parents calculate the solution applicable to the smallest
      # total FS progeny and reject that.
      for (grp in seq_along(FSgrouping$FSgroups)) {
        FSconflicts <- TRUE
        while (FSconflicts) {
          FSconflicts <- FALSE
          FSgrpstruct <-
            getFSgroupStructure(FSgrouping=FSgrouping, grp=grp, FSdata=FSdata,
                                pcos=solstats[[cycle]]$pcos, parents=parents)
          # It is possible that no FSs are left in a group, if none of the FSs
          # has at least one solution. In that case we skip the rest of the loop
          # and leave FSconflicts FALSE to exit the loop.
          if (length(FSgrpstruct$availFSs) > 1) {
            # find and resolve conflicts between FSs:
            # where both FS have one or more solutions but the shared parent
            # cannot have the same hac.
            # For each parental solution, see for how many total FS progeny (over
            # multiple FSs) it is acceptable.
            # Discard the one with the smallest total FS progeny;
            # check for which FSs this was the only solution for that parent and
            # remove these FSs from further consideration (i.e. treat as unrelated).
            # Continue until no further conflicts: all remaining parental solutions
            # fit all remaining FSs where they are a parent
            miscount <- vector("list", length(FSgrpstruct$availparents)) # for each
            #           availparent a vector with for each hac the nr of FS progeny
            #           it does not match
            names(miscount) <- FSgrpstruct$availparents
            for (prow in seq_along(FSgrpstruct$availparents)) {
              # prow indicates the parent, row number in the FSgrpstruct$FSmatlist's
              # (and in allparcombs, but that doesn't exist at this stage)
              miscount[[prow]] <- rep(0, length(FSgrpstruct$parhapcombs[[prow]]))
              for (cix in seq_along(FSgrpstruct$parhapcombs[[prow]])) {
                # cix is the index of the current parental ahccol in
                # FSgrpstruct$parhapcombs[[prow]]
                phcmb <- FSgrpstruct$parhapcombs[[prow]][cix]
                # phcmb is the actual parental ahccol indicated by cix
                # now, for all FSs in grp:
                for (avfs in seq_along(FSgrpstruct$availFSs)) {
                  fs <- FSgrpstruct$availFSs[avfs]
                  if (FSgrpstruct$availparents[prow] %in% parents[fs,]) {
                    # if the prow parent is a parent of this FS ...
                    if(!(phcmb %in% FSgrpstruct$FSmatlist[[avfs]][prow,])) {
                      # ... and this ahccol for this parent is not in one of the
                      #     solutions of this FS, then this parental solution
                      #     "misses" all GENOTYPED indivs of the entire FS
                      miscount[[prow]][cix] <-
                        miscount[[prow]][cix] + FSdata[[fs]]$mrksize
                    }
                  }
                } # for avfs
              } # for cix
            } # for prow
            #now for each parental combination we know how many fully genotyped
            # FS progeny it misses
            #find the (any) maximum:
            maxmissedprogcount <- sapply(miscount, max)
            FSconflicts <- any(maxmissedprogcount > 0)
            if (!FSconflicts) break #ends the "while FSconflicts" loop
            maxprow <- which.max(maxmissedprogcount) #the parent with the most
            #                                missed progeny in its worst solution
            maxcix <- which.max(miscount[[maxprow]]) #the index to the solution itself
            maxphcmb <- FSgrpstruct$parhapcombs[[maxprow]][maxcix] #the ahccol
            #now we find in which availFSs this maxhcmb occurs for this parent;
            #the corresponding column is deleted from the FSmatlist matrices,
            #and for each FS in which this ahccol is a solution:
            #the index to the corresponding row in FSdata[[fs]]$parcombok
            #(and to the corresponding $Pseg) is deleted from
            #solstats[[cycle]]$pcos[[fs]]
            for (avfs in seq_along(FSgrpstruct$availFSs)) {
              cix <- which(FSgrpstruct$FSmatlist[[avfs]][maxprow,] == maxphcmb)
              if (length(cix) > 0) {
                # the maxprow parent occurs in this FS and the hapcomb to be deleted
                # is one of its possibilities in this FS
                fs <- FSgrpstruct$availFSs[avfs]
                solstats[[cycle]]$pcos[[fs]] <-
                  solstats[[cycle]]$pcos[[fs]][-cix]
              }
            } # for avfs
          } else if(length(FSgrpstruct$availFSs) == 1) {
            # exactly 1 FS in this group, just apply its best solution:
            fs <- FSgrpstruct$availFSs
            solstats[[cycle]]$pcos[[fs]] <- which.max(FSdata[[fs]]$Pseg)
          } else {
            # length(FSgrpstruct$availFSs) == 0, no FS's left in this group
            # nothing to do ?
          }
        } # while FSconflicts
        # If there are parents with multiple solutions left (these fit all
        # remaining FSs of which they are parents, so these FSs also have
        # multiple solutions):
        # For each combination of remaining parental solutions over all
        # remaining parents calculate an overall P-value (by multiplying the
        # P-values of all remaining FSs for that combination).
        # (TODO: perhaps use selcrit instead of Pval? Then selcrit must also
        #        be saved in FSresults)
        # If any single one of these parental haplotype combinations has the
        # largest P-value, select that and apply it to all these remaining FSs;
        # else treat all remaining FSs as unrelated.

        FSgrpstruct <- getFSgroupStructure(FSgrouping, grp, FSdata,
                                           solstats[[cycle]]$pcos, parents)
        allparcombs <- getAllParCombs(parents, FSdata, FSgrpstruct)
        combinedPval <- allparcombs$Pval
        allparcombs <- allparcombs$allparcombs

        # debug check:
        if ((length(FSgrpstruct$availFSs) == 0 && nrow(allparcombs) > 0) ||
            (length(FSgrpstruct$availFSs) > 0 && nrow(allparcombs) == 0))
          stop("mismatch availFSs and allparcombs")
        # All FSs in a group that have 0 rows in FSdata[[fs]]$parcombok are not
        # among the availFSs anymore. This may be because there was no solution
        # in the first place (conflict parent/offspring, missing data parents,
        # too many missing data offspring, skipped because of too many parcombs)
        # or because the solutions all conflicted with other FSs in the group.
        # These FSs will not be solved but their FSdata$messages must be
        # filled (if that hasn't been done yet).
        # If there are FSs left over (still available), allparcombs has 1 valid
        # combination per row with ahccols of all available parents
        #
        # first treat all FSs in grp that are not in availFSs:
        # (including the case that there are no availFSs left over)
        for (fs in FSgrouping$FSgroups[[grp]]) {
          if (!(fs %in% FSgrpstruct$availFSs)) {
            #check for debugging:
            #if (nrow(FSsof[[fs]]$parcombok) > 0)
            if (length(solstats[[cycle]]$pcos[[fs]]) > 0)
              stop("error: length(pcos) should be 0")
          }
        }
        if (nrow(allparcombs) == 0) {
          # no solutions remaining for this group
          solstats[[cycle]]$knownHap <- knownHap # the original, specified knownhap
          # we leave solstats[[cycle]]$Pahccols and $solvedFSindiv
          # unchanged; $pcos should now all be integer(0) for the group FSs
        } else {
          # at least one solution for this group, and FSgrpstruct$availFSs
          # not empty (but perhaps not including all FSs)
          # keep the combination(s) of parental ahccols that has/have the
          # highest combinedPval:
          maxP <- max(combinedPval)
          bestcomb <- which(combinedPval > 0.99 * maxP) #1% tolerance for rounding errors
          # Now we first do the things that use all best solutions:
          # - mark the haplotypes that occur in all solutions of each parent as
          #   knownHap for this cycle
          # - mark the hacs of the parent(s) with the same hac over all
          #   solutions as "known"
          # - mark all FSs in group as FSfit

          tmpapc <- allparcombs[bestcomb,, drop=FALSE]
          # limit the solstats[[cycle]]$pcos to those allowed in any
          # of the bestcomb group solutions:
          for (fs in FSgrpstruct$availFSs) {
            p1ahccols <- tmpapc[, colnames(tmpapc)==parents[fs, 1]]
            p2ahccols <- tmpapc[, colnames(tmpapc)==parents[fs, 2]]
            pcos <- integer(0)
            for (i in seq_along(p1ahccols)) {
              pco <- which(FSdata[[fs]]$parcombok[, 1] == p1ahccols[i] &
                           FSdata[[fs]]$parcombok[, 2] == p2ahccols[i])
              if (length(pco) == 1 && !(pco %in% pcos)) {
                pcos <- c(pcos, pco)
              }
            }
            names(pcos) <- NULL # all added pco have name P1col, confusing
            # debug check:
            if (!all(pcos %in% solstats[[cycle]]$pcos[[fs]]))
              stop("mismatch in pcos")
            solstats[[cycle]]$pcos[[fs]] <- pcos
          } # for fs

          parmrkdids <- allmrkdids[match(colnames(tmpapc), names(allmrkdids))]
          # for each parent we determine the haplotypes common to all selected
          # solutions and add these to the knownHap for this cycle:
          for (p in seq_len(ncol(tmpapc))) {
            phacs <- getAllHapcomb(mrkdid=parmrkdids[p], nmrk=ncol(ihl$allhap),
                                   ahcinfo=ahcinfo)[, tmpapc[, p], drop=FALSE]
            for (i in seq_len(ncol(phacs))) {
              if (i == 1) {
                phaps <- unique(phacs[, 1])
              } else phaps <- intersect(phaps, phacs[, i])
            }
            solstats[[cycle]]$knownHap <-
              union(solstats[[cycle]]$knownHap, phaps)
          }
          solstats[[cycle]]$knownHap <- sort(solstats[[cycle]]$knownHap)

          # for the parents with the same solution (same ahccol in all group
          # solutions), set that solution in solstats[[cycle]]$Pahscols
          uniqPar <- which(apply(tmpapc, MARGIN=2, max) ==
                           apply(tmpapc, MARGIN=2, min) )
          uniqPar <- colnames(tmpapc)[uniqPar] # names of the parents with a
          # unique solution
          for (fs in FSgrpstruct$availFSs) {
            for (p in 1:2) {
              up <- which(uniqPar == parents[fs, p])
              if (length(up) == 1) {
                solstats[[cycle]]$Pahccols[fs, p] <-
                  tmpapc[1, colnames(tmpapc) == uniqPar[up]]
              }
            }
            # add the number of genotyped FS indivs if this FS has
            # a unique solution:
            if (!anyNA(solstats[[cycle]]$Pahccols[fs,])) {
              solstats[[cycle]]$solvedFSindiv <-
                solstats[[cycle]]$solvedFSindiv + FSdata[[fs]]$mrksize
            }
          }
          # debug check: the number of rows for the FSs with one single
          # Pahccol for both parents should be 1:
          for (fs in FSgrpstruct$availFSs) {
            if (!anyNA(solstats[[cycle]]$Pahccols[fs,])) {
              if (length(solstats[[cycle]]$pcos[[fs]]) != 1) {
                stop("both parents one solution but multiple parcombok rows")
              }
            }
          }


          # Next we decide which one (or no) solution to apply
          selectAmongBest <- 2
          # selectAmongBest: if there are <= selectAmongBest equally good
          # solutions we select the first one. If > selectAmongBest we don't
          # select any (but we already marked these FSs as "fitted" so they
          # also won't be solved later as unrelated material, which would
          # be an even worse solution than a random one among the best).
          #TODO?: we could then still check which solution requires the smallest
          #number of haplotypes; or we could check for each solution how many of
          #the "rest" individuals can be assigned and select the maximum.
          if (length(bestcomb) == 1) {
            selbestcomb <- bestcomb
          } else if (length(bestcomb) <= selectAmongBest) {
            selbestcomb <- bestcomb[1] #so we take the first of a small number
            # (2) of equally good solutions
          } else {
            selbestcomb <- NA # too many equally good solutions, cannot select
          }
          # now length(selbestcomb) == 1

          if (is.na(selbestcomb)) {
            # we can't assign more parental genotypes than those that are the
            # same in all solutions:
            solstats[[cycle]]$selPahccols[FSgrpstruct$availFSs,] <-
              solstats[[cycle]]$Pahccols[FSgrpstruct$availFSs,]
          } else {
            # we selected one solution (possible the first among several) and
            # assign the selPahccol for that solution:
            # note that we keep all solstats[[cycle]]$pcos[fs], i.e. also
            # for the non-selected bestcomb solutions
            tmpapc <- allparcombs[selbestcomb,, drop=FALSE] # 1 row
            for (fs in FSgrpstruct$availFSs) {
              for (p in 1:2) {
                pcol <- which(colnames(tmpapc) == parents[fs, p])
                #if (!is.na(pcol) && (pcol %in% uniqParAhccol))
                solstats[[cycle]]$selPahccols[fs, p] <- tmpapc[1, pcol]
                #debug check: Pahccol should be all equal to selPahccol or NA
                if (!is.na(solstats[[cycle]]$Pahccol[fs, p]) &&
                    (is.na(solstats[[cycle]]$selPahccol[fs, p]) ||
                     solstats[[cycle]]$selPahccol[fs, p] !=
                     solstats[[cycle]]$Pahccol[fs, p]))
                  stop("Pahccol and selPahccol don't match")
              } # for p
            } # for fs
          }
          for (fs in FSgrpstruct$availFSs) {
            selpco <-
              which(FSdata[[fs]]$parcombok[, 1] ==
                      solstats[[cycle]]$selPahccol[fs, 1] &
                      FSdata[[fs]]$parcombok[, 2] ==
                      solstats[[cycle]]$selPahccol[fs, 2])
            # may be integer(0) if one or both selPahccol[fs,] are NA
            # debug checkS:
            if (!(length(selpco) %in% 0:1))
              stop("selpco incorrect length")
            if (xor(length(selpco) == 0,
                    anyNA(solstats[[cycle]]$selPahccol[fs,])))
              stop("selpco and selPahccol don't match NA")
            if (length(selpco) == 1) {
              # more debug checks:
              if (!all(FSdata[[fs]]$parcombok[selpco,] ==
                       solstats[[cycle]]$selPahccol[fs,]))
                stop("selpco and selPahccol don't match non-NA")
              if (!(selpco %in% solstats[[cycle]]$pcos[[fs]]))
                stop("selpco and pcos don't match")
              # all ok:
              solstats[[cycle]]$selpco[fs] <- selpco
            } else {
              # length(selpco)==0, no remaining solution for this fs
              solstats[[cycle]]$selpco[fs] <- NA
            }
          } # for fs
        } # nrow(allparcombs) > 0, at least one group solution
      } # for grp
    } # while(cycle)

    # hapOneBlock after cycles ####

    # now we know which cycle (bestcycle) has the optimum solution,
    # and with FSdata and solstats[[bestcycle]] we will now get all the results.
    # We assign the parental hacs from solstats[[bestcycle]]$selPahccols
    # we assign the FS individuals from the FS calculations based on these Phacs
    # we assign the rest (including the non-fitted FSs and their parents) based
    # on the known haplotypes in solstats[[bestcycle]]$knownhap
    hac <- matrix(integer(0), nrow=ahcinfo$ploidy, ncol=ncol(mrkDosage),
                  dimnames=list(NULL, colnames(mrkDosage)))
    FSfit <- logical(length(FS)) # all FALSE
    FSpval <- rep(NA_real_, length(FS)) # all NA
    rest <- setdiff(indiv, c(parents, unlist(FS)))
    imputedGeno <- matrix(integer(0), nrow=nrow(mrkDosage),
                          dimnames=list(rownames(mrkDosage), NULL))
    for (fs in seq_along(FS)) {
      pcorows <- solstats[[bestcycle]]$pcos[[fs]] # even if one
      # group solution (bestcomb) was chosen among several (2) equally good
      # ones, solstats[[cycle]]$pcos[[fs]] is not reduced to 1 but keeps all
      # rows allowed by bestcomb (the one selected solution, if any, is in
      # solstats[[cycle]]$selpco[fs] ). Similarly, solstats[[cycle]]$Pahccol
      # contains only non-NA ahccols for parents that have the same solution
      # across all bestcomb, while solstats[[cycle]]$selPahccol has also
      # ahccols for parents based on the selected selbestcomb group solution.
      if (length(pcorows) == 0) {
        # should never happen?
        if ((nrow(FSdata[[fs]]$parcombok) > 0)) {
          FSdata[[fs]]$message <- "parental conflict with other FS(s)"
        } else if (FSdata[[fs]]$message == "") {
          FSdata[[fs]]$message <- "no solution found"
        }
        # else: length(pcorows)==0, nrow(parcombok)==0 and message not '';
        # message already contains "no solution found" or something more
        # informative(?), we don't overwrite it.
        FSpval[fs] <- FSdata[[fs]]$maxPseg
        rest <- union(rest, c(parents[fs,], FS[[fs]])) # note that some parents
        # may still be filled in from other FSs
      } else {
        # length(pcorows) > 0 :
        # there are 1 or more solutions (pcorows) for this fs in bestcycle
        # It is possible that there is just one parental combination valid
        # (if length(pcorows)==1; in that case
        # solstats[[bestcycle]]$selPahccols[fs, ] are both known), or that
        # there are multiple parental combinations possible; in that case
        # we don't want to assign the parent(s) with multiple solutions but
        # still we want to assign the FS indivs with the same solution in all
        # possible parental combinations.
        # we assign the parent(s) with the same solution over all, if any:
        for (p in 1:2) {
          if (!is.na(solstats[[bestcycle]]$selPahccols[fs, p])) {
            # one single selected hac for this parent, assign:
            Pdid <- allmrkdids[names(allmrkdids)==parents[fs, p]]
            hac[, colnames(hac)==parents[fs, p]] <-
              getAllHapcomb(mrkdid=Pdid, nmrk=ncol(ihl$allhap),
                            ahcinfo=ahcinfo)[, solstats[[bestcycle]]$selPahccols[fs, p]]
          }
        }
        if (!is.na(solstats[[bestcycle]]$selpco[fs])) {
          # one solution for this FS selected (but there might be more
          # than one solstats[[bestcycle]]$pcos[[fs]])
          FSfit[[fs]] <- TRUE
          FSpval[fs] <- FSdata[[fs]]$Pseg[solstats[[bestcycle]]$selpco[fs]]
        } else {
          # no selected pco for this FS but >=1 FS solutions (pcorows) found
          FSdata[[fs]]$message <- "multiple FS solutions"
          FSpval[fs] <- max(FSdata[[fs]]$Pseg[pcorows])
        }
        # Next we assign all FS indivs with a mrkdid that has a unique solution
        # over all (1 or more) possible parental combinations:
        if (!is.na(solstats[[bestcycle]]$selpco[fs]))
          pcorows <- solstats[[bestcycle]]$selpco[fs]
        unexp_dids <- integer(0)
        multisol_dids <- integer(0)
        for (pcor in seq_along(pcorows)) {
          pcr <- pcorows[pcor]
          sof <- solveOneFS(mrkdids=allmrkdids, ihl=ihl, FS=FS[[fs]],
                            P=parents[fs,],
                            knownPhac=FSdata[[fs]]$parcombok[pcr,], #known Phac both parents
                            triedPhacs=list(integer(0), integer(0)), #to have sof recomputed
                            maxparcombs=maxparcombs,
                            ahcinfo=ahcinfo, minfrac=minfrac, errfrac=errfrac,
                            DRrate=DRrate, minPseg=minPseg)
          # This sof has only one row in its parcombok because only one, valid,
          # combination of parents was tested! So sof$parcombok is different
          # from FSdata[[fs]]$parcombok, which may have more than 1 row
          if (pcor == 1) {
            expdid <- unique(sof$didNhapcomb[[1]]$mrkdid) #original may contain (consecutive) duplicates
            mdhacs <- selOneMrkdidHac(sof, pco=1, expdid=expdid) # mdhacs is matrix of hacs, not vector of ahccols
            if (length(pcorows) == 1) {
              #TODO remove debug check
              if (FSdata[[fs]]$message != "") stop("FS message not empty")
              # this FS has a unique solution, we can test for unexpected mrkdids
              fsmrkdids <- allmrkdids[names(allmrkdids) %in% FS[[fs]]]
              fsdidstable <- table(fsmrkdids)
              unexp_dids <- setdiff(names(fsdidstable), expdid) #colnames(somh$didproblems)[somh$didproblems[1,]]
              sum_unexp_dids <- sum(fsmrkdids %in% unexp_dids) #sof$didNhapcomb[[1]]$mrkdid %in% unexp_dids)
              if (sum_unexp_dids > 0) {
                FSdata[[fs]]$message <-
                  paste0("unexpected mrkdid(s) ",
                         paste(unexp_dids, collapse=", "),
                         ", n=", sum_unexp_dids)
              }
              multisol_dids <-
                intersect(colnames(mdhacs)[is.na(mdhacs[1,])], names(fsdidstable))
              sum_multisol_dids <- sum(fsmrkdids %in% multisol_dids)
              if (sum_multisol_dids > 0) {
                if (FSdata[[fs]]$message != "")
                  FSdata[[fs]]$message <- paste0(FSdata[[fs]]$message, "; ")
                FSdata[[fs]]$message <-
                  paste0(FSdata[[fs]]$message,
                         "multiple solutions for mrkdid(s) ",
                         paste(multisol_dids, collapse=", "),
                         ", n=", sum_multisol_dids)
              }
              # length(pcorows)==1, so a unique parental combination
              # we can check if there are partially genotyped FS individuals
              # that can have only one expected mrkdid; these can be imputed,
              # and if the mrkdid has only one allowed hac it can also be
              # haplotyped
              NAind <- which(is.na(colSums(mrkDosage)) &
                               colnames(mrkDosage) %in% FS[[fs]])
              if (length(NAind) > 0 && length(NAind) < 0.5 * length(FS[[fs]])) {
                expmrkdos <- mrkdid2mrkdos(dosageIDs=expdid, nmrk=nrow(mrkDosage),
                                           ploidy=ahcinfo$ploidy)
                impdids <-
                  matrix(c(NAind, rep(NA, length(NAind))), nrow=2, byrow=TRUE,
                         dimnames=list(c("NAind", "mrkdid"), NULL))
                for (i in seq_along(NAind)) {
                  nai <- NAind[i]
                  if (!all(is.na(mrkDosage[,nai]))) {
                    md <- 1; matches <- integer(0)
                    while (md <= ncol(expmrkdos) && length(matches) < 2) {
                      if (all(is.na(mrkDosage[,nai]) |
                              mrkDosage[,nai] == expmrkdos[, md]))
                        matches <- c(matches, md)
                      md <- md + 1
                    }
                    if (length(matches) == 1)
                      impdids[2, i] <- as.integer(colnames(expmrkdos)[matches])
                  }
                }
                # impdids (imputed dids) calculated; check if the
                # segregation still fits parents:
                impdids <- impdids[, !is.na(impdids[2,]), drop=FALSE]
                if (ncol(impdids) > 0) {
                  parhac <- hac[, colnames(hac) %in% parents[fs, ], drop=FALSE]
                  # perhaps in reverse order; this doesn't matter
                  # debug check:
                  if (ncol(parhac) != length(unique(parents[fs,])) ||
                      anyNA(parhac)) stop("error in parhac")
                  totFSmrkdids <- c(allmrkdids[names(allmrkdids) %in% FS[[fs]]],
                                    impdids[2,]) # includes NAs
                  newFSdidstable <- table(totFSmrkdids, useNA="no")
                  impstats <-
                    testFSseg(parhac=parhac, DRrate=DRrate, errfrac=errfrac,
                               FSmrkdidsTable=newFSdidstable, allhap=ihl$allhap)
                  #TODO: remove debug check:
                  # the imputed mrkdids should be among the expected mrkdids, check:
                  if (anyNA(impdids[2,]) || !all(impdids[2,] %in% expdid))
                    stop("error in impdids")
                  if (impstats$stats$P > 0.1 * FSpval[fs]) {
                    # with the extra imputed FS mrkdids the fit of the FS family
                    # is not too much worse, so add the imputed mrkdids
                    allmrkdids[impdids[1,]] <- impdids[2,]
                    # and we include the imputed mrk genotypes as a component
                    # of the return value:
                    FSimpGeno <-
                      mrkdid2mrkdos(dosageIDs=impdids[2,], nmrk=nrow(mrkDosage),
                                    ploidy=ahcinfo$ploidy)
                    rownames(FSimpGeno) <- rownames(mrkDosage)
                    colnames(FSimpGeno) <- colnames(mrkDosage)[impdids[1,]]
                    imputedGeno <- cbind(imputedGeno, FSimpGeno)
                  }
                }
              }
            } # if (length(pcorows) == 1)
            mdhacs <- mdhacs[, !is.na(mdhacs[1,]), drop=FALSE]
            # end of pcor == 1
          } else if (ncol(mdhacs) > 0) {
            # already looked at pcorows[1 ... pcor-1]; now see if there are common
            # mdhacs between those and the current pcorow.
            # each solution has a different set of expected mrkdid;
            # delete the ones that are not common to all.
            # get the hacs from the new pcorow that also occurred in the earlier pcorows:
            tmphacs <- selOneMrkdidHac(sof, pco=1, #because this sof has only one row in parcombok
                                       expdid=colnames(mdhacs)) #$hacmatrix #the earlier mdhacs
            tmphacs <- tmphacs[, !is.na(tmphacs[1,])]
            # keep only the mdhacs from the earlier rows that also occur
            # in the current row:
            mdhacs <- mdhacs[, colnames(mdhacs) %in% colnames(tmphacs),
                             drop=FALSE]
            # keep only the mdhacs for mrkdids that HAVE THE SAME SELECTED HACs
            # in all solutions:
            id <- logical(ncol(mdhacs))
            for (md in seq_along(id)) {
              did <- colnames(mdhacs)[md]
              id[md] <- identical(mdhacs[, md], tmphacs[, colnames(tmphacs)==did])
            }
            mdhacs <- mdhacs[, id, drop=FALSE]
          }
        } # for pcor
        # mdhacs now has one ahccolnr for all FS mrkdids that
        # have the same selected hac over all pcorows.
        # we now assign all FS individuals that have one of these mrkdids:
        for (mdhcol in seq_len(ncol(mdhacs))) {
          did <- as.numeric(colnames(mdhacs)[mdhcol])
          ind <- intersect(FS[[fs]],
                           names(allmrkdids)[!is.na(allmrkdids) & allmrkdids==did])
          hac[, colnames(hac) %in% ind] <-
            mdhacs[, mdhcol]
          #TODO remove DEBUG CHECK:
          # if (length(ind) > 0) {
          #   nwdid <-
          #     mrkdos2mrkdid(
          #       mrkDosage=hapcomb2mrkdos(hapcomb=hac[, colnames(hac) %in% ind],
          #                                allhap=ihl$allhap),
          #       ploidy=ahcinfo$ploidy,
          #       check=FALSE)
          # }
        } # for mdhcol
      } #length(pcorows) > 0
    } # for fs

    # finally add the rest = the original rest plus the FSs for which no
    # solution or parental conflicts (but not the FSs that were left unsolved
    # because of too many equally good solutions)
    rest_ind <- intersect(rest, names(allmrkdids)[!is.na(allmrkdids)])
    rest_ihl <- inferHaps_noFS(mrkdids=allmrkdids[names(allmrkdids) %in% rest_ind],
                                ahcinfo=ahcinfo,
                                markers=rownames(mrkDosage),
                                minfrac=minfrac,
                                knownHap=solstats[[bestcycle]]$knownHap,
                                useKnownDosage=TRUE)
    for (md in seq_along(rest_ihl$mrkdidsTable)) {
      did <- names(rest_ihl$mrkdidsTable)[md]
      if (length(rest_ihl$hclist[[md]]) == 1) {
        ind <- names(allmrkdids)[(names(allmrkdids) %in% rest_ind) &
                                  allmrkdids==did]
        if (length(ind) > 0) {
          hac[, colnames(hac) %in% ind] <-
          getAllHapcomb(mrkdid=did, nmrk=ncol(ihl$allhap),
                        ahcinfo=ahcinfo)[, rest_ihl$hclist[[md]]]
          #TODO remove DEBUG CHECK:
          # nwdid <-
          #   mrkdos2mrkdid(
          #     mrkDosage=hapcomb2mrkdos(hapcomb=hac[, colnames(hac) %in% ind],
          #                              allhap=ihl$allhap),
          #     ploidy=ahcinfo$ploidy,
          #     check=FALSE)
        }
      }
    } # for md
    #for all indiv with a unique combination of haplotypes we now have
    #stored the solution (haplotype combination) in hac
    #
    #Finally, for parents that are among the rest (i.e. parents of FSs for which
    #no solution FS was found, and that were therefore solved as unrelated
    #material) we check if their solution conflicts with too many
    #of their offspring - the idea is that the parent may be mis-genotyped
    #which results in a failed FS solution.
    # There are 2 situations:
    # - both parents are in rest: then both parents must be checked separately
    #   against the FS progeny, and if both pass then the combination must be
    #   checked
    # - one of the parents is also parent of a non-rejected FS. If that parent
    #   has a final hac, the "rest" parent must be checked in combination with
    #   this; if the other-FS parent did not get a hac in the end, the "rest"
    #   parent is checked separately
    for (fs in seq_along(FS)) {
      restpar <- which(parents[fs,] %in% rest)
      FShac <- hac[, !is.na(hac[1,]) & colnames(hac) %in% FS[[fs]], drop=FALSE]
      if (length(restpar) > 0 && !all(is.na(hac[1, parents[fs, restpar]])) &&
          ncol(FShac[]) > 4) {
        #only in this case (some rest-parents with assigned hacs, and there
        #are "sufficient" (>4) haplotyped FS progeny to test against) we need
        #to do something.
        # first get a list of parents; gametes for each parent, with DR:
        pargamlist <- vector(mode="list", 2)
        for (p in 1:2) {
          if (is.na(hac[1, parents[fs, p]])) {
            pargamlist[[p]] <- "no data"
          } else {
            pargamlist[[p]] <- list()
            pargamlist[[p]]$DR <-
              getGameteFreqs(parhac=hac[, parents[fs, p]], DRrate=0.04)
            #  we only need the set of possible gametes; their freq (and so the exact
            #  DRrate) is not used. Only with DRrate=0 we would exclude all DR gametes
          }
        }
        # next get the FShac ordered so that identical columns are adjacent:
        o <- do.call("order", lapply(1:nrow(FShac), function(i) FShac[i,]))
        FShac <- FShac[, o, drop=FALSE]
        # is one of the parents also parent of a non-rest FS?
        if (length(restpar) == 1) {
          #one of the parents in rest (and we know it has a non-missing hac)
          if (parents[fs, 1] %in% rest) p <- 1 else p <- 2
          # p is the "rest" parent, 3-p the non-rest; we test only p separately
          if (is.na(hac[1, parents[3-p]])) {
            # the non-rest parent has missing hac
            chk <- checkFS_OneParent(pargamlist[[p]], FShac=FShac, nhap=nhap)
            if (nrow(chk) != 1) stop("error")
            nok <- rowSums(chk)
            if (nok/ncol(chk) < 0.9 && ncol(chk)-nok > 4) {
              # this parent has too many conflicts with its progeny: more than 10%
              # and more than 4 of the haplotyped progeny conflict.
              hac[, parents[fs, p]] <- NA
            }
          } else {
            # the non-rest parent has hac, we test p in combination with the other
            chk <- checkFS_TwoParents(pargamlist, FShac)
            if (nrow(chk) != 1) stop("error")
            nok <- rowSums(chk)
            if (nok/ncol(chk) < 0.9 && ncol(chk)-nok > 4) {
              hac[, parents[fs, p]] <- NA
              # parent 3-p might also not fit this FS but it is supported by
              # other FSs
            }
          }
        } else {
          # both parents in rest; first test separately:
          for (p in 1:2) if (!is.na(hac[1, parents[fs, p]])) {
            chk <- checkFS_OneParent(pargamlist[[p]], FShac=FShac, nhap=nhap)
            if (nrow(chk) != 1) stop("error")
            nok <- rowSums(chk)
            if (nok/ncol(chk) < 0.9 && ncol(chk)-nok > 4) {
              # this parent has too many conflicts with its progeny: more than 10%
              # and more than 4 of the haplotyped progeny conflict.
              hac[, parents[fs, p]] <- NA
            }
          }
          # next we test the parental combination if both still not NA
          if (!anyNA(hac[1, parents[fs,]])) {
            chk <- checkFS_TwoParents(pargamlist, FShac)
            if (nrow(chk) != 1) stop("error")
            nok <- rowSums(chk)
            if (nok/ncol(chk) < 0.9 && ncol(chk)-nok > 4) {
              hac[, parents[fs,]] <- NA
              # both parent hacs set to NA
            }
          }
        }
      }
    } # for fs

    hapdos <- hapcomb2hapdos(hapcomb=hac, nhap=nrow(ihl$allhap))
    rownames(hapdos) <-
      getHaplotypeNames(haploblock=hbname, hapcount=nrow(ihl$allhap))
    result <-
      list(hapdos=hapdos, mrkdids=allmrkdids, markers=colnames(ihl$allhap),
           FSfit=FSfit,
           FSmessages=sapply(FSdata, `[[`, "message"),
           FSpval=FSpval, imputedGeno=imputedGeno)
    # TODO?:(if rest_ihl with two minfrac values results in additional haplotypes,
    # see if these can help solve the remaining FSs? And then re-iterate?
  } # one or more FSs
  if (dropUnused) {
    result$hapdos <-
      result$hapdos[rowSums(result$hapdos, na.rm=TRUE) > 0, , drop=FALSE]
    #has 0 rows if no indiv has a unique solution
  }
  result
} # hapOneBlock

getLinkedFSgroups <- function(parents) {
  #parents: a matrix with 2 columns and one row for each FS family,
  #         with the names of the mother and father of each family.
  #Result: a list with 4 items:
  # $FSgroups: a list with one item for each group of linked FSs;
  #            each of the items is an integer vector with the row numbers
  #            of the FSs in the group
  # $FSgroupParents: a list with the same number of items as FSgroups,
  #            each item is a character vector with the names of all
  #            parents in that group (without duplicates)
  # $FSgroupnrs: an integer vector with for each FS the group it belongs to
  # $ParentGroupnrs: a named integer vector with for each parent the group
  #                  it belongs to
  # test with:
  # parents <- matrix(c("A","B","C","D","E","A","G","H","H","E"), ncol=2, byrow=TRUE)
  #
  #initially each FS gets its own item in both lists:
  FSgroups <- as.list(seq_len(nrow(parents)))
  FSgroupParents <- lapply(as.list(data.frame(t(parents), stringsAsFactors=FALSE)), sort)
  #merge groups based on shared parents:
  changes <- TRUE
  while (changes) {
    changes <- FALSE
    g <- 2
    while (!changes && g <= length(FSgroupParents)) {
      h <- 1
      while (!changes && h<g) {
        if (length(intersect(FSgroupParents[[g]], FSgroupParents[[h]])) > 0) {
          changes <- TRUE
          FSgroupParents[[h]] <- sort(unique(c(FSgroupParents[[g]], FSgroupParents[[h]])))
          FSgroups[[h]] <- sort(c(FSgroups[[g]], FSgroups[[h]]))
          FSgroupParents <- FSgroupParents[-g]
          FSgroups <- FSgroups[-g]
        }
        h <- h+1
      }
      g <- g+1
    }
  }
  FSgroupnrs <- integer(nrow(parents))
  uniparents <- sort(unique(as.vector(parents)))
  Parentgroupnrs <- integer(length(uniparents))
  names(Parentgroupnrs) <- uniparents
  for (g in seq_along(FSgroups)) {
    FSgroupnrs[FSgroups[[g]]] <- g
    Parentgroupnrs[names(Parentgroupnrs) %in% FSgroupParents[[g]]] <- g
  }
  list(FSgroups=FSgroups, FSgroupParents=FSgroupParents,
       FSgroupnrs=FSgroupnrs, Parentgroupnrs=Parentgroupnrs)
} #getLinkedFSgroups

getNewMrkdids <- function(mrkDosages, haploblock, maxmrk, ahcinfo) {
  # non-exported function only for use by function inferHaplotypes.
  # Parameters as the corresponding ones of inferHaplotypes but already
  # checked, and mrkDosages as matrix not data.frame with only the selected
  # indiv.
  # output: a list with one element for each haploblock size (nr of markers),
  # from 1 to the maximum haploblock size in haploblock (or to maxmrk if that
  # is smaller);
  # each element has all mrkdid (as integers) that occur over all haploblocks
  # of that size and for which no allhapcomb is yet available in
  # ahclist or ahccompletelist
  hbsizes <- sapply(haploblock, length)
  maxmrk <- min(maxmrk, max(hbsizes))
  result <- vector(mode="list", maxmrk)
  if (maxmrk <= ahcinfo$complete_nmrk) return(result)
  # we need to check all mrkdids for new ones
  knownmrkdid <- result # list of same length
  for (hbsize in (ahcinfo$complete_nmrk+1):maxmrk) {
    if (length(ahcinfo$ahclist) < hbsize) {
      knownmrkdid[[hbsize]] <- integer(0)
    } else knownmrkdid[[hbsize]] <- as.integer(names(ahcinfo$ahclist[[hbsize]]))
  }
  for (hb in seq_along(haploblock)) {
    if (hbsizes[[hb]] %in% (ahcinfo$complete_nmrk+1):maxmrk) {
      # mrkdids <-
      #   unique(mrkdos2mrkdid(mrkDosages[rownames(mrkDosages) %in%
      #                                     haploblock[[hb]],,drop=FALSE],
      #                        ploidy=ploidy, check=FALSE))
      uniqdos <-
        unique(mrkDosages[match(haploblock[[hb]], rownames(mrkDosages)),,
                          drop=FALSE], MARGIN=2)
      mrkdids <- mrkdos2mrkdid(uniqdos, ploidy=ahcinfo$ploidy, check=FALSE)
      result[[hbsizes[hb]]] <- c(result[[hbsizes[hb]]],
                                 setdiff(mrkdids, knownmrkdid[[hbsizes[hb]]]))
    }
  }
  for (hbsize in (ahcinfo$complete_nmrk+1):maxmrk) {
    if (!is.null(result[[hbsize]]))
        result[[hbsize]] <- setdiff(unique(result[[hbsize]]), NA)
  }
  result
} #getNewMrkdids

addNewMrkdids <- function(newMrkdids, ahcinfo, progress) {
  # non-exported function only for use by function inferHaplotypes.
  # newMrkdids: list as produced by getNewMrkdids
  # Parameters as the corresponding ones of inferHaplotypes but already
  # checked, and mrkDosages as matrix not data.frame with only the selected
  # indiv.
  # output: a list with one element for each haploblock size (nr of markers),
  # from 1 to the maximum haploblock size in haploblock (or to maxmrk if that
  # is smaller);
  # each element has all mrkdid (as integers) that occur over all haploblocks
  # of that size and for which no allhapcomb is yet available in
  # ahclist or ahccompletelist
  newdidcounts <- sapply(newMrkdids, length)
  totnew <- sum(newdidcounts)
  if (any(newdidcounts > 0)) {
    maxnmrk <- max(which(newdidcounts > 0))
    # add the haplotype combinations for the newdids to ahcinfo$ahclist
    # and save the new ahclist in ahcdir (overwriting any earlier version)
    if (progress && (maxnmrk + ahcinfo$ploidy >= 12 ||
                     sum(newdidcounts) > 100)) {
      cat(paste(sum(newdidcounts),
                "sets of haplotype combinations need to be calculated;\n"))
      cat("this may take quite a long time\n")
    }
    cum <- 0
    for (m in which(newdidcounts > 0)) {
      if (length(ahcinfo$ahclist) < m || !is.list(ahcinfo$ahclist[[m]])) {
        ahcinfo$ahclist[[m]] <- list()
        # initialize also all lower nmrk:
        for (i in seq_len(m)) {
          if (!is.list(ahcinfo$ahclist[[i]]))
            ahcinfo$ahclist[[i]] <- list()
        }
      }
      allhap <- allHaplotypes(1:m)
      for (i in seq_along(newMrkdids[[m]])) {
        did <- newMrkdids[[m]][i]
        if (progress) {
          cat(paste0("\rcalculation of new mrkdids: ", cum + i, " / ", totnew))
          flush.console()
        }
        ahcinfo$ahclist[[m]][[as.character(did)]] <-
          calcAllhapcomb4mrkdid(mrkdid=did, ploidy=ahcinfo$ploidy, allhap=allhap)
        # save after every 500 mrkdids so we can resume after crash:
        if (i %% 500 == 0) {
          save(ahcinfo$ahclist,
               file=paste0(ahcinfo$ahcdir, "ahclist_", ahcinfo$ploidy, "x.RData"),
               version=2)
        }
      }
      cum <- cum + newdidcounts[m]
    }
    cat("\n"); flush.console()
    ahclist <- ahcinfo$ahclist # because save will only save an entire R object, not a list element
    save(ahclist,
         file=paste0(ahcinfo$ahcdir, "ahclist_", ahcinfo$ploidy, "x.RData"),
         version=2)
  }
  ahcinfo
}

#'@title infer haplotypes for one or more haploblocks
#'@description infer haplotypes for one or more haploblocks, for all individuals,
#'using FS family(s) (with parents) if present, and infer haplotypes for
#'non-FS material as well
#'@usage inferHaplotypes(mrkDosage, indiv=NULL, ploidy, haploblock,
#'parents=NULL, FS=NULL, minfrac=c(0.1, 0.01), errfrac=0.025, DRrate=0.025,
#'maxmrk=0, dropUnused=TRUE, maxparcombs=150000, minPseg=1e-8,
#'knownHap=integer(0), progress=TRUE, printtimes=FALSE, ahcdir)
#'@param mrkDosage matrix or data.frame. Markers are in rows, individuals in
#'columns, each cell has a marker dosage. Names of individuals are the column
#'names, marker names are the row names or (if a data.frame) in a column named
#'MarkerNames. All marker dosages must be in 0:ploidy or NA.
#'@param indiv NULL (default) or a character vector with names of all individuals
#'to be considered. If NULL, all columns of mrkDosage are selected.\cr
#'All indivs that are not in parents or FS vectors (see below) are considered
#'unrelated, i.e. we have no implementation for pedigrees (yet).
#'@param ploidy all marker dosages should be in 0:ploidy or NA
#'@param haploblock a list of character vectors. The names are the names of the
#'haploblocks, the character vectors have the names of the markers in each
#'haploblock. Haplotype names are constructed from the haploblock names, which
#'are used as prefixes to which the (zero-padded) haplotype numbers are are
#'appended with separator '_'.
#'@param parents a matrix with one row for each FS family and two columns
#'for the two parents, containing the names of the female and male parent
#' of each family.
#'@param FS a list of character vectors. Each character vector has the names of
#'the individuals of one FS family. The items of the list should correspond
#'to the rows of the parents matrix, in the same order.
#'@param minfrac vector of two fractions, default 0.1 and 0.01. A haplotype is
#'considered to be certainly present if it must occur in at least a fraction
#'minfrac[1] of all individuals; in the final stage for the "other"
#'individuals (those that do not belong to the FS or its parents) this fraction
#'is lowered to minfrac[2]; see also inferHaps_noFS
#'@param errfrac the assumed fraction marker genotypes with an error (over all
#'markers in the haploblock). The errors are assumed to be uniformly distributed
#'over all except the original marker dosage combinations (mrkdids)
#'@param DRrate default 0.025. The rate of double reduction per meiosis (NOT
#'per allele!); e.g. with a DRrate of 0.04, a tetraploid parent with
#'genotype ABCD will produce a fraction of 0.04 of DR gametes AA, BB, CC and DD
#'(each with a frequency of 0.01), and a fraction of 0.96 of the non-DR gametes
#'AB, AC, AD, BC, BD, CD (each with a frequency of 0.16)
#'@param maxmrk Haploblocks with more than maxmrk markers will be skipped.
#'Default 0: no haploblocks are skipped
#'@param dropUnused TRUE (default) if the returned matrix should only contain
#'rows for haplotypes that are present; if FALSE matrix contains rows for all
#'possible haplotypes
#'@param maxparcombs Parent 1 and 2 both may have multiple possible haplotype
#'combinations. For each pair of haplotype combinations (one from P1 and one
#'from P2) the expected FS segregation must be checked against the observed.
#'This may take a long time if many such combinations need to be checked.
#'This parameter sets a limit to the number of allowed combinations per
#'haploblock; default 150000 takes about 45 min.
#'@param minPseg The minimum P-value of a chisquared test for segregation in
#'FS families. The best solution for an FS family is selected based on a
#'combination of P-value and number of required haplotypes, among all
#'candidate solutions with a P-value of at least minPseg. If no such solution
#'is found the FS and its parents are treated as unrelated material
#'@param knownHap integer vector with haplotype numbers (haplotypes that must be
#'present according to prior inference or knowledge, numbers refer to rows of
#'matrix produced by allHaplotypes); default integer(0), i.e. no known
#'haplotypes
#'@param progress if TRUE, and new haplotype combinations need to be calculated,
#'and the number of markers and the ploidy are both >= 6, progress is indicated
#'by printed messages
#'@param printtimes if TRUE, the time needed to process each haploblock is
#'printed
#'@param ahcdir a single directory, or not specified.
#'inferHaplotypes uses lists that for each combination of marker
#'dosages give all possible combinations of haplotype dosages. These lists
#'(ahclist and ahccompletelist) are loaded and saved at the directory
#'specified by ahcdir. If no ahcdir is specified it is set to the current
#'working directory.\cr
#'If an ahclist or ahccompletelist for the correct ploidy is already in
#'GlobalEnv this is used and no new list is loaded, even if ahcdir is specified.
#'@details
#'First we consider the case where one or more FS families and their
#'parents are present in the set of samples. In that case, initially the
#'possible haplotype configurations of the parents are determined.
#'From that, all their possible gametes (assuming polysomic
#'inheritance) are calculated and all possible FS haplotype configurations.
#'Comparing this with the observed FS marker dosages the most likely parental
#'and FS configurations are found.\cr
#'It is possible that multiple parental combinations can explain the observed
#'marker dosages in the FS. In that case, if one is clearly more likely and/or
#'needs less haplotypes, that one is chosen. If there is no clear best solution
#'still the parents and FS individuals that have the same haplotype
#'configuration over all likely solutions are assigned that configuration.\cr
#'For FS where no good solution is found (because of an error in the marker
#'dosages of a parent, or because the correct solution was not considered) the
#'parents and individuals will be considered as unrelated material.\cr
#'If several FS families share common parents they are treated as a group,
#'and only solutions are considered that are acceptable for all families
#'in the group.\cr
#'Finally (or if no FS families are present, immediately) the other samples
#'are haplotyped, which are considered as unrelated material. If FS families
#'have been solved the haplotypes in their parents are considered "known",
#'and known haplotypes can also be supplied (parameter knownHap). For these
#'samples we consecutively add haplotypes that must be present in a minimum
#'number of individuals, always trying to minimize the number of needed
#'haplotypes.\cr
#'InferHaplotypes uses tables that, for each combination of dosages of the
#'markers in the haploblock, list all haplotype combinations (ahc) that
#'result in these marker dosages. In principle inferHaplotypes uses a list
#'(ahccompletelist) that, for a given ploidy, has all the haplotype combinations
#'for haploblocks from 1 up to some maximum number of markers. This list can be
#'computed with function build_ahccompletelist. If this list is not available
#'(or is some haploblocks contain more markers than the list), the ahc for
#'the (extra) marker.\cr
#'See the PolyHaplotyper vignette for an illustrated explanation.
#'@return a list with for each haploblock one item that itself is a list
#'with items:\cr
#'message; if this is "" the haploblock is processed and further
#'         elements are present; else this message says why the haploblock was
#'         skipped (currently only if it contains too many markers)\cr
#'hapdos: a matrix with the dosages of each haplotype (in rows) for each
#'        individual (in columns). For each individual the haplotype dosages
#'        sum to the ploidy. If dropUnused is TRUE Only the haplotypes that
#'        occur in the population are shown, else all haplotypes\cr
#'mrkdids: a vector of the mrkdid (marker dosage ID) for each individual
#'         (each combination of marker dosages has its own ID; if any of the
#'         markers has an NA dosage the corresponding mrkdid is also NA).\cr
#'         The mrkdids can be converted to the marker dosages with function
#'         mrkdid2mrkdos.\cr
#'markers: a vector with the names of the markers in the haploblock\cr
#'imputedGeno: a matrix in the same format as param mrkDosage, with one row
#'        for each marker in the haploblock and one column per imputed
#'        individual, with the dosages of the markers. These are the individuals
#'        that have incomplete data in mrkDosage but where the available marker
#'        dosages match only one of the expected marker genotypes in the FS
#'        family (only individuals in FS families are imputed). It is possible
#'        that an individual with imputed marker dosages is not haplotyped (as
#'        is the case for individuals with complete marker data) if the
#'        marker dosages match different possible haplotype combinations.
#'The next elements are only present if one or more FS families were
#'specified:\cr
#'FSfit: a logical vector with one element per FS family; TRUE if a (or
#'        more than one) acceptable solution for the FS is found (although
#'        if multiple solution are found they might not be used if unclear
#'        which one is the best solution). (Even if no
#'        solution was found for an FS, still its individuals may have a
#'        haplotype combination assigned ignoring their pedigree)\cr
#'FSmessages: a character vector with one item per FS family: any
#'        message relating to the fitting of a model for that FS,
#'        not necessarily an error\cr
#'FSpval: a vector of the chi-squared P-value associated with the selected
#'        FS model for each FS family, or the maximum P value over all
#'        models in case none was selected\cr
#'If for new combinations of marker dosages the possible haplotype combinations
#'have to be calculated, an ahclist file is written to ahcdir
#'@examples
#'\donttest{
#'# this example takes about 1 minute to run:
#'data(PolyHaplotyper_small)
#'results <- inferHaplotypes(mrkDosage=phdos, ploidy=6,
#'haploblock=phblocks, parents=phpar, FS=phFS)
#'names(results)
#'names(results[[1]])
#'}
#'@export
inferHaplotypes <- function(mrkDosage, indiv=NULL, ploidy,
                            haploblock, parents=NULL, FS=NULL,
                            minfrac=c(0.1, 0.01), errfrac=0.025, DRrate=0.025,
                            maxmrk=0,
                            dropUnused=TRUE, maxparcombs=150000,
                            minPseg=1e-8,
                            knownHap=integer(0), progress=TRUE,
                            printtimes=FALSE,
                            ahcdir) {

  mrkDosage <- checkmrkDosage(mrkDosage, ploidy=ploidy, indiv=indiv)
  tmp <- checkpops(parents, FS, mrkDosage)
  if (tmp$message != "") stop(tmp$message)
  parents <- tmp$parents; FS <- tmp$FS; rm(tmp) #a.o. all as character
  haploblock <- checkHaploblock(haploblock, mrkDosage)
  nmrk <- max(sapply(haploblock, length))
  if (maxmrk == 0) maxmrk <- nmrk else maxmrk <- min(maxmrk, nmrk)
  if (missing(ahcdir)) ahcdir <- NULL
  ahcinfo <- loadAllHapcombLists(ploidy=ploidy, nmrk=maxmrk, ahcdir)
  # see if there are new mrkdids for which all haplotype combinations
  # must be calculated:
  newdids <- getNewMrkdids(mrkDosage, haploblock, maxmrk, ahcinfo)
  ahcinfo <- addNewMrkdids(newdids, ahcinfo, progress)
  nmrk <- sapply(haploblock, length)
  result <- vector(mode="list", length(haploblock))
  for (hb in seq_along(haploblock)) {
    starttime <- proc.time()
    cat(paste0("haploblock ", hb, " of ", length(haploblock), ": ",
               names(haploblock)[hb])); flush.console() # no \n !
    result[[hb]] <- list()
    if (nmrk[hb] > maxmrk) {
      result[[hb]]$message <- paste("skipped:", nmrk[hb],
                                    "markers while the limit is", maxmrk)
    } else {
      result[[hb]] <-
        hapOneBlock(mrkDosage[match(haploblock[[hb]], rownames(mrkDosage)),,
                              drop=FALSE],
                    hbname=names(haploblock)[hb],
                    ahcinfo=ahcinfo,
                    parents=parents, FS=FS, minfrac=minfrac,
                    errfrac=errfrac, DRrate=DRrate,
                    dropUnused=dropUnused,
                    maxparcombs=maxparcombs,
                    minPseg=minPseg,
                    knownHap=knownHap, progress=progress)
    }
    if (is.null(result[[hb]]$message)) result[[hb]]$message <- ""
    if (printtimes) {
      cat(paste(" -", format((proc.time()-starttime)[3], nsmall=2), "sec\n"))
    } else cat("\n")
    flush.console()
  } # for hb
  names(result) <- names(haploblock)
  result
} #inferHaplotypes

#'@title show marker and haplotype dosages for one FS family
#'@description show marker and haplotype dosages for one FS family and its
#'parents
#'@usage showOneFS(FSnr, hbresults, mrkDosage, FS, parents)
#'@param FSnr the number of the FS family (indexes the FS list and parents)
#'@param hbresults a list with the haplotyping results for one haploblock:
#'one element of a list as returned by inferHaplotypes
#'@param mrkDosage a matrix of marker dosages, may contain rows for more
#'markers than only those in the current haploblock (the relevant markers are
#'specified in the hbresults list)
#'@param FS a list of which each item is a (character) vector with the
#'names of the individuals in that FS family
#'@param parents a (character) matrix with 2 columns and one row for each
#'FS family in FS, with the names of the two parents of each family
#'@return a list with 3 elements:\cr
#'$mrkdat: a matrix with info on the marker dosages distribution in the FS. The
#'first two columns are for the parents: the parent name is between (brackets)
#'in the column name, their mrkdid (marker dosage ID) in row 1 and their marker
#'dosages below that. The remaining columns are for the different mrkdids
#'observed in the FS: the mrkdid itself in the column name, its frequency in
#'row 1 and its marker dosages below that. The final column gives the frequency
#'of individuals with one or more missing marker dosages.\cr
#'$hapdat: a matrix with similar layout as mrkdat, but now with the haplotype
#'dosages rather than the marker dosages. Some mrkdids (columns) may not have
#'a haplotype dosage combination assigned (if multiple possible haplotype
#'combinations result in the same marker dosages)\cr
#'usedhap: a matrix with the dosage (0 or 1) of each marker in each of the
#'used haplotypes; haplotype nrs in columns, markers in rows
#'@examples
#'data(PolyHaplotyper_small)
#'# show the results of the first FS family in the first haploblock:
#'showOneFS(FSnr=1, hbresults=phresults[[1]], mrkDosage=phdos,
#'          FS=phFS, parents=phpar)
#'
#'@export
showOneFS <- function(FSnr, hbresults, mrkDosage, FS, parents) {
  if (!all(colnames(hbresults$hapdos) %in% colnames(mrkDosage)))
    stop("colnames of hbresults$hapdos and mrkDosage don't match")
  if (!identical(colnames(hbresults$hapdos), names(hbresults$mrkdids)))
    stop("colnames of hbresults$hapdos and hbresults$mrkdids don't match")
  markers <- hbresults$markers
  if (!all(markers %in% rownames(mrkDosage)))
    stop("not all markers from hbresults in mrkDosage")
  mrkDosage <- mrkDosage[rownames(mrkDosage) %in% markers,]
  #check parents, make sure FS and parents are character:
  if (!is.matrix(parents) || nrow(parents) != length(FS))
    stop("parents must be a matrix matching FS")
  parents <- matrix(as.character(parents), nrow=nrow(parents))
  if (length(FSnr) != 1 || !(FSnr %in% (1:length(FS))))
    stop("invalid FSnr")
  FS <- as.character(FS[[FSnr]])
  FSdidstable <-
    table(hbresults$mrkdids[names(hbresults$mrkdids) %in% FS],
          useNA="always")
  res <- list()

  # first matrix mrkdat: showing the frq of all mrkdids and their marker dosages
  res$mrkdat <- matrix(NA_integer_, nrow=nrow(mrkDosage)+1,
                       ncol=length(FSdidstable)+2)
  rownames(res$mrkdat) <- c("frq", rownames(mrkDosage))
  colnames(res$mrkdat) <- c(paste0("(",parents[FSnr, 1], ")"),
                            paste0("(",parents[FSnr, 2], ")"),
                            names(FSdidstable))
  res$mrkdat[1, 1:2] <- hbresults$mrkdids[match(parents[FSnr,],
                                            colnames(mrkDosage))]
  res$mrkdat[1, 3:ncol(res$mrkdat)] <- FSdidstable
  res$mrkdat[2:nrow(res$mrkdat), 1:2] <-
    mrkDosage[, match(parents[FSnr,], colnames(mrkDosage))]
  for (i in seq_along(FSdidstable)) {
    mrkDosagecol <- which.max(hbresults$mrkdids == names(FSdidstable)[i])
    if (length(mrkDosagecol) == 1) #i.e. not for mrkdid <NA>
      res$mrkdat[2:nrow(res$mrkdat), i+2] <- mrkDosage[,mrkDosagecol]
  }

  # second matrix hapdat: same, but with haplotype dosages
  res$hapdat <- matrix(NA_integer_, nrow=nrow(hbresults$hapdos)+1,
                       ncol=length(FSdidstable)+2)
  rownames(res$hapdat) <- c("frq", rownames(hbresults$hapdos))
  colnames(res$hapdat) <- c(paste0("(",parents[FSnr, 1], ")"),
                            paste0("(",parents[FSnr, 2], ")"),
                            names(FSdidstable))
  res$hapdat[1, 1:2] <- hbresults$mrkdids[match(parents[FSnr,],
                                            colnames(hbresults$hapdos))]
  res$hapdat[1, 3:ncol(res$hapdat)] <- FSdidstable
  res$hapdat[2:nrow(res$hapdat), 1:2] <-
    hbresults$hapdos[, match(parents[FSnr,], colnames(hbresults$hapdos))]
  for (i in seq_along(FSdidstable)) {
    matcol <- which.max(hbresults$mrkdids == names(FSdidstable)[i])
    if (length(matcol) == 1) #i.e. not for mrkdid <NA>
      res$hapdat[2:nrow(res$hapdat), i+2] <- hbresults$hapdos[,matcol]
  }
  res$hapdat <- res$hapdat[rowSums(res$hapdat, na.rm=TRUE) > 0,]

  # third matrix hapcomp: composition of all haplotypes in res$hapdat
  hapnrs <- split_hapnames(hapnames=rownames(res$hapdat)[-1])$hapnrs
  r <- t(allHaplotypes(markers))
  colnames(r) <- seq_len(ncol(r))
  res$usedhap <- r[, hapnrs, drop=FALSE]
  res
} #showOneFS

#TODO:
# - impute missing marker dosages (1) in FS (2) in other?

getHaploblockResults <- function(hapresults, haploblock) {
  # Internal function used by allhap and usedhap.
  # if hapresults is already the results of a single haploblock and
  # haploblock not specified, hapresults itself is resturned.
  # Else hapresults is the list of results for all haploblocks, and the results
  # of the one indicated by haploblock is returned
  if (missing(haploblock) || is.null(haploblock)) {
    if (names(hapresults)[1] != "hapdos") {
      stop("haploblock must be specified")
    } else {
      # hapresults is already the results of a single haploblock
      return(hapresults)
    }
  }
  #haploblock specified, hapresults must be a list of results for all haploblocks
  if (length(haploblock) != 1) {
    stop("only one haploblock allowed")
  }
  if (names(hapresults)[1] == "hapdos") {
    stop("hapresults must be a list as returned by inferHaplotypes")
  }
  if (is.numeric(haploblock)) {
    if (!(haploblock %in% seq_along(hapresults)))
      stop(paste("invalid haploblock number:", haploblock))
    hbnr <- haploblock
  } else {
    if (!(haploblock %in% names(hapresults)))
      stop(paste("haploblock", haploblock, "not found"))
    hbnr <- which(names(hapresults) == haploblock)
  }
  hapresults[[hbnr]]
} #getHaploblockResult


#'@title Find all possible haplotypes
#'@description Find all possible haplotypes for a haploblock from the
#'haplotyping result
#'@usage allhap(hapresults, haploblock)
#'@param hapresults list as returned by inferhaplotypes, or one element
#'of such a list (i.e. the results for one haploblock)
#'@param haploblock if hapresults is one element of the return value of
#'inferHaplotypes, haploblock should be missing of NULL; else haploblock is
#'a single value indicating the haploblock: either its name of its index
#'in hapresults
#'@details This function works with the results of inferHaplotypes; the setting
#'of dropUnused does not affect this function
#'@return an array with all possible haplotypes. The haplotypes are in columns,
#'with the haplotype numbers as colnames; the markers are in rows.
#'@examples
#'data(PolyHaplotyper_small)
#'# show the composition of all possible haplotypes with the markers
#'# in the first haploblock:
#'allhap(hapresults=phresults, haploblock=1)
#'@export
allhap <- function(hapresults, haploblock) {
  hapresults <- getHaploblockResults(hapresults, haploblock)
  r <- t(allHaplotypes(hapresults$markers))
  colnames(r) <- seq_len(ncol(r))
  r
} # allhap

#'@title Find all used (inferred) haplotypes
#'@description Find all haplotypes for a haploblock that were inferred to be
#'present in the population (i.e. all haplotypes used for haplotyping any of
#'the individuals)
#'@usage usedhap(hapresults, haploblock)
#'@param hapresults list as returned by inferhaplotypes, or one element
#'of such a list (i.e. the results for one haploblock)
#'@param haploblock if hapresults is one element of the return value of
#'inferHaplotypes, haploblock should be missing of NULL; else haploblock is
#'a single value indicating the haploblock: either its name of its index
#'in hapresults
#'@details This function works with the results of inferHaplotypes; the setting
#'of dropUnused does not affect this function
#'@return an array with the haplotypes that are used in the
#'population. The haplotypes are in columns, with the haplotype numbers as
#'colnames; the markers are in rows.
#'@examples
#'data(PolyHaplotyper_small)
#'# show the composition of haplotypes inferred to be present
#'# in the first haploblock:
#'usedhap(hapresults=phresults, haploblock=1)
#'@export
usedhap <- function(hapresults, haploblock) {
  hapresults <- getHaploblockResults(hapresults, haploblock)
  hapnrs <-
    split_hapnames(hapnames=rownames(hapresults$hapdos)[
                       rowSums(hapresults$hapdos, na.rm=TRUE) > 0])$hapnrs
  # this takes care of the situation with dropUnused FALSE
  allhap(hapresults)[, hapnrs]
} # usedhap

#'@title merge replicate samples in dosage matrix
#'@description merge replicate samples in dosage matrix
#'@usage mergeReplicates(mrkDosage, replist, solveConflicts=TRUE)
#'@param mrkDosage a dosage matrix with markers in rows and individuals in columns.
#'row names are marker names, column names are individual names.
#'@param replist a list of character vectors, each of which has the sample names
#'of a set of replicates
#'@param solveConflicts if TRUE (default) and there are conflicting dosage
#'assignments between replicates for the same marker, the one with highest
#'frequency is used, provided the total freq of other dosages = 1 OR
#'<= 10\% of the frequency of the most frequent dosage. If
#'solveConflicts is FALSE and there are conflicting dosages, the consensus for
#'that marker will be NA
#'@details This function merges all sets of replicates, each to one column.
#'The column name of the one retained column is the first one for that set
#'in replist.\cr
#'For each set of replicates it calls getConsensusmrkDosage
#'@return a version of mrkDosage in which only one column of each set of replicates
#'is retained; this column (the first in its set as specified in replist)
#'now has the consensus scores over all replicates. Also, if mrkDosage was a
#'data.frame, it is converted into a matrix.
#'@examples
#'# construct a dosage matrix with some missing data:
#'dosmat <-
#'  matrix(c(rep(c(3,0,1),3), rep(c(1,1,2),4)), nrow=3,
#'         dimnames=list(c("mrk1","mrk2","mrk3"),
#'                       c("a1","a2","a3", "b1","b2","b3","b4")))
#'ix <- matrix(c(1,1, 3,1, 2,2, 1,3, 2,4, 1,5, 2,6), ncol=2, byrow=TRUE)
#'dosmat[ix] <- NA
#'dosmat
#'# define 2 sets of replicates:
#'reps <- list(c("a1","a2","a3"), c("b1","b2","b3","b4"))
#'# merge:
#'mergeReplicates(mrkDosage=dosmat, replist=reps)
#'# introduce a conflicting dosage:
#'dosmat[3,2] <- 2
#'# merge:
#'mergeReplicates(mrkDosage=dosmat, replist=reps)
#'@export
mergeReplicates <- function(mrkDosage, replist, solveConflicts=TRUE) {
  mrkDosage <- checkmrkDosage(mrkDosage)
  if (is.null(colnames(mrkDosage)))
    colnames(mrkDosage) <- as.character(seq_len(ncol(mrkDosage)))
  if (!is.list(replist)) replist <- list(replist)
  for (p in seq_along(replist)) {
    replist[[p]] <- as.character(replist[[p]])
    consdos <- getConsensusmrkDosage(mrkDosage=mrkDosage,
                                      indiv=replist[[p]],
                                      solveConflicts=solveConflicts)
    mrkDosage[, replist[[p]][1]] <- consdos #1st sample gets the consensus
    delcol <- colnames(mrkDosage) %in% replist[[p]][-1]
    mrkDosage <- mrkDosage[, !delcol] #delete cols of all other samples in this set
  }
  mrkDosage
} #mergeReplicates

# not exported: user function is mergeReplicates
#_#'@title get consensus marker dosages from one or more samples
#_#'@description get consensus marker dosages from one or more samples
#_#'@usage getConsensusmrkDosage(mrkDosage, indiv, solveConflicts=TRUE)
#_#'@param mrkDosage a dosage matrix with markers in rows and individuals in columns.
#_#'row names are marker names, column names are individual names.
#_#'@param indiv character vector with the names of the samples from which to
#_#'obtain consensus scores
#_#'@param solveConflicts if TRUE (default) and there are conflicting dosage
#_#'assignments between the samples for the same marker, the one with highest
#_#'frequency is used, provided the total freq of other dosages  <= 1 OR
#_#'<= 10\% of the frequency of the maximum dosage. If
#_#'solveConflicts is FALSE and there are conflicting dosages, the consensus for
#_#'that marker will be NA
#_#'@details getConsensusmrkDosage is primarily used by mergeReplicates
#_#'but can also be used by itself
#_#'@return a vector with one consensus dosage for each row of mrkDosage, and the
#_#'row names of mrkDosage as names
getConsensusmrkDosage <- function(mrkDosage, indiv, solveConflicts=TRUE) {
  indiv <- as.character(indiv)
  if (is.null(colnames(mrkDosage)))
    colnames(mrkDosage) <- as.character(seq_len(ncol(mrkDosage)))
  if (is.vector(mrkDosage)) mrkDosage <- matrix(mrkDosage, ncol=1)
  if (!is.matrix(mrkDosage)) stop("mrkDosage must be a matrix")
  if (!all(indiv %in% colnames(mrkDosage)))
    stop("not all indiv in mrkDosage")
  mrkDosage <- mrkDosage[,colnames(mrkDosage) %in% indiv, drop=FALSE]
  suppressWarnings({ #suppress warnings about all dosages NA
    mindos <- apply(mrkDosage, MARGIN=1, FUN=min, na.rm=TRUE)
    maxdos <- apply(mrkDosage, MARGIN=1, FUN=max, na.rm=TRUE)
  })
  if (solveConflicts) {
    for (i in which(mindos != maxdos)) {
      tb <- table(mrkDosage[i,]) #integer(0) is all scores NA
      suppressWarnings({ maxtb <- which(tb == max(tb)) })
      if (length(maxtb) == 1) { # 0 if all scores NA
        othercount <- sum (!is.na(mrkDosage[i,])) - tb[maxtb]
        if (othercount > 1 && othercount > 0.1 * tb[maxtb]) {
          mindos[i] <- NA
        } else mindos[i] <- as.integer(names(tb)[maxtb])
      } else mindos[i] <- NA
    }
  } else {
    #not solve conflicts, set any marker with conflicting dosages to NA:
    mindos[mindos != maxdos] <- NA #includes rows with all NA because -Inf != Inf
  }
  names(mindos) <- row.names(mrkDosage)
  mindos
} #getConsensusmrkDosage

# not exported
#_#'@title get all FS haplotype combinations expected from two parental haplotype
#_#'combinations
#_#'@description get all FS haplotype combinations expected from two parental
#_#'haplotype combinations, no considering DR, and without their frequencies
#_#'@usage getFScombs(parcomb)
#_#'@param parcomb matrix with one column for each parent and <ploidy> rows,
#_#'giving the haplotype combinations of the parents
#_#'@return a matrix with <ploidy> rows and one column for each
#_#'haplotype combination that can be generated from these two parents,
#_#'assuming polysomic inheritance with no double reduction.\cr
#_#'The colnames are NULL: haplotype combinations are not named.
getFScombs <- function(parcomb) {
  gamP1 <-   unique(getAllGametes(parcomb[, 1]), MARGIN=2)
  gamP2 <-   unique(getAllGametes(parcomb[, 2]), MARGIN=2)
  FScomb <- matrix(NA_integer_, ncol=ncol(gamP1) * ncol(gamP2),
                   nrow=nrow(parcomb))
  colnames(FScomb) <- 1:ncol(FScomb)
  for (g1 in seq_len(ncol(gamP1))) for (g2 in seq_len(ncol(gamP2))) {
    FScol <- (g1-1) * ncol(gamP2) + g2
    FScomb[, FScol] <-
      sort(c(gamP1[,g1], gamP2[,g2]))
  }
  unique(FScomb, MARGIN=2)
} #getFScombs

getAllGametes <- function (hapcomb) {
  #hapcomb: vector of the haplotype combination of a parent with <ploidy>
  #        elements
  #result: matrix with (ploidy/2) rows and one column per gamete,
  #        each element is the number (ID) of a haplotype). Within columns the
  #        haplotypes are sorted from low to high; the columns are ordered
  #        from left to right, first on row 1, then on row 2 etc (therefore
  #        duplicated columns are consecutive)
  # get from hapcomb a vector of length ploidy with all haplotypes present:
  #get the matrix with gametes:
  ploidy <- length(hapcomb)
  combn(hapcomb, ploidy/2)
} #getAllGametes

getHapcombFreq <- function(HapcombMatrix) {
  #takes a matrix, sorts it such that all columns that are identical are
  #consecutive, and returns a list with 2 elements:
  #hapcomb: a matrix like HapcombMatrix, with only unique columns
  #hcfreq: the number of times each column occurs in HapcombMatrix
  #(utility function, works with all types of matrices)
  o <- do.call("order",
               lapply(1:nrow(HapcombMatrix), function(i) HapcombMatrix[i,]))
  HapcombMatrix <- HapcombMatrix[, o, drop=FALSE]
  du <- duplicated(HapcombMatrix, MARGIN=2)
  uniqcomb <- which(!du)
  hapcounts <- c(uniqcomb[-1], length(du)+1) - uniqcomb
  list(hapcomb=HapcombMatrix[, !du, drop=FALSE], freq=hapcounts)
} #getHapcombFreq

#'@title get all gametes and their frequencies for a parental
#'haplotype combination
#'@description get all gametes and their frequencies for a parental
#'haplotype combination
#'@usage getGameteFreqs(parhac, DRrate)
#'@param parhac vector of the parental haplotype combination of length <ploidy>,
#'giving the <ploidy> haplotype numbers present
#'per haplotype giving the dosage of that haplotype, summing to ploidy
#'@param DRrate the rate of double reduction per meiosis (NOT per allele!); e.g.
#'with a DRrate of 0.04, a tetraploid parent with genotype ABCD will produce
#'a fraction of 0.04 of DR gametes AA, BB, CC and DD (each with a frequency of
#'0.01), and a fraction of 0.96 of the non-DR gametes AB, AC, AD, BC, BD, CD
#'(each with a frequency of 0.16)
#'@details for hexaploids the DR gametes consist of a duplication of one of the
#'6 parental alleles, combined with one copy of one of the other 5 alleles.\cr
#'Calculation is faster if DRrate is 0.0
#'@return a list of 2 elements:\cr
#'hapcomb: a matrix with one column per unique gamete and ploidy/2 rows.
#'Each element is the number (ID) of a haplotype). Within columns the
#'haplotypes are sorted from low to high; the columns are ordered from left
#'to right, first on row 1, then on row 2 etc\cr
#'freq: a vector of length ncol(hapcomb), with for each gamete in hapcomb its
#'frequency (0.0 ... 1.0)
#'@examples
#'# specify combination of haplotypes in a tetraploid parent:
#'hapcomb <- c(2,2,5,6) # 2 copies of haplotype 2, 1 each of 5 and 6
#'# gamete frequencies without double reduction:
#'getGameteFreqs(parhac=hapcomb, DRrate=0)
#'# gamete frequencies with 5\% double reduction:
#'getGameteFreqs(parhac=hapcomb, DRrate=0.05)
#'@export
getGameteFreqs <- function(parhac, DRrate) {
  noDR <- getHapcombFreq(getAllGametes(parhac))
  noDR$freq <- (1.0-DRrate) * noDR$freq / sum(noDR$freq)
  #get DR gametes:
  if (DRrate == 0.0) noDR else {
    # calculate the gamete freqs in case of DR:
    # make a matrix in which each of the columns will be a gamete
    ploidy <- length(parhac)
    if (ploidy==4) {
      # all DR gametes are just 1 of the 4 alleles duplicated
      haps <- unique(parhac)
      hapfrq <- table(parhac)
      DR <- list(hapcomb=rbind(haps, haps), freq=hapfrq)
      rownames(DR$hapcomb) <- NULL
    } else {
      # In a hexaploid a DR gamete contains one duplicated allele
      # combined with any of the other 5 alleles.
      # This does not generalize to ploidy 8 and higher because then there
      # may be more than one quadrivalent, and each quadrivalent can produce
      # a DR, i.e. each gamete may have more than one duplicated haplotype.
      # However, in an octo- or decaploid the chances of having two
      # quadrivalents, each producing a DR are small (about DRrate squared)
      # so we ignore that.
      x <- choose((ploidy-1), (ploidy/2-2)) #nr of gametes with one specific
      #                                      allele DR
      HMx <- matrix(integer(0), nrow=ploidy/2,
                    ncol=ploidy*x)
      for (d in 1:ploidy) {
        cols <- ((d-1) * x + 1) : (d * x)
        HMx[1:2, cols] <- parhac[d]
        HMx[3:(ploidy/2), cols] <- combn(parhac[-d], ploidy/2-2)
      }
      HMx <- apply(HMx, 2, sort) #sort each column separately
      HMx <- HMx[, order(HMx[1,], HMx[2,], HMx[3,])] #reorder all columns
      # HMx is now sorted, and comparable to that from getAllGametes
      DR <- getHapcombFreq(HMx)
    }
    # combine the noDR and DR gametes:
    DR$freq <- DRrate * DR$freq / sum(DR$freq)
    freq <- c(noDR$freq, DR$freq)
    hapcomb <- cbind(noDR$hapcomb, DR$hapcomb)
    hapcomb <- hapcomb[, freq>0, drop=FALSE]
    freq <- freq[freq>0]
    #find the order to sort the hapcomb columns (on all 2 or 3 rows):
    o <- do.call("order", lapply(1:nrow(hapcomb), function(i) hapcomb[i,]))
    hapcomb <- hapcomb[, o]
    freq <- freq [o]
    du <- duplicated(hapcomb, MARGIN=2)
    uniqcomb <- which(!du)
    #sum the freqs per set of duplicates:
    uniqfreq <- numeric(length(uniqcomb))
    uhx <- c(uniqcomb, length(du)+1)
    for (u in seq_along(uniqfreq))
      uniqfreq[u] <- sum(freq[uhx[u]:(uhx[u+1]-1)])
    list(hapcomb=hapcomb[, uniqcomb, drop=FALSE], freq=uniqfreq)
  }
} #getGameteFreqs

#'@title get all FS haplotype combinations expected from two parental haplotype
#'combinations, with their frequencies
#'@description get all FS haplotype combinations expected from two parental
#'haplotype combinations, with their frequencies. Different from getFScombs in
#'that it returns the unique FS combinations as well as their frequencies, and
#'that it also considers Double Reduction
#'@usage getFSfreqs(parhac, DRrate)
#'@param parhac matrix with one column for each parent and <ploidy> rows,
#'giving the haplotype combinations for each parent
#'@param DRrate the rate of double reduction per meiosis (NOT per allele!); e.g.
#'with a DRrate of 0.04, a tetraploid parent with genotype ABCD will produce
#'a fraction of 0.04 of DR gametes AA, BB, CC and DD (each with a frequency of
#'0.01), and a fraction of 0.96 of the non-DR gametes AB, AC, AD, BC, BD, CD
#'(each with a frequency of 0.16)
#'@return a list of 2 elements:\cr
#'FShac: a matrix with one column per unique FS haplotype combination and
#'the same row count as parhac, giving the FS haplotype combinations. There are
#'no duplicated columns but several (not necessarily adjacent) columns
#'may correspond to the same mrkdid\cr
#'freq: a vector of length ncol(hapcomb), with for each FS haplotype combination
#'in hapcomb its frequency (0.0 ... 1.0). The colnames are NULL: haplotype
#'combinations are not named.
#'@examples
#'# specify combinations of haplotypes in two tetraploid parents:
#'hapcomb <- matrix(c(2,2,5,6, 1,2,5,5), ncol=2)
#'# FS frequencies without double reduction:
#'getFSfreqs(parhac=hapcomb, DRrate=0)
#'# FS frequencies with 5\% double reduction:
#'getFSfreqs(parhac=hapcomb, DRrate=0.05)
#'@export
getFSfreqs <- function(parhac, DRrate) {
  gamP1 <-   getGameteFreqs(parhac[, 1], DRrate)
  gamP2 <-   getGameteFreqs(parhac[, 2], DRrate)
  getFSfreqs_gamfrq(gamP1, gamP2)
}

getFSfreqs_gamfrq <- function(gamP1, gamP2) {
  #function that does most of the actual work of getFSfreqs;
  #also called by checkPedigree
  #gamP1, gamP2: result of getGameteFrqs for two parents and specified DRrate
  #return value: see getFSFreqs
  FScomb <- matrix(NA_integer_, ncol=ncol(gamP1$hapcomb) * ncol(gamP2$hapcomb),
                   nrow=nrow(gamP1$hapcomb) + nrow(gamP2$hapcomb))
  FSfreq <- numeric(ncol(FScomb))
  colnames(FScomb) <- 1:ncol(FScomb)
  for (g1 in seq_len(ncol(gamP1$hapcomb)))
    for (g2 in seq_len(ncol(gamP2$hapcomb))) {
    FScol <- (g1-1) * ncol(gamP2$hapcomb) + g2
    FScomb[, FScol] <-
      sort(c(gamP1$hapcomb[,g1], gamP2$hapcomb[,g2]))
    FSfreq[FScol] <- gamP1$freq[g1] * gamP2$freq[g2]
  }
  #find the order to sort the FScomb columns (on all rows):
  o <- do.call("order", lapply(1:nrow(FScomb), function(i) FScomb[i,]))
  FScomb <- FScomb[, o, drop=FALSE]
  FSfreq <- FSfreq [o]
  du <- duplicated(FScomb, MARGIN=2)
  uniqcomb <- which(!du)
  #sum the freqs per set of duplicates:
  uniqfreq <- numeric(length(uniqcomb))
  uhx <- c(uniqcomb, length(du)+1)
  for (u in seq_along(uniqfreq))
    uniqfreq[u] <- sum(FSfreq[uhx[u]:(uhx[u+1]-1)])
  #remove the duplicated columns and the column names:
  FScomb <- FScomb[, uniqcomb, drop=FALSE]
  colnames(FScomb) <- NULL
  list(FShac=FScomb, freq=uniqfreq)
} #getFSfreqs_gamfrq

checkFS_TwoParents <- function(pargamlist, FShac) {
  # pargamlist: a list with two elements, one for each parent. Each element
  #             is itseld a list with two elements named $noDR and $DR; each of
  #             these is a list as returned by getGameteFreqs, the noDR obtained
  #             with DRrate=0 and the other with DRrate>0 (say 0.04, we only
  #             need the set of possible gametes with DR; their freq (and so
  #             the exact DRrate) is not used).
  # FShac: matrix with ploidy rows and one column for each FS indiv, does
  #        not contain NAs: the haplotype combinations of the FS indiv.
  #        FShac is sorted such that identical columns are adjacent
  nrows <- sapply(pargamlist, length)
  if (max(nrows) != min(nrows)) stop("nrows differ")
  matchparents <- matrix(FALSE, nrow=nrows[1], ncol=ncol(FShac),
                         dimnames=list(NULL, colnames(FShac)))
  if (!is.null(names(pargamlist[[1]])))
    rownames(matchparents) <- names(pargamlist[[1]])
  # matchparents is the return value: for each FS indiv it is TRUE if it
  # is a possible progeny of the parents, without and with assuming DR,
  # and else FALSE. Cannot be NA as all FS indiv with missing hac were already
  # excluded from FShac
  du <- duplicated(FShac, MARGIN=2)
  uniqcomb <- which(!du)
  #sum the freqs per set of duplicates:
  uhx <- c(uniqcomb, length(du)+1)
  uniqfreq <- diff(uhx)
  #now we have the inferred haplotype combinations in FShac[, uniqcomb] and
  #their frequencies in uniqfreq
  #next we determine how many of these are expected without and with DR:
  for (m in seq_len(nrow(matchparents))) { #modes: 1 = no DR, 2 = with DR
    FSexp <- getFSfreqs_gamfrq(pargamlist[[1]][[m]],
                               pargamlist[[2]][[m]])
    for (hc in seq_along(uniqcomb)) {
      indiv <- colnames(FShac)[uhx[hc]:(uhx[hc+1]-1)]
      alleq <- function(x, y) { isTRUE(all.equal(x, y)) }
      #        identical returns FALSE, all.equal TRUE for a matching column,
      #        but all.equal returns a character for a non-matching column
      matchcols <- apply(FSexp$FShac, MARGIN=2,
                         FUN=alleq, FShac[, uniqcomb[hc]])
      # so, for each expected hac (each column of FSexp$FShac) we get
      # a TRUE or FALSE result on whether it matches the currently
      # studied observed hapcomb (FShac[, uniqcomb[hc]])
      if (any(matchcols)) {
        matchparents[m, indiv] <- TRUE
      }
    } #for hc
  } # for m
  matchparents
} # checkFS_TwoParents

checkFS_OneParent <- function(pargamlist, FShac, nhap) {
  # pargamlist: a list with two elements named $noDR and $DR; each a list
  #             as returned by getGameteFreqs, the noDR obtained with DRrate=0
  #             and the other with DRrate>0 (say 0.04, we only need the set of
  #             possible gametes with DR; their freq (and so the exact
  #             DRrate) is not used). Or a list with only one such element.
  # FShac: matrix with ploidy rows and one column for each FS indiv, does
  #        not contain NAs: the haplotype combinations of the FS indiv.
  #        FShac is sorted such that identical columns are adjacent
  # nhap: the total number of possible haplotypes in this haploblock
  matchparent <- matrix(FALSE, nrow=length(pargamlist), ncol=ncol(FShac),
                        dimnames=list(NULL, colnames(FShac)))
  if (!is.null(names(pargamlist)))
    rownames(matchparent) <- names(pargamlist)
  # matchparents is the return value: for each FS indiv it is TRUE if it
  # is a possible progeny of the parent, without and with assuming DR,
  # and else FALSE. Cannot be NA as all FS indiv with missing hac were already
  # excluded from FShac
  du <- duplicated(FShac, MARGIN=2)
  uniqcomb <- which(!du)
  uniqhad <- hapcomb2hapdos(FShac[, uniqcomb, drop=FALSE], nhap)
  #sum the freqs per set of duplicates:
  uhx <- c(uniqcomb, length(du)+1)
  uniqfreq <- diff(uhx) #uhx[-1] - uhx[-length(uhx)]
  #now we have the inferred haplotype combinations and their frequencies;
  #next we determine how many of these are expected without and with DR.
  #we have only the expected gametes of par1; a hapcomb in HS can only
  #occur if there is at least one gamete of par1 that has less than or
  #equal the nr of copies of all haplotypes
  for (m in seq_along(pargamlist)) { #modes: 1 = no DR, 2 = with DR
    expgamhad <- hapcomb2hapdos(pargamlist[[m]]$hapcomb, nhap)
    for (hc in seq_along(uniqcomb)) {
      indiv <- colnames(FShac)[uhx[hc]:(uhx[hc+1]-1)]
      smeq <- function(x, y) { all(x <= y) }
      matchcols <- apply(expgamhad, MARGIN=2,
                         FUN=smeq, uniqhad[, hc]) #0, 1 or more TRUE
      # so, for each expected gamete (each column of expgamhad) we get
      # a TRUE or FALSE result on whether it could have participated
      # in currently studied observed FShad (uniqhad[, hc])
      if (any(matchcols)) {
        #indivs with this uniqhad match with the parent
        matchparent[m, indiv] <- TRUE
      }
    } # for hc
  } # for m
  matchparent
} #checkFS_OneParent


##!title check the consistency of a single haploblock genotypes over a pedigree
##!description For one haploblock, check whether the inheritance of inferred
##!haplotypes over the pedigree is consistent (without or with DR, and assuming
##!polysomic inheritance)
##!usage pedigreeHapCheck_OneHaploblock(ped, parentnames, mrkDosage, hapcomb,
##!nhap)
##!param ped data.frame or matrix listing the pedigree. Column 1 has names of
##!all individuals (no duplicates, no NA), column 2 and 3 have the parents
##!(duplicates and NA allowed, also individuals with one known parent or
##!with two identical parents allowed). No sorting of the pedigree is needed but
##!all parents must also be listed as individuals
##!param parentnames names of all individuals that occur as parents in ped
##!(must be correct, not checked, serves to avoid double work)
##!param mrkDosage a matrix with one row for each marker in the haploblock
##!and one column for each individual, with the observed dosages. Individuals
##!with one or more NA dosages are considered to have no marker data. Individual
##!names are the colnames, no duplicates allowed; the set of individuals
##! need not be the same as those in ped; individual names common to mrkDosage and
##! ped must be identical
##!param hapcomb a matrix with <ploidy> rows and one column for each individual,
##!with the inferred haplotype combinations per individual. Individual names are
##!the colnames, the set of individuals need not be the same as those in ped;
##!individual names common to hapcomb and ped must be identical. hapcomb has
##!no missing values: all individuals with no assigned hapcomb were omitted.
##!return a list with two items:\cr
##! $ped: the same pedigree as the input, with two additional columns
##! matchParents_noDR and matchParents_withDR. FALSE indicates a conflict with
##! (one of) the parents, FALSE no conflict, NA if individual or both its
##! parents NA.
##! $parents: a data frame with one row for each individual that occurs as a
##! parent, with some statistics about all its progeny and its match with all
##! its progeny
# don't export
pedigreeHapCheck_OneHaploblock <- function(ped, parentnames, mrkDosage,
                                           hapcomb, nhap, nohapdata,
                                           imputedGeno) {
  # - for each parent in the pedigree generate all its gametes
  #   * without and with DR
  # - sort the pedigree into FS families
  # - for each (reciprocal) cross generate all its possible progeny
  #   * without and with DR
  # - tabulate the inferred haplotype combinations for each FS family
  #   and check for each if it is expected without and with DR
  # if nohapdata is TRUE the arrays in the result are filled with NAs or zeroes
  # as appropriate
  # Result: a list with 2 elements:
  # $ped: a logical matrix with one row per individual and 5 columns:
  #       mrk, imp, hap, noDR, withDR: is it genotyped, imputed, haplotyped,
  #       does it conflict with its parents without or with allowing DR
  # $parents: an integer matrix with one row per individuals that is a parent in
  #       ped, and columns par_mrkdata, par_hapdata, totprogeny, prog_mrkdata,
  #       prog_hapdata, nonDRmatch, DRmatch: first 2 are 0/1 (is parent
  #       genotyped and haplotyped), then the nr of direct progeny and the nr
  #       of progeny that is genotyped and haplotyped, and finally the nr of
  #       progeny that matches (can be produced) without or with allowing DR
  peddat <- matrix(NA, nrow=nrow(ped), ncol=5,
                   dimnames=list(ped[,1], c("mrk", "imp", "hap", "noDR", "withDR")))
  mode0_ped <- which(colnames(peddat) == "noDR") - 1
  peddat[, "mrk"] <- #mrk: complete marker dosage data for this indiv in this haploblock?
    !is.na(colSums(mrkDosage[, match(rownames(peddat), colnames(mrkDosage)), drop=FALSE]))
  if (is.null(imputedGeno)) peddat[, "imp"] <- FALSE else
  peddat[, "imp"] <- #imp: was the marker genotype of this indiv imputed for this block?
    rownames(peddat) %in% colnames(imputedGeno)
  parcols <- c("par_mrkdata", "par_hapdata", "totprogeny",
               "prog_mrkdata", "prog_hapdata", "nonDRmatch", "DRmatch")
  parents <- matrix(as.integer(0), nrow=length(parentnames), ncol=length(parcols),
                    dimnames=list(parentnames, parcols))
  mode0_par <- which(parcols == "nonDRmatch") - 1 # base column number for modes
  parents[, "par_mrkdata"] <-
    !is.na(colSums(mrkDosage[, match(parentnames, colnames(mrkDosage)), drop=FALSE]))

  if (nohapdata) {
    peddat[, "hap"] <- 0
    parents[, c("par_hapdata", "prog_hapdata")] <- 0
  } else {
    if (any(is.na(hapcomb[1,])))
        stop("hapcomb should not be empty and not contain missing data")
    # sort hapcomb such that identical columns are adjacent:
    # (hapcomb only contains individuals that have a haplotype combination
    # assigned, i.e. no NAs)
    o <- do.call("order", lapply(1:nrow(hapcomb), function(i) hapcomb[i,]))
    hapcomb <- hapcomb[, o, drop=FALSE]
    peddat[, "hap"] <- #hap: assigned haplotype combi for this indiv?
      rownames(peddat) %in% colnames(hapcomb) #no test for NAs needed: not present
    parentsNodata <- setdiff(parentnames, colnames(hapcomb))
    parents[, "par_hapdata"] <-
      parentnames %in% colnames(hapcomb) #hapcomb has only haplotyped individuals
    # first create a list with the possible gametes from each parent:
    pargamlist <- list()
    for (par in parentnames) {
      if (par %in% parentsNodata) {
        pargamlist[[par]] <- "nodata"
      } else {
        pargamlist[[par]] <- list()
        pargamlist[[par]]$noDR <-
          getGameteFreqs(parhac=hapcomb[, par], DRrate=0)
        pargamlist[[par]]$DR <-
          getGameteFreqs(parhac=hapcomb[, par], DRrate=0.04)
        #  we only need the set of possible gametes; their freq (and so the exact
        #  DRrate) is not used
      }
    }
    #for each fullsib family, get all possible progeny haplotype combinations
    #and check the inferred ones against these:
    for (p1 in seq_along(parentnames)) {
      par1 <- parentnames[p1]
      for (p2 in p1:nrow(parents)) {
        par2 <- parentnames[p2]
        FS <- ped[!is.na(ped[, 2]) & ped[, 2] %in% c(par1, par2) &
                    !is.na(ped[, 3]) & ped[, 3] %in% c(par1, par2),, drop=FALSE]
        if (par1 != par2) FS <- FS[FS[, 2] != FS[, 3],, drop=FALSE]
        if (nrow(FS) > 0) {
          # this parental combination occurs in ped
          parents[c(p1, p2), "totprogeny"] <-
            parents[c(p1, p2), "totprogeny"] + nrow(FS)
          parents[c(p1, p2), "prog_mrkdata"] <-
            parents[c(p1, p2), "prog_mrkdata"] +
            sum(!is.na(colSums(mrkDosage[, colnames(mrkDosage) %in% FS[, 1],
                                         drop=FALSE])))
          FShac <- hapcomb[, colnames(hapcomb) %in% FS[, 1], drop=FALSE]
          parents[c(p1, p2), "prog_hapdata"] <-
            parents[c(p1, p2), "prog_hapdata"] + ncol(FShac)
          if (all(parents[c(p1, p2), "par_hapdata"] > 0) && ncol(FShac) > 0) {
            #both parents have haplotypes and at least some FS indiv in hapcomb
            chk <- checkFS_TwoParents(pargamlist=pargamlist[c(p1,p2)],
                                      FShac=FShac)
            pedrows <- match(colnames(chk), rownames(peddat))
            peddat[pedrows, mode0_ped+(1:2)] <- t(chk)
            parents[c(p1, p2), mode0_par+(1:2)] <-
              parents[c(p1, p2), mode0_par+(1:2)] +
              rep(rowSums(chk, na.rm=TRUE), each=2)
          }
        } # nrow(FS) > 0
      } # for p2
      #and for the individuals of which only one parent is known in ped and has
      #been haplotyped:
      HS <- ped[(!is.na(ped[, 2]) & ped[, 2] == par1 &
                   (is.na(ped[, 3]) | ped[, 3] %in% parentsNodata)) |
                (!is.na(ped[, 3]) & ped[, 3] == par1 &
                   (is.na(ped[, 2]) | ped[, 2] %in% parentsNodata)),, drop=FALSE]
      # includes all indiv with par1 as one parent and the other parent
      # unknown or with missing hapcomb
      if (nrow(HS) > 0) {
        # there are individuals in ped with p1 as only parent
        #
        # only the HS indivs without known second parent must be added to
        # the totals, the others (where the second parent is known but
        # has no assigned hapcomb) were already counted above (because the
        # FS were defined on the ped, which also includes parents without data)
        OneParent <- HS[rowSums(is.na(HS[, 2:3, drop=FALSE])) == 1,, drop=FALSE]
        parents[p1, "totprogeny"] <- parents[p1, "totprogeny"] + nrow(OneParent)
        parents[p1, "prog_mrkdata"] <-
          parents[p1, "prog_mrkdata"] +
          sum(!is.na(colSums(mrkDosage[, colnames(mrkDosage) %in% OneParent[, 1],
                                       drop=FALSE])))
        parents[p1, "prog_hapdata"] <-
          parents[p1, "prog_hapdata"] +
          sum(!is.na(hapcomb[1, colnames(hapcomb) %in% OneParent[, 1]]))
        # for the whole HS (all indiv with only one known parent plus
        # all indiv with two parents where the second parent has missing haplodata)
        # we calculate the nr that match with this parent (which was not done
        # earlier)
        HShac <- hapcomb[, colnames(hapcomb) %in% HS[, 1], drop=FALSE]
        if (parents[p1, 2] && ncol(HShac) > 0 ) {
          chk <- checkFS_OneParent(pargamlist=pargamlist[[p1]],
                                    FShac=HShac, nhap=nhap)
          pedrows <- match(colnames(chk), rownames(peddat))
          peddat[pedrows, mode0_ped+(1:2)] <- t(chk)
          parents[p1, mode0_par+(1:2)] <-
            parents[p1, mode0_par+(1:2)] + rowSums(chk, na.rm=TRUE)
        }
      } # nrow(HS) > 0
    } # for p1
  } # at least some indiv with haplotypes
  list(ped=peddat, parents=parents)
} #pedigreeHapCheck_OneHaploblock

#'@title check the consistency of haploblock genotypes over a pedigree
#'@description For all haploblocks, check whether the inheritance of inferred
#'haplotypes over the pedigree is consistent without or with allowing for
#'double reduction (DR), and assuming polysomic inheritance.
#'@usage pedigreeHapCheck(ped, mrkDosage, haploblock, hapresults)
#'@param ped data.frame or matrix listing the pedigree. Column 1 has names of
#'all individuals (no duplicates, no NA), column 2 and 3 have the parents
#'(duplicates and NA allowed, also individuals with one known parent or
#'with two identical parents allowed). No sorting of the pedigree is needed.
#'All parents should also be listed as individuals; if that is not the case
#'the missing parents will be added as founders and a warning will be issued.
#'@param mrkDosage a matrix with one row for each marker
#'and one column for each individual, with the observed dosages. Individuals
#'with one or more NA dosages are considered to have no marker data. Individual
#'names are the colnames, no duplicates allowed; the set of individuals
#'needs not be the same as those in ped; individual names common to mrkDosage
#'and ped must be exactly identical (upper/lower case and whitespace included)
#'@param haploblock a list of character vectors. The names are the names of the
#'haploblocks, the character vectors have the names of the markers in each
#'haploblock
#'@param hapresults a list as returned by inferHaplotypes, with one
#'item per haploblock, containing at least all those in the list
#'specified by haploblock
#'@return a list with two items:\cr
#'- ped_arr: a 3D logical array with dimensions individuals, diagnostics and
#'haploblocks. For each individual and each haploblock there are 4
#'diagnostics:\cr
#'  * mrk: does the individual have complete marker dosage data?\cr
#'  * imp: were the marker dosages for this individual imputed?\cr
#'  * hap: is there a haplotype combination assigned to the individual?\cr
#'  * noDR: does the haplotype genotype of the individual match that of its
#'  parents, assuming polysomic inheritance but no Double Reduction? NA if the
#'  individual or both its parents do not have a haplotype genotype assigned\cr
#'  * withDR: as noDR, but allowing Double Reduction\cr
#'- parents_arr: a 3D integer array with dimensions parents (all individuals
#'that occur as parents in the pedigree), diagnostics and haploblocks.
#'For each parent and each haploblock there are 7 diagnostics:\cr
#' * par_mrkdata: 0=FALSE, 1=TRUE, does this parent have complete marker data?\cr
#' * par_hapdata: 0=FALSE, 1=TRUE, does this parent have a haplotype
#' genotype assigned?\cr
#' * totprogeny: how many first-generation progeny (children) does this parent
#' have (combined over all its matings, both as mother and as father)\cr
#' * prog_mrkdata: how many progeny have complete marker dosage data?\cr
#' * prog_hapdata: how many progeny have a haplotype combination assigned?\cr
#' * nonDRmatch: how many progeny have a haplotype combination that is
#' compatible with their parent's haplotype combinations, assuming polysomic
#' inheritance but no Double Reduction\cr
#' * DRmatch: as nonDRmatch, but also allowing DR\cr
#'Both ped_arr and parents_arr contain all haploblocks in haploblock, also
#'those skipped because of too many markers and those without any haplotyped
#'individuals. These can be excluded by excluding them from the haploblock list.
#'@examples
#'data(PolyHaplotyper_small)
#'phchk <- pedigreeHapCheck(ped=phped, mrkDosage=phdos, haploblock=phblocks,
#'                          hapresults=phresults)
#'# show the top of the ped_arr for haploblock 1:
#'phchk$ped_arr[1:6,,1]
#'# show the top of the parents_arr for haploblock 1:
#'phchk$parents_arr[1:6,,1]

#'@export
pedigreeHapCheck <- function(ped, mrkDosage, haploblock,
                                hapresults) {
  hbnames <- names(haploblock)
  if (!all(hbnames %in% names(hapresults)))
    stop("not all haploblocks specified by haploblock occur in hapresults")
  hapresults <- hapresults[hbnames] # same length and same order as haploblock
  #which haploblocks were done (i.e. have a non-error result):
  hbislist <- sapply(hapresults, FUN=is.list)
  hbnomess <- sapply(hapresults, FUN=function(x) (is.null(x$message) || x$message==""))
  hbrows <- sapply(hapresults, FUN=function(x) (!is.null(x$hapdos) && nrow(x$hapdos) > 0)) # this last condition added by Giorgio
  hbdone <- hbislist & hbnomess & hbrows
  #check if all parents are also listed as individuals in pedigree:
  ped <- cbind(as.character(ped[,1]), as.character(ped[,2]),
                 as.character(ped[,3])) #works for data.frame and matrix
  parentnames <- sort(setdiff(c(ped[,2], ped[,3]), NA))
  missingpar <- !(parentnames %in% ped[,1])
  if (sum(missingpar) > 0) {
    message(paste(sum(missingpar),
                  "parents are not listed as individuals in ped;",
                  "added as founders"))
    addparents <- parentnames[missingpar]
    ped <- rbind(ped, cbind(addparents, NA, NA))
  }
  pedarr <- array(NA, dim=c(nrow(ped), 5, length(haploblock)),
                  dimnames=list(as.character(ped[,1]),
                                c("mrk", "imp", "hap", "noDR", "withDR"),
                                hbnames))
  pararr <- array(NA_integer_, dim=c(length(parentnames), 7, length(haploblock)),
                  dimnames=list(parentnames,
                                c("par_mrk", "par_hap", "prog_tot",
                                  "prog_mrk", "prog_hap",
                                  "nonDRmatch", "DRmatch"),
                                hbnames))
  for (h in seq_along(haploblock)) {
    hbname <- hbnames[h]
    cat(paste0("checking haploblock ", h, " of ", length(hbdone), ": ",
               hbname,"\n")); flush.console()
    if (hbdone[h]) {
      nhap <- 2 ^ length(hapresults[[hbname]]$markers)
      tmp <-
        hapresults[[hbname]]$hapdos[, !is.na(hapresults[[hbname]]$hapdos[1,]),
                                    drop=FALSE]
    } else {
      # not hbdone[h]
      tmp <- matrix(integer(0), ncol=0) # to signal that no individuals haplotyped
    }
    if (ncol(tmp) == 0) {
        cat(paste("haploblock ", h, ": no individuals without missing values\n"))
    } else {
        cs <- colSums(tmp)
        if (anyNA(cs)) stop(paste("haploblock ", h, ": inconsistent missing values"))
        # so tmp has 1 or more rows and cols and no missing values
        if (min(cs) != max(cs)) stop(paste("haploblock ", h, ": inconsistent haplotype dosages"))
        ploidy <- min(cs)
        tmp <- hapdos2hapcomb(hapdos=tmp, ploidy=ploidy)
    }
    if ("imputedGeno" %in% names(hapresults[[h]]))
      imputedGeno <- hapresults[[h]]$imputedGeno else imputedGeno <- NULL
    pedchk <- pedigreeHapCheck_OneHaploblock(ped=ped,
      parentnames=parentnames,
      mrkDosage=mrkDosage[rownames(mrkDosage) %in% haploblock[[hbname]],, drop=FALSE],
      hapcomb=tmp, nhap=nhap, nohapdata=!hbdone[h], imputedGeno=imputedGeno)
    pedarr[,, h] <- pedchk$ped
    pararr[,, h] <- pedchk$parents
  } # for h
  list(ped_arr=pedarr, parents_arr=pararr)
} #pedigreeHapCheck

#'@title generate an overview of the results by haploblock and by FS family
#'@description generate an overview of the results by haploblock and by
#'FS family
#'@usage overviewByFS(haploblock, parents, FS, hapresults)
#'@param haploblock a list of character vectors. The names are the names of the
#'haploblocks, the character vectors have the names of the markers in each
#'haploblock.
#'@param parents a matrix of two colums and one row per FS, containing the names
#'of the two parents of each FS
#'@param FS a list of character vectors with the names of the samples for
#'each FS
#'@param hapresults a list as returned by inferHaplotypes, with one
#'item per haploblock, containing at least all those in the list
#'specified by haploblock
#'@return a list with two items:\cr
#'- ovw: an integer matrix with one row per haploblock and the following
#'columns:\cr
#'  * nmrk: the number of markers in the haploblock\cr
#'  * nhap: the number of different haplotypes assigned over all individuals
#'    (NA if no solution was found for this haplotype; the reason for that is
#'    listed in the first column of item messages of the return value)\cr
#'  * for each FS family a set of 6 columns:\cr
#'  + parmrk (0, 1 or 2: the number of parents with complete marker data)\cr
#'  + fit (0=FALSE or 1=TRUE), indicating if a solution for the FS was found
#'    based on polysomic inheritance)\cr
#'  + mrk: the number of FS progeny with complete marker data\cr
#'  + imp: the number of FS progeny where complete marker data were imputed\cr
#'  + hap: the number of FS progeny with assigned haplotype combinations.
#'    hap will be less than the mrk value if the same FS marker genotype can be
#'    produced with different combinations of haplotypes that are all compatible
#'    with the parental haplotype combinations) or if some FS marker genotypes
#'    cannot be produced by the inferred parental haplotype combinations\cr
#'  + P: the chi-squared P-value of the best fitting solution, even if this is
#'    discarded because of lack of fit).\cr
#'  * For "rest" (all individuals that are not part of the FS's or their
#'  parents) and "all" (all individuals) there are also columns mrk and hap,
#'  and for "all" there is also a column imp, similar to those for the
#'  FS families. The numbers for "all" are the sums
#'  of those for the FS families, the FS parents (some FS may share a parent
#'  but shared parents are counted only one) and the "rest".\cr
#'- messages  : a character matrix with one row per haploblock and the following
#'columns:\cr
#'  * haploblock: the reason why there is no solution for the
#'  haploblock ("" if there is a solution)
#'  * one column for each FS family with a possible message or "". A message
#'  can indicate a failure to find a solution for the FS family but may also
#'  describe less significant problems, such as some progeny with unexpected
#'  marker dosages etc.
#'@examples
#'data(PolyHaplotyper_small)
#'overviewByFS(haploblock=phblocks, parents=phpar, FS=phFS,
#'             hapresults=phresults)
#'@export
overviewByFS <- function(haploblock, parents, FS, hapresults) {
  hbnames <- names(haploblock)
  if (!all(hbnames %in% names(hapresults)))
    stop("not all haploblocks specified by haploblock occur in hapresults")
  hapresults <- hapresults[hbnames] # same length and same order as haploblock
  nFS <- length(FS)
  if (nFS == 0) stop("no FS's")
  FS <- lapply(FS, as.character)
  parents <- matrix(as.character(parents), ncol=2)
  perFScols <- c("parmrk", "fit", "P", "mrk", "imp", "hap")
  lastFScol <- 2 + length(perFScols)*nFS
  FSsok <- TRUE
  ovw <- matrix(NA_integer_, nrow=length(hapresults),
                ncol=lastFScol + 5)
  rownames(ovw) <- names(hapresults)
  colnames(ovw) <- c("nmrk", "nhap",
                     paste0(perFScols,
                            rep(1:nFS, each=length(perFScols))),
                     "mrkrest", "haprest", "mrkall", "impall", "hapall")
  ovw[, 1] <- sapply(haploblock, length) # nmrk
  messages <- matrix("", nrow=length(hapresults), ncol=1 + nFS)
  rownames(messages) <- names(hapresults)
  colnames(messages) <- c("haploblock", paste0("FS", padded(1:nFS)))
  for (hb in seq_along(hapresults)) {
    if (hapresults[[hb]]$message != "") {
      messages[hb, 1] <- hapresults[[hb]]$message
    } else {
      FSsok <- FSsok && "FSfit" %in% names(hapresults[[hb]])
      if (FSsok)
        messages[hb, 2:(nFS+1)] <- hapresults[[hb]]$FSmessages
      ovw[hb, 2] <-
        sum(rowSums(hapresults[[hb]]$hapdos, na.rm=TRUE) > 0) #nhap
      totimputed <- 0
      for (fs in 1:nFS) {
        # parmrk: how many parents have complete marker data?
        basecol <- 3 + length(perFScols)*(fs-1)
        ovw[hb, basecol] <-
          sum(!is.na(hapresults[[hb]]$mrkdids[
            names(hapresults[[hb]]$mrkdids) %in% parents[fs,] ]))
        # imp: how many FS indiv have imputed mrk data?
        if ("imputedGeno" %in% names(hapresults[[hb]])) {
          ovw[hb, basecol+4] <-
            sum(colnames(hapresults[[hb]]$imputedGeno) %in% FS[[fs]])
          totimputed <- totimputed + ovw[hb, basecol+4]
        } else ovw[hb, basecol+4] <- 0
        # mrk: how many FS indiv have complete marker data?
        ovw[hb, basecol+3] <-
          sum(!is.na(hapresults[[hb]]$mrkdids[
            names(hapresults[[hb]]$mrkdids) %in% FS[[fs]] ])) -
          ovw[hb, basecol+4] # mrkdids present for original mrk and imputed
        # hap: how many FS indiv have inferred haplotype dosages?
        if (nrow(hapresults[[hb]]$hapdos) == 0) {
          ovw[hb, basecol+5] <- 0
        } else {
          hap <- which(!is.na(hapresults[[hb]]$hapdos[1,]))
          hap <- colnames(hapresults[[hb]]$hapdos)[hap]
          hap <- hap[hap %in% FS[[fs]]]
          ovw[hb, basecol+5] <-
            sum(!is.na(hapresults[[hb]]$hapdos[1,
              colnames(hapresults[[hb]]$hapdos) %in% FS[[fs]] ]))
        }
        if (FSsok) {
          # fit: was an FS solution found?
          ovw[hb, basecol+1] <- hapresults[[hb]]$FSfit[fs]
          # P: P-value for best segregation, even if rejected:
          ovw[hb, basecol+2] <- hapresults[[hb]]$FSpval[fs]
        }
      } # for fs
      #data for "rest": the non-FS, non-FS-parents indiv
      rest <- setdiff(names(hapresults[[hb]]$mrkdids),
                      c(unlist(FS), parents))
      # mrkrest: how many rest indiv with complete marker data
      ovw[hb, lastFScol + 1] <-
        sum(!is.na(hapresults[[hb]]$mrkdids[
          names(hapresults[[hb]]$mrkdids) %in% rest]))
      # haprest: how many rest indiv haplotyped
      if (nrow(hapresults[[hb]]$hapdos) == 0) {
        ovw[hb, lastFScol + 2] <- 0
      } else {
        ovw[hb, lastFScol + 2] <-
          sum(!is.na(hapresults[[hb]]$hapdos[1,
            colnames(hapresults[[hb]]$hapdos) %in% rest]))
      }
      #data for all indiv together:
      # mrkall: how many indiv with complete marker data
      ovw[hb, lastFScol + 3] <-
        sum(!is.na(hapresults[[hb]]$mrkdids)) - totimputed
      # impall: how many indiv with imputed mrk genotypes
      ovw[hb, lastFScol + 4] <- totimputed
      # hapall: how many indiv haplotyped
      if (nrow(hapresults[[hb]]$hapdos) == 0) {
        ovw[hb, lastFScol + 5] <- 0
      } else {
        ovw[hb, lastFScol + 5] <-
        sum(!is.na(hapresults[[hb]]$hapdos[1,]))
      }
    } # message == ""
  } # for hb
  if (!FSsok)
    warning("hapresults were calculated without specifying FS families")
  list(ovw=ovw, messages=messages)
} # overviewByFS

calcPedstats <- function(pedchk, indiv, haploblocks) {
  # function used only in calcStatistics; for input and output see there
  if (missing(indiv) || is.null(indiv))
    indiv <- dimnames(pedchk$ped_arr)[[1]]
  indix <- match(indiv, dimnames(pedchk$ped_arr)[[1]])
  if (length(indix) == 0 || anyNA(indix))
    stop("names in indiv not found in pedchk")
  if (missing(haploblocks) || is.null(haploblocks)) {
    haploblocks <- dimnames(pedchk$ped_arr)[[3]]
  } else if (is.numeric(haploblocks)) {
    haploblocks <- dimnames(pedchk$ped_arr)[[3]][haploblocks]
  } else if (is.list(haploblocks)) {
    haploblocks <- names(haploblocks)
  }
  if (!all(haploblocks %in% dimnames(pedchk$ped_arr)[[3]]))
    stop("not all haploblocks occur in pedchk")
  hbix <- match(haploblocks, dimnames(pedchk$ped_arr)[[3]])
  if (length(hbix) == 0 || anyNA(hbix))
    stop("no selected haploblocks present in pedchk")
  pedstats <- data.frame(
    haploblock = haploblocks,
    stringsAsFactors = FALSE
  )
  pedstats$mrk <- apply(pedchk$ped_arr[indix, "mrk", hbix, drop=FALSE],
                        MARGIN=3, FUN=sum, na.rm=TRUE)
  pedstats$imp <- apply(pedchk$ped_arr[indix, "imp", hbix, drop=FALSE],
                        MARGIN=3, FUN=sum, na.rm=TRUE)
  pedstats$hap <- apply(pedchk$ped_arr[indix, "hap", hbix, drop=FALSE],
                        MARGIN=3, FUN=sum, na.rm=TRUE)
  pedstats$match.NA <- apply(is.na(pedchk$ped_arr[indix, "noDR", hbix, drop=FALSE]),
                             MARGIN=3, FUN=sum, na.rm=TRUE)
  pedstats$noDR.TRUE <- apply(pedchk$ped_arr[indix, "noDR", hbix, drop=FALSE],
                              MARGIN=3, FUN=sum, na.rm=TRUE)
  pedstats$noDR.FALSE <- apply(!(pedchk$ped_arr[indix, "noDR", hbix, drop=FALSE]),
                               MARGIN=3, FUN=sum, na.rm=TRUE)
  pedstats$withDR.TRUE <- apply(pedchk$ped_arr[indix, "withDR", hbix, drop=FALSE],
                                MARGIN=3, FUN=sum, na.rm=TRUE)
  pedstats$withDR.FALSE <- apply(!(pedchk$ped_arr[indix, "withDR", hbix, drop=FALSE]),
                                 MARGIN=3, FUN=sum, na.rm=TRUE)
  pedstats
}

calcFSstats <- function(ovwFS, haploblocks) {
  # function used only in calcStatistics; for input and output see there
  if (missing(haploblocks) || is.null(haploblocks)) {
    haploblocks <- rownames(ovwFS$ovw)
  } else if (is.numeric(haploblocks)) {
    haploblocks <- rownames(ovwFS$ovw)[haploblocks]
  } else if (is.list(haploblocks)) {
    haploblocks <- names(haploblocks)
  }
  if (!all(haploblocks %in% rownames(ovwFS$ovw)))
    stop("not all haploblocks occur in ovwFS")
  hbix <- match(haploblocks, rownames(ovwFS$ovw))
  if (length(hbix) == 0 || anyNA(hbix))
    stop("no selected haploblocks present in ovwFS")
  startcols <- which(substr(colnames(ovwFS$ovw), 1, 6) == "parmrk")
  FSstats <- data.frame(FSfam = 1:length(startcols))
  FSstats$bothparmrk <- colSums(ovwFS$ovw[hbix, startcols, drop=FALSE] == 2) #na.rm should not be needed!
  FSstats$FSFit <- colSums(ovwFS$ovw[hbix, startcols+1, drop=FALSE] == 1)
  FSstats$mrk <- colSums(ovwFS$ovw[hbix, startcols+3, drop=FALSE]) / length(hbix)
  FSstats$imp <- colSums(ovwFS$ovw[hbix, startcols+4, drop=FALSE]) / length(hbix)
  FSstats$hap <- colSums(ovwFS$ovw[hbix, startcols+5, drop=FALSE]) / length(hbix)
  startcols <- which(colnames(ovwFS$ovw) %in% c("mrkrest", "mrkall"))
  reststats <-
    data.frame(FSfam=c("rest", "all"), bothparmrk=c(NA, NA),
               FSFit=c(NA, NA),
               mrk=colSums(ovwFS$ovw[hbix, startcols, drop=FALSE]) / length(hbix),
               imp=c(0, sum(FSstats$imp)), #only imputation in FS families
               hap=c(sum(ovwFS$ovw[hbix, startcols[1]+1]),
                     sum(ovwFS$ovw[hbix, startcols[2]+2])) / length(hbix)
    )
  FSstats <- rbind(FSstats, reststats)
  rownames(FSstats) <- NULL
  list(FSstats=FSstats, haploblocks=haploblocks)
}

#'@title calculate some statistics of the solutions of all haploblocks
#'@description calculate some statistics of the solutions of all haploblocks
#'@usage calcStatistics(pedchk, ovwFS, indiv, haploblocks)
#'@param pedchk a list as generated by pedigreeHapCheck
#'@param ovwFS a list as produced by overviewFS, or NA if no FS families
#'available. If ovwFS is specified, pedchk and ovwFS should have
#'been calculated from the same underlying data
#'@param indiv missing, NULL or a vector of names of individuals. If present,
#'the component pedstats is calculated only over the individuals in the vector,
#'else over all individuals. Component FSstats is not affected by indiv
#'@param haploblocks missing, NULL or a vector of names or numbers of
#'haploblocks, or a list of haploblocks as in the inferHaplotypes input.
#'If present, components pedstats and FSstats are calculated only over the
#'specified haploblocks, else over all haploblocks.
#'@return a list with one or two items:\cr
#'pedstats: a data.frame with one row per haploblock, showing numbers
#'of individuals in columns\cr
#'- mrk: individuals with complete marker data\cr
#'- hap: individuals with an assigned haplotype combination\cr
#'- match.NA: individuals where a match between the haplotypes of the individual
#'  and its parents could not be verified (either the individual itself or both
#'  parents have no haplotype combination assigned)\cr
#'- noDR.TRUE: individuals whose haplotype combination matches that of its
#'  parents, assuming polysomic inheritance without Double Reduction\cr
#'- noDR.FALSE: individuals whose haplotype combination doesn't match that of its
#'  parents, assuming polysomic inheritance without Double Reduction\cr
#'- withDR.TRUE: individuals whose haplotype combination matches that of its
#'  parents, assuming polysomic inheritance, allowing Double Reduction\cr
#'- withDR.FALSE: individuals whose haplotype combination doesn't match that of
#'  its parents, assuming polysomic inheritance, allowing Double Reduction\cr
#'The total number of individuals is \cr
#'match.NA + noDR.TRUE + noDR.FALSE == match.NA + withDR.TRUE + withDR.FALSE;\cr
#'the number of individuals with missing marker data or without assigned
#'haplotype dosages can be calculated by subtracting mrk or hap from that
#'total.\cr
#'FSstats (only if ovwFS specified): a data.frame with for each FS family,
#'"rest" and "all" one row, with columns\cr
#'- pop: the FS family numbers, "rest" (all individuals not belonging to the
#'  FSs or their parents) and "all" (all individuals)\cr
#'- bothparmrk: for how many haploblocks have both parents complete marker data?\cr
#'- FSdone: for how many haploblocks has a solution based on polysomic
#'  inheritance been found for this FS family?\cr
#'- mrk: the average number of (FS) individuals with complete marker data\cr
#'- imp: the average number of (FS) individuals with imputed marker data\cr
#'- hap: the average number of (FS) individuals with an assigned haplotype
#'  combination\cr
#'The numbers for "all" include also the FS parents (each parent counted only
#'once, even if some parents produced more than one FS family), which are not
#'included in any of the other rows; therefore "all" is not equal to the sum of
#'the rows above.
#'@examples
#'data(PolyHaplotyper_small)
#'phchk <- pedigreeHapCheck(ped=phped, mrkDosage=phdos, haploblock=phblocks,
#'                           hapresults=phresults)
#'phovw <- overviewByFS(haploblock=phblocks, parents=phpar, FS=phFS,
#'                      hapresults=phresults)
#'calcStatistics(pedchk=phchk, ovwFS=phovw)
#'@export
calcStatistics <- function(pedchk, ovwFS, indiv, haploblocks) {
  if (missing(indiv)) indiv <- NULL
  if (missing(haploblocks)) haploblocks <- NULL
  if (missing(pedchk)) pedchk <- NULL else
  if (!is.list(pedchk) || !identical(names(pedchk), c("ped_arr", "parents_arr")))
    stop ("invalid pedchk")
  if (missing(ovwFS)) ovwFS <- NULL else
    if (!is.list(ovwFS) || !identical(names(ovwFS), c("ovw", "messages")))
      stop("invalid ovwFS")
  if (is.null(pedchk) && is.null(ovwFS))
    stop("no data")

  result <- list()
  if (!is.null(pedchk))
    result$pedstats <- calcPedstats(pedchk, indiv, haploblocks)
  if (!is.null(ovwFS)) {
    FSstats <- calcFSstats(ovwFS, haploblocks)
    result$FSstats <- FSstats$FSstats
  }
  if (!is.null(pedchk) && !is.null(ovwFS) &&
      !setequal(result$pedstats$haploblock, FSstats$haploblocks)) {
    cat("The haploblocks in pedchk and ovwFS do not match\n"); flush.console()
  }
  result
} #calcStatistics

#'@title produce a table of nr of markers vs nr of haplotypes
#'@description produce a table of nr of markers vs nr of haplotypes
#'@usage calcMrkHaptable(ovwFS)
#'@param ovwFS a list as produced by overviewFS. Only the first 2 columns of
#'list item ovw are used
#'@return a frequency table with the numbers of haploblocks, with
#'all combinations of marker counts and inferred haplotype counts per haploblock.
#'The column with haplotype count NA (if any) shows the haploblocks for which
#'no haplotype solution was found (the reason for that would usually be found
#'in column 1 of ovwFS$messages)
#'@examples
#'data(PolyHaplotyper_small)
#'phovw <- overviewByFS(haploblock=phblocks, parents=phpar, FS=phFS,
#'                      hapresults=phresults)
#'# in this small dataset there are only 2 haploblocks, each with 4 markers:
#'calcMrkHaptable(ovwFS=phovw)
#'# in both haploblocks 5 haplotypes are inferred
#'@export
calcMrkHaptable <- function(ovwFS){
  table(ovwFS$ovw[ ,1], ovwFS$ovw[, 2], useNA="ifany",
                     dnn=c("mrk", "haplotypes"))
} # calcMrkHaptable

#'@title convert a haploblock-defining data frame to a list
#'@description convert a haploblock-defining data frame to a list as needed
#'by inferHaplotypes
#'@usage haploblock_df2list(df, mrkcol, hbcol, sorted=TRUE)
#'@param df a data.frame with at least the following two columns
#'@param mrkcol the name or number of the column with the markers
#'@param hbcol the name or number of the column with the haploblocks
#'@param sorted if TRUE (default) the haploblock list is sorted in alphabetical
#'order of haploblock name (the markers within haploblocks are not sorted);
#'if FALSE the haploblocks will be in order of first occurrence in df
#'@details function inferHaplotypes needs a list where each item is a
#'vector of all markers in one haploblock. This function produces such a list
#'from a data.frame where the markers are in one column and the haploblocks
#'in another column. The markers and haploblocks columns may be character,
#'factor or numeric, and the columns may be indicated by name or number.
#'@return the desired list
#'@examples
#'df1 <- data.frame(
#'  marker=paste0("mrk",1:9),
#'  block=LETTERS[c(1,2,3,2,3,1,1,2,2)],
#'  extracol=runif(9))
#'haploblock_df2list(df=df1, mrkcol="marker", hbcol=2)
#'@export
haploblock_df2list <- function(df, mrkcol, hbcol, sorted=TRUE) {
  uhb <- unique(df[, hbcol])
  if (is.factor(df[, mrkcol])) df[, mrkcol] <- as.character(df[, mrkcol])
  lst <- vector(length(uhb), mode="list")
  for (hb in seq_along(uhb)) {
    lst[[hb]] <- df[df[, hbcol]==uhb[hb], mrkcol]
  }
  names(lst) <- as.character(uhb)
  if (sorted) lst <- lst[order(names(lst))]
  lst
  # split(df[,mrkcol], df[,hbcol]) would also work but:
  # - the haploblocks might be reordered if they are not alphabetically sorted
  # - if df[,mrkcol] is a factor, the returned list will also consist of
  #   factors rather than character vectors
} # haploblock_df2list

# .........................................................
# Functions to compare with other haplotyping software ####
# .........................................................

#'@title convert PolyHaplotyper marker data to ShesisPlus format
#'@description convert PolyHaplotyper marker data (input) to ShesisPlus
#'format for a single haploblock.
#'@usage make.ShesisPlus.input(mrkDosage, indiv=NULL, markers,
#'ploidy, phenotype=0, fname="")
#'@param mrkDosage matrix or data.frame of allele dosages; same as input for
#'inferHaplotypes. Markers are in rows, individuals in columns, each cell has
#'a marker dosage. All marker dosages must be in 0:ploidy or NA
#'@param indiv the names of the individuals to include in the ShesisPlus
#'input data. Default NULL includes all individuals
#'@param markers character vector with the names of the markers in the
#'haploblock; all must occur in mrkDosage
#'@param ploidy single integer: the ploidy level
#'@param phenotype vector with the phenotypes of all individuals, in order
#'of the columns of dosmat; default 0
#'@param fname filename of the output file: this will contain the data
#'in the ShesisPlus format (the saved data.frame is also the return value).
#'If "" (default) no file is written
#'@details ShesisPlus needs the data formatted as: 1 column with names of
#'individuals, 1 column with the phenotypes (with missing values represented
#'as "NA"), and for each of the selected markers <ploidy> columns,
#'with the (unphased) marker alleles. Here we use only biallelic markers;
#'their alleles are indicated by 1 and 2 for the ref and alt allele,
#'and 0 for missing data. \cr
#'The contents of file fname can be pasted into the input data box of
#'the ShesisPlus web interface at http://shesisplus.bio-x.cn/SHEsis.html
#'@return a data.frame in the described ShesisPlus input format
#'@examples
#'data(PolyHaplotyper_small)
#'SSPin <- make.ShesisPlus.input(mrkDosage=phdos, markers=phblocks[[1]],
#'                               ploidy=6)
#'SSPin[1:6,1:8]
#'@export
make.ShesisPlus.input <- function (mrkDosage, indiv=NULL, markers,
                                   ploidy, phenotype=0, fname="") {
  dm <- checkmrkDosage(mrkDosage, ploidy=ploidy, markers=markers) # also converts to matrix
  if (length(phenotype) == 1) phenotype <- rep(phenotype, ncol(dm))
  if (length(phenotype) != ncol(dm)) stop("Incorrect number of phenotypes")
  if (anyNA(phenotype)) phenotype[is.na(phenotype)] <- "NA"
  result <- matrix(numeric(0), nrow=ncol(dm), ncol=ploidy*nrow(dm),
                   dimnames=list(colnames(dm), rep(rownames(dm), each=ploidy)))
  for (m in seq_along(markers)) {
    # in order to use the hapdos2hapcomb function of PolyHaplotyper
    # we need two rows: one for the 1 and one for the 2 allele:
    mdat <- rbind(ploidy-dm[m,], dm[m,])
    result[, (ploidy*(m-1)+1):(ploidy*m)] <- t(hapdos2hapcomb(mdat, ploidy))
  }
  result[is.na(result)] <- 0 # ShesisPlus code for missing allele
  result <- cbind(phenotype, result)
  if (fname != "") write.table(result, file=fname, row.names=TRUE, col.names=FALSE,
                               quote=FALSE, sep=" ")
  result
}

#'@title Read the haplotyping results from the ShesisPlus output
#'@description Read the haplotyping results from the ShesisPlus output
#'@usage read.ShesisPlus.output(SSPout, order.by=c("", "hapnr", "count")[2])
#'@param SSPout filename of a text file with ShesisPlus output (as copied from
#'the web page produced by running ShesisPlus web interface
#'at http://shesisplus.bio-x.cn/SHEsis.html); or a character vector with the
#'same output
#'@param order.by how to order the rows of the hapstat data.frame. "hapnr" means
#'ordering by haplotype number, "count" means ordering by decreasing Total.count,
#'anything else results in no reordering. By default order by hapnr.
#'@details If present, the markernames and the haplotype statistics are read from
#'the file. ShesisPlus does not provide haplotype combinations for individuals
#'@return a list with 2 items: $markernames has the marker names in the markers in
#'the haploblock (if present in the file); hapstat contains the haplotype statistics
#'as read from the file with an additional (first) column hapnr: the haplotype
#'numbers as defined in PolyHaplotyper. The other columns are Haplotype (as sequences
#'of marker alleles), Total.count (of the haplotype over the whole population) and,
#'according to an email from Zhiqiang Li of 13-04-2020) BETA: Regression coefficient,
#'SE: Standard error, R2: Regression r-squared, T: t-distribution statistics,
#'P: p-value
#'@examples
#'# we give a typical SSP output as character vector; instead we could also
#'# give the name of a text file
#'SSPout <- c(
#'  " Please cite:",
#'  "",
#'  "    Shen, J. et al. SHEsisPlus, a toolset for genetic studies ...",
#'  "    Shi, Y. et al. SHEsis, a powerful software platform ...",
#'  "    Li, Z. et al. A partition-ligation-combination-subdivision ...",
#'  "",
#'  "if you find this tool useful in your research. Thanks!",
#'  "",
#'  "Haplotype Analysis:",
#'  "",
#'  "Haplotypes with frequency <0.03 are ignored.",
#'  "Loci chosen for haplotype analysis: m1, m2, m3, m4",
#'  "Haplotype 	Total count 	Beta 	SE 	R2 	T 	p",
#'  "1122 	2232 	0.002 	0.01 	1.01e-04 	0.252 	0.8",
#'  "1222 	230 	-0.019 	0.017 	0.002 	-1.15 	0.25",
#'  "1221 	152 	-0.01 	0.024 	2.94e-04 	-0.43 	0.667",
#'  "2222 	288 	0.008 	0.02 	2.93e-04 	0.429 	0.667",
#'  "1121 	142 	-0.009 	0.022 	2.81e-04 	-0.42 	0.674"
#')
#'read.ShesisPlus.output(SSPout)
#'@export
read.ShesisPlus.output <- function(SSPout, order.by=c("", "hapnr", "count")[2]) {
  if (missing(order.by) || is.null(order.by) || is.na(order.by)) order.by <- ""
  order.by <- order.by[1]
  fromfile <- length(SSPout)==1 && file.exists(SSPout)
  if (!fromfile && !(is.character(SSPout) && length(SSPout)>1))
    stop("SSPout must be a text file name or a character vector")
  if (fromfile) con <- file(SSPout) else con=textConnection(SSPout)
  txt <- readLines(con)
  i <- grep("Loci chosen for haplotype analysis:", txt, value=FALSE)
  markernames <- character(0)
  if (length(i) == 1) {
    # read marker names from line, after ": " and separated by ", "
    colon <- regexpr(": ", txt[i])
    commas <- gregexpr(", ", txt[i])[[1]]
    commas <- c(colon, commas, nchar(txt)[i]+1)
    markernames <- substring(txt[i], commas[1:(length(commas)-1)]+2,
                             commas[2:length(commas)]-1)
  }
  # next, get the haplotype statistics:
  i <- grep("Haplotype \tTotal count \tBeta \tSE \tR2 \tT \tp", txt, value=FALSE)
  if (length(i) != 1) stop("header line of haplotype statistics not found")
  # find out if there is an LD section after the statistics:
  LD <- grep("Linkage Disequilibrium Analysis:", txt, value=FALSE)
  if (length(LD) > 0) nrows <- LD[1] - i - 1 else nrows <- -1
  if (fromfile) con <- file(SSPout) else con <- textConnection(SSPout)
  hapstat <- read.table(con, header=TRUE, sep="\t", skip=i-1, nrows=nrows,
                        na.strings=c("NA", "-NA", "nan", "-nan", ""))
  close(con)
  if (order.by == "count") hapstat <- hapstat[order(hapstat$Total.count, decreasing=TRUE),]
  # check length of haplotypes:
  nc <- nchar(hapstat$Haplotype)
  if (min(nc) != max(nc)) {
    warning("Haplotypes have different lengths")
    hapstat <- cbind(rep(NA, nrow(hapstat)), hapstat)
    } else {
      if (length(markernames) == 0) {
        markernames <- paste0("m", seq_len(min(nc)))
      } else if (length(markernames) != min(nc)) {
        warning("markernames don't match length of haplotypes")
      }
      # translate haplotypes to haplotype numbers
      allhap <- apply(allHaplotypes(character(min(nc))) + 1, # in SSP we use alleles 1&2, not 0&1
                      MARGIN=1, FUN=paste, collapse="")
      hapstat <- cbind(match(hapstat$Haplotype, allhap), hapstat)
      if (order.by == "hapnr") hapstat <- hapstat[order(hapstat[,1]),]
    }
  names(hapstat)[1] <- "hapnr"
  list(hapstat=hapstat, markernames=markernames)
} #read.ShesisPlus.output

#'@title convert PolyHaplotyper marker data to SATlotyper format
#'@description convert PolyHaplotyper marker data (input) to SATlotyper
#'format for a single haploblock.
#'@usage make.SATlotyper.input(mrkDosage, indiv=NULL, markers,
#'ploidy, phenotype=0, fname)
#'@param mrkDosage matrix or data.frame of allele dosages; same as input for
#'inferHaplotypes. Markers are in rows, individuals in columns, each cell has
#'a marker dosage. All marker dosages must be in 0:ploidy or NA
#'@param indiv the names of the individuals to include in the ShesisPlus
#'input data. Default NULL includes all individuals
#'@param markers character vector with the names of the markers in the
#'haploblock; all must occur in mrkDosage
#'@param ploidy single integer: the ploidy level
#'@param phenotype vector with the phenotypes of all individuals, in order
#'of the columns of dosmat; default 0
#'@param fname filename of a *.csv output file: this will contain the data
#'in the SATlotyper format (the saved data.frame is also the return value).
#'If "" no file is written
#'@return a data.frame in the SATlotyper input format: a header row with
#'"Genotype" and the marker names, and one row per individual with the
#'individual name plus for each marker the genotype as a sorted
#'string of <ploidy> A's and B's, or <ploidy> N's
#'@examples
#'data(PolyHaplotyper_small)
#'SATin <- make.SATlotyper.input(mrkDosage=phdos, markers=phblocks[[1]],
#'                               ploidy=6, fname="")
#'head(SATin)
#'@export
make.SATlotyper.input <- function (mrkDosage, indiv=NULL, markers,
                                   ploidy, phenotype=0, fname) {
  sep <- "\t" # separator of the columns in the output file
  dm <- checkmrkDosage(mrkDosage, ploidy=ploidy, markers=markers) # also converts to matrix
  result <- matrix(NA_character_, nrow=ncol(dm), ncol=nrow(dm),
                   dimnames=list(colnames(dm), rownames(dm)))
  dos2str <- function(dosage, ploidy) {
    # e.g. (ploidy=4): dosage 3 -> ABBB, dosage 0 -> AAAA, dosage NA -> NNNN
    alleles <- c("A", "B")
    miss <- "N"
    if (is.na(dosage)) paste(rep(miss, ploidy), collapse="") else {
      str0 <- paste(rep(alleles[1], ploidy-dosage), collapse="")
      str1 <- paste(rep(alleles[2], dosage), collapse="")
      paste0(str0, str1)
    }
  }
  for (i in seq_len(ncol(result))) {
    result[, i] <- sapply(dm[i,], FUN=dos2str, ploidy=ploidy)
  }
  if (length(fname)==1 && !is.na(fname) && fname != "") {
    con <- file(fname)
    writeLines(paste("Genotype", paste(colnames(result), collapse=sep), sep=sep),
               con=con, sep="\n")
    close(con)
    write.table(result, file=fname, col.names=FALSE, row.names=TRUE, quote=FALSE,
                sep=sep, append=TRUE)
  }
  result
} #make.SATlotyper.input

#'@title read the haplotyping results from the SATlotyper output
#'@description read the haplotyping results from the SATlotyper output
#'@usage read.SATlotyper.output(fname, output="", allelecodes=c("A", "B"),
#'sep, haploblockname="")
#'@param fname filename of an xml file produced by SATlotyper
#'@param output character vector, console output of SATlotyper;
#'this contains a table of calculated haplotypes with more info
#'than the xml file
#'@param allelecodes the codes used for the ref and alt SNP alleles in
#'make.SATlotyper.input, default "A" and "B"
#'@param sep the separator used in make.SATlotyper.input (SATlotyper cleverly
#'identifies this separator and uses it for its output), default tab
#'@param haploblockname if not "" (default) the haplotype IDs are given as
#'<haploblockname>_<hapnr>, else just as <hapnr>, left_padded with zeroes
#'@details The xml file is parsed using the package XML. The resulting list has
#'3 items named source, bootstrapping and haplotyping. In this function source is
#'ignored and the results are obtained from haplotypings[[1]] (i.e. even if multiple
#'haplotypings were done only the first is extracted), with the haplotype
#'score from the bootstrapping item added to the haplotype.info.
#'@return a list with 3 items: \cr
#'$hapdos is a matrix with individuals in columns and haplotypes in rows,
#'giving the dosages of the haplotypes in each individual (summing to ploidy).
#'This is the same format as the hapdos components of the inferHaplotypes
#'results of PolyHaplotyper \cr
#'$haplotype.info is a data.frame that is a combination of information from two
#'or three tables:\cr
#'columns Haplotype and necessity are from element
#'haplotypings[[1]]$haplotypes in the xml file, \cr
#'column hapnr gives the PolyHaplotyper equivalents of the Haplotypes, \cr
#'column score is from element bootstrapping in the xml file, \cr
#'columns id, number, frequency, homozygous, necessary, distance and
#'neighbours are from (the console) output if present. \cr
#'The haplotype.info data.frame is ordered by hapnr. \cr
#'$HaplotypingScore is the single number in $haplotyping[[1]]$score.
#'
#'@export
read.SATlotyper.output <- function(fname, output="",
                                   allelecodes=c("A", "B"), sep,
                                   haploblockname="") {
  if (missing(sep)) sep <- "\t" #separator in input file
  getAllhapStrings <- function(nmrk, allelecodes){
    result <- character(2^nmrk)
    allhap <- allHaplotypes(seq_len(nmrk))
    result <- allhap #dim correct
    result[allhap==0] <- allelecodes[1]
    result[allhap==1] <- allelecodes[2]
    apply(result, MARGIN=1, paste, collapse="")
  } #getAllhapStrings in read.SATlotyper.result

  xpl <- XML::xmlToList(XML::xmlParse(file=fname))
  markernames <- strsplit(xpl$source[[1]], sep)[[1]][-1]
  names(markernames) <- NULL # all names are "genotype"
  nmrk <- length(markernames)
  allHapstrings <- getAllhapStrings(nmrk, allelecodes)
  HaplotypingScore <- xpl$haplotypings$haplotyping$score
  if (!anyNA(as.numeric(HaplotypingScore)))
    HaplotypingScore <- as.numeric(HaplotypingScore)
  usedhap2 <- t(sapply(xpl$haplotypings[[1]]$haplotypes, FUN=function(x) unlist(x)))
  usedhaplo <- data.frame(Haplotype=usedhap2[,1])
  usedhaplo$hapnr <- match(usedhaplo$Haplotype, allHapstrings)
  usedhaplo$necessity=usedhap2[,2] != "false"
  usedhaplo <- usedhaplo[order(usedhaplo$hapnr),]
  bootstrapping <- t(sapply(xpl$bootstrapping, FUN=function(x) unlist(x)))
  if (anyNA(as.numeric(bootstrapping[, 2]))) score <- bootstrapping[, 2] else
    score <- as.numeric(bootstrapping[, 2])
  usedhaplo$score <- score[match(bootstrapping[,1], usedhaplo$Haplotype)]
  hapnames <- padded(usedhaplo$hapnr, max(usedhaplo$hapnr))
  if (haploblockname[1] != "")  hapnames <- paste(haploblockname[1], hapnames, sep="_")
  strgeno <- sapply(xpl$haplotypings[[1]]$genotypes, FUN=function(x) unlist(x))
  hapgeno <- matrix(match(strgeno[-1,], allHapstrings), ncol=ncol(strgeno))
  if (any(is.na(hapgeno) & !is.na(strgeno[-1,])))
    stop("some assigned haplotypes not among used haplotypes")
  if (length(setdiff(usedhaplo$hapnr, unique(hapgeno))) > 0)
    stop("some used haplotypes not among assigned haplotypes")
  hapdos <- matrix(0, nrow=nrow(usedhaplo), ncol=ncol(hapgeno),
                   dimnames=list(hapnames, strgeno[1,]))
  for (i in seq_len(nrow(usedhaplo))) {
    hn <- usedhaplo$hapnr[i]
    hapdos[i,] <- apply(hapgeno, MARGIN=2, FUN=function(x) sum(x==hn, na.rm=TRUE))
  }
  names(colnames(hapdos)) <- NULL
  # finally we check if the table with extended haplotype info is
  # available from the console output:
  header <- grep("id | sequence |", output, fixed=TRUE)
  end <- grep("calculated phased genotypes", output, fixed=TRUE)
  if (length(header) == 1 && length(end) == 1) {
    if (end - header != nrow(usedhaplo) + 3) {
      warning("table calculated.haplotypes in output doesn't match haplotypes in xml file")
    } else {
      headers <- trimws(strsplit(output[header], " | ", fixed=TRUE)[[1]])
      valu <- trimws(t(simplify2array(strsplit(output[(header+2):(end-2)], " | ", fixed=TRUE))))
      # we assume here that the columns are always the same:
      if (length(headers) != ncol(valu) || length(headers) != 8) {
        warning("in output table calculated.genotypes has unexpected format")
      } else {
        numconv <- function(x)
          if (anyNA(suppressWarnings(as.numeric(x)))) x else as.numeric(x)
        xhapinfo <- data.frame(
          id=numconv(valu[,1]),
          sequence=valu[,2],
          number=numconv(valu[,3]),
          frequency=numconv(valu[,4]),
          homozygous=valu[,5]=="true",
          necessary=valu[,6]=="true",
          distance=numconv(valu[,7]),
          neighbours=numconv(valu[,8])
        )
        names(xhapinfo) <- headers # just in case they are different
        o <- match(xhapinfo$sequence, usedhaplo$Haplotype)
        if (anyNA(o)) {
          warning("table calculated.haplotypes in output doesn't match haplotypes in xml file")
        } else {
          usedhaplo <- cbind(usedhaplo, xhapinfo[o, -2])
        }
      }

    }
  }
  list(hapdos=hapdos, haplotype.info=usedhaplo, HaplotypingScore=HaplotypingScore,
       markers=markernames)
} #read.SATlotyper.output

#'@title A simple interface to run SATlotyper
#'@description A simple interface to run SATlotyper
#'@usage run.SATlotyper(path_to_SATlotyper, infile, outfile,
#'SAT_solver="sat4j.conf")
#'@param path_to_SATlotyper path to the folder where SATlotyper.jar and
#'the *.conf files of the SAT-solvers are located
#'@param infile name (and path) to the SATlotyper input file
#'@param outfile name (and path) for the SATlotyper output file (an xml file)
#'@param SAT_solver name of the *.conf file for the SAT-solver to use.
#'Default "sat4j.conf" because this works under both Windows and Linux.
#'@details This function issues a system command to invoke SATlotyper.
#'Java and SATlotyper must be installed. This is just a simple interface for
#'convenience; for more control run SATlotyper directly from the command window.
#'@return The return is a list with elements:\cr
#'$cmd : the command passed to the system\cr
#'$result: screen output from SATlotyper; this
#'contains some extra info not present in the output xml file)\cr
#'The main result is the outfile
#'@export
run.SATlotyper <- function(path_to_SATlotyper, infile, outfile,
                           SAT_solver="sat4j.conf") {
  cmd <- paste0('java -jar "', file.path(path_to_SATlotyper, 'SATlotyper.jar"'),
                ' -all -i "', infile, '" -o "', outfile, '"',
                ' -sat "', file.path(path_to_SATlotyper, SAT_solver), '"')
  # in R.Version 4.X this will be improved using "raw strings" so that
  # no \" sequences appear which prevent copy-pasting the cmd into
  # a Windows system prompt window
  args <- c(paste0('-jar "', file.path(path_to_SATlotyper, 'SATlotyper.jar"'),
            ' -all',
            paste0(' -i "', infile, '"'),
            paste0(' -o "', outfile, '"'),
            paste0(' -sat "', file.path(path_to_SATlotyper, SAT_solver), '"')))
  result <- system2(command="java", args=args, stdout=TRUE, stderr=TRUE)
  invisible(c(cmd, result))
}

#'@title convert PolyHaplotyper input data to Happy-inf format
#'@description convert PolyHaplotyper input data to Happy-inf format
#'for a single haploblock
#'@usage make.Happyinf.input(mrkDosage, indiv=NULL, haploblock,
#'ploidy, fname)
#'@param mrkDosage matrix or data.frame of allele dosages; same as input for
#'inferHaplotypes. Markers are in rows, individuals in columns, each cell has
#'a marker dosage. All marker dosages must be in 0:ploidy or NA
#'@param indiv the names of the individuals to include in the Happy-inf
#'input data. Default NULL includes all individuals
#'@param haploblock a list of character vectors. The names are the names of the
#'haploblocks, the character vectors have the names of the markers in each
#'haploblock. Only the markers in haploblock will be included in the
#'Happy-inf input data
#'@param ploidy single integer: the ploidy level
#'@param fname filename of a tab-separated output file: this will contain the
#'data in Happy-inf format (the saved data.frame is also the return value).
#'If "" no file is written
#'@return a data.frame in the Happy-inf input format: a header row with
#'"SNPID", "block" and names of the individuals and one row per marker, with
#'only the individuals, markers and blocks as specified.
#'#'SNPID has the marker names, "Block" the haploblock names. All markers
#'of a haploblock are in contiguous rows. Missing dosages
#'are represented by "NA"
#'@examples
#'data(PolyHaplotyper_small)
#'HAPin <- make.Happyinf.input(mrkDosage=phdos, haploblock=phblocks,
#'                             ploidy=6, fname="")
#'HAPin[,1:8]
#'@export
make.Happyinf.input <- function (mrkDosage, indiv=NULL, haploblock,
                                   ploidy, fname) {
  result <- checkmrkDosage(mrkDosage, ploidy=ploidy, indiv=indiv,
                           markers=unlist(haploblock))
  result <- data.frame(SNPID=rownames(result),
                       Block=rep(names(haploblock), times=sapply(haploblock, length)),
                       result,
                       check.names=FALSE,
                       stringsAsFactors=FALSE)
  rownames(result) <- NULL
  if (length(fname)==1 && !is.na(fname) && fname != "") {
    write.table(result, file=fname,
                row.names=FALSE, col.names=TRUE,
                sep="\t", na="NA", quote=FALSE)
  }
  invisible(result)
} #make.Happyinf.input

#'@title read the haplotyping results from the Happy-inf output
#'@description read the haplotyping results from the Happy-inf output
#'@usage read.Happyinf.output(file_prefix, dropUnused=TRUE)
#'@param file_prefix the prefix for the output files generated by
#'Happy-inf (the value of the -o parameter)
#'@param dropUnused TRUE (default) if the returned matrix should only contain
#'rows for haplotypes that are present; if FALSE matrix contains rows for all
#'possible haplotypes
#'@details This function reads the <file_prefix>.stat.dat and
#'<file_prefix>.stats.dat files. The first contains one row per SNP marker,
#'grouped per haploblock and <ploidy> columns per individual with the SNP
#'haplotypes, where SNP alleles are represented as 0 or 1, or NA for unknown.
#'For one SNP, the <ploidy> alleles in an individual are all known or all
#'unknown It is possible that of the different SNPs in a haploblock
#'some are known and some are not; in the conversion to PolyHaplotyper format
#'all these partially known haplotypes are made unknown.
#'The latter file has a statistics "mismatch" and "ratio" for each individual /
#'haploblock combination; for their meaning see the Happyinf readme file.
#'@return a list with 2 items: \cr
#'$hapdos is a matrix with individuals in columns and haplotypes in rows,
#'giving the dosages of the haplotypes in each individual (summing to ploidy).
#'This is the same format as the hapdos components of the inferHaplotypes
#'results of PolyHaplotyper except that \cr
#'$stats is a matrix with individuals in rows (matching dosmat) and two rows
#'named mismatch and ratio.
#'
#'@export
read.Happyinf.output <- function(file_prefix, dropUnused=TRUE) {
  dat <- read.table(paste0(file_prefix, ".haplotypes.dat"), header=TRUE,
                    sep="\t", stringsAsFactors=FALSE, check.names=FALSE,
                    na.strings="NA")
  stat <- read.table(paste0(file_prefix, ".stat.dat"), header=TRUE,
                     sep="\t", stringsAsFactors=FALSE, check.names=FALSE,
                     na.strings=c("NA", "nan"))
  block <- unique(dat$Block)
  indiv <- names(stat)[-(1:2)]
  ploidy <- (ncol(dat)-2) / length(indiv)
  if (names(dat)[1]!="Block" || names(stat)[1]!="Block" ||
      !identical(block, unique(stat$Block)) ||
      !identical(indiv, names(stat)[-(1:2)]))
    stop("haplotypes and stat files don't match")
  results <- vector(length=length(block), mode="list")
  names(results) <- block
  for (b in seq_along(results)) {
    # first hapdos:
    brows <- dat$Block==block[b]
    markers <- dat$SNP_ID[brows]
    nhap <- 2^length(markers)
    datb <- as.matrix(dat[brows, -(1:2)])
    # convert haplotypes to haplotype numbers:
    NAcols <- is.na(colSums(datb))
    datb[, NAcols] <- NA
    print(paste(b, block[b])); flush.console()
    hapnrs <- mrkdos2mrkdid(mrkDosage=datb, indiv=NULL, ploidy=1, check=FALSE)
    #hapnrs is now a vector of hapnrs, of length ploidy * length(indiv)
    hapdos <- hapcomb2hapdos(hapcomb=matrix(hapnrs, nrow=ploidy, byrow=FALSE),
                             nhap=nhap)
    rownames(hapdos) <- paste0(block[b], "_", padded(1:nhap))
    colnames(hapdos) <- indiv
    if (dropUnused) {
      hapdos <-
        hapdos[rowSums(hapdos, na.rm=TRUE) > 0, , drop=FALSE]
    }
    results[[b]] <- list(hapdos=hapdos, markers=markers)
    # then stats:
    brows <- stat$Block==block[b]
    results[[b]]$stats <- as.matrix(stat[brows, -(1:2)])
    rownames(results[[b]]$stats) <- stat$type[brows] # mismatch and ratio
  }
  results
} # read.Happyinf.output

#'@title compare two haplotyping results
#'@description compare two haplotyping results, e.g. PolyHaplotyper and
#'SATlotyper
#'@usage compareHapresults(haploblock, hapresultsA, hapresultsB)
#'@param haploblock a list of character vectors. The names are the names of the
#'haploblocks, the character vectors have the names of the markers in each
#'haploblock.
#'@param hapresultsA and
#'@param hapresultsB two list as returned by inferHaplotypes,
#'with one item (itself a list) per haploblock with at least a matrix hapdos
#'and a character vector markers. All
#'haploblocks in param haploblock must occur in hapresultsA and in hapresultsB.
#'The individual names (colnames of the hapdos items for each haploblock)
#'must be identical and in the same order in hapresultsA and hapresultsB
#'@return a list with one element per haploblock in param haploblock.
#'Each element is itself a list with elements:
#'$identical: TRUE or FALSE\cr
#'$message: a single string, "" if the comparison is possible, else the reason
#'why not (if $message is not "", $identical is always FALSE). The next
#'elements are only present if $message is "":\cr
#'$compindiv: a matrix comparing the two hapdos, with one column per individual
#'and 5 rows: Both_NA, A_NA, B_NA, Equal, Uneq. The last 2 have NA values if A
#'and/or B is NA
#'$haplofreq: a matrix with one row per haplotype occurring in A and/or B,
#'and columns A and B, with the total frequency of each haplotype in
#'hapdos A or hapdos B
#'@export
compareHapresults <- function(haploblock, hapresultsA, hapresultsB) {
  if (!all(names(haploblock) %in% names(hapresultsA)) ||
      !all(names(haploblock) %in% names(hapresultsB)))
    stop("haploblocks in hapresultsA or hapresultsB don't match haploblock")
  if (ncol(hapresultsA[[1]]$hapdos) != ncol(hapresultsB[[1]]$hapdos) ||
      !all(colnames(hapresultsA[[1]]$hapdos) == colnames(hapresultsB[[1]]$hapdos)))
    stop("individuals in hapresultsA and hapresultsB not identical")
  compresults <- vector("list", length(haploblock))
  names(compresults) <- names(haploblock)
  for (hb in names(haploblock)) {
    compresults[[hb]]$identical <- FALSE
    if (!all(haploblock[[hb]] == hapresultsA$markers) ||
        !all(haploblock[[hb]] == hapresultsB$markers)) {
      compresults[[hb]]$message <-
        "markers in hapresultsA or hapresultsB different from haploblock"
      break
    }
    A <- hapresultsA[[hb]]$hapdos
    B <- hapresultsB[[hb]]$hapdos
    ploidy <- setdiff(unique(c(colSums(A, na.rm=TRUE), colSums(B, na.rm=TRUE))),
                      0)
    if (length(ploidy) != 1) {
      compresults[[hb]]$message <- "not all haplotype dosages sum to ploidy"
      break()
    }
    compresults[[hb]]$message <- ""
    A <- A[rowSums(A, na.rm=TRUE) > 0,, drop=FALSE]
    Ahap <- split_hapnames(rownames(A))$hapnrs
    B <- B[rowSums(B, na.rm=TRUE) > 0,, drop=FALSE]
    Bhap <- split_hapnames(rownames(B))$hapnrs
    # make Anw and Bnw as A and B, but containing all haplotypes in either
    x <- setdiff(Bhap, Ahap)
    xm <- matrix(rep(ploidy-colSums(A), length(x)), byrow=TRUE, ncol=ncol(A), nrow=length(x),
                 dimnames=list(rownames(B)[Bhap %in% x], colnames(A)))
    Anw <- rbind(A, xm)
    Anw <- Anw[order(c(Ahap, x)),,drop=FALSE]
    x <- setdiff(Ahap, Bhap)
    xm <- matrix(rep(ploidy-colSums(B), length(x)), byrow=TRUE, ncol=ncol(B), nrow=length(x),
                 dimnames=list(rownames(A)[Ahap %in% x], colnames(B)))
    Bnw <- rbind(B, xm)
    Bnw <- Bnw[order(c(Bhap, x)),,drop=FALSE]
    compresults[[hb]]$compIndiv <-
      matrix(NA, nrow=5, ncol=ncol(Anw),
             dimnames=list(c("Both_NA", "A_NA", "B_NA", "Equal", "Uneq"),
                           colnames(Anw)))
    compresults[[hb]]$compIndiv["Both_NA",] <- is.na(colSums(Anw)) & is.na(colSums(Bnw))
    compresults[[hb]]$compIndiv["A_NA",] <- is.na(colSums(Anw)) & !is.na(colSums(Bnw))
    compresults[[hb]]$compIndiv["B_NA",] <- !is.na(colSums(Anw)) & is.na(colSums(Bnw))
    compresults[[hb]]$compIndiv["Equal",] <- colSums(Anw != Bnw) == 0
    compresults[[hb]]$compIndiv["Uneq",] <- colSums(Anw != Bnw) > 0
    compresults[[hb]]$identical <-
      sum(compresults[[hb]]$compIndiv["Both_NA",]) +
      sum(compresults[[hb]]$compIndiv["Equal",], na.rm=TRUE) == ncol(Anw)
    compresults[[hb]]$haplofreq <-
      matrix(NA_integer_, nrow=nrow(Anw), ncol=2,
             dimnames=list(rownames(Anw), c("A", "B")))
    compresults[[hb]]$haplofreq[,"A"] <- rowSums(Anw, na.rm=TRUE)
    compresults[[hb]]$haplofreq[,"B"] <- rowSums(Bnw, na.rm=TRUE)
  } # for hb
  compresults
} # compareHapresults

#'@title compare haplotyping results with observed markers dosages
#'@description compare haplotyping results with observed markers dosages
#'@usage compareHapMrkDosages(mrkDosage, hapresults)
#'@param mrkDosage a data.frame or matrix with the marker dosages, in
#'the format of inferHaplotypes()
#'@param hapresults a list as returned by inferHaplotypes,
#'with one item (itself a list) per haploblock with at least a matrix hapdos
#'and a character vector markers. All markers in all haploblocks must also
#'occur in mrkDosages, and all individuals in the hapdos matrices must also
#'occur in mrkDosages
#'@return a 3-D array with dimensions haploblock, individual, and chkresult:
#'mrkNA (ALL markers in the haploblock have missing data for the individual),
#'hapNA (the haplotype dosages do not sum to ploidy and/or are missing),
#'match (TRUE if the non-missing marker dosages match the haplotype dosages,
#'FALSE is there is a conflict, NA if mrkNA and/or hapNA are TRUE)
#'Each element is itself a list with elements:
#'@examples
#'data(PolyHaplotyper_small)
#'chmd <- compareHapMrkDosages(mrkDosage=phdos, hapresults=phresults)
#'# show results for first haploblock, first 8 individuals:
#'chmd[1, 1:8,]
#'@export
compareHapMrkDosages <- function(mrkDosage, hapresults) {
  if (!all(sapply(hapresults,
                  FUN=function(x) all(c("hapdos", "markers") %in% names(x)))))
    stop("all elements of hapresults must at least have components hapdos and markers")
  hblist <- lapply(hapresults, FUN=function(x) x$markers)
  ploidy <- max(sapply(hapresults,
                       FUN=function(x) max(colSums(x$hapdos), na.rm=TRUE)))
  mrkDosages <- checkmrkDosage(mrkDosage=mrkDosage, ploidy=ploidy,
                               indiv=colnames(hapresults[[1]]$hapdos),
                               markers=unlist(hblist))
  # now with columns in same order as the hapdos columns
  result <- array(NA, dim=c(length(hapresults), ncol(mrkDosage), 3),
                  dimnames=list(names(hapresults), colnames(mrkDosage),
                                c("mrkNA", "hapNA", "match")))
  for (hb in names(hapresults)) {
    mrkdosobs <- mrkDosage[hblist[[hb]],, drop=FALSE]
    nomrk <- colSums(!is.na(mrkdosobs)) == 0
    mrkdoshap <- hapdos2mrkdos(hapdos=hapresults[[hb]]$hapdos,
                               allhap=allHaplotypes(hblist[[hb]]))
    nohap <- colSums(hapresults[[hb]]$hapdos, na.rm=TRUE) != ploidy
    conflict <- mrkdoshap != mrkdosobs
    result[hb,, "mrkNA"] <- nomrk
    result[hb,, "hapNA"] <- nohap
    m <- !nomrk & !nohap & colSums(!is.na(conflict)) > 0
    result[hb, m, "match"] <-
      colSums(conflict[, m, drop=FALSE], na.rm=TRUE) == 0
  }
  result
} #compareHapMrkDosages

#'@title convert the PedigreeSim true haplotypes to marker dosages and
#'haplotype dosages
#'@usage pedigreeSim2PH(ps_geno, haploblock, indiv=NULL, dropUnused=TRUE)
#'@param ps_geno the filename of a *_genotypes.dat file as produced by
#'PedigreeSim, or a data.frame read from such a file (with a column marker,
#'followed by <ploidy> columns per individual with the marker alleles for each
#'homolog; the names of these columns must be <indivname>_<homolog number>).
#'ps_geno should not contain missing data (this will be true unless the
#'PedigreeSim output is midified)
#'@param haploblock a list of character vectors. The names are the names of the
#'haploblocks, the character vectors have the names of the markers in each
#'haploblock. Haplotype names are constructed from the haploblock names, which
#'are used as prefixes to which the (zero-padded) haplotype numbers are are
#'appended with separator '_'.\cr
#'All markers in haploblock must be present in ps_geno. If haploblock is NULL
#'only the marker dosages are generated, not the haplotype dosages.
#'@param indiv the names of the individuals to be extracted; default NULL: all
#'individuals in ps_geno. If not NULL, all indiv must be present in ps_geno
#'@param dropUnused TRUE (default) if the returned matrix should only contain
#'rows for haplotypes that are present; if FALSE matrix contains rows for all
#'possible haplotypes
#'@details if all alleles are in 0/1, these are kept. If only 2 allele symbols
#'occur in ps_geno all are converted (alphabetically) to 0/1 in
#'the same way (e.g. if the alleles are A and B, A -> 0 and B -> 1, even in
#'markers with only B's). If different allele symbols are used between
#'markers, per marker the (alphabetically or numerically) lowest -> 0 and
#'the highest -> 1.
#'So if more than two different allele symbols occur in ps_geno, and
#'one marker has only A and another only B alleles, both are converted
#'to 0's.\cr
#'in mrkDosage, the dosage of the alleles (converted to) 1 is reported
#'@return a list with 2 items:\cr
#'$mrkDosage: a matrix of marker dosages in the input format of inferHaplotypes,
#'with all markers that occur in haploblock, sorted according to haploblock,
#'and for the selected indiv in the specified order\cr
#'$haplist: (only if haploblock is not NULL) a list with one element for each
#'haploblock, similar to the
#'inferHaplotypes output). Each element is itself a list, with two components:\cr
#'$hapdos is a matrix with the haplotype dosages for that haploblock for each
#'individual\cr
#'$markers:  a vector with the names of
#'the markers in the haploblock in the output of inferHaplotypes
pedigreeSim2PH <- function(ps_geno, haploblock, indiv=NULL, dropUnused=TRUE) {
  if (is.character(ps_geno))
    ps_geno <- read.table(ps_geno[1], header=TRUE, sep="\t",
                               check.names=FALSE, stringsAsFactors=FALSE)
  if (is.data.frame(ps_geno)) {
    mrk <- ps_geno$marker
    ps_geno <- as.matrix(ps_geno[, -1])
    rownames(ps_geno) <- mrk
  }
  # now ps_geno is a matrix; check against non-biallelic markers
  # and get the allele symbols per marker:
  # (we do this over all indiv in ps_geno, not only for the selected indiv
  #  because we must know about all allele symbols per marker to be consistent
  #  between different selctions of individuals)
  alleles <- apply(ps_geno, MARGIN=1, FUN=function(x) sort(unique(x)))
  if (is.list(alleles) && any(sapply(alleles, length) > 2))
    stop("ps_geno contains markers with more than 2 alleles")
  if (is.matrix(alleles))
    alleles <- as.data.frame(alleles, stringsAsFactors=FALSE)
  all_alleles <- sort(unique(unlist(alleles)))
  # reduce ps_geno to only the requested markers and individuals:
  if (!(is.null(haploblock))) {
    if (!all(unlist(haploblock) %in% rownames(ps_geno)))
      stop("not all markers in haploblock occur in ps_geno")
    ps_geno <- ps_geno[unlist(haploblock),, drop=FALSE]
  }
  # get indiv names and ploidy and select specified indiv:
  us <- gregexpr("_", colnames(ps_geno))
  last_ <- sapply(us, FUN=function(x) x[length(x)]) # the last "_" in each indiv name
  PSindiv <- unique(substr(colnames(ps_geno), 1, last_-1))
  ploidy <- ncol(ps_geno) / length(PSindiv)
  if (is.null(indiv)) indiv <- PSindiv
  if (!all(indiv %in% PSindiv))
    stop("Not all indiv found in PedigreeSim genotypes file")
  # keep only the specified indiv in the specified order:
  indix <- match(indiv, PSindiv)
  colix <- (indix-1) * ploidy + 1 # the 1st column in ps_geno for each indiv
  colix <- rep(colix, each=ploidy) + (0:(ploidy-1))
  ps_geno <- ps_geno[, colix, drop=FALSE]
  # convert the allele symbols to 0/1:
  if (!all(all_alleles %in% 0:1)) {
    if (length(all_alleles) == 1) ps_geno[] <- 0 else {
      # all markers have the same 2 allele symbols, these are translated
      # to 0/1 identically for all markers, even in markers that have only one
      # of them
      if (length(all_alleles) ==2) {
        newalleles <- (0:1)[match(ps_geno, all_alleles)]
        newalleles <- matrix(newalleles, nrow=nrow(ps_geno),
                             dimnames=dimnames(ps_geno))
      } else {
        # allele symbols differ between markers but all are mono- or bi-allelic
        newalleles <- matrix(NA_integer_, nrow=nrow(ps_geno),
                             ncol=ncol(ps_geno),
                             dimnames=dimnames(ps_geno))
        for (m in seq_len(nrow(newalleles)))
          newalleles[m,] <- (0:1)[match(ps_geno[m,], alleles[[m]])]
      }
      ps_geno <- newalleles; rm(newalleles)
    }
  }
  rm(alleles, all_alleles)
  # get mrkDosage:
  mrkDosage <- ps_geno
  dim(mrkDosage) <- c(nrow(mrkDosage), ploidy, length(indiv))
  mrkDosage <- apply(mrkDosage, MARGIN=c(1,3), sum)
  dimnames(mrkDosage) <- list(rownames(ps_geno), indiv)
  result <- list(mrkDosage=mrkDosage)
  # get list of hapdos:
  if (!is.null(haploblock)) {
    haplist <- vector("list", length(haploblock))
    names(haplist) <- names(haploblock)
    for (hb in seq_along(haploblock)) {
      nhap <- 2^length(haploblock[[hb]])
      dat <- ps_geno[haploblock[[hb]],, drop=FALSE]
      hapnrs <- mrkdos2mrkdid(mrkDosage=dat, indiv=NULL, ploidy=1, check=FALSE)
      #hapnrs is now a vector of hapnrs, of length ploidy * length(indiv)
      hapdos <- hapcomb2hapdos(hapcomb=matrix(hapnrs, nrow=ploidy, byrow=FALSE),
                               nhap=nhap)
      rownames(hapdos) <- paste0(names(haploblock)[hb], "_", padded(1:nhap))
      colnames(hapdos) <- indiv
      if (dropUnused) {
        hapdos <-
          hapdos[rowSums(hapdos, na.rm=TRUE) > 0, , drop=FALSE]
      }
      haplist[[hb]] <- list(hapdos=hapdos, markers=haploblock[[hb]])
    }
    result$haplist <- haplist
  }
  result
}
