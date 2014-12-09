
## ----vignetteSetup, echo=FALSE, message=FALSE, warning = FALSE-----------
## Track time spent on making the vignette
startTime <- Sys.time()

## Bib setup
library('knitcitations')

## Load knitcitations with a clean bibliography
cleanbib()
cite_options(hyperlink = 'to.doc', citation_format = 'text', style = 'html')
# Note links won't show for now due to the following issue
# https://github.com/cboettig/knitcitations/issues/63

## Write bibliography information
bibs <- c(knitcitations = citation('knitcitations'),
    derfinder = citation('derfinder')[1], 
    knitrBootstrap = citation('knitrBootstrap'), 
    knitr = citation('knitr')[3],
    rmarkdown = citation('rmarkdown'),
    brainspan = RefManageR::BibEntry(bibtype = 'Unpublished', key = 'brainspan', title = 'Atlas of the Developing Human Brain [Internet]. Funded by ARRA Awards 1RC2MH089921-01, 1RC2MH090047-01, and 1RC2MH089929-01.', author = 'BrainSpan', year = 2011, url = 'http://developinghumanbrain.org'),
    originalder = citation('derfinder')[2],
    R = citation(),
    IRanges = citation('IRanges'),
    devtools = citation('devtools'),
    testthat = citation('testthat'),
    GenomeInfoDb = citation('GenomeInfoDb'),
    GenomicRanges = citation('GenomicRanges'),
    ggplot2 = citation('ggplot2'),
    biovizBase = citation('biovizBase'),
    bumphunter = citation('bumphunter'),
    TxDb.Hsapiens.UCSC.hg19.knownGene = citation('TxDb.Hsapiens.UCSC.hg19.knownGene'),
    AnnotationDbi = citation('AnnotationDbi'),
    BiocParallel = citation('BiocParallel'),
    derfinderHelper = citation('derfinderHelper'),
    GenomicAlignments = citation('GenomicAlignments'),
    GenomicFeatures = citation('GenomicFeatures'),
    GenomicFiles = citation('GenomicFiles'),
    Hmisc = citation('Hmisc'),
    qvalue = citation('qvalue'),
    Rsamtools = citation('Rsamtools'),
    rtracklayer = citation('rtracklayer'),
    S4Vectors = citation('S4Vectors'),
    bumphunterPaper = RefManageR::BibEntry(bibtype = 'article', key = 'bumphunterPaper', title = 'Bump hunting to identify differentially methylated regions in epigenetic epidemiology studies', author = 'Jaffe, Andrew E and Murakami, Peter and Lee, Hwajin and Leek, Jeffrey T and Fallin, M Daniele and Feinberg, Andrew P and Irizarry, Rafael A', year = 2012, journal = 'International Journal of Epidemiology'),
    derfinderData = citation('derfinderData'),
    polyester = citation('polyester')
)

write.bibtex(bibs,
    file = 'derTutorRef.bib')
bib <- read.bibtex('derTutorRef.bib')

## Assign short names
names(bib) <- names(bibs)

## Working on Windows?
windowsFlag <- .Platform$OS.type == 'windows'

evalDER <- TRUE


## ----'install', eval = FALSE---------------------------------------------
## ## Install derfinder from BioC
## source('http://bioconductor.org/biocLite.R')
## biocLite('derfinder')


## ----'install2', eval = FALSE--------------------------------------------
## biocLite('TxDb.Hsapiens.UCSC.hg19.knownGene')


## ----'rhelp', eval = FALSE-----------------------------------------------
## help(package = 'derfinder')


## ----'loadDer', message = FALSE, eval = evalDER--------------------------
library('derfinder')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')

## Files
files <- paste0('http://lcolladotor.github.io/derTutor/data/sample', 
    1:30, '.bw')
names(files) <- paste0('sample', 1:30)

## Load data
system.time( fullCov <- fullCoverage(files, 'chr22', verbose = FALSE) )

## You are using Windows?
# Use http://lcolladotor.github.io/derTutor/data/fullCov.Rdata


## ----'adjustments', eval = evalDER---------------------------------------
## Calculate library size adjustment
system.time( collapsedFull <- collapseFullCoverage(fullCov) )
lapply(collapsedFull[[1]], head)
sampleDepths <- sampleDepth(collapsedFull, probs = 1)



## ----'models', eval = evalDER--------------------------------------------
## Create models
groupInfo <- factor(rep(c('A', 'B', 'C'), each = 10))
models <- makeModels(sampleDepths = sampleDepths, testvars = groupInfo)


## ----'analysis', eval = evalDER------------------------------------------
## Filter Data
covData <- filterData(fullCov$chr22, cutoff = 0)

## Run analysis for chr22
dir.create('chr22', showWarnings = FALSE)
system.time(
    res <- analyzeChr(chr = 'chr22', coverageInfo = covData, models = models, 
        cutoffFstat = 1e-03, cutoffPre = 0,
        nPermute = 100, seeds = seq_len(100) + 20141202, maxClusterGap = 3000,
        groupInfo = groupInfo, mc.cores = 1,
        lowMemDir = file.path(tempdir(), 'chr22', 'chunksDir'),
        writeOutput = TRUE, returnOutput = TRUE)
)


## ----'exploreRes', eval = evalDER----------------------------------------
names(res)
dir('chr22')
names(res$coveragePrep)
names(res$regions)


## ----'merge', eval = evalDER---------------------------------------------
## Genomic state
system.time(gs <- makeGenomicState(txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, chrs='chr22'))

## Merge results from different chrs
mergeResults('chr22', genomicState = gs$fullGenome, optionsStats = res$optionsStats)


## ----'final', eval = evalDER---------------------------------------------
dir(pattern = 'Rdata')
load('fullRegions.Rdata')
class(fullRegions)
length(fullRegions)


## ----'final2', eval = evalDER--------------------------------------------
colnames(mcols(fullRegions))
table(fullRegions$significantFWER)


## ----'citation'----------------------------------------------------------
## Citation info
citation('derfinder')


## ----createVignette, eval=FALSE------------------------------------------
## ## Create this page
## library('rmarkdown')
## render('index.Rmd')
## 
## ## Clean up
## file.remove('derTutorRef.bib')
## 
## ## Extract the R code
## library('knitr')
## knit('index.Rmd', tangle = TRUE)


## ----reproducibility1, echo=FALSE----------------------------------------
## Date the vignette was generated
Sys.time()


## ----reproducibility2, echo=FALSE----------------------------------------
## Processing time in seconds
totalTime <- diff(c(startTime, Sys.time()))
round(totalTime, digits=3)


## ----reproducibility3, echo=FALSE----------------------------------------
## Session info
library('devtools')
options(width = 120)
session_info()$platform


## ----reproducibility4, echo=FALSE----------------------------------------
## Session info packages
subset(session_info()$packages, package %in% c('derfinder', unique(names(getNamespaceImports('derfinder'))),  'TxDb.Hsapiens.UCSC.hg19.knownGene', 'rmarkdown', 'knitcitations', 'Matrix'))


## ----vignetteBiblio, results = 'asis', echo = FALSE, warning = FALSE-----
## Print bibliography
bibliography()


