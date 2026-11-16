# EPIC-array preprocessing with ChAMP: import, filtering, BMIQ and ComBat.
# This optional upstream workflow requires raw IDAT files and SampleSheet.csv
# under rawdata/. It does not derive the clinical analysis data frame or compute
# DNAmTL/clocks; those upstream steps were not included in the supplied scripts.
source("functions.R", local = TRUE)
require_packages("ChAMP")

raw_directory <- "rawdata"
sample_sheet <- file.path(raw_directory, "SampleSheet.csv")
if (!file.exists(sample_sheet)) stop("Place IDAT files and SampleSheet.csv in rawdata/.")
samples <- read.csv(sample_sheet, stringsAsFactors = FALSE, check.names = FALSE)
# Inspect/rename documented columns only; never drop a column by position.
names(samples)[names(samples) == "Slide"] <- "Sentrix_ID"
names(samples)[names(samples) == "Array"] <- "Sentrix_Position"
assert_columns(samples, c("Sample_Name", "Sentrix_ID", "Sentrix_Position"))
if (anyDuplicated(samples$Sample_Name)) stop("Sample_Name must be unique.")
# ChAMP reads the sheet in the IDAT directory. If names need updating, preserve
# an untouched backup before writing the standardized local input sheet.
original_names <- names(read.csv(sample_sheet, check.names = FALSE))
if (!identical(original_names, names(samples))) {
  backup <- paste0(sample_sheet, ".original")
  if (!file.exists(backup)) file.copy(sample_sheet, backup, overwrite = FALSE)
  write.csv(samples, sample_sheet, row.names = FALSE)
}

myLoad <- ChAMP::champ.import(directory = raw_directory, arraytype = "EPIC")
save(myLoad, file = result_file("methylation", "myLoad.RData"))
myLoad_filter <- ChAMP::champ.filter(
  beta = myLoad$beta, pd = myLoad$pd,
  filterDetP = TRUE, detPcut = 0.01, detP = myLoad$detP,
  ProbeCutoff = 0.1, SampleCutoff = 0.2,
  beadcount = myLoad$beadcount, filterBeads = TRUE, beadCutoff = 1,
  filterNoCG = TRUE, filterSNPs = FALSE, filterMultiHit = FALSE,
  filterXY = FALSE, fixOutlier = TRUE, arraytype = "EPIC"
)
save(myLoad_filter, file = result_file("methylation", "myLoad_filter.RData"))

# Explicitly use the BMIQ method described in the manuscript.
myNorm <- ChAMP::champ.norm(beta = myLoad_filter$beta, arraytype = "EPIC", method = "BMIQ",
                            cores = getOption("dnamtl.cores", 1L))
save(myNorm, file = result_file("methylation", "myNorm.RData"))

# Align metadata to retained samples by ID. This replaces a hard-coded sample
# deletion and adapts to whichever samples pass filtering in the current run.
pd <- myLoad_filter$pd
assert_columns(pd, c("Sample_Name", "Pool_ID", "age"))
if (anyDuplicated(pd$Sample_Name)) stop("Filtered sample identifiers are not unique.")
idx <- match(colnames(myNorm), pd$Sample_Name)
if (anyNA(idx)) stop("Normalized methylation samples do not match the sample sheet.")
pd <- pd[idx, , drop = FALSE]
stopifnot(identical(as.character(pd$Sample_Name), colnames(myNorm)))
# Preserve the source's three-group batch-protection design, correcting age 50
# so it belongs to the middle group. Review notes explain this upstream change.
pd$Sample_Group <- factor(ifelse(pd$age < 50, "EOCRC", ifelse(pd$age < 55, "IOCRC", "LOCRC")))
pd$Sample_Group2 <- factor(ifelse(pd$age < 50, "EOCRC", ifelse(pd$age < 69, "IOCRC", "LOCRC")))
save(pd, file = result_file("methylation", "pdmyload.RData"))

myCombat <- ChAMP::champ.runCombat(beta = myNorm, pd = pd, batchname = "Pool_ID")
myCombat <- as.matrix(myCombat)
save(myCombat, file = result_file("methylation", "myCombat.RData"))
ChAMP::champ.SVD(beta = myCombat, pd = pd,
                resultsDir = result_file("methylation", "SVD"))
beta_all <- data.frame(ProbeID = rownames(myCombat), myCombat, check.names = FALSE)
save(beta_all, file = result_file("methylation", "beta_alldf.RData"))
save_session_info("01_methylation")
