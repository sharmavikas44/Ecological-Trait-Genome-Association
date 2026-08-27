## ---------------------------
## 0. Libraries and file paths
## ---------------------------
library(vcfR)
library(GAPIT)
library(EMMREML)  # for kinship (used internally by GAPIT)
library(rrBLUP)   # for alternative kinship if needed

## ------------------------------------
## 1. Read VCF and convert to HapMap
## ------------------------------------
# 1.1 Read VCF
vcf <- read.vcfR("biallelic_GTonly.vcf.gz")

fix_section <- vcf@fix                         # CHROM, POS, ID, REF, ALT, ...
gt_section  <- extract.gt(vcf, element = "GT") # genotype calls (0/0, 0/1, etc.)

dim(fix_section)
dim(gt_section)

# 1.2 Make sure we have a SNP ID column
if (!"ID" %in% colnames(fix_section)) {
  fix_section <- cbind(
    fix_section,
    ID = paste(fix_section[, "CHROM"], fix_section[, "POS"], sep = "_")
  )
}
head(fix_section[, "ID"])

# 1.3 Convert GT calls (0/0, 0/1, 1/1) → AA, AT, TT HapMap codes
geno_hapmap <- apply(gt_section, 2, function(x) {
  sapply(seq_along(x), function(i) {
    geno <- x[i]
    ref  <- fix_section[i, "REF"]
    alt  <- fix_section[i, "ALT"]
    if (geno %in% c("0/0", "0|0")) {
      paste0(ref, ref)           # homozygous REF
    } else if (geno %in% c("0/1", "1/0", "0|1", "1|0")) {
      paste0(ref, alt)           # heterozygous
    } else if (geno %in% c("1/1", "1|1")) {
      paste0(alt, alt)           # homozygous ALT
    } else {
      "NA"
    }
  })
})

# Confirm the structure of geno_hapmap
dim(geno_hapmap)
head(geno_hapmap[, 1:5])

# 1.4 Build HapMap data.frame (standard columns + genotypes)
hapmap_df <- data.frame(
  rs       = fix_section[, "ID"],                              # SNP ID
  alleles  = paste(fix_section[, "REF"], fix_section[, "ALT"], sep = "/"),
  chrom    = fix_section[, "CHROM"],
  pos      = fix_section[, "POS"],
  strand   = "+",
  assembly = "AGPv1",
  center   = "WHBetula",
  protLSID = "NA",
  assayLSID= "NA",
  panel    = "Betula",
  QCcode   = "NA",
  geno_hapmap,
  check.names = FALSE
)

dim(hapmap_df)
head(hapmap_df[, 1:10])

# 1.5 Save HapMap in text format for GAPIT
write.table(
  hapmap_df,
  file = "hapmap_file",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

## OPTIONAL: numeric coding 0/1/2 and map file
convert_to_numeric <- function(geno, ref, alt) {
  if (geno == paste0(ref, ref)) 0
  else if (geno %in% c(paste0(ref, alt), paste0(alt, ref))) 1
  else if (geno == paste0(alt, alt)) 2
  else NA
}

hapmap_numeric <- hapmap_df
geno_idx <- 12:ncol(hapmap_df)  # genotype columns

hapmap_numeric[, geno_idx] <- apply(hapmap_df[, geno_idx], 2, function(col) {
  sapply(seq_len(nrow(hapmap_df)), function(i) {
    alle <- hapmap_df$alleles[i]
    ref  <- substr(alle, 1, 1)
    alt  <- substr(alle, 3, 3)
    convert_to_numeric(col[i], ref, alt)
  })
})

write.table(
  hapmap_numeric,
  file = "hapmap_numeric.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Map file (for GWASpoly etc.)
map_file <- hapmap_df[, c("rs", "chrom", "pos")]
colnames(map_file) <- c("Marker", "Chromosome", "Position")
write.csv(
  map_file,
  "map_file.csv",
  row.names = FALSE,
  quote = FALSE
)

## ------------------------------------
## 2. Phenotype & environmental data
## ------------------------------------

# Phenotype (make sure this is individual-wise if using GAPIT)
myY <- read.table("Traits_ind.txt", header = TRUE, sep = "\t")
head(myY)

# Environmental covariates (soil)
myev <- read.csv("Soil_ind_data.csv")
head(myev)

## ------------------------------------
## 3. Genotype for GAPIT (HapMap)
## ------------------------------------

myG <- read.delim("hapmap_file", header = FALSE)
head(myG[, 1:10])

## ------------------------------------
## 4. Run GAPIT
## ------------------------------------

gapit_results <- GAPIT(
  Y = myY,                     # phenotype
  G = myG,                     # HapMap genotype
  CV = myev,                   # covariates (optional, can be NULL)
  PCA.total = 3,               # number of PCs
  SNP.MAF = 0.05,              # MAF filter
  SNP.fraction = 0.8,          # subset SNPs if needed
  kinship.algorithm = "VanRaden",
  model = c("FarmCPU", "Blink"),
  kinship.cluster = "average",
  kinship.group   = "Mean",
  SNP.test = TRUE,
  SNP.P3D  = TRUE,
  cutOff   = 0.01,             # significance threshold for GWAS
  PCA.View.output  = TRUE,
  Geno.View.output = TRUE,
  file.output      = TRUE      # write all results to working directory
)
