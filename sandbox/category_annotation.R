library(rtracklayer)

# 1. Import your GTF (the gzipped file)
gtf <- import("data/flair_filter_transcripts.gtf.gz", format="gtf")

# 2. Subset to the transcript records
tx   <- gtf[gtf$type == "transcript"]
library(ensembldb)
library(EnsDb.Hsapiens.v86)

# 1. Pull all transcripts + biotype
tx_ens <- transcripts(
  EnsDb.Hsapiens.v86,
  columns = c("tx_id", "tx_name", "tx_biotype")
)

# 2. Subset to your transcripts
df_ens <- as.data.frame(tx_ens)
sub    <- df_ens[df_ens$tx_name %in% mcols(tx)$transcript_id, ]

# 3. Attach biotypes onto your tx GRanges
mcols(tx)$transcript_biotype <- 
  sub$tx_biotype[match(mcols(tx)$transcript_id, sub$tx_name)]

# 4. Check
table(mcols(tx)$transcript_biotype)