library(plyranges)
library(dplyr)
library(tidyr)

# set up a toy example

df <- data.frame(
  seqnames = "chr1",
  start = c(1, 11, 21, 31),
  width = 5,
  rank = 1:4,
  gene = rep(1, 4),
  txp = rep(1, 4),
  coef = rep("+", 4),
  candidate = c(0,0,1,0)
)


gr <- df |> as_granges()

gr <- bind_ranges(gr, gr) |>
  mutate(
    txp = rep(1:2, each=4),
    txp = rep(1:2, each=4),
    coef = rep(c("+","-"), each=4),
    candidate = case_when(coef == "-" ~ 0, TRUE ~ candidate) #only + exons can be candidates
  
  )

gr <- gr |>
  shift_right(rep(c(0,10), c(6,2))) #shift right last 2 exons

# add an identifier
gr <- gr |>
  mutate(txp_rank = paste0(txp, "-", rank))

# candidate skip exons
skip <- gr |>  
  filter(candidate == 1)

left_exons <- gr |>
  slice(match(paste0(skip$txp, "-", skip$rank-1), txp_rank))

right_exons <- gr |>
  slice(match(paste0(skip$txp, "-", skip$rank+1), txp_rank))

neg_exons <- gr |>
  filter(coef == "-")

# check if left and right exons
# of candidate skip exons are both
# present (technically, just overlap)
# the exons of negative coefficient txps
skip |>
  mutate(left_and_right =
           left_exons %over% neg_exons &
           right_exons %over% neg_exons
  )

# this is just a little demo, maybe needs
# testing or assumes also that the candidates
# are not initial or terminal exons within
# their transcripts...
