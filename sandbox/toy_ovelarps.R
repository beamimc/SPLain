library(plyranges)

gr <- data.frame(
  seqnames = 1,
  start=c(1,11,21,1,21,31,41,51),
  width=5,
  txp=rep(letters[1:3], c(3,2,3))
) |>
  as_granges()

gr <- bind_ranges(c(gr, gr)) |>
  mutate(
    gene = rep(1:2, each=8),
    coef = ifelse(txp == "b", "-", "+"),
    unit = 1
  )

gr <- flat_sig_exons|>
  mutate(unit=1)

gr |>
  as.data.frame()

gr |>
  # filter(txp != "c") |>
  group_by(gene) |>
  reduce_ranges(cover = sum(unit))

gr |>
  # filter(txp != "c") |>
  group_by(gene) |>
  reduce_ranges(cover = sum(unit)) |>
  group_by(gene) |>
  summarize(
    total_width = sum(width),
    intersection_width = sum(width[cover > 1]),
    overlap_score = sum(width[cover > 1]) /sum(width)
  )|> tibble::as_tibble() |> print(n=25)




table <- gr |>
  # filter(txp != "c") |>
  group_by(gene,isoform) |>
  summarize(
    total_width = sum(width),
    range = max(end) - min(start),
    start = min(start),
    end = max(end)
    
    )|>
  tibble::as_tibble() |> print(n=25)



library(dplyr)

df2 <- table %>%
  group_by(gene) %>%
  # grab the start/end of the widest interval in each gene
  mutate(
    ref_start = start[which.max(range)],
    ref_end   = end  [which.max(range)]
  ) %>%
  ungroup() %>%
  # compute raw overlap length and then the score
  mutate(
    overlap = pmax(0, pmin(end, ref_end) - pmax(start, ref_start)),
    score   = overlap / range
  ) %>%
  select(gene, ids, range, start, end, score)


df2 |> filter(gene == "ENSG00000145730.20")

