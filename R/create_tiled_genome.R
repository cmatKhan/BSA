# tiled_genome_df = tileGenome(seqlengths(kn99_genome), tilewidth = 10000, cut.last.tile.in.chrom = TRUE)
#
#
# tiled_genome_df = tiled_genome_df %>% as_tibble() %>% left_join(chr_map)
#
# chr_map = tibble(seqnames = seqnames(kn99_genome), chr=c(seq(1:14), "M"))

# x = filtered_allPoolsInOneBSA$oneMouse$brain_15 %>%
#   filter(CHROM == 1) %>%
#   mutate(tile = cut(POS, filter(tiled_genome_df, chr == 1)$start, include.lowest = TRUE, right = FALSE))
