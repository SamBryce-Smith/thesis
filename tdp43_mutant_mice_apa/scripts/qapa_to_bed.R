library(tidyverse)

# Extract last exon coordinates from QAPA results to BED file



dexseq_m323k_mef <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/M323K_MEFs/dexseq_apa.HOMvsWT.results.processed.tsv")
dexseq_f210i_mef <- read_tsv("data/QAPA_mouse_data_gdrive/QAPA_mouse_data/F210I_MEFs/dexseq_apa.HOMvsWT.results.processed.tsv")

m323k_mef_bed <- dexseq_m323k_mef %>%
  group_by(Gene_Name, Chr) %>%
  summarise(start = min(LastExon.Start),
            end = max(LastExon.End),
            strand = unique(Strand),
            score = min(padj)) %>%
  ungroup() %>%
  rename(chr = Chr, name = Gene_Name) %>%
  mutate(sig = if_else(score < 0.05, "sig", "bg"),
         name = paste(name, sig, sep = "|")) %>%
  select(chr, start, end, name, score, strand) %>%
  arrange(chr, start, end)


f210i_mef_bed <- dexseq_f210i_mef %>%
  group_by(Gene_Name, Chr) %>%
  summarise(start = min(LastExon.Start),
            end = max(LastExon.End),
            strand = unique(Strand),
            score = min(padj)) %>%
  ungroup() %>%
  rename(chr = Chr, name = Gene_Name) %>%
  mutate(sig = if_else(score < 0.05, "sig", "bg"),
         name = paste(name, sig, sep = "|")) %>%
  select(chr, start, end, name, score, strand) %>%
  arrange(chr, start, end)

write_tsv(m323k_mef_bed, "processed/2024-02-29/M323K_MEFs/m323k_mefs.dexseq_last_exons.bed", col_names = FALSE)
write_tsv(m323k_mef_bed, "processed/2024-02-29/F210I_MEFs/f210i_mefs.dexseq_last_exons.bed", col_names = FALSE)
