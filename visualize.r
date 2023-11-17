naive_strand_counts <- naive_reads %>% 
  group_by(sample_name, haplotype) %>% 
  summarize(num_strands = n() + sum(is.na(strand))) %>% 
  ungroup(haplotype) %>% 
  mutate(
    proportion = num_strands / sum(num_strands), # normalize so that all facets have the same scale
    haplotype = ordered(haplotype, levels = c("ac", "ag", "gc", "gg"))
  ) 

pbaa_strand_counts <- pbaa_reads %>% 
  group_by(sample_name, haplotype) %>% 
  summarize(num_strands = n() + sum(is.na(strand))) %>% 
  ungroup(haplotype) %>% 
  mutate(
    proportion = num_strands / sum(num_strands), # normalize so that all facets have the same scale
    haplotype = ordered(haplotype, levels = c("ac", "ag", "gc", "gg"))
  ) 

haplotype_colors = c(
  "ac" = "#FA8072", # red
  "ag" = "#00FF7F", # green
  "gc" = "#87CEEB", # blue
  "gg" = "#F0E68C", # yellow
  "NA" = "grey"
)

# method, sample_name, num_strands
slope_df <- bind_rows(
    naive_strand_counts,
    pbaa_strand_counts,
    .id = "method"
  ) %>% 
  mutate(method = factor(ifelse(method == 1, "naive", "pbaa")))
  
  
library(CGPfunctions)

# compare method of determining of haplotypes 
ht <- "ac"; newggslopegraph(
  dataframe = filter(
    slope_df, 
    haplotype == ht
  ),
  Measurement = num_strands,
  Times = method,
  Grouping = sample_name,
  Title = paste("Haplotype", ht, "assignment"),
  SubTitle = "Naive vs. pbAA"
)
  
# view haplotype assignments
df <- naive_strand_counts; ggplot(df, aes(x = "", y = proportion, fill = haplotype)) +
  geom_col() +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = num_strands), 
    position = position_stack(vjust = 0.5)
  ) +
  facet_wrap(~sample_name) +
  scale_fill_manual(values = haplotype_colors) +
  theme_void() # remove background, grid, labels
