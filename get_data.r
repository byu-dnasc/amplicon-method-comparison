library(readr) # read_csv
library(tidyr) # separate
library(purrr) # list_rbind
library(dplyr) # mutate

naive_reads <- list_rbind(
  map(
    list.files(
      file.path("reads_by_haplotype"), 
      recursive = TRUE,
      full.names = TRUE,
      pattern = "reads_by_haplotype.csv"
    ), 
    function(reads_file_path) {
      read_csv(
          reads_file_path,
          col_names = c("haplotype", "read_name"),
          col_types = "fc"
        ) %>% 
        separate(
          read_name, 
          into = c("movie_id", "zmw_id", "read_type + strand"), 
          sep = "/", 
          extra = "merge",
          convert = TRUE
        ) %>% 
        separate(
          `read_type + strand`, 
          into = c(NA, "strand"), 
          fill = "right"
        ) %>% 
        mutate(
          sample_name = factor(basename(dirname(reads_file_path))), 
          strand = factor(strand),
          movie_id = factor(movie_id)
        )
    }
  )
)
  
pbaa_reads_by_cluster <- list_rbind(
  map(
    list.files(
      file.path("execution"), 
      recursive = TRUE,
      full.names = TRUE,
      pattern = "*read_info.txt"
    ),
    function(read_info_file_path) {
      read_delim(
          read_info_file_path,
          delim = " ", 
          col_names = c("read_name", "guide_name", "strand", "SecondBestGuideName", "Score", "FirstHighest/SecondHighest/UniqueHitSum", "Sample, input fastq", "Sequence length", "average_read_quality", "cluster_id", "cluster_size"),
          show_col_types = FALSE
        ) %>%
        separate(
          read_name, 
          sep = "/", 
          into = c(NA, "zmw_id", NA, "strand"), 
          fill = "right",
          convert = TRUE
        ) %>% 
        mutate(
          strand = factor(strand),
          sample_name = factor(basename(dirname(read_info_file_path)))
        ) %>% 
        select(
          sample_name, 
          zmw_id, 
          strand, 
          average_read_quality, 
          cluster_size, 
          cluster_id
        )
    }
  )
)

clusters_by_haplotype <- read_csv(
    file.path(file.choose()),
    show_col_types = FALSE
  ) %>%
  mutate(haplotype = factor(haplotype)) %>%
  mutate(sample_name = factor(sample_name))

pbaa_reads <- left_join(
  pbaa_reads_by_cluster,
  clusters_by_haplotype, 
  by = c("sample_name", "cluster_id")
); rm(pbaa_reads_by_cluster); rm(clusters_by_haplotype)
