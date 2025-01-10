# This script is used to generate a comparison between time and memory used 
# for Cell Ranger and Alevin-fry 

renv::load()
library(ggplot2)
library(ggforce)

theme_set(
  theme_classic() + 
    theme(
      aspect.ratio = 1,
      text = element_text(size=14),
      axis.text.x = element_text(angle = 90, vjust = 0.5), 
      legend.position = "top",
      legend.text = element_text(size = rel(1.175)),
      strip.background = element_rect(linewidth = rel(0.75)), 
      axis.line = element_line(linewidth = rel(0.5)),
      axis.ticks = element_line(linewidth = rel(0.5))
      )
)

# Set up -----------------------------------------------------------------------

# nextflow logs directory 
logs_dir <- here::here("nextflow_logs")
af_log_file <- file.path(logs_dir, "alevin-fry-trace.txt")
cellranger_log_file <- file.path(logs_dir, "cellranger-trace.txt")

# create a named list of log files
log_file_paths <- list("Alevin-fry" = af_log_file,
                      "Cell Ranger" = cellranger_log_file)

# method palette
palette_file <- here::here("palettes", "method-palette.tsv")

# path to output plot
pdf_dir <- here::here("figures", "pdfs")
output_pdf_file <- file.path(pdf_dir, "FigS1A_time-memory-benchmarking.pdf")

# define single-cell and single-nuc samples
single_cell <- c("SCPCR000003", "SCPCR000126", "SCPCR000127")
single_nuclei <- c("SCPCR000220", "SCPCR000221", "SCPCR000495")

# Prep data for plotting -------------------------------------------------------

# read in log files and combine into one data frame
log_df <- log_file_paths |> 
  purrr::map(readr::read_tsv) |> 
  dplyr::bind_rows(.id = "method")


modified_log_df <- log_df |>
  # pull out process tag into process, run_id, and index
  tidyr::separate(name, into = c("process", "run_id"), sep = " ") |> 
  tidyr::separate(run_id, into = c("run_id", "index_name"), sep = "-") |> 
  # remove any extra () hanging around 
  dplyr::mutate(
    run_id = stringr::str_remove(run_id, "[()]")
  ) |> 
  # separate peak rss into rss and units
  tidyr::separate(peak_rss, into = c("peak_rss", "rss_units"), sep = " ") |>
  ## convert all MB to GB 
  dplyr::mutate(peak_rss = as.numeric(peak_rss),
                peak_rss = dplyr::case_when(rss_units == "MB" ~ (peak_rss/1000),
                              TRUE ~ peak_rss),
                rss_units = "GB") |>
  # separate hours, minutes, seconds
  dplyr::mutate(
    hours = stringr::str_extract(duration, "\\d+h") |> 
      stringr::str_remove("h"),
    minutes = stringr::str_extract(duration, "\\d+m") |> 
      stringr::str_remove("m"),
    seconds = stringr::str_extract(duration, "\\d+s") |>
      stringr::str_remove("s")
  ) |> 
  # make sure everything is numeric
  dplyr::mutate_at(c("hours", "minutes", "seconds"), as.numeric) |> 
  # replace any NA's with 0
  dplyr::mutate_if(is.numeric, ~tidyr::replace_na(., 0)) |>
  # get total time in minutes 
  dplyr::mutate(total_time = hours*60 + minutes + seconds/60)
  
# total time and memory across all processes
grouped_log_df <- modified_log_df |> 
  dplyr::group_by(method, run_id) |> 
  # sum time across all processes
  # get max peak rss used 
  dplyr::summarise(total_time = sum(total_time),
                   total_memory = max(peak_rss)) |>
  # add in seq unit 
  dplyr::mutate(
    seq_unit = dplyr::case_when(
      run_id %in% single_cell ~ "Single-cell",
      run_id %in% single_nuclei ~ "Single-nuclei"
    )
  ) 

# Plot -------------------------------------------------------------------------

# read in palette colors 
palette <- readr::read_tsv(palette_file)

# get list of all colors 
method_colors <- palette$color |> 
  purrr::set_names(palette$method)

# order by time
time_ordered_df <- grouped_log_df |> 
  dplyr::mutate(run_id = forcats::fct_reorder(run_id, total_time))

# compare run time in minutes
time_plot <- ggplot(time_ordered_df, aes(x = run_id, y = total_time, fill = method)) + 
                      geom_col(position="dodge") +
  facet_wrap(vars(seq_unit), 
             scales = "free") + 
  labs(x = "",
       y = "Total run time (minutes)",
       fill = "") +
  scale_fill_manual(values = method_colors)

# order by memory
mem_ordered_df <- grouped_log_df |> 
  dplyr::mutate(run_id = forcats::fct_reorder(run_id, total_memory))


# compare peak memory in GB
memory_plot <- ggplot(mem_ordered_df, aes(x = run_id, y = total_memory, fill = method)) + 
  geom_col(position="dodge") +
  facet_wrap(vars(seq_unit), 
             scales = "free") + 
  labs(x = "",
       y = "Peak memory (GB)", 
       fill = "") +
  scale_fill_manual(values = method_colors)

# combine into one side by side plot
combined_plot <- patchwork::wrap_plots(list(time_plot, memory_plot), ncol = 1, guides = "collect") & theme(legend.position = "top")
# export as png
# using width and height that were exported when width and height weren't specified
ggsave(output_pdf_file, plot = combined_plot, width = 10, height = 10)

