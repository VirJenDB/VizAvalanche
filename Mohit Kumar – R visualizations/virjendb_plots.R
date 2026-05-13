Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")

# always save outputs next to this script
local({
  args <- commandArgs(trailingOnly = FALSE)
  f    <- sub("--file=", "", args[startsWith(args, "--file=")])
  if (length(f)) setwd(dirname(normalizePath(f)))
})

suppressPackageStartupMessages({
  library(httr2)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(plotly)
  library(lubridate)
  library(scales)
  library(htmlwidgets)
})

API_URL      <- "https://api2.virjendb.org/v2/search"
BRAND_COLORS <- c("#E37B40", "#58A4B0")
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# ── API ────────────────────────────────────────────────────────────────────────

.fetch_page <- function(query, from = 0L, size = 50L) {
  request(API_URL) |>
    req_body_json(list(query = query, size = size, from = from)) |>
    req_retry(max_tries = 3L, backoff = \(i) 2^i) |>
    req_perform() |>
    resp_body_json()
}

count_q <- function(query) {
  tryCatch(
    { Sys.sleep(0.05); .fetch_page(query, size = 1L)$total },
    error = function(e) { message("  skip: ", query); NA_integer_ }
  )
}

.flatten_source <- function(src) {
  map_chr(src, \(v) {
    v <- v %||% NA_character_
    if (is.list(v) || length(v) > 1L) paste(v, collapse = "|") else as.character(v)
  })
}

fetch_diverse <- function(query, n_positions = 8L) {
  first <- .fetch_page(query, from = 0L)
  total <- first$total
  rows  <- map(first$results, \(r) .flatten_source(r$source))
  if (total > 50L && n_positions > 1L) {
    offsets <- as.integer(seq(50L, max(51L, total - 50L), length.out = n_positions - 1L))
    for (from in offsets) {
      tryCatch(
        rows <- c(rows, map(.fetch_page(query, from = from)$results,
                            \(r) .flatten_source(r$source))),
        error = function(e) invisible(NULL)
      )
      Sys.sleep(0.05)
    }
  }
  bind_rows(map(rows, as.list))
}

diverse_sample <- function() {
  message("  sampling diverse records …")
  bind_rows(
    fetch_diverse('"Genome Completeness":complete',    n_positions = 5L),
    fetch_diverse('"Genome Completeness":partial',     n_positions = 5L),
    fetch_diverse('"Family NCBI":"Orthomyxoviridae"', n_positions = 3L), # influenza
    fetch_diverse('"Family NCBI":"Flaviviridae"',     n_positions = 3L), # dengue/HCV
    fetch_diverse('"Family NCBI":"Herpesviridae"',    n_positions = 3L), # dsDNA, high-GC
    fetch_diverse('"Molecule Type":dsDNA',            n_positions = 3L)  # DNA viruses
  ) |> distinct(`VirJenDB Accession`, .keep_all = TRUE)
}

# ── Plot 1: GC Content ─────────────────────────────────────────────────────────

plot_gc_clusters <- function() {
  message("Plot 1: GC content …")
  raw <- diverse_sample()

  df <- raw |>
    transmute(
      gc    = suppressWarnings(as.numeric(`GC Content BV-BRC`)),
      acc   = `VirJenDB Accession`,
      clref = `Cluster Reference`
    ) |>
    filter(!is.na(gc), gc > 0, gc < 100) |>
    mutate(
      type   = if_else(!is.na(clref) & clref == acc,
                       "Cluster Representative", "Cluster Member"),
      gc_bin = cut(gc, breaks = seq(0, 100, by = 5),
                   labels = paste0(seq(0, 95, by = 5), "–", seq(5, 100, by = 5)),
                   include.lowest = TRUE)
    ) |>
    filter(!is.na(gc_bin))

  counts <- df |>
    count(gc_bin, type) |>
    mutate(type = factor(type, levels = c("Cluster Representative", "Cluster Member")))

  p <- ggplot(counts, aes(gc_bin, n, fill = type)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = setNames(BRAND_COLORS,
                                        c("Cluster Representative", "Cluster Member"))) +
    scale_y_log10(labels = label_comma()) +
    labs(title = "Sequences by GC Content Bin",
         subtitle = "Cluster role · log-scale Y",
         x = "GC Content (%)", y = "Count (log₁₀)", fill = NULL) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")

  list(
    html = ggplotly(p) |> layout(legend = list(orientation = "h", x = 0, y = 1.12)),
    pdf  = p
  )
}

# ── Plot 2: Sequence Length by Family ─────────────────────────────────────────

plot_seq_len_by_family <- function() {
  message("Plot 2: Seq length by family …")
  sample_raw <- diverse_sample()

  top_families <- sample_raw |>
    filter(!is.na(`Family NCBI`), `Family NCBI` != "") |>
    count(`Family NCBI`) |>
    filter(n >= 2L) |>
    slice_max(n, n = 15L) |>
    pull(`Family NCBI`)

  message(sprintf("  fetching records for %d families …", length(top_families)))

  df <- map_dfr(top_families, \(fam) {
    tryCatch(
      fetch_diverse(sprintf('"Family NCBI":"%s"', fam), n_positions = 4L) |>
        transmute(family = `Family NCBI`,
                  len    = suppressWarnings(as.numeric(`Sequence Length`))) |>
        filter(!is.na(len), len > 0),
      error = function(e) tibble(family = character(), len = numeric())
    )
  })

  family_order <- df |>
    group_by(family) |>
    summarise(med = median(len), .groups = "drop") |>
    arrange(med) |>
    pull(family)

  df <- df |> mutate(family = factor(family, levels = family_order))

  p <- ggplot(df, aes(family, len)) +
    geom_boxplot(fill = BRAND_COLORS[2], color = "black",
                 outlier.size = 0.7, outlier.alpha = 0.3) +
    scale_y_log10(labels = label_comma()) +
    coord_flip() +
    labs(title = "Sequence Length by Viral Family",
         subtitle = "Dedicated sample per family · log-scale Y",
         x = "NCBI Family", y = "Sequence Length (bp, log₁₀)") +
    theme_classic(base_size = 14)

  list(html = ggplotly(p), pdf = p)
}

# ── Plot 3: Submissions Over Time ──────────────────────────────────────────────

plot_submissions_over_time <- function() {
  message("Plot 3: Submissions over time …")
  years <- 1990L:as.integer(format(Sys.Date(), "%Y"))
  message(sprintf("  querying %d years …", length(years)))

  yearly <- map_dfr(years, \(y) tibble(
    year = y,
    n    = count_q(sprintf('"Release Date":%d', y))
  )) |> filter(!is.na(n), n > 0L)

  p <- ggplot(yearly, aes(year, n)) +
    geom_col(fill = BRAND_COLORS[1]) +
    scale_x_continuous(breaks = seq(1990, max(yearly$year), by = 5)) +
    scale_y_log10(labels = label_comma()) +
    labs(title = "VirJenDB Sequences Released per Year",
         subtitle = "Exact full-DB counts · log-scale Y",
         x = "Year", y = "Sequences (log₁₀)") +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  list(html = ggplotly(p), pdf = p)
}

# ── Plot 4: World Map ──────────────────────────────────────────────────────────

plot_world_map <- function() {
  message("Plot 4: World map …")
  sample_raw <- diverse_sample()

  countries <- sample_raw |>
    filter(!is.na(Country), Country != "") |>
    distinct(Country) |>
    pull(Country)

  message(sprintf("  querying %d countries …", length(countries)))

  df <- map_dfr(countries, \(ctry) tibble(
    Country   = ctry,
    sequences = count_q(sprintf('"Country":"%s"', ctry))
  )) |> filter(!is.na(sequences), sequences > 0L) |>
    arrange(desc(sequences))

  p_html <- plot_geo(df,
    locationmode = "country names",
    locations    = ~Country,
    z            = ~log10(sequences),
    text         = ~paste0(Country, "<br>", format(sequences, big.mark = ","), " sequences"),
    hoverinfo    = "text",
    colorscale   = list(list(0, "#ffffff"), list(0.25, "#c8e6f0"),
                        list(0.5, "#7bbfd4"), list(0.75, "#3890b8"),
                        list(1, BRAND_COLORS[2])),
    colorbar     = list(title = "Sequences\n(log₁₀)",
                        tickvals = 0:6,
                        ticktext = c("1","10","100","1K","10K","100K","1M"))
  ) |>
    layout(
      title = list(text = "VirJenDB Sequences by Collection Country"),
      geo   = list(showframe = FALSE, showcoastlines = TRUE,
                   coastlinecolor = "#aaaaaa", showland = TRUE,
                   landcolor = "#f5f5f5", projection = list(type = "natural earth"))
    )

  world <- map_data("world")
  p_pdf <- ggplot() +
    geom_map(data = world, map = world,
             aes(long, lat, map_id = region),
             fill = "#f0f0f0", color = "#aaaaaa", linewidth = 0.2) +
    geom_map(data = df, map = world,
             aes(fill = log10(sequences), map_id = Country)) +
    scale_fill_gradient(low = "#c8e6f0", high = BRAND_COLORS[2],
                        name = "Sequences\n(log₁₀)") +
    coord_fixed(1.3) +
    labs(title = "VirJenDB Sequences by Collection Country") +
    theme_void(base_size = 12)

  list(html = p_html, pdf = p_pdf)
}

# ── Plot 5: Molecule Type ──────────────────────────────────────────────────────

plot_molecule_type <- function() {
  message("Plot 5: Molecule type …")
  sample_raw <- diverse_sample()

  mol_types <- sample_raw |>
    filter(!is.na(`Molecule Type`), `Molecule Type` != "") |>
    distinct(`Molecule Type`) |>
    pull(`Molecule Type`)

  message(sprintf("  querying %d molecule types …", length(mol_types)))

  df <- map_dfr(mol_types, \(mt) tibble(
    `Molecule Type` = mt,
    n = count_q(sprintf('"Molecule Type":"%s"', mt))
  )) |>
    filter(!is.na(n), n > 0L) |>
    arrange(n) |>
    mutate(`Molecule Type` = factor(`Molecule Type`, levels = `Molecule Type`))

  p_html <- plot_ly(df,
    x = ~n, y = ~`Molecule Type`, type = "bar", orientation = "h",
    marker        = list(color = rep_len(BRAND_COLORS, nrow(df))),
    text          = ~format(n, big.mark = ","),
    textposition  = "outside",
    hovertemplate = "<b>%{y}</b><br>Count: %{x:,}<extra></extra>"
  ) |>
    layout(title = list(text = "Sequences by Molecule Type"),
           xaxis = list(title = "<b>Count (log scale)</b>", type = "log"),
           yaxis = list(title = ""),
           margin = list(l = 130, r = 90),
           font = list(size = 14, color = "black"),
           plot_bgcolor = "white", paper_bgcolor = "white")

  p_pdf <- ggplot(df, aes(n, `Molecule Type`, fill = `Molecule Type`)) +
    geom_col() +
    scale_x_log10(labels = label_comma()) +
    scale_fill_manual(values = rep_len(BRAND_COLORS, nrow(df))) +
    labs(title = "Sequences by Molecule Type",
         x = "Count (log scale)", y = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")

  list(html = p_html, pdf = p_pdf)
}


message("\n=== VirJenDB R Visualizations ===\n")

p1 <- plot_gc_clusters()
p2 <- plot_seq_len_by_family()
p3 <- plot_submissions_over_time()
p4 <- plot_world_map()
p5 <- plot_molecule_type()

plots      <- list(p1, p2, p3, p4, p5)
html_files <- c("gc_content_clusters.html", "seq_len_by_family.html",
                "submissions_over_time.html", "world_map.html",
                "molecule_type_log.html")

walk2(plots, html_files, \(p, f) {
  saveWidget(p$html, f, selfcontained = TRUE)
  message("  Saved: ", f)
})

pdf("virjendb_plots.pdf", width = 12, height = 7)
walk(plots, \(p) print(p$pdf))
invisible(dev.off())
message("  Saved: virjendb_plots.pdf")

message("\nDone.")
