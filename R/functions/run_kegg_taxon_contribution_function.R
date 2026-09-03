################################################################################
# Script: run_kegg_taxon_contribution_function.R
# Purpose:
#   Shared, documented helpers for Step 5's MetaCyc/KEGG contribution-table
#   filtering, pathway annotation, zero-complete condition summaries, taxonomy-
#   rank aggregation, interactive heatmaps, and KEGG workbook export.
#
# Statistical convention:
#   A taxon absent from a sample contributes zero; it must not be dropped from
#   a condition mean. complete_contribution_grid() explicitly inserts those
#   zeros before means, sample standard deviations, and heatmap percentages are
#   calculated.
################################################################################

# Resolve a PICRUSt2 table regardless of whether it is stored as plain TSV or
# gzip-compressed TSV.
resolve_picrust2_table_path <- function(path_without_gz) {
  if (file.exists(path_without_gz)) return(path_without_gz)
  compressed_path <- paste0(path_without_gz, ".gz")
  if (file.exists(compressed_path)) return(compressed_path)
  NA_character_
}

# Load one documented ggpicrust2 package-data object without depending on it
# already being attached to the caller's workspace or package namespace.
load_ggpicrust2_reference <- function(data_name) {
  reference_environment <- new.env(parent = emptyenv())
  suppressWarnings(utils::data(
    list = data_name, package = "ggpicrust2",
    envir = reference_environment
  ))
  if (!exists(data_name, envir = reference_environment, inherits = FALSE)) {
    stop("ggpicrust2's bundled ", data_name, " dataset could not be loaded.")
  }
  as.data.frame(get(data_name, envir = reference_environment, inherits = FALSE))
}

# Build a memory-safe shell command that retains the contribution-table header
# and rows whose function column matches one of the requested identifiers.
# Compressed inputs are streamed through gzip and remain compressed on disk.
resolve_required_executable <- function(command_name) {
  executable <- unname(Sys.which(command_name))
  if (!nzchar(executable)) {
    stop("Required executable '", command_name, "' was not found on PATH.")
  }
  executable
}

build_contribution_filter_command <- function(contribution_path, function_ids) {
  function_ids <- unique(as.character(function_ids))
  function_ids <- function_ids[!is.na(function_ids) & nzchar(function_ids)]
  if (length(function_ids) == 0) stop("At least one function ID is required.")
  if (any(grepl(",", function_ids, fixed = TRUE))) {
    stop("Function IDs must not contain commas.")
  }
  awk_executable <- resolve_required_executable("awk")
  awk_filter <- sprintf(
    paste0(
      "%s -F'\\t' -v requested_ids=%s ",
      "'BEGIN{n=split(requested_ids, requested, \",\"); ",
      "for(i=1; i<=n; i++) ids[requested[i]]=1} ",
      "NR==1 || ($2 in ids)'"
    ),
    shQuote(awk_executable),
    shQuote(paste(function_ids, collapse = ","))
  )
  if (grepl("\\.gz$", contribution_path, ignore.case = TRUE)) {
    gzip_executable <- resolve_required_executable("gzip")
    paste(shQuote(gzip_executable), "-cd", shQuote(contribution_path), "|", awk_filter)
  } else {
    paste(awk_filter, shQuote(contribution_path))
  }
}

# Build one standardized pathway-ID-to-name lookup for tables and plots.
# Local references are used first so routine execution does not require a
# network connection. Step 4's annotated workbook is an optional second local
# source. KEGG REST is queried only for KEGG IDs still missing a name.
build_pathway_name_lookup <- function(pathway_ids, pathway_level,
                                      step4_workbook = NULL) {
  pathway_ids <- unique(as.character(pathway_ids))
  pathway_ids <- pathway_ids[!is.na(pathway_ids) & nzchar(pathway_ids)]
  lookup <- data.frame(
    function_id = pathway_ids,
    pathway_name = rep(NA_character_, length(pathway_ids)),
    stringsAsFactors = FALSE
  )

  fill_missing_names <- function(candidate_lookup) {
    if (is.null(candidate_lookup) || nrow(candidate_lookup) == 0) return(invisible(NULL))
    candidate_lookup <- candidate_lookup %>%
      dplyr::transmute(
        function_id = as.character(function_id),
        pathway_name = as.character(pathway_name)
      ) %>%
      dplyr::filter(!is.na(pathway_name), nzchar(pathway_name)) %>%
      dplyr::distinct(function_id, .keep_all = TRUE)
    matched <- match(lookup$function_id, candidate_lookup$function_id)
    replace <- is.na(lookup$pathway_name) & !is.na(matched)
    lookup$pathway_name[replace] <<- candidate_lookup$pathway_name[matched[replace]]
    invisible(NULL)
  }

  if (identical(pathway_level, "MetaCyc")) {
    # ggpicrust2 2.5.x ships this reference as the documented package-data
    # object `metacyc_reference`; it is not an internal namespace object.
    # Load it into a private environment so the helper does not depend on a
    # caller having attached the dataset earlier in the R session.
    metacyc_reference <- load_ggpicrust2_reference("metacyc_reference")
    fill_missing_names(data.frame(
      function_id = as.character(metacyc_reference[[1]]),
      pathway_name = as.character(metacyc_reference[[2]]),
      stringsAsFactors = FALSE
    ))

  }

  # Reuse Step 4's already annotated results when available. This is optional:
  # Step 5 still runs from Step 2 alone if the workbook is absent.
  if (!is.null(step4_workbook) && file.exists(step4_workbook) &&
      requireNamespace("openxlsx", quietly = TRUE)) {
    sheet_names <- openxlsx::getSheetNames(step4_workbook)
    annotation_sheets <- intersect(
      c("Method_Agreement", "LinDA", "Maaslin2"),
      sheet_names
    )
    for (annotation_sheet in annotation_sheets) {
      annotated <- openxlsx::read.xlsx(step4_workbook, sheet = annotation_sheet)
      id_column <- intersect(c("Feature_ID", "feature", "function_id"), names(annotated))
      name_column <- intersect(c("pathway_name", "description"), names(annotated))
      if (length(id_column) > 0 && length(name_column) > 0) {
        fill_missing_names(data.frame(
          function_id = annotated[[id_column[[1]]]],
          pathway_name = annotated[[name_column[[1]]]],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  if (identical(pathway_level, "KEGG")) {
    # Current ggpicrust2 releases store the KEGG reference in long format,
    # including pathway IDs and names on every pathway-to-KO row.
    kegg_reference <- load_ggpicrust2_reference("ko_to_kegg_reference")
    if (all(c("pathway_id", "pathway_name") %in% names(kegg_reference))) {
      fill_missing_names(data.frame(
        function_id = kegg_reference$pathway_id,
        pathway_name = kegg_reference$pathway_name,
        stringsAsFactors = FALSE
      ))
    }

    missing_ids <- lookup$function_id[is.na(lookup$pathway_name)]
    if (length(missing_ids) > 0 && requireNamespace("KEGGREST", quietly = TRUE)) {
      live_names <- tryCatch({
        listed <- KEGGREST::keggList("pathway")
        data.frame(
          function_id = sub("^map", "ko", sub("^path:", "", names(listed))),
          pathway_name = sub(" - [^;]+$", "", as.character(listed)),
          stringsAsFactors = FALSE
        )
      }, error = function(error) {
        warning(
          "KEGG pathway names could not be refreshed from KEGG REST: ",
          conditionMessage(error),
          ". IDs without a local annotation will remain explicitly labeled as unavailable."
        )
        NULL
      })
      fill_missing_names(live_names)
    }
  }

  lookup$pathway_name[is.na(lookup$pathway_name) |
                        !nzchar(lookup$pathway_name)] <- "Name unavailable"
  lookup
}

# Insert explicit zero contributions for every sample/function/taxon
# combination supported by a function's observed taxa. PICRUSt2's stratified
# tables generally omit zero rows; without this completion, mean() and sd()
# would summarize only samples where a taxon contributed and would therefore
# overstate condition-level contribution.
complete_contribution_grid <- function(contributions, sample_ids) {
  required <- c("sample", "function_id", "taxon_label", "contribution")
  missing <- setdiff(required, names(contributions))
  if (length(missing) > 0L) {
    stop("Contribution grid input is missing columns: ",
         paste(missing, collapse = ", "))
  }
  sample_ids <- unique(as.character(sample_ids))
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  if (length(sample_ids) == 0L) {
    stop("At least one non-missing sample ID is required to complete contributions.")
  }

  contributions %>%
    dplyr::filter(sample %in% sample_ids) %>%
    dplyr::group_by(function_id) %>%
    tidyr::complete(
      sample = sample_ids,
      taxon_label,
      fill = list(contribution = 0)
    ) %>%
    dplyr::ungroup()
}

# Aggregate raw per-ASV contributions to a requested DADA2 taxonomy rank.
# Taxa beyond the top-N contributors for each pathway are retained as "Other"
# in the returned table for complete accounting; the heatmap helper removes
# that aggregate before calculating displayed relative percentages.
aggregate_contributions_by_taxonomic_rank <- function(contributions,
                                                      taxonomy_table,
                                                      taxonomic_level,
                                                      top_n_taxa,
                                                      sample_ids) {
  required <- c("sample", "function_id", "taxon", "contribution")
  missing <- setdiff(required, names(contributions))
  if (length(missing) > 0) {
    stop("Taxonomic-rank aggregation input is missing columns: ",
         paste(missing, collapse = ", "))
  }
  if (is.null(taxonomy_table)) {
    stop(
      "The ", taxonomic_level, " heatmap requires data/asv_taxonomy.csv. ",
      "Supply a DADA2-compatible taxonomy table before running Step 5."
    )
  }
  if (!("ASV_ID" %in% names(taxonomy_table))) {
    stop("data/asv_taxonomy.csv must contain an ASV_ID column.")
  }
  if (!(taxonomic_level %in% names(taxonomy_table))) {
    stop("Taxonomy level '", taxonomic_level,
         "' is not present in data/asv_taxonomy.csv.")
  }

  rank_lookup <- taxonomy_table %>%
    dplyr::transmute(
      taxon = as.character(ASV_ID),
      taxon_label = as.character(.data[[taxonomic_level]])
    ) %>%
    dplyr::distinct(taxon, .keep_all = TRUE)

  ranked <- contributions %>%
    dplyr::left_join(rank_lookup, by = "taxon") %>%
    dplyr::mutate(
      taxon_label = ifelse(
        is.na(taxon_label) | !nzchar(taxon_label),
        paste("Unclassified", taxonomic_level),
        taxon_label
      )
    ) %>%
    dplyr::group_by(sample, function_id, taxon_label) %>%
    dplyr::summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop")

  ranked <- complete_contribution_grid(ranked, sample_ids)

  top_taxa <- ranked %>%
    dplyr::group_by(function_id, taxon_label) %>%
    dplyr::summarise(total = sum(contribution, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(function_id) %>%
    dplyr::slice_max(total, n = top_n_taxa, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(function_id, taxon_label) %>%
    dplyr::mutate(is_top = TRUE)

  ranked %>%
    dplyr::left_join(top_taxa, by = c("function_id", "taxon_label")) %>%
    dplyr::mutate(taxon_label = ifelse(is.na(is_top), "Other", taxon_label)) %>%
    dplyr::group_by(sample, function_id, taxon_label) %>%
    dplyr::summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(taxonomic_level = taxonomic_level)
}

# Convert per-sample contributions to group-wise relative contribution values
# suitable for a condition-split heatmap. "Other" is excluded before the
# displayed taxonomic groups are renormalized within each pathway and condition.
build_condition_heatmap_data <- function(contributions_with_group, group_column,
                                         pathway_lookup) {
  required <- c("function_id", "taxon_label", "contribution", group_column)
  missing <- setdiff(required, names(contributions_with_group))
  if (length(missing) > 0) {
    stop("Condition heatmap input is missing columns: ", paste(missing, collapse = ", "))
  }

  heatmap_data <- contributions_with_group %>%
    dplyr::filter(taxon_label != "Other") %>%
    dplyr::group_by(function_id, taxon_label, .data[[group_column]]) %>%
    dplyr::summarise(
      mean_contribution = mean(contribution, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(function_id, .data[[group_column]]) %>%
    dplyr::mutate(
      displayed_total = sum(mean_contribution, na.rm = TRUE),
      relative_contribution = ifelse(
        displayed_total > 0,
        100 * mean_contribution / displayed_total,
        0
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(pathway_lookup, by = "function_id") %>%
    dplyr::mutate(
      pathway_name = dplyr::coalesce(pathway_name, "Name unavailable"),
      pathway_label = ifelse(
        pathway_name == "Name unavailable",
        function_id,
        paste0(pathway_name, " (", function_id, ")")
      )
    )

  if (nrow(heatmap_data) == 0) {
    stop("No displayed taxon contributions were available for the interactive heatmap.")
  }
  heatmap_data
}

# Create one interactive heatmap with a side-by-side panel for every condition.
# A shared white-yellow-red continuous scale makes percentages comparable
# between panels. Pathway names remain above the cells, while condition labels,
# the main plot title, and a plain-language explanation of the top-taxa cutoff
# are placed below the heatmap. Stable IDs stay available in hover text and
# result tables.
build_condition_heatmap_widget <- function(heatmap_data, group_column,
                                           title, taxonomic_level,
                                           top_n_taxa) {
  if (!is.numeric(top_n_taxa) || length(top_n_taxa) != 1L ||
      is.na(top_n_taxa) || !is.finite(top_n_taxa) || top_n_taxa < 1) {
    stop("top_n_taxa must be one positive number.")
  }

  group_levels <- unique(as.character(heatmap_data[[group_column]]))
  group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels)]
  if (length(group_levels) == 0) stop("No non-missing condition levels were available.")

  pathway_order <- unique(heatmap_data$function_id)
  pathway_lookup <- heatmap_data %>%
    dplyr::distinct(function_id, pathway_name)
  pathway_labels <- pathway_lookup$pathway_name[
    match(pathway_order, pathway_lookup$function_id)
  ]
  pathway_labels <- make.unique(pathway_labels, sep = " ")
  taxon_order <- heatmap_data %>%
    dplyr::group_by(taxon_label) %>%
    dplyr::summarise(total = sum(relative_contribution, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(total) %>%
    dplyr::pull(taxon_label)
  shared_max <- max(heatmap_data$relative_contribution, na.rm = TRUE)
  if (!is.finite(shared_max) || shared_max <= 0) shared_max <- 100

  panel_plots <- lapply(seq_along(group_levels), function(index) {
    group_value <- group_levels[[index]]
    panel_data <- heatmap_data[
      as.character(heatmap_data[[group_column]]) == group_value,
      , drop = FALSE
    ]
    grid <- expand.grid(
      taxon_label = taxon_order,
      function_id = pathway_order,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::left_join(
        panel_data %>%
          dplyr::select(function_id, taxon_label, mean_contribution,
                        relative_contribution, pathway_name),
        by = c("function_id", "taxon_label")
      ) %>%
      dplyr::mutate(
        mean_contribution = dplyr::coalesce(mean_contribution, 0),
        relative_contribution = dplyr::coalesce(relative_contribution, 0),
        pathway_name = dplyr::coalesce(pathway_name, "Name unavailable"),
        pathway_label = pathway_labels[match(function_id, pathway_order)],
          hover_text = paste0(
            "<b>", taxon_label, "</b><br>",
            "Pathway: ", pathway_name, "<br>",
            group_column, ": ", group_value, "<br>",
            "Relative contribution: ", sprintf("%.1f%%", relative_contribution), "<br>",
            "Mean predicted contribution: ", format(mean_contribution, digits = 5)
        )
      )

    z_matrix <- xtabs(
      relative_contribution ~ taxon_label + pathway_label,
      data = grid,
      drop.unused.levels = FALSE
    )[taxon_order, pathway_labels, drop = FALSE]
    text_matrix <- matrix(
      grid$hover_text[
        match(
          paste(rep(taxon_order, times = length(pathway_order)),
                rep(pathway_order, each = length(taxon_order)), sep = "\r"),
          paste(grid$taxon_label, grid$function_id, sep = "\r")
        )
      ],
      nrow = length(taxon_order),
      ncol = length(pathway_order),
      dimnames = list(taxon_order, pathway_labels)
    )

    plotly::plot_ly(
      x = pathway_labels,
      y = taxon_order,
      z = unclass(z_matrix),
      text = text_matrix,
      type = "heatmap",
      hovertemplate = "%{text}<extra></extra>",
      colors = c("#FFFFFF", "#FFD54F", "#D73027"),
      zmin = 0,
      zmax = shared_max,
      xgap = 1,
      ygap = 1,
      showscale = index == length(group_levels),
      colorbar = list(title = list(text = "Relative<br>contribution (%)"))
    ) %>%
      plotly::layout(
        xaxis = list(
          title = list(text = ""),
          side = "top",
          tickangle = -35,
          automargin = TRUE
        ),
        yaxis = list(
          title = list(
            text = if (index == 1) taxonomic_level else "",
            standoff = if (index == 1) 50 else 0
          ),
          automargin = TRUE
        )
      )
  })

  combined_widget <- plotly::subplot(
    panel_plots,
    nrows = 1,
    shareY = TRUE,
    titleX = TRUE,
    titleY = TRUE,
    margin = 0.005
  )

  # Compact rows keep pathway-specific top-taxa unions readable without pushing
  # the title/subtitle footer below a typical browser viewport.
  cell_height_pixels <- 22
  plot_height <- max(850, cell_height_pixels * length(taxon_order) + 540)
  plot_margins <- list(l = 230, r = 120, b = 60, t = 300)
  plot_inner_height <- plot_height - plot_margins$t - plot_margins$b
  footer_height_pixels <- 150
  footer_fraction <- footer_height_pixels / plot_inner_height
  combined_widget$x$layout$yaxis$domain <- c(footer_fraction, 1)

  shared_y_domain <- combined_widget$x$layout$yaxis$domain
  panel_outline_shapes <- lapply(seq_along(group_levels), function(index) {
    xaxis_name <- if (index == 1) "xaxis" else paste0("xaxis", index)
    x_domain <- combined_widget$x$layout[[xaxis_name]]$domain
    list(
      type = "rect",
      x0 = x_domain[[1]],
      x1 = x_domain[[2]],
      y0 = shared_y_domain[[1]],
      y1 = shared_y_domain[[2]],
      xref = "paper",
      yref = "paper",
      line = list(color = "#D3D3D3", width = 1.5),
      fillcolor = "rgba(0,0,0,0)",
      layer = "above"
    )
  })
  # subplot() initializes an empty shapes list that otherwise overrides shapes
  # supplied by the final layout() call during plotly_build().
  combined_widget$x$layout$shapes <- NULL

  # Reserve a fixed-height footer inside the Plotly paper area so condition
  # labels, titles, and subtitles remain visible for both short Phylum and tall
  # Genus heatmaps.
  panel_centers <- (seq_along(group_levels) - 0.5) / length(group_levels)
  condition_annotations <- lapply(seq_along(group_levels), function(index) {
    list(
      text = paste0("<b>", group_levels[[index]], "</b>"),
      x = panel_centers[[index]],
      y = 116 / plot_inner_height,
      xref = "paper",
      yref = "paper",
      xanchor = "center",
      yanchor = "bottom",
      showarrow = FALSE
    )
  })
  below_plot_annotations <- c(
    condition_annotations,
    list(list(
      text = paste0("<b>", title, "</b>"),
      x = 0.5,
      y = 82 / plot_inner_height,
      xref = "paper",
      yref = "paper",
      xanchor = "center",
      yanchor = "bottom",
      showarrow = FALSE,
      font = list(size = 16)
    )),
    list(list(
      text = paste0(
        "Top ", top_n_taxa,
        " taxa are selected separately for each pathway.<br>",
        "Within each pathway, taxa are ranked highest to lowest by predicted contribution summed across all samples.<br>",
        "The combined plot may therefore contain more than ", top_n_taxa,
        " taxa in total."
      ),
      x = 0.5,
      y = 10 / plot_inner_height,
      xref = "paper",
      yref = "paper",
      xanchor = "center",
      yanchor = "bottom",
      align = "center",
      showarrow = FALSE,
      font = list(size = 12, color = "#555555")
    ))
  )

  suppressWarnings(
    plotly::layout(
      combined_widget,
      title = list(text = ""),
      annotations = below_plot_annotations,
      shapes = panel_outline_shapes,
      plot_bgcolor = "#D3D3D3",
      paper_bgcolor = "#FFFFFF",
      height = plot_height,
      margin = plot_margins
    )
  )
}

# Remove obsolete Step 5 plot files so the output folder contains only the
# heatmaps created by the current run.
remove_stale_taxon_plot_files <- function(plots_folder) {
  stale <- list.files(
    plots_folder,
    pattern = "^taxon_contribution_.*\\.(png|html)$",
    full.names = TRUE
  )
  if (length(stale) > 0) unlink(stale, recursive = TRUE, force = TRUE)
  invisible(stale)
}

# Keep the raw ASV identity visible even when a higher taxonomy rank is used.
# If the selected taxonomy label already contains the ASV ID (for example the
# DADA2 Unique_Tax column), do not append it a second time.
build_taxon_display_labels <- function(taxonomy_labels, asv_ids) {
  mapply(function(label, asv_id) {
    label <- as.character(label)
    asv_id <- as.character(asv_id)
    if (is.na(label) || !nzchar(label)) return(asv_id)
    if (grepl(asv_id, label, fixed = TRUE)) return(label)
    paste0(label, " (", asv_id, ")")
  }, taxonomy_labels, asv_ids, USE.NAMES = FALSE)
}

# Confirm that every Step 5 heatmap includes both required annotations before
# it is written. This guards MetaCyc and KEGG plots equally because they share
# the same save path.
validate_taxon_heatmap_annotations <- function(widget) {
  built_widget <- plotly::plotly_build(widget)
  annotations <- built_widget$x$layout$annotations
  if (is.null(annotations)) annotations <- list()
  annotation_text <- vapply(
    annotations,
    function(annotation) {
      if (is.null(annotation$text)) "" else as.character(annotation$text)
    },
    character(1)
  )

  has_main_title <- any(grepl(
    "Relative Taxon Contribution",
    annotation_text,
    fixed = TRUE
  ))
  top_taxa_subtitle_match <-
    grepl(
      "taxa are selected separately for each pathway",
      annotation_text,
      fixed = TRUE
    ) &
      grepl(
        "Within each pathway, taxa are ranked highest to lowest by predicted contribution summed across all samples",
        annotation_text,
        fixed = TRUE
      ) &
      grepl(
        "The combined plot may therefore contain more than",
        annotation_text,
        fixed = TRUE
      )
  has_top_taxa_subtitle <- any(top_taxa_subtitle_match)
  heatmap_bottom <- built_widget$x$layout$yaxis$domain[[1]]
  subtitle_inside_plotting_area <- has_top_taxa_subtitle &&
    any(vapply(
      annotations[top_taxa_subtitle_match],
      function(annotation) {
        is.numeric(annotation$y) && length(annotation$y) == 1L &&
          annotation$y > 0 && annotation$y < heatmap_bottom &&
          identical(annotation$yref, "paper") &&
          identical(annotation$yanchor, "bottom")
      },
      logical(1)
    ))

  if (!has_main_title || !has_top_taxa_subtitle ||
      !subtitle_inside_plotting_area) {
    missing_annotations <- c(
      if (!has_main_title) "main title",
      if (!has_top_taxa_subtitle) "top-taxa subtitle",
      if (has_top_taxa_subtitle && !subtitle_inside_plotting_area) {
        "top-taxa subtitle placement inside the plotting area"
      }
    )
    stop(
      "Step 5 heatmap is missing its required ",
      paste(missing_annotations, collapse = " and "),
      "."
    )
  }
  invisible(widget)
}

# Save a standalone widget when Pandoc is available. On systems where Pandoc
# is not on PATH, save a normal HTML file plus its companion dependency folder
# instead of failing the whole notebook after the plots have been computed.
save_interactive_taxon_widget <- function(widget, path) {
  validate_taxon_heatmap_annotations(widget)
  pandoc_available <- requireNamespace("rmarkdown", quietly = TRUE) &&
    rmarkdown::pandoc_available()
  if (!pandoc_available) {
    warning(
      "Pandoc was not detected; saving the interactive chart with a companion ",
      "dependency folder rather than as one self-contained HTML file."
    )
  }
  htmlwidgets::saveWidget(widget, path, selfcontained = pandoc_available)
  invisible(path)
}

# Run the KEGG branch of Step 5 end to end.
#
# The function selects or validates KEGG pathways, maps them to KOs with the
# bundled ggpicrust2 reference, stream-filters the compressed PICRUSt2 KO
# contribution table, aggregates KO contributions to pathways, completes
# absent sample/taxon combinations as zeros, creates rank-level heatmaps, and
# writes a self-documenting KEGG workbook. `group_column` names the metadata
# condition field; `top_n_taxa` is applied independently within each pathway.
# `step4_workbook` is annotation-only and optional: Step 5 does not depend on a
# significant Step 4 result and does not filter contributions by p-value.
run_kegg_taxon_contribution <- function(sample_metadata, group_column, top_n_taxa,
                                        community_run_dir, contribution_run_dir,
                                        output_root_folder,
                                        pathway_ids = character(0),
                                        auto_select_n_features = 10,
                                        taxonomy_table = NULL, tax_level = "Genus",
                                        heatmap_tax_levels = c("Genus", "Phylum"),
                                        step4_workbook = NULL,
                                        contribution_run_summary = NULL) {
  output_folder <- here::here(output_root_folder, "KEGG")
  plots_folder <- here::here(output_folder, "plots")
  dir.create(plots_folder, recursive = TRUE, showWarnings = FALSE)
  remove_stale_taxon_plot_files(plots_folder)

  # Accept the common KEGG pathway ID forms used by KEGG and ggpicrust2.
  pathway_ids <- unique(as.character(pathway_ids))
  pathway_ids <- pathway_ids[!is.na(pathway_ids) & nzchar(pathway_ids)]
  pathway_ids <- sub("^map", "ko", pathway_ids)
  pathway_ids <- ifelse(grepl("^[0-9]{5}$", pathway_ids),
                        paste0("ko", pathway_ids), pathway_ids)

  # If no explicit pathways were requested, build the same KEGG pathway
  # abundance table used elsewhere in the PICRUSt2 ecosystem and select the
  # most abundant pathways. This keeps the contribution
  # notebook self-contained: its automatic selection depends only on Step 2.
  if (length(pathway_ids) == 0) {
    ko_abundance_path <- resolve_picrust2_table_path(here::here(
      community_run_dir, "KO_metagenome_out", "pred_metagenome_unstrat.tsv"
    ))
    if (is.na(ko_abundance_path)) {
      stop(
        "No KEGG pathway IDs were supplied and the Step 2 KO abundance table ",
        "was not found: ", ko_abundance_path
      )
    }

    ko_abundance <- readr::read_tsv(ko_abundance_path, show_col_types = FALSE)
    ko_abundance[[1]] <- sub("^ko:", "", ko_abundance[[1]])
    ko_input_path <- tempfile(fileext = ".tsv")
    on.exit(unlink(ko_input_path), add = TRUE)
    readr::write_tsv(ko_abundance, ko_input_path)
    kegg_abundance <- ggpicrust2::ko2kegg_abundance(file = ko_input_path)

    pathway_totals <- rowSums(as.matrix(kegg_abundance), na.rm = TRUE)
    pathway_ids <- head(
      names(sort(pathway_totals, decreasing = TRUE)),
      auto_select_n_features
    )
  }

  if (length(pathway_ids) == 0) {
    stop(
      "No KEGG pathways were available for taxon contribution analysis. Supply one or ",
      "more kegg_pathway_ids_of_interest values in the notebook Configuration."
    )
  }

  pathway_lookup <- build_pathway_name_lookup(
    pathway_ids = pathway_ids,
    pathway_level = "KEGG",
    step4_workbook = step4_workbook
  )

  reference <- load_ggpicrust2_reference("ko_to_kegg_reference")
  if (all(c("pathway_id", "ko_id") %in% names(reference))) {
    # ggpicrust2 2.5.x long format: one pathway/KO pair per row.
    pathway_to_ko <- reference %>%
      dplyr::filter(pathway_id %in% pathway_ids) %>%
      dplyr::transmute(
        function_id = as.character(pathway_id),
        ko_id = as.character(ko_id)
      ) %>%
      dplyr::filter(!is.na(ko_id), nzchar(ko_id)) %>%
      dplyr::distinct()
    pathway_to_ko$ko_id <- ifelse(
      grepl("^ko:", pathway_to_ko$ko_id),
      pathway_to_ko$ko_id,
      paste0("ko:", pathway_to_ko$ko_id)
    )
  } else {
    # Backward compatibility for older wide references where the first
    # column is a pathway ID and the remaining cells contain KO IDs.
    pathway_to_ko <- dplyr::bind_rows(lapply(pathway_ids, function(pathway_id) {
      row <- reference[reference[[1]] == pathway_id, -1, drop = FALSE]
      kos <- unique(as.character(unlist(row, use.names = FALSE)))
      kos <- kos[!is.na(kos) & nzchar(kos)]
      data.frame(
        function_id = pathway_id,
        ko_id = ifelse(grepl("^ko:", kos), kos, paste0("ko:", kos)),
        stringsAsFactors = FALSE
      )
    }))
  }
  if (nrow(pathway_to_ko) == 0)
    stop("No KO mappings were found for the selected KEGG pathways: ", paste(pathway_ids, collapse = ", "))
  pathway_to_ko <- pathway_to_ko %>%
    dplyr::left_join(pathway_lookup, by = "function_id") %>%
    dplyr::select(function_id, pathway_name, ko_id)

  ko_ids <- unique(pathway_to_ko$ko_id)
  if (any(grepl(",", ko_ids, fixed = TRUE))) stop("Unexpected comma in a KO identifier.")
  contribution_path <- here::here(
    contribution_run_dir, "KO_metagenome_out", "pred_metagenome_contrib.tsv.gz"
  )
  if (!file.exists(contribution_path)) {
    stop("The Step 5 KO contribution table was not found: ", contribution_path)
  }
  awk_command <- build_contribution_filter_command(contribution_path, ko_ids)
  raw <- as.data.frame(data.table::fread(cmd = awk_command, sep = "\t", header = TRUE))
  if (nrow(raw) == 0) stop("No stratified KO rows matched the selected KEGG pathways.")

  contributions <- raw %>%
    dplyr::transmute(sample = as.character(sample), ko_id = as.character(`function`),
                     taxon = as.character(taxon),
                     contribution = as.numeric(taxon_function_abun)) %>%
    dplyr::inner_join(pathway_to_ko, by = "ko_id") %>%
    dplyr::select(sample, function_id, taxon, contribution) %>%
    dplyr::group_by(sample, function_id, taxon) %>%
    dplyr::summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop")

  rank_source_contributions <- contributions
  contributions <- contributions %>%
    dplyr::mutate(taxon_label = taxon)

  if (!is.null(taxonomy_table)) {
    if (!("ASV_ID" %in% names(taxonomy_table))) {
      stop("data/asv_taxonomy.csv must contain an ASV_ID column.")
    }
    taxonomy_id_column <- "ASV_ID"
    if (!(tax_level %in% names(taxonomy_table)))
      stop("tax_level ('", tax_level, "') is not present in data/asv_taxonomy.csv.")
    taxonomy_lookup <- taxonomy_table %>%
      dplyr::transmute(taxon = as.character(.data[[taxonomy_id_column]]),
                       named_taxon = as.character(.data[[tax_level]]))
    contributions <- contributions %>%
      dplyr::left_join(taxonomy_lookup, by = "taxon") %>%
      dplyr::mutate(taxon_label = build_taxon_display_labels(named_taxon, taxon)) %>%
      dplyr::select(-named_taxon, -taxon) %>%
      dplyr::group_by(sample, function_id, taxon_label) %>%
      dplyr::summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop")
  } else {
    contributions <- contributions %>% dplyr::select(-taxon)
  }

  contributions <- complete_contribution_grid(
    contributions,
    sample_ids = sample_metadata$SampleID
  )

  top_taxa <- contributions %>%
    dplyr::group_by(function_id, taxon_label) %>%
    dplyr::summarise(total = sum(contribution, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(function_id) %>%
    dplyr::slice_max(total, n = top_n_taxa, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(function_id, taxon_label) %>%
    dplyr::mutate(is_top = TRUE)
  aggregated <- contributions %>%
    dplyr::left_join(top_taxa, by = c("function_id", "taxon_label")) %>%
    dplyr::mutate(taxon_label = ifelse(is.na(is_top), "Other", taxon_label)) %>%
    dplyr::group_by(sample, function_id, taxon_label) %>%
    dplyr::summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop") %>%
    dplyr::left_join(pathway_lookup, by = "function_id") %>%
    dplyr::select(sample, function_id, pathway_name, taxon_label, contribution)

  with_group <- aggregated %>%
    dplyr::left_join(sample_metadata %>% dplyr::select(SampleID, dplyr::all_of(group_column)),
                     by = c("sample" = "SampleID"))
  group_summary <- with_group %>%
    dplyr::group_by(function_id, pathway_name, taxon_label, .data[[group_column]]) %>%
    dplyr::summarise(mean_contribution = mean(contribution),
                     sd_contribution = stats::sd(contribution), n_samples = dplyr::n(),
                     .groups = "drop")

  heatmap_results <- stats::setNames(lapply(heatmap_tax_levels, function(level) {
    ranked <- aggregate_contributions_by_taxonomic_rank(
      contributions = rank_source_contributions,
      taxonomy_table = taxonomy_table,
      taxonomic_level = level,
      top_n_taxa = top_n_taxa,
      sample_ids = sample_metadata$SampleID
    ) %>%
      dplyr::left_join(
        sample_metadata %>% dplyr::select(SampleID, dplyr::all_of(group_column)),
        by = c("sample" = "SampleID")
      )
    heatmap_data <- build_condition_heatmap_data(
      contributions_with_group = ranked,
      group_column = group_column,
      pathway_lookup = pathway_lookup
    )
    heatmap_widget <- build_condition_heatmap_widget(
      heatmap_data = heatmap_data,
      group_column = group_column,
      title = paste0(level, "-level Relative Taxon Contribution — KEGG Pathways"),
      taxonomic_level = level,
      top_n_taxa = top_n_taxa
    )
    heatmap_path <- here::here(
      plots_folder,
      paste0("taxon_contribution_heatmap_kegg_", tolower(level), ".html")
    )
    save_interactive_taxon_widget(heatmap_widget, heatmap_path)
    list(path = heatmap_path, widget = heatmap_widget)
  }), heatmap_tax_levels)

  workbook <- here::here(output_folder, "taxon_contribution_summary.xlsx")
  if (file.exists(workbook)) unlink(workbook)
  if (!is.null(contribution_run_summary)) {
    add_sheet_to_excel(workbook, "PICRUSt2_Contribution_Run",
                       as.data.frame(contribution_run_summary),
                       rownames = FALSE, overwrite = TRUE)
  }
  add_sheet_to_excel(workbook, "Metadata", as.data.frame(sample_metadata), rownames = FALSE, overwrite = TRUE)
  add_sheet_to_excel(workbook, "Pathway_KO_Mapping", pathway_to_ko, rownames = FALSE, overwrite = TRUE)
  add_sheet_to_excel(workbook, "Taxon_Contributions_per_Sample", with_group, rownames = FALSE, overwrite = TRUE)
  add_sheet_to_excel(workbook, "Contribution_by_Condition", group_summary, rownames = FALSE, overwrite = TRUE)

  per_sample_descriptions <- c(
    sample = "Sample identifier.",
    function_id = "KEGG pathway identifier.",
    pathway_name = "Human-readable KEGG pathway name.",
    taxon_label = "Contributing ASV with its selected taxonomy label, or Other in the exported table.",
    contribution = "Summed contribution of this taxon's constituent KOs to the pathway."
  )
  per_sample_descriptions[group_column] <- paste0(
    "Sample's ", group_column, " level, joined from data/metadata.tsv."
  )
  group_summary_descriptions <- c(
    function_id = "KEGG pathway identifier.",
    pathway_name = "Human-readable KEGG pathway name.",
    taxon_label = "Contributing ASV with its selected taxonomy label, or Other in the exported table.",
    mean_contribution = "Mean predicted contribution across all samples in this condition, including explicit zeros when the taxon did not contribute.",
    sd_contribution = "Sample standard deviation of predicted contribution across all samples in this condition, including explicit zeros.",
    n_samples = "Number of samples summarized in this condition."
  )
  group_summary_descriptions[group_column] <- paste0(group_column, " level being summarized.")

  dictionary_parts <- list()
  if (!is.null(contribution_run_summary)) {
    dictionary_parts[[length(dictionary_parts) + 1]] <- build_column_dictionary(
      "PICRUSt2_Contribution_Run", as.data.frame(contribution_run_summary),
      c(
        Stage = "PICRUSt2 contribution-generation stage.",
        `Run Status` = "Whether Step 5 executed this stage or reused an already complete, verified contribution output.",
        Command = "Fully resolved env/conda command configured for the stage; executed when Run Status says the stage ran.",
        `Exit Code` = "Exit code reported by the command.",
        `Verified Output` = "TRUE when the expected output file was found after a zero exit code.",
        `Duration (minutes)` = "Wall-clock duration in minutes.",
        `Output File` = "Expected output file produced by the stage.",
        `PICRUSt2 Version` = "PICRUSt2 version used for contribution generation."
      ),
      workbook_path = workbook, rownames = FALSE
    )
  }
  dictionary_parts <- c(dictionary_parts, list(
    build_column_dictionary(
      "Metadata", as.data.frame(sample_metadata),
      stats::setNames(paste("Column from data/metadata.tsv:", names(sample_metadata)),
                      names(sample_metadata)), workbook_path = workbook, rownames = FALSE
    ),
    build_column_dictionary(
      "Pathway_KO_Mapping", pathway_to_ko,
      c(function_id = "Selected KEGG pathway identifier.",
        pathway_name = "Human-readable KEGG pathway name.",
        ko_id = "Constituent KO identifier used to aggregate contributions."),
      workbook_path = workbook, rownames = FALSE
    ),
    build_column_dictionary(
      "Taxon_Contributions_per_Sample", with_group,
      per_sample_descriptions,
      workbook_path = workbook, rownames = FALSE
    ),
    build_column_dictionary(
      "Contribution_by_Condition", group_summary,
      group_summary_descriptions,
      workbook_path = workbook, rownames = FALSE
    )
  ))
  dictionary <- dplyr::bind_rows(dictionary_parts)
  add_sheet_to_excel(workbook, "Column_Dictionary", dictionary,
                     rownames = FALSE, overwrite = TRUE)
  list(pathway_ids = pathway_ids, pathway_lookup = pathway_lookup,
       workbook = workbook,
       heatmaps = vapply(heatmap_results, `[[`, character(1), "path"),
       heatmap_widgets = lapply(heatmap_results, `[[`, "widget"),
       n_rows = nrow(aggregated))
}
