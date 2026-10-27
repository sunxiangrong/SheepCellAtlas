#!/usr/bin/env Rscript
# Grouped, flow-conserving GO alluvial plot. Each displayed gene-GO edge occupies
# one equally sized slot at both ends; BP, CC and MF remain in one coordinate system.
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: plot_AD_GO_alluvial_enrichment_V6.R <cross_species_base>")
base <- normalizePath(args[[1]], mustWork = TRUE)
out <- file.path(base, "plots", "strict_matched_scatter_20260721", "atopic_dermatitis_skin_softhorn", "GO_sankey_publication_V7")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

ratio_num <- function(x) as.numeric(sub("/.*", "", x)) / as.numeric(sub(".*/", "", x))
res_dir <- file.path(base, "results", "strict_matched_scatter_20260721", "atopic_dermatitis_skin_softhorn", "enrichment_conserved_disease_deg")
read_go <- function(component, ontology) {
  read_tsv(file.path(res_dir, paste0(component, ".GO_", ontology, ".tsv")), show_col_types = FALSE) |>
    mutate(component = component, ontology = ontology)
}
terms_all <- bind_rows(
  read_go("Spinous_cells", "BP"), read_go("Granular_cells", "BP"),
  read_go("Spinous_cells", "CC"), read_go("Granular_cells", "CC"),
  read_go("Follicular_cells", "MF"), read_go("Spinous_cells", "MF")
)

# Fixed order only. Vertical centres are calculated from real displayed edge counts below.
wanted <- tibble(
  display_order = seq_len(11),
  component = c("Spinous_cells", "Granular_cells", "Granular_cells",
                "Spinous_cells", "Spinous_cells", "Spinous_cells", "Spinous_cells",
                "Spinous_cells", "Spinous_cells", "Spinous_cells", "Follicular_cells"),
  description = c("chromatin remodeling", "mRNA processing", "regulation of mRNA splicing, via spliceosome",
                  "nuclear ribonucleoprotein granule", "nuclear speck", "proteasome complex", "U12-type spliceosomal complex",
                  "transcription coregulator activity", "chromatin binding", "cadherin binding", "mRNA binding"),
  ontology = c("BP", "BP", "BP", "CC", "CC", "CC", "CC", "MF", "MF", "MF", "MF"),
  # Focus terms use the strongest shade; background terms use muted shades.
  term_colour = c("#3F7F8F", "#5D9FA2", "#A8C9CA", "#77629B", "#AD9BC2", "#C1B4CF", "#D6CDE0",
                  "#A46C42", "#BE8B60", "#C9823C", "#DFB08B")
)
focus_terms <- c("cadherin binding", "chromatin remodeling", "mRNA processing", "nuclear ribonucleoprotein granule")
row_map <- wanted |>
  left_join(terms_all |> select(component, description = Description, ontology, ID, GeneRatio, Count, p.adjust, geneID), by = c("component", "description", "ontology")) |>
  mutate(
    GeneRatio_num = ratio_num(GeneRatio),
    logFDR = -log10(p.adjust),
    is_focus = description %in% focus_terms,
    ontology = factor(ontology, levels = c("BP", "CC", "MF"))
  ) |>
  # Keep ontology blocks together and enforce a reproducible GeneRatio order.
  arrange(ontology, desc(GeneRatio_num), desc(Count), description) |>
  mutate(
    display_order = row_number(),
    ontology = as.character(ontology),
    display_label = if_else(
      description == "regulation of mRNA splicing, via spliceosome",
      "regulation of mRNA splicing",
      description
    )
  )
if (any(is.na(row_map$ID))) stop("Required GO terms are absent from the source table.")
ratio_order_ok <- row_map |>
  group_by(ontology) |>
  summarise(ok = all(diff(GeneRatio_num) <= 1e-12), .groups = "drop")
if (any(!ratio_order_ok$ok)) stop("GO rows are not ordered by decreasing GeneRatio within ontology.")

pair_dir <- file.path(base, "results", "strict_matched_scatter_20260721", "atopic_dermatitis_skin_softhorn")
meta <- bind_rows(lapply(c("Spinous_cells", "Granular_cells"), function(component) {
  read_tsv(file.path(pair_dir, paste0("GSE197023_AD_Keratinocytes__sheep_", component, ".conserved_disease_deg.tsv")), show_col_types = FALSE) |>
    transmute(gene = as.character(gene), disease_FDR = as.numeric(disease_FDR), disease_log2FC = as.numeric(disease_log2FC))
})) |>
  filter(!is.na(disease_FDR), !is.na(disease_log2FC)) |>
  group_by(gene) |>
  summarise(disease_FDR = min(disease_FDR), disease_log2FC = disease_log2FC[which.min(disease_FDR)][1], .groups = "drop")

all_focus_edges <- row_map |>
  filter(is_focus) |>
  select(description, display_order, term_colour, ontology, geneID) |>
  separate_rows(geneID, sep = "/") |>
  rename(gene = geneID) |>
  inner_join(meta, by = "gene") |>
  distinct(gene, description, display_order, term_colour, ontology, disease_FDR, disease_log2FC)

# Stratified, reproducible gene selection: preserve six representative genes per focus term,
# then add multi-term genes and fill to 28 using real FDR/effect-size ordering.
term_priority <- all_focus_edges |>
  arrange(description, disease_FDR, desc(abs(disease_log2FC)), gene) |>
  group_by(description) |>
  slice_head(n = 6) |>
  ungroup() |>
  pull(gene) |> unique()
gene_summary_all <- all_focus_edges |>
  group_by(gene) |>
  summarise(n_terms = n(), best_fdr = min(disease_FDR), best_lfc = disease_log2FC[which.min(disease_FDR)][1], .groups = "drop") |>
  arrange(desc(n_terms), best_fdr, desc(abs(best_lfc)), gene)
selected_genes <- unique(c(term_priority, gene_summary_all$gene)) |> head(28)
focus_edges <- all_focus_edges |> filter(gene %in% selected_genes)
coverage_check <- focus_edges |> count(description, name = "n")
if (n_distinct(focus_edges$gene) != 28L) stop("Expected exactly 28 displayed genes.")
if (any(!focus_terms %in% coverage_check$description) ||
    any(coverage_check$n[match(focus_terms, coverage_check$description)] < 6L)) {
  stop("Each focus GO term must retain at least six displayed genes.")
}

# One common vertical unit controls source nodes, pathway nodes and ribbon width.
# The larger unit makes the flow field as visually full as the reference schematic.
unit_h <- .165
min_term_h <- .36
term_edge_counts <- focus_edges |> count(description, name = "n_edges")
row_map <- row_map |> left_join(term_edge_counts, by = "description") |>
  mutate(n_edges = coalesce(n_edges, 0L), node_height = pmax(n_edges * unit_h, min_term_h)) |>
  arrange(display_order)

# Dynamic pathway node bounds: fixed order, heights determined from the actual displayed edges.
cursor <- 9.0
term_bounds <- vector("list", nrow(row_map))
for (i in seq_len(nrow(row_map))) {
  h <- row_map$node_height[i]
  ymax <- cursor; ymin <- cursor - h
  term_bounds[[i]] <- tibble(description = row_map$description[i], ymin = ymin, ymax = ymax, y = (ymin + ymax) / 2)
  gap <- if (i < nrow(row_map) && row_map$ontology[i] != row_map$ontology[i + 1]) .30 else .12
  cursor <- ymin - gap
}
term_bounds <- bind_rows(term_bounds)
row_map <- row_map |> left_join(term_bounds, by = "description")

# Target pathway position determines ordering only. The gene stack spans above the
# first and below the last focus pathway, producing the broad reference-style fan.
focus_edges <- focus_edges |>
  left_join(row_map |> select(description, display_order, ymin, ymax, y),
            by = c("description", "display_order"))
gene_meta <- focus_edges |>
  group_by(gene) |>
  summarise(
    target_y = mean(y), n_edges = n(), best_fdr = min(disease_FDR),
    best_lfc = disease_log2FC[which.min(disease_FDR)][1],
    dominant_order = display_order[which.min(abs(y - mean(y)))][1],
    .groups = "drop"
  ) |>
  arrange(desc(target_y), dominant_order, desc(n_edges), best_fdr, desc(abs(best_lfc)), gene)

focus_top <- max(row_map$ymax[row_map$is_focus])
focus_bottom <- min(row_map$ymin[row_map$is_focus])
gene_top <- focus_top + .90
gene_bottom <- focus_bottom - .20
gene_span <- gene_top - gene_bottom
node_height_total <- sum(gene_meta$n_edges * unit_h)
module_break_after <- vapply(seq_len(nrow(gene_meta) - 1L), function(i) {
  gene_meta$dominant_order[i] != gene_meta$dominant_order[i + 1L]
}, logical(1))
module_gap_ratio <- 1.35
gap_weight_total <- sum(ifelse(module_break_after, module_gap_ratio, 1))
available_gap_height <- gene_span - node_height_total
if (available_gap_height <= 0) stop("Gene nodes are taller than the focus-pathway span.")
base_gene_gap <- available_gap_height / gap_weight_total
module_gene_gap <- base_gene_gap * module_gap_ratio

g_cursor <- gene_top
gene_bounds <- vector("list", nrow(gene_meta))
for (i in seq_len(nrow(gene_meta))) {
  h <- gene_meta$n_edges[i] * unit_h
  ymax <- g_cursor
  ymin <- ymax - h
  centre <- (ymin + ymax) / 2
  gene_bounds[[i]] <- tibble(
    gene = gene_meta$gene[i], gene_ymin = ymin,
    gene_ymax = ymax, gene_y = centre
  )
  if (i < nrow(gene_meta)) {
    gap_i <- if (module_break_after[i]) module_gene_gap else base_gene_gap
    g_cursor <- ymin - gap_i
  }
}
gene_meta <- gene_meta |> left_join(bind_rows(gene_bounds), by = "gene")
if (abs(max(gene_meta$gene_ymax) - gene_top) > 1e-8 ||
    abs(min(gene_meta$gene_ymin) - gene_bottom) > 1e-8) {
  stop("The expanded gene stack does not match its requested bounds.")
}

# Low target terms occupy low source slots and high target terms occupy high slots.
# This fixes the reversed slot ordering in v8 for multi-term genes.
focus_edges <- focus_edges |>
  left_join(gene_meta |> select(gene, gene_ymin, gene_ymax, gene_y), by = "gene") |>
  arrange(gene, y, description) |>
  group_by(gene) |>
  mutate(gene_slot_y = gene_ymin + (row_number() - .5) * unit_h) |>
  ungroup() |>
  # Match incoming slot order to source order to reduce crossings at term nodes.
  arrange(description, gene_y, gene_slot_y, gene) |>
  group_by(description) |>
  mutate(term_slot_y = ymin + (row_number() - .5) * unit_h) |>
  ungroup()

if (any(focus_edges$gene_slot_y <= focus_edges$gene_ymin) ||
    any(focus_edges$gene_slot_y >= focus_edges$gene_ymax)) {
  stop("A gene ribbon slot falls outside its source node.")
}
if (any(focus_edges$term_slot_y <= focus_edges$ymin) ||
    any(focus_edges$term_slot_y >= focus_edges$ymax)) {
  stop("A ribbon slot falls outside its pathway node.")
}
if (nrow(gene_meta) > 1L && any(head(gene_meta$gene_ymin, -1) <= tail(gene_meta$gene_ymax, -1))) {
  stop("Collision-aware gene placement produced overlapping source nodes.")
}

# Each source slot is coloured by its destination term; multi-term genes therefore
# become meaningful stacked nodes rather than arbitrary decorative colours.
gene_slots <- focus_edges |>
  transmute(
    gene, xmin = .190, xmax = .214,
    ymin = gene_slot_y - unit_h * .46,
    ymax = gene_slot_y + unit_h * .46,
    slot_colour = term_colour
  )

# Cubic Bezier boundaries leave/enter nodes horizontally and bend in the middle.
make_ribbon <- function(x0, x1, y0, y1, id, half_width = unit_h * .47, curve_pull = .42, n = 120) {
  t <- seq(0, 1, length.out = n)
  bezier <- function(p0, p1, p2, p3, tt) {
    (1 - tt)^3 * p0 + 3 * (1 - tt)^2 * tt * p1 + 3 * (1 - tt) * tt^2 * p2 + tt^3 * p3
  }
  dx <- x1 - x0
  xx <- bezier(x0, x0 + curve_pull * dx, x1 - curve_pull * dx, x1, t)
  y_top <- bezier(y0 + half_width, y0 + half_width, y1 + half_width, y1 + half_width, t)
  y_bottom <- bezier(y0 - half_width, y0 - half_width, y1 - half_width, y1 - half_width, t)
  tibble(x = c(xx, rev(xx)), y = c(y_top, rev(y_bottom)), ribbon = id)
}
ribbons <- bind_rows(lapply(seq_len(nrow(focus_edges)), function(i) {
  z <- focus_edges[i, ]
  make_ribbon(.214, .695, z$gene_slot_y, z$term_slot_y, paste(z$gene, z$description, sep = "__"))
}))

# A fixed real GeneRatio scale keeps zero exactly at the left edge of each dot frame.
dot_xmin <- .748; dot_xmax <- 1.000; gr_min <- 0; gr_max <- .085
row_map <- row_map |> mutate(dot_x = dot_xmin + (GeneRatio_num - gr_min) / (gr_max - gr_min) * (dot_xmax - dot_xmin))
fdr_min <- 1; fdr_max <- 3.2 # display cap only; raw FDR values remain in row_map/source data.
fdr_pal <- colorRampPalette(c("#F4B2A6", "#B42025"))(160)
row_map <- row_map |> mutate(fdr_index = pmax(1, pmin(160, round((logFDR - fdr_min) / (fdr_max - fdr_min) * 159 + 1))), fdr_col = fdr_pal[fdr_index])
ticks <- tibble(v = c(0, .02, .04, .06, .08)) |> mutate(x = dot_xmin + v / gr_max * (dot_xmax - dot_xmin))
count_min <- min(row_map$Count); count_max <- max(row_map$Count)
count_to_size <- function(x) 4.05 + (sqrt(x) - sqrt(count_min)) / (sqrt(count_max) - sqrt(count_min)) * (5.55 - 4.05)
row_map <- row_map |> mutate(bubble_size = count_to_size(Count))
count_legend <- tibble(Count = c(10, 30, 50), x = c(.765, .810, .855)) |>
  mutate(size = count_to_size(Count))
fdrbar_xmin <- .885; fdrbar_xmax <- 1.000
bar_dat <- tibble(v = seq(fdr_min, fdr_max, length.out = 160), x = seq(fdrbar_xmin, fdrbar_xmax, length.out = 160), fdr_col = fdr_pal)
bar_ticks <- tibble(v = c(1, 2, 3.2), label = c("1", "2", "3.2+")) |>
  mutate(x = fdrbar_xmin + (v - fdr_min) / (fdr_max - fdr_min) * (fdrbar_xmax - fdrbar_xmin))
group_ranges <- row_map |>
  group_by(ontology) |>
  summarise(group_ymin = min(ymin), group_ymax = max(ymax), .groups = "drop") |>
  mutate(ontology = factor(ontology, levels = c("BP", "CC", "MF"))) |>
  arrange(ontology)
group_band <- group_ranges |>
  mutate(
    xmid = 1.017,
    ymin_band = group_ymin - .035,
    ymax_band = group_ymax + .035,
    band_colour = c(BP = "#3F7F8F", CC = "#77629B", MF = "#A46C42")[as.character(ontology)]
  )
separator_dat <- tibble(
  y = c(
    mean(c(group_ranges$group_ymin[group_ranges$ontology == "BP"], group_ranges$group_ymax[group_ranges$ontology == "CC"])),
    mean(c(group_ranges$group_ymin[group_ranges$ontology == "CC"], group_ranges$group_ymax[group_ranges$ontology == "MF"]))
  )
)

p <- ggplot() +
  geom_polygon(data = ribbons, aes(x, y, group = ribbon), fill = "#D2D2D2", colour = NA, alpha = .46) +
  geom_rect(data = gene_slots, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = slot_colour),
            colour = "white", linewidth = .16) +
  geom_text(data = gene_meta, aes(x = .183, y = gene_y, label = gene), hjust = 1, size = 2.68,
            family = "Arial", fontface = "plain", colour = "#202020") +
  # A subtle white halo keeps pathway names readable over the ribbon field.
  geom_text(data = row_map, aes(x = .685, y = y, label = str_wrap(display_label, 27)), hjust = 1,
            size = 2.48, lineheight = .86, family = "Arial", fontface = "bold", colour = "white") +
  geom_text(data = row_map, aes(x = .685, y = y, label = str_wrap(display_label, 27), colour = is_focus), hjust = 1,
            size = 2.25, lineheight = .86, family = "Arial") +
  # Exactly one pathway stratum per GO term; focus ribbons enter this same rectangle.
  geom_rect(data = row_map, aes(xmin = .695, xmax = .735, ymin = ymin, ymax = ymax, fill = term_colour), colour = NA) +
  annotate("rect", xmin = dot_xmin, xmax = dot_xmax, ymin = min(row_map$ymin) - .16,
           ymax = max(row_map$ymax) + .16, fill = NA, colour = "black", linewidth = .36) +
  geom_segment(data = separator_dat, aes(x = dot_xmin, xend = dot_xmax, y = y, yend = y),
               colour = "#C8C8C8", linewidth = .25) +
  geom_point(data = row_map, aes(x = dot_x, y = y, fill = fdr_col, size = bubble_size), shape = 21,
             stroke = .20, colour = "#5C3030") +
  geom_rect(data = group_band, aes(xmin = 1.010, xmax = 1.024, ymin = ymin_band, ymax = ymax_band, fill = band_colour), colour = NA, alpha = .88) +
  geom_text(data = group_band, aes(x = xmid, y = (ymin_band + ymax_band) / 2, label = ontology), angle = 90,
            size = 2.05, family = "Arial", fontface = "bold", colour = "white") +
  geom_segment(data = ticks, aes(x = x, xend = x, y = min(row_map$ymin) - .16,
                                 yend = min(row_map$ymin) - .26), linewidth = .30) +
  geom_text(data = ticks, aes(x = x, y = min(row_map$ymin) - .45, label = format(v, nsmall = 2)),
            size = 1.95, family = "Arial") +
  annotate("text", x = mean(c(dot_xmin, dot_xmax)), y = min(row_map$ymin) - .76,
           label = "Gene ratio", size = 2.55, family = "Arial") +
  geom_point(data = count_legend, aes(x = x, y = max(row_map$ymax) + .78, size = size),
             shape = 21, fill = "white", colour = "#555555", stroke = .24) +
  geom_text(data = count_legend, aes(x = x, y = max(row_map$ymax) + .48, label = Count),
            size = 1.72, family = "Arial", colour = "#303030") +
  annotate("text", x = mean(count_legend$x), y = max(row_map$ymax) + 1.10,
           label = "Gene count", size = 2.25, family = "Arial") +
  geom_rect(data = bar_dat, aes(xmin = x, xmax = x + .0009, ymin = max(row_map$ymax) + .69,
                                ymax = max(row_map$ymax) + .87, fill = fdr_col), colour = NA) +
  geom_segment(data = bar_ticks, aes(x = x, xend = x, y = max(row_map$ymax) + .68,
                                     yend = max(row_map$ymax) + .60), linewidth = .26) +
  geom_text(data = bar_ticks, aes(x = x, y = max(row_map$ymax) + .47, label = label),
            size = 1.78, family = "Arial") +
  annotate("text", x = mean(c(fdrbar_xmin, fdrbar_xmax)), y = max(row_map$ymax) + 1.10,
           label = "-log[10](FDR)", parse = TRUE, size = 2.25, family = "Arial") +
  scale_fill_identity() +
  scale_size_identity() +
  scale_colour_manual(values = c(`TRUE` = "#303030", `FALSE` = "#717171"), guide = "none") +
  coord_cartesian(xlim = c(.02, 1.040), ylim = c(min(row_map$ymin) - .86, max(row_map$ymax) + 1.16), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(2, 8, 7, 8), plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA))

prefix <- file.path(out, "AD_conserved_disease_GO_alluvial_enrichment_V7")
ggsave(paste0(prefix, ".pdf"), p, width = 180, height = 145, units = "mm", device = cairo_pdf, bg = "white")
ggsave(paste0(prefix, ".svg"), p, width = 180, height = 145, units = "mm", device = grDevices::svg, bg = "white")
ggsave(paste0(prefix, ".png"), p, width = 180, height = 145, units = "mm", dpi = 450, bg = "white")
write_tsv(row_map, file.path(out, "row_map_V7.tsv"))
write_tsv(focus_edges, file.path(out, "displayed_gene_GO_edges_V7.tsv"))
write_tsv(gene_meta, file.path(out, "displayed_genes_V7.tsv"))
script_arg <- commandArgs(FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
if (!is.na(script_path) && file.exists(script_path)) file.copy(script_path, file.path(out, "plot_AD_GO_alluvial_enrichment_V7.R"), overwrite = TRUE)
