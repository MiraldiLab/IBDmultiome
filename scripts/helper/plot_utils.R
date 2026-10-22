# plot_utils.R
# Plotting helpers for data frames and Seurat objects.
#
# Dependencies: ggplot2 and Seurat.
# Both helpers return plot objects for further customization or saving.
# No local file paths require configuration.

# Load packages --------------------------------------------------------------
suppressPackageStartupMessages({
    library(ggplot2)
    library(Seurat)
})

# gg_vln_bar -----------------------------------------------------------------
# Create a violin or identity bar plot with optional box plots and QC lines.
# Column mappings (x, y, fill) are supplied as column-name strings.
# Only 'vln' and 'bar' add a primary geometry; other type values add neither.
# For a box-only plot, use type = 'box' with flag_box = TRUE.
# Returns a ggplot object.
gg_vln_bar <- function(
    dfp, # Input data frame
    type = 'vln', # Primary geometry: 'vln' or 'bar'; 'box' uses flag_box
    x = NULL, # x column
    y = NULL, # y column
    fill = NULL, # fill column
    fill_col = NULL, # fill colors
    lab_x = '', # x axis label
    lab_y = '', # y axis label
    lab_fill = NULL, # fill label
    hline = NULL, # y intercept of horizontal line(s)
    hline_type = 'dashed', # horizontal line type
    hline_col = 'red', # horizontal line color
    bar_width = 0.8, # bar width
    flag_leg = TRUE, # plot legend
    flag_box = FALSE, # add box plot
    box_width = 0.3, # box width
    box_col = 'white', # box fill color
    size_text = 12, # axis text size
    size_title = 14, # axis title size
    face_text = 'bold', # axis text face
    face_title = 'bold', # axis title face
    angle_x = 45, # x axis text angle
    pos_x = 'bottom', # x axis position
    x_hjust = 1, # x axis text justification
    flag_log10 = FALSE, # Apply a log10 scale to the y axis
    flag_grid_v = TRUE # vertical grid lines
) {
    p <- ggplot(dfp, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]]))
    if (type == 'vln') {
        p <- p + geom_violin()
    } else if (type == 'bar') {
        p <- p + geom_bar(stat = 'identity', width = bar_width)
    }
    p <- p + theme_minimal() +
            labs(x = lab_x, y = lab_y, fill = lab_fill) +
            theme(axis.text = element_text(size = size_text, face = face_text)) +
            theme(axis.title = element_text(size = size_title, face = face_title)) +
            theme(axis.text.x = element_text(angle = angle_x)) +
            geom_hline(yintercept = hline, linetype = hline_type, color = hline_col)
    if (pos_x == 'top') {
        p <- p + scale_x_discrete(position = 'top')
        x_hjust <- 0
    }
    p <- p + theme(axis.text.x = element_text(hjust = x_hjust))
    if (!(flag_leg)) {
        p <- p + theme(legend.position = 'none')
    }
    if (!(is.null(fill_col))) {
        p <- p + scale_fill_manual(values = fill_col)
    }
    if (flag_box) {
        p <- p + geom_boxplot(width = box_width, fill = box_col)
    }
    if (flag_log10) {
        p <- p + scale_y_continuous(trans = 'log10') + annotation_logticks(sides = 'l')
                            #theme(axis.line.y=element_line(color='grey40',size=1))
    }
    if (!(flag_grid_v)) {
        p <- p + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
    }
    p <- p + theme(panel.background = element_rect(color = 'grey40', size = 1))
    return(p)
}

# gg_vln_seurat_single --------------------------------------------------------
# Plot one gene or metadata feature across all cells as a single identity group.
# The identity assignment applies to the function's local object.
# Returns the object produced by VlnPlot() with the requested plot layers.
gg_vln_seurat_single <- function(
    obj, # Input Seurat object
    feature = NULL, # metadata or gene feature
    size_pt = 0, # point size
    title = '', # plot title
    title_size = 12, # title font size
    title_face = 'bold', # title face
    hline = NULL, # y intercept of horizontal line(s)
    hline_type = 'dashed', # horizontal line type
    hline_col = 'red', # horizontal line color
    fill_col = 'white', # Violin fill color
    flag_box = FALSE, # add box plot
    box_width = 0.3, # box width
    flag_log10 = FALSE # Apply a log10 scale to the y axis
) {
    if (is.null(feature)) {
        stop('ERROR: SELECT FEATURE TO PLOT')
    }
    Idents(obj) <- 'IDENT'
    p <- VlnPlot(obj, features = feature, pt.size = size_pt, cols = fill_col) + labs(x = '', title = title) +
                        geom_hline(yintercept = hline, linetype = hline_type, color = hline_col) + NoLegend() +
                        theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                                axis.ticks.x = element_blank(), axis.line.x = element_blank(),
                                title = element_text(size = title_size, face = title_face))
    if (flag_box) {
        p <- p + geom_boxplot(width = box_width)
    }
    if (flag_log10) {
        p <- p + scale_y_continuous(trans = 'log10') + annotation_logticks(sides = 'l')
    }
    return(p)
}
