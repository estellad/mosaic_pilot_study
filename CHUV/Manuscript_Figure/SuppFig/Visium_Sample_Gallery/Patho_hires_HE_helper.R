plotVisium_new <- function(spe, spots = TRUE, annotate = NULL, image = TRUE, 
                           zoom = FALSE, show_axis = FALSE, assay = "counts", 
                           trans = "identity", legend.position = "right", 
                           x_coord = NULL, y_coord = NULL, y_reverse = TRUE, sample_ids = NULL, 
                           image_ids = NULL, palette = NULL, pt.size = 1, guide.pt.size = 4){
  stopifnot(is(spe, "SpatialExperiment"), is.logical(spots), 
            length(spots) == 1, is.logical(image), length(image) == 
              1, is.logical(y_reverse), length(y_reverse) == 1)
  stopifnot(legend.position %in% c("left", "right", "top", 
                                   "bottom", "none"))
  if(is.null(x_coord)) 
    x_coord <- spatialCoordsNames(spe)[1]
  if(is.null(y_coord)) 
    y_coord <- spatialCoordsNames(spe)[2]
  plt_df <- data.frame(colData(spe), spatialCoords(spe))
  if(!is.null(annotate)){
    stopifnot(is.character(annotate), length(annotate) == 1)
    if(!annotate %in% c(names(plt_df), rownames(spe))){
      stop("'annotate' should be in rownames(spe) or names(colData(spe))")
    }
    if(annotate %in% rownames(spe)){
      stopifnot(is.character(assay))
      plt_df[[annotate]] <- assay(spe, assay)[annotate, ]
    }
    if(is.numeric(plt_df[[annotate]]) & is.null(palette)){
      palette <- "seuratlike"
    }
    palette <- .get_pal(palette, plt_df[[annotate]])
  }else{
    annotate <- "foo"
    plt_df[[annotate]] <- "black"
  }
  if(is.null(sample_ids)){
    sample_ids <- unique(spe$sample_id)
  }else{
    spe <- spe[, spe$sample_id %in% sample_ids]
  }
  img_df <- .sub_imgData(spe, sample_ids, image_ids)
  rownames(img_df) <- img_df$sample_id
  if (image) {
    images <- lapply(sample_ids, function(s) {
      spi <- img_df[s, "data"]
      img <- imgRaster(spi[[1]])
      layer(data = data.frame(sample_id = s), inherit.aes = FALSE, 
            stat = "identity", position = "identity", geom = ggplot2::GeomCustomAnn, 
            params = list(grob = rasterGrob(img), xmin = 0, 
                          xmax = ncol(img), ymin = 0, ymax = nrow(img)))
    })
    img <- img_df$data[[1]]
    xlim <- c(0, ncol(img))
    ylim <- c(0, nrow(img))
    if (zoom) {
      xlim <- ylim <- NULL
    }
  }else{
    img <- NULL
    images <- xlim <- ylim <- NULL
  }
  for (s in sample_ids) {
    ix <- plt_df$sample_id == s
    xy <- c(x_coord, y_coord)
    sf <- img_df[s, "scaleFactor"]
    plt_df[ix, xy] <- sf * plt_df[ix, xy]
    if (y_reverse) 
      plt_df <- .y_reverse(plt_df, ix, y_coord, img)
  }
  if (spots){
    guide <- ifelse(is.numeric(plt_df[[annotate]]), guide_colorbar, guide_legend)
    points <- list(guides(fill = guide(title = annotate, order = 1, override.aes = list(col = NA, size = guide.pt.size)),
                          color = "none"), 
                   geom_point(shape = 21, size = pt.size, stroke = 0, 
                              alpha = 1))
  }else{
    points <- geom_point(col = "transparent")
  }
  scale <- if(annotate != "foo"){
    if(is.numeric(plt_df[[annotate]])){
      if(length(palette) == 1 && palette == "viridis"){
        scale_fill_viridis_c(trans = trans)
      }else if(length(palette) == 1 && palette == "seuratlike"){
        scale_fill_gradientn(colors = colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100), 
                             trans = trans, 
                             limits = c(min(plt_df[[annotate]]), max(plt_df[[annotate]])))
      }else{
        scale_fill_gradient(low = palette[1], high = palette[2], trans = trans)
      }
    }else if(is.factor(plt_df[[annotate]]) | is.character(plt_df[[annotate]])){
      plt_df[[annotate]] <- as.factor(plt_df[[annotate]])
      if (is.null(palette)) {
        scale_fill_manual(name = annotate, values = (scales::hue_pal())(length(unique(plt_df[[annotate]]))))
      }else if(!is.null(palette)){
        scale_fill_manual(values = palette)
      }
    }
  }else{
    scale_fill_identity()
  }
  p <- ggplot(plt_df, aes_string(x_coord, y_coord, fill = annotate, color = annotate)) + images + points +
    scale_fill_manual(values = palette) + scale_color_manual(values = palette) + coord_fixed(xlim = xlim, ylim = ylim)
  if (show_axis) {
    p <- p + theme_bw() + 
      theme(strip.background = element_blank(), 
            legend.text = element_text(size=10),
            legend.position = legend.position) + 
      labs(x = paste0("pxl_col_in_", img_df[s, "image_id"]), y = paste0("pxl_col_in_", img_df[s, "image_id"])) 
  }else{
    p <- p + theme_void() + 
      theme(strip.text = element_text(margin = margin(0, 0, 0.5, 0, "lines")), legend.position = legend.position) 
  }
  return(p)
}


# -------------------------------------------------------------------------
.get_pal <- function(pal, val) {
  
  if (length(pal) == 1) {
    pal <- switch(pal, 
                  "libd_layer_colors" = c(
                    "#F0027F", "#377EB8", "#4DAF4A", "#984EA3", 
                    "#FFD700", "#FF7F00", "#1A1A1A", "#666666"), 
                  "Okabe-Ito" = c(
                    "#000000", "#E69F00", "#56B4E9", "#009E73", 
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 
                  # use 'scale_fill_viridis_c' for the following options
                  "viridis" = pal, 
                  "magma" = pal, 
                  "inferno" = pal, 
                  "plasma" = pal, 
                  "viridis" = pal, 
                  "cividis" = pal, 
                  "rocket" = pal, 
                  "mako" = pal, 
                  "turbo" = pal, 
                  "seuratlike" = pal, 
                  # for a single color name, combine with "gray95" for continuous color scale
                  c("gray95", pal)
    )
  }
  
  # if length(pal) == 0 (i.e. 'pal' is NULL), leave 'pal' unchanged and the
  # plotting functions will select default palettes instead
  
  # if length(pal) > 1, use 'pal' as provided (e.g. multiple colors for discrete
  # labels, or length 2 for continuous gradient)
  
  return(pal)
}


#' @importFrom SpatialExperiment imgData
.sub_imgData <- function(spe, sample_ids, image_ids) {
  .get_img_idx <- SpatialExperiment:::.get_img_idx
  if (is.null(image_ids)) {
    # default to first available image for each sample
    idx <- .get_img_idx(spe, TRUE, NULL)
  } else {
    if (length(image_ids) == 1) {
      idx <- .get_img_idx(spe, TRUE, image_ids)
    } else {
      stopifnot(length(image_ids) == length(sample_ids))
      idx <- mapply(s = sample_ids, i = image_ids,
                    function(s, i) .get_img_idx(spe, s, i))
    }
  }
  imgData(spe)[idx, ]
}


.y_reverse <- function(df, ix, y, img) {
  y_tmp <- df[ix, y]
  if (!is.null(img)) {
    y_tmp <- nrow(img) - y_tmp
  } else {
    y_tmp <- max(y_tmp) - y_tmp
  }
  df[ix, y] <- y_tmp
  return(df)
}

