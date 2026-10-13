#' @importFrom EBImage normalize otsu bwlabel
.calcTissueMask <- function(image, tissue_index, size_selection = 10) {
  if (is.null(tissue_index)) {
    tissue_index <- seq_len(dim(image)[3])
  }

  ## add intensities of markers to try and highlight the tissue structure
  tissue <- rowMeans(image[, , tissue_index, drop = FALSE], dims = 2)

  tissue <- EBImage::normalize(log10(tissue + 0.01),
    ft = c(0, 1)
  )

  tissue_otsu <- EBImage::otsu(tissue,
    range = c(0, 1)
  )

  tissue <- tissue > tissue_otsu

  ## Do a size selection
  tissue[.smallObjects(EBImage::bwlabel(tissue), size_selection)] <- 0

  tissue
}
