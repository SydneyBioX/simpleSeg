# helper function to load data for testing.
load_data <- function() {
  # Get path to image directory
  path_to_images <- system.file("extdata", package = "simpleSeg")

  # Get directories of images

  image_dirs <- dir(path_to_images, "Point", full.names = TRUE)
  names(image_dirs) <- dir(path_to_images, "Point", full.names = FALSE)

  # Get files in each directory
  files <- lapply(
    image_dirs,
    list.files,
    pattern = "tif",
    full.names = TRUE
  )

  # Read files with readImage from EBImage
  images <- lapply(files, EBImage::readImage, as.is = TRUE)

  # Convert to cytoImageList
  images <- cytomapper::CytoImageList(images)
  mcols(images)$imageID <- names(images)

  return(images)
}
