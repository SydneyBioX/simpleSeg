test_that("Test masks are the sames as the saved ones.", {
    # load saved masks
    saved_masks <- readRDS("saved_masks.rds")

    # compute masks on example data
    # Get path to image directory
    pathToImages <- system.file("extdata", package = "simpleSeg")

    # Get directories of images

    imageDirs <- dir(pathToImages, "Point", full.names = TRUE)
    names(imageDirs) <- dir(pathToImages, "Point", full.names = FALSE)





    # Get files in each directory
    files <- sapply(imageDirs, list.files, pattern = "tif", full.names = TRUE, simplify = FALSE)

    # Read files with readImage from EBImage
    images <- lapply(files, EBImage::readImage, as.is = TRUE)

    # Convert to cytoImageList
    images <- cytomapper::CytoImageList(images)
    mcols(images)$imageID <- names(images)

    masks <- simpleSeg::simpleSeg(images,
                                 nucleus = "HH3",
                                 transform = "sqrt")

    expect_equal(saved_masks, masks)
})