library(simpleSeg)
library(cytomapper)

pathToImages <- "/home/nick/Downloads/images/"
# Get directories of images
imageDirs <- dir(pathToImages, full.names = TRUE)
names(imageDirs) <- dir(pathToImages, full.names = FALSE)
# Get files in each directory
files <- sapply(imageDirs, list.files, pattern = "tif", full.names = TRUE, simplify = FALSE)
# Read files with readImage from EBImage
images <- lapply(files, EBImage::readImage, as.is = TRUE)
images <- cytomapper::CytoImageList(images)
image = images[[1]]
discSize = 3
smooth = 1
nucleusIndex <- "191Ir_DNA1"
sizeSelection <- 40
ext <- 1
# Mask the nuclei
nuc <- image[,, "191Ir_DNA1"]
nth <- EBImage::otsu(nuc, range = range(nuc))
nMask <- nuc > nth
# Dilate out from the nuclei
kern <- EBImage::makeBrush(discSize, shape = "disc")
cell <- EBImage::dilate(nMask, kern)
# Bit of smoothing
image <- apply(image, 3, function(x) {
  x <- (x)
  EBImage::gblur(x, smooth)
}, simplify = FALSE)
# Do PCA only on the dilated nuclei
image <- EBImage::abind(image, along = 3)
image.long <- apply(image, 3, as.numeric)
use <- as.vector(cell)
pca <- prcomp(image.long[use, apply(image.long, 2, sd) >
                           0], scale = TRUE)
usePC <- 1
if (any(nucleusIndex %in% colnames(image.long))) {
  ind <- intersect(nucleusIndex, colnames(image.long))
  usePC <- which.max(abs(apply(pca$x, 2, cor, image.long[use,
                                                         nucleusIndex[nucleusIndex != "PCA"][1]])))
  PC <- pca$x[, usePC]
  PC <- PC * sign(cor(PC, image.long[use, nucleusIndex[nucleusIndex !=
                                                      "PCA"][1]]))
} else {
  PC <- pca$x[, usePC]
}
imagePC <- as.matrix(image[, , 1])
imagePC[] <- 0
imagePC[cell] <- PC -min(PC)
# Chuck out all the gunk from the dilation
nuc <- imagePC
nth <- EBImage::otsu(nuc, range = range(nuc))
nuc[nuc < nth] <- 0


#Estimate tolerance for watershed
y <- EBImage::distmap(nuc>0)
fit <- lm(as.numeric(nuc[y > 0 & y < 9]) ~ as.numeric(y[y > 0 & y < 9]))
tolerance <- coef(fit)[2]

# Do watershed
wMask <- EBImage::watershed(nuc, tolerance = tolerance, ext = ext)
tabNuc <- table(wMask)
wMask[wMask %in% names(which(tabNuc <= sizeSelection))] <- 0
EBImage::display(colorLabels(wMask))
cytImage <- cytomapper::CytoImageList(images[[1]])
cytMask <- cytomapper::CytoImageList(Image(wMask))
mcols(cytMask) <- mcols(cytImage) <- DataFrame(imageID = "a")
cytomapper::plotPixels(image = cytImage[1],
                       mask = cytMask[1],
                       img_id = "imageID",
                       colour_by = c("89Y_aSMA", "145Nd_ECadherin", "175Lu_CK5", "191Ir_DNA1", "171Yb_ER"),
                       display = "single",
                       colour = list(`175Lu_CK5` = c("black","blue"),
                                     `191Ir_DNA1` = c("black","purple"),
                                     `171Yb_ER` = c("black","yellow"),
                                     `145Nd_ECadherin` = c("black", "red"),
                                     `89Y_aSMA` = c("black", "green")),
                       bcg = list(`175Lu_CK5` = c(0, 1, 1.5),
                                  `191Ir_DNA1` = c(0, 1, 2),
                                  `171Yb_ER` = c(0, 1, 1.5),
                                  `145Nd_ECadherin` = c(0, 1, 1.5),
                                  `89Y_aSMA` = c(0, 1, 1.5)),
                       legend = NULL, return_plot= TRUE)