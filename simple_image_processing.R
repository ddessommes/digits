#' ---
#' title: "Simple Image Processing and Clustering"
#' author: "Michael Hahsler"
#' output:
#'  html_document:
#'    toc: true
#' ---

#' ![CC](https://i.creativecommons.org/l/by/4.0/88x31.png)
#' This work is licensed under the
#' [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/). For questions please contact
#' [Michael Hahsler](http://michael.hahsler.net).
#'

#' # Numbers data set
#'
#' * data is available at http://michael.hahsler.net/SMU/EMIS7332/data/numbers/
#' * data contains hand-written digits as 28x28 pixels
#' * data was extracted from http://yann.lecun.com/exdb/mnist/

library(seriation) # used for `pimage()`
set.seed(1234)

numbers <- read.csv("numbers.csv", header=TRUE)
dim(numbers)

#' # Helper functions
#' Take a row of pixels and make it into a 28x28 matrix
toMatrix <- function(x) matrix(as.numeric(x), nrow=28, byrow=TRUE)

#' Convert a matrix back to a vector
toVector <- function(x) as.vector(t(x))

#' Test if it works
all(toVector(toMatrix(numbers[2,])) == numbers[2,])

toMatrix(numbers[2,])
pimage(toMatrix(numbers[2,]))

#' # Some basic image processing

#' # 2D convolution
#' see http://en.wikipedia.org/wiki/Kernel_%28image_processing%29
#' and http://graphics.stanford.edu/courses/cs178/applets/convolution.html
#'
#' __Note:__ This is very slow $O(n^2 k^2)$ ($n$=image size, $k$=kernel size).
#' So we do it in C++ (Source: http://permalink.gmane.org/gmane.comp.lang.r.rcpp/2926)
#'
#' For Windows you will need to install Rtools from https://cran.r-project.org/bin/windows/Rtools/ to be able to compile code.
library(Rcpp)
library(inline)
convolve_2d <- cxxfunction(signature(sampleS = "numeric",
  kernelS = "numeric"),
  plugin = "Rcpp", '
    Rcpp::NumericMatrix sample(sampleS), kernel(kernelS);
    int x_s = sample.nrow(), x_k = kernel.nrow();
    int y_s = sample.ncol(), y_k = kernel.ncol();

    Rcpp::NumericMatrix output(x_s + x_k - 1, y_s + y_k - 1);
    for (int row = 0; row < x_s; row++) {
     for (int col = 0; col < y_s; col++) {
       for (int i = 0; i < x_k; i++) {
         for (int j = 0; j < y_k; j++) {
           output(row + i, col + j) += sample(row, col) * kernel(i, j);
         }
       }
     }
    }
    return output;
  ')


#' ## Blurring images (with a Gaussian Kernel)
blurr <- rbind(
  c(0  , .5, 0 ),
  c(0.5, 1 , .5),
  c(0  , .5, 0 )
)

blurr <- blurr/sum(blurr) #' normalize to sum to 1
pimage(blurr)

#m <- toMatrix(numbers[1,])
#m <- toMatrix(numbers[2,])
#m <- toMatrix(numbers[3,])
m <- toMatrix(numbers[4,])

pimage(m)
mb <- convolve_2d(m, blurr) # blurr
pimage(mb)
mbb <- convolve_2d(mb, blurr) # blurr some more
pimage(mbb)

#' ## Keep only 40% of the darkest (non-white) pixels
mbb2 <- mbb>=quantile(mbb[mbb>0], 1-.4)
pimage(mbb2)


#' ## Extract patterns
#' vertical pattern of 2 dark pixels
conv1 <- rbind(
  c(0,0,1,1,0,0),
  c(0,0,1,1,0,0),
  c(0,0,1,1,0,0),
  c(0,0,1,1,0,0),
  c(0,0,1,1,0,0),
  c(0,0,1,1,0,0)
)
#' make the mean of the kernel 0
conv1 <- conv1-mean(conv1)

pimage(conv1)

#' the horizontal pattern
conv2 <- t(conv1)
pimage(conv2)

#' we start with a b/w image
mc1 <- convolve_2d(mbb2, conv1)
pimage(mc1)
# use 5% of the highest values
pimage(mc1>quantile(mc1, .95))


#' find horizontal lines
mc2 <- convolve_2d(mbb2, conv2)
pimage(mc2)
pimage(mc2>quantile(mc1, .95))

#' ## Edge detection
conv3 <- rbind(
  c(0 ,-2 ,0 ),
  c(-2, 8 ,-2),
  c(0 ,-2 ,0 )
)

mc3 <- convolve_2d(mbb2, conv3)
pimage(mc3)
pimage(mc3<0)

#' ## End point detection
#'
#' We use thinning and then the edge detection kernel.
#' _Note:_ There are better ways to do this!
source("thinning.R")
pimage(m>100)
m_thin <- thinImage(m>100)
pimage(m_thin)

mc4 <- convolve_2d(m_thin, conv3)
pimage(mc4)
pimage(mc4>4)
sum(mc4>4)

#' # Clustering
#' _Simple approach:_ Take a small sample and
#' cluster the pixels directly. I use 10 clusters, but
#' different writing styles mean that we should use more clusters.
#' You should probably use some preprocessing (e.g., blurring, rotation) and
#' feature engineering instead of the raw pixels.
s <- numbers[sample(1:nrow(numbers), 1000),]
km <- kmeans(s, c=10, nstart=5)
#km

#' How many are in each cluster?
km$size

#' Within Sum of Squares of the clusters
km$withinss

#' ## Look at some members of cluster 1
s1 <- s[km$cluster==1,]
pimage(toMatrix(s1[1,]))
pimage(toMatrix(s1[2,]))
pimage(toMatrix(s1[3,]))
pimage(toMatrix(s1[4,]))
pimage(toMatrix(s1[5,]))
pimage(toMatrix(s1[6,]))

#' ## Look at centroids
pimage(toMatrix(km$centers[1,]))
pimage(toMatrix(km$centers[2,]))
pimage(toMatrix(km$centers[3,]))
pimage(toMatrix(km$centers[4,]))
pimage(toMatrix(km$centers[5,]))
pimage(toMatrix(km$centers[6,]))
pimage(toMatrix(km$centers[7,]))
pimage(toMatrix(km$centers[8,]))
pimage(toMatrix(km$centers[9,]))
pimage(toMatrix(km$centers[10,]))

#' ## How similar are the cluster centers?
#'
#' Measured as Euclidean distance
plot(hclust(dist(km$centers)))

#' Measured as Pearson Correlation
plot(hclust(as.dist(1-cor(t(km$centers)))))

#' Measured as amount of ink on the paper
plot(hclust(as.dist(
  abs(outer(rowSums(km$centers), rowSums(km$centers), FUN = "-"))
  )))

#' # Feature Reduction with PCA
pc <- prcomp(s)

# How important are the first principal components?
plot(pc)

str(pc)
#' `$x` contains the data projected on the PCs. Let's look at the
#' first two PCs.

plot(pc$x[,1:2])

#' show the row number in s for the image
plot(pc$x[,1:2], col = 0)
text(pc$x[,1], pc$x[,2], 1:nrow(s), cex = .6)

#' look at some images to the far left
pimage(toMatrix(s[986,]))
pimage(toMatrix(s[432,]))
pimage(toMatrix(s[110,]))
pimage(toMatrix(s[315,]))

#' cluster the images using only the first 10 PCs
data_subspace <- pc$x[,1:10]
k <- 12
km <- kmeans(data_subspace, centers = k)
plot(data_subspace, col = km$cluster)

#' show average image for all clusters
for(i in 1:k)
pimage(toMatrix(colMeans(s[km$cluster == i,])))

#' # 25 Best and Worst Written Digits using LOF
library(dbscan)

l <- lof(s)

# ' Best handwriting
for(i in order(l)[1:25]) pimage(toMatrix(s[i,]))

#' Worst handwriting
for(i in order(l, decreasing = TRUE)[1:25]) pimage(toMatrix(s[i,]))
