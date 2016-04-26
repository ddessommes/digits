#A large part of this Code is taken from Michael Hahsler (michael.hahsler.net)

#LIBRARIES
library(dplyr)
library(Rcpp)
library(inline)
library(seriation)
library(dbscan)

#CREATE CLUSTERING ENTROPY AND PURITY FUNCTIONS 
entropy <- function(cluster, truth) {
  k <- max(cluster, truth)
  cluster <- factor(cluster, levels = 1:k)
  truth <- factor(truth, levels = 1:k)
  m <- length(cluster)
  mi <- table(cluster)
  
  cnts <- split(truth, cluster)
  cnts <- sapply(cnts, FUN = function(n) table(n))
  p <- sweep(cnts, 1, rowSums(cnts), "/")
  p[is.nan(p)] <- 0
  e <- -p * log(p, 2)
  sum(rowSums(e, na.rm = TRUE) * mi/m)
}

purity <- function(cluster, truth) {
  k <- max(cluster, truth)
  cluster <- factor(cluster, levels = 1:k)
  truth <- factor(truth, levels = 1:k)
  m <- length(cluster)
  mi <- table(cluster)
  
  cnts <- split(truth, cluster)
  cnts <- sapply(cnts, FUN = function(n) table(n))
  p <- sweep(cnts, 1, rowSums(cnts), "/")
  p[is.nan(p)] <- 0
  
  sum(apply(p, 1, max) * mi/m)
}


#LOAD/SAVE DATA
#numbers <- read.csv("numbers.csv", header=TRUE)
#num_labels <- read.csv("numbers_labels.csv", header=TRUE)
#save(numbers, file="numbers.rda")
#save(num_labels, file="num_labels.rda")
#load("numbers.rda")
#load("num_labels.rda")

#CREATE HELPER FUNCTIONS
toMatrix <- function(x) matrix(as.numeric(x), nrow=28, byrow=TRUE)
toVector <- function(x) as.vector(t(x))

#Validate Functions
all(toVector(toMatrix(numbers[2,])) == numbers[2,])

#VISUALIZE SAMPLE DATA
pimage(toMatrix(numbers[2,]))

#CREATE C++ 2D Convolution FUNCTION
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


#IMAGE BLURR with a Gaussian Kernel
blurr <- rbind(
  c(0.0, 0.5, 0.0),
  c(0.5, 1.0, 0.5),
  c(0.0, 0.5, 0.0)
)
blurr <- blurr/sum(blurr) #normalize to sum to 1
pimage(blurr)

#Sample operations to visualize a number
x <- toMatrix(numbers[4,])
pimage(x)
x_blur <- convolve_2d(x, blurr)
pimage(x_blur)
x_bblur <- convolve_2d(x_blur, blurr) # blurr some more
pimage(x_bblur)
x_bb40 <- x_bblur>=quantile(x_bblur[x_bblur>0], 1-.4) #Keep only 40% darkest
pimage(x_bb40)

#PATTERN EXTRACTION
#Example of vertical pattern with 2 dark pixels - i.e. NUMBER ONE
conv_vc <- rbind(c(0,0,1,1,0,0),
                 c(0,0,1,1,0,0),
                 c(0,0,1,1,0,0),
                 c(0,0,1,1,0,0),
                 c(0,0,1,1,0,0),
                 c(0,0,1,1,0,0)
)

conv_vc <- conv_vc-mean(conv_vc) #make the mean of the kernel = 0
pimage(conv_vc)

#Make the horizontal
conv_hc <- t(conv_vc)
pimage(conv_hc)

#Start with a B/W image (i.e. the 40% darkest)
c1 <- convolve_2d(x_bb40, conv_vc) #this attempts to find vertical lines
pimage(c1)
pimage(c1>quantile(c1, .95)) # use 5% of the highest values

c2 <- convolve_2d(x_bb40, conv_hc) #this attempts to find horizontal lines
pimage(c2)
pimage(c2>quantile(c1, .95)) # use 5% of the highest values


#IMAGE EDGE DETECT with a Gaussian Kernel
edge <- rbind(
  c(0.0,-2.0, 0.0),
 c(-2.0, 8.0,-2.0),
  c(0.0,-2.0, 0.0)
)

c3 <- convolve_2d(x_bb40, edge)
pimage(c3)
pimage(c3<0)

#CLUSTERING A SAMPLE
s <- numbers[sample(1:nrow(numbers), 1000),]
kms <- kmeans(s, c=20, nstart=5)
#kms
kms$size
kms$withinss

#View Cluster Sample Images (i.e. Cluster1)
s1 <- s[kms$cluster==1,]
pimage(toMatrix(s1[1,]))
pimage(toMatrix(s1[2,]))

#View the Centroid Images
pimage(toMatrix(kms$centers[1,]))
pimage(toMatrix(kms$centers[2,]))

#View the similarity of Cluster Centroids
plot(hclust(dist(kms$centers)))              #Euclidean Distance
plot(hclust(as.dist(1-cor(t(kms$centers))))) #Pearson Correlation
plot(hclust(as.dist(
  abs(outer(rowSums(km$centers), rowSums(kms$centers), FUN = "-"))
)))
#As the amount of ink on image
