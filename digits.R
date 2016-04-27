#A large part of this Code is taken from Michael Hahsler (michael.hahsler.net)

#LIBRARIES
library(dplyr)
library(Rcpp)
library(inline)
library(seriation)
library(dbscan)

#LOAD/SAVE DATA
#numbers <- read.csv("numbers.csv", header=TRUE)
#num_labels <- read.csv("numbers_labels.csv", header=TRUE)
#save(numbers, file="numbers.rda")
#save(num_labels, file="num_labels.rda")
#save(gnumbers, file="gnumbers.rda")
#save(gfeatures, file="gfeatures.rda")
load("numbers.rda")
#write.csv(gnumbers, file = "gnumbers.csv", row.names = FALSE)
#write.csv(gfeatures, file = "gnumbers.csv", row.names = FALSE)


#DATA OPERATIONS dplyr
#Reference
#http://stat545.com/block010_dplyr-end-single-table.html for more dplyr single DS operations

gnumbers <- numbers %>% tbl_df
gnumbers %>% glimpse

gnumbers <- gnumbers %>% mutate(pixelSUM_ALL = rowSums(gnumbers))
#NOTE: If pixel columns have been rearranged prior to doing this IT may NOT work!
gnumbers$pixelSUM_H1 <- select(gnumbers, pixel0:pixel391) %>% rowSums
gnumbers$pixelSUM_H2 <- select(gnumbers, pixel392:pixel783) %>% rowSums
gnumbers$pixelSUM_Q1 <- select(gnumbers, pixel0:pixel195) %>% rowSums
gnumbers$pixelSUM_Q2 <- select(gnumbers, pixel196:pixel391) %>% rowSums
gnumbers$pixelSUM_Q3 <- select(gnumbers, pixel392:pixel587) %>% rowSums
gnumbers$pixelSUM_Q4 <- select(gnumbers, pixel588:pixel783) %>% rowSums
#Get Ratios
gnumbers <- gnumbers %>% mutate(pixelPCT_H1 = (pixelSUM_H1/pixelSUM_ALL))
gnumbers <- gnumbers %>% mutate(pixelPCT_H2 = (pixelSUM_H2/pixelSUM_ALL))
gnumbers <- gnumbers %>% mutate(pixelPCT_Q1 = (pixelSUM_Q1/pixelSUM_ALL))
gnumbers <- gnumbers %>% mutate(pixelPCT_Q2 = (pixelSUM_Q2/pixelSUM_ALL))
gnumbers <- gnumbers %>% mutate(pixelPCT_Q3 = (pixelSUM_Q3/pixelSUM_ALL))
gnumbers <- gnumbers %>% mutate(pixelPCT_Q4 = (pixelSUM_Q4/pixelSUM_ALL))
#Get Averages
gnumbers <- gnumbers %>% mutate(pixelAVG_ALL = rowMeans(gnumbers))
gnumbers$pixelAVG_H1 <- select(gnumbers, pixel0:pixel391) %>% rowMeans
gnumbers$pixelAVG_H2 <- select(gnumbers, pixel392:pixel783) %>% rowMeans
gnumbers$pixelAVG_Q1 <- select(gnumbers, pixel0:pixel195) %>% rowMeans
gnumbers$pixelAVG_Q2 <- select(gnumbers, pixel196:pixel391) %>% rowMeans
gnumbers$pixelAVG_Q3 <- select(gnumbers, pixel392:pixel587) %>% rowMeans
gnumbers$pixelAVG_Q4 <- select(gnumbers, pixel588:pixel783) %>% rowMeans
gnumbers %>% glimpse

gfeatures <- select(gnumbers, pixelSUM_ALL:pixelAVG_Q4)

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
#pimage(blurr)

#Sample operations to visualize a number
#x1 <- toMatrix(numbers[420,])
#pimage(x1)
#x1_blur <- convolve_2d(x1, blurr)
#pimage(x1_blur)
#x1_bblur <- convolve_2d(x1_blur, blurr) # blurr some more
#pimage(x1_bblur)
#x1_bb40 <- x1_bblur>=quantile(x1_bblur[x1_bblur>0], 1-.4) #Keep only 40% darkest
#pimage(x1_bb40)


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


#CONVOLVE THE DATA SET
l <- 0
for(i in 1:5) 
{
  l <- l+i
  }
l

Zoo2 <- Zoo
for(i in 1:ncol(Zoo2)) Zoo2[[i]] <- as.factor(Zoo2[[i]])
sapply(Zoo2, class)


x <- toMatrix(numbers[3,])
x_25 <- x>=quantile(x[x>0], 1-.25) #Keep only 25% darkest

xh <- convolve_2d(x_25, conv_vc) #this attempts to find vertical lines
xv <- convolve_2d(x_25, conv_hc) #this attempts to find horizontal lines

xhq <- (xh>quantile(xh, .97)) # use (1-x)% of the highest values
xvq <- (xv>quantile(xv, .97)) # use (1-x)% of the highest values

(hq <-sum(toVector(xhq)))
(vq <-sum(toVector(xvq)))

pimage(x)
pimage(x_25)
pimage(xhq) 
pimage(xvq) 


#CONFIRM CLUSTERING TENDENCY
#Distance Matrix
#gfeatures %>% select(pixelSUM_H1,pixelSUM_H2) %>% scale
sg <- gfeatures[sample(1:nrow(gfeatures), 1000),]
plot(sg$pixelSUM_ALL)
sgs <- scale(sg)
dist_sgs <- dist(sgs)
VAT(dist_sgs)
iVAT(dist_sgs)

#CLUSTERING A SAMPLE
s <- numbers[sample(1:nrow(numbers), 1000),]
kms <- kmeans(s, c=20, nstart=5)
#kms
kms$size
kms$withinss

#View Cluster Sample Images (i.e. Cluster1)
s3 <- s[kms$cluster==3,]
pimage(toMatrix(s3[10,]))
pimage(toMatrix(s3[20,]))

#View the Centroid Images
pimage(toMatrix(kms$centers[1,]))
pimage(toMatrix(kms$centers[2,]))

#View the similarity of Cluster Centroids
plot(hclust(dist(kms$centers)))              
#Euclidean Distance
plot(hclust(as.dist(1-cor(t(kms$centers))))) 
#Pearson Correlation
plot(hclust(as.dist(
  abs(outer(rowSums(kms$centers), rowSums(kms$centers), FUN = "-"))
)))
#As the amount of ink on image
