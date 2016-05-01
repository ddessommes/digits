#A large part of this Code is taken from Michael Hahsler (michael.hahsler.net)

#LIBRARIES
library(plyr)
library(dplyr)
library(Rcpp)
library(inline)
library(seriation)
library(dbscan)
library(IM)
library(spatialfil)
#library(EBImage)
#library(imager)
#library(adimpro)
#library(ripa)
#library(smoothie)



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
#all(toVector(toMatrix(numbers[2,])) == numbers[2,])

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


#VISUALIZE SAMPLE DATA
pimage(toMatrix(numbers[2,]))


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
#c1 <- convolve_2d(x_bb40, conv_vc) #this attempts to find vertical lines
#pimage(c1)
#pimage(c1>quantile(c1, .95)) # use 5% of the highest values

#c2 <- convolve_2d(x_bb40, conv_hc) #this attempts to find horizontal lines
#pimage(c2)
#pimage(c2>quantile(c1, .95)) # use 5% of the highest values


#IMAGE EDGE DETECT with a Gaussian Kernel
#edge <- rbind(
#  c(0.0,-2.0, 0.0),
# c(-2.0, 8.0,-2.0),
#  c(0.0,-2.0, 0.0)
#)
#
#c3 <- convolve_2d(x_bb40, edge)
#pimage(c3)
#pimage(c3<0)

#spatialfil() easily generates different Gaussian kernel types
cG <- convKernel(sigma = 1.5, k = "gaussian")   #Gaussian Filter
G <- cG$matrix
cS <- convKernel(sigma = 1.5, k = "sobel")      #Sobel Filter
S <- cS$matrix
cL <- convKernel(sigma = 1.5, k = "laplacian")  #Laplacian Filter
L <- cL$matrix
cK <- convKernel(sigma = 1.5, k = "LoG")        #Laplacian of Gaussian Filter
K <- cK$matrix
pimage(G)


#CONVOLVE THE DATA SET FOR HORIZONTAL AND VERTICAL LINE PATTERNS
#Initialize Variables
x <- 0
obj = matrix(nrow = 42000, ncol = 2, byrow = TRUE) #create empty matrix to dump results into
#Execute For Loop
for(i in 1:42000){
x <- toMatrix(numbers[i,])
x_25 <- x>=quantile(x[x>0], 1-.25) #Keep only 25% darkest
xh <- convolve_2d(x_25, conv_vc) #this attempts to find vertical lines
xv <- convolve_2d(x_25, conv_hc) #this attempts to find horizontal lines
xhq <- (xh>quantile(xh, .97)) # use (1-x)% of the highest values
xvq <- (xv>quantile(xv, .97)) # use (1-x)% of the highest values
obj[i,1] <- sum(toVector(xhq))
obj[i,2] <- sum(toVector(xvq))
}
colnames(obj) <- c("pixelSUM_Horiz","pixelSUM_Vert")
#pimage(x)
#pimage(x_25)
#pimage(xhq) 
#pimage(xvq) 


#CONVOLVE THE DATA SET FOR SOBEL FILTER from spatialfil() AND CENTROID (x,y)
#Initialize Variables
z <- 0
obj2 = matrix(nrow = 42000, ncol = 3, byrow = TRUE) #create empty matrix to dump results into
#Execute For Loop
for(i in 1:42000){
  z <- toMatrix(numbers[i,])
  #z_25 <- z>=quantile(z[z>0], 1-.25) #Keep only 25% darkest
  zS <- convolve_2d(z, S) #this attempts to convolve with Sobel
  zSq <- (zS>quantile(zS, .97)) # use (1-x)% of the highest values
  cI <- calcCentroid(z)
  obj2[i,1] <- sum(toVector(zSq))
  obj2[i,2] <- cI[1]
  obj2[i,3] <- cI[2]
}
colnames(obj2) <- c("pixelSUM_Sobel","pixelAVG_CentX", "pixelAVG_CentY")
#pimage(z)
#pimage(z_25)
#pimage(zG) 
#pimage(zGq)

mutate(gfeatures, 
       pixelSUM_Horiz = obj[,"pixelSUM_Horiz"],
       pixelSUM_Vert = obj[,"pixelSUM_Vert"],
       pixelSUM_Sobel = obj2[,"pixelSUM_Sobel"],
       pixelAVG_CentX = obj2[,"pixelAVG_CentX"],
       pixelAVG_CentY = obj2[,"pixelAVG_CentY"])

scaled_features <- scale(gfeatures)
#scaled_features <- gfeatures %>% scale

#END POINT DETECTION

source("thinning.R")
pimage(m>100)
m_thin <- thinImage(m>100)
pimage(m_thin)
mc4 <- convolve_2d(m_thin, edge)
pimage(mc4)
pimage(mc4>4)
sum(mc4>4)

#CONFIRM CLUSTERING TENDENCY
#Distance Matrix
sample_scaled <- scaled_features[sample(1:nrow(scaled_features), 1000),]
sample_scaled <- as.data.frame(sample_scaled)
plot(sample_scaled$pixelSUM_ALL)
dist_sample_scl <- dist(sample_scaled)
#VAT(dist_sample_scl)
iVAT(dist_sample_scl)

#CLUSTERING WITH SAMPLES OF 5000 AND CHANGING # OF CENTERS FROM 10 TO 20 BY 2
s10 <- scaled_features[sample(1:nrow(scaled_features), 5000),]
kms10 <- kmeans(s10, centers = 10, nstart = 5)

s12 <- scaled_features[sample(1:nrow(scaled_features), 5000),]
kms12 <- kmeans(s12, centers = 12, nstart = 5)

s14 <- scaled_features[sample(1:nrow(scaled_features), 5000),]
kms14 <- kmeans(s14, centers = 14, nstart = 5)

s16 <- scaled_features[sample(1:nrow(scaled_features), 5000),]
kms16 <- kmeans(s16, centers = 16, nstart = 5)

s18 <- scaled_features[sample(1:nrow(scaled_features), 5000),]
kms18 <- kmeans(s18, centers = 18, nstart = 5)

s20 <- scaled_features[sample(1:nrow(scaled_features), 5000),]
kms20 <- kmeans(s20, centers = 20, nstart = 5)


str(kms10)
str(kms12)
str(kms14)
str(kms16)
str(kms18)
str(kms20)

kms10$size
kms12$size
kms14$size
kms16$size
kms18$size
kms20$size

100*(kms10$betweenss / kms10$totss)
100*(kms12$betweenss / kms12$totss)
100*(kms14$betweenss / kms14$totss)
100*(kms16$betweenss / kms16$totss)
100*(kms18$betweenss / kms18$totss)
100*(kms20$betweenss / kms20$totss)

100*(kms10$tot.withinss / kms10$totss)
100*(kms12$tot.withinss / kms12$totss)
100*(kms14$tot.withinss / kms14$totss)
100*(kms16$tot.withinss / kms16$totss)
100*(kms18$tot.withinss / kms18$totss)
100*(kms20$tot.withinss / kms20$totss)

plot(hclust(dist(kms10$centers)))
plot(hclust(dist(kms12$centers)))
plot(hclust(dist(kms14$centers)))
plot(hclust(dist(kms16$centers)))
plot(hclust(dist(kms18$centers)))
plot(hclust(dist(kms20$centers)))

plot(hclust(as.dist(1-cor(t(kms10$centers)))))
plot(hclust(as.dist(1-cor(t(kms12$centers)))))
plot(hclust(as.dist(1-cor(t(kms14$centers)))))
plot(hclust(as.dist(1-cor(t(kms16$centers)))))
plot(hclust(as.dist(1-cor(t(kms18$centers)))))
plot(hclust(as.dist(1-cor(t(kms20$centers)))))

#CLUSTERING ON THE RAW DATA
#View the Centroid Images
#pimage(toMatrix(kms$centers[1,]))
#pimage(toMatrix(kms$centers[2,]))

#View Cluster Sample Images (i.e. Cluster1)
#s3 <- s[kms1$cluster==3,]
#pimage(toMatrix(s3[10,]))
#pimage(toMatrix(s3[20,]))

#View the similarity of Cluster Centroids
#plot(hclust(dist(kms1$centers)))             
#Euclidean Distance
#plot(hclust(as.dist(1-cor(t(kms1$centers)))))
#Pearson Correlation
#plot(hclust(as.dist(
#  abs(outer(rowSums(kms$centers), rowSums(kms$centers), FUN = "-"))
#)))
#As the amount of ink on image
