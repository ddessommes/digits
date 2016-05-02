#A large part of this Code is taken from Michael Hahsler (michael.hahsler.net)

#LIBRARIES ####
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



#LOAD/SAVE DATA ####
#gnumbers <- read.csv("gnumbers.csv", header=TRUE)
#num_labels <- read.csv("numbers_labels.csv", header=TRUE)
#save(numbers, file="numbers.rda")
#save(num_labels, file="num_labels.rda")
#load("numbers.rda")
#write.csv(gnumbers, file = "gnumbers.csv", row.names = FALSE)
#write.csv(gfeatures, file = "gfeatures.csv", row.names = FALSE)


#FUNCTIONS ####
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

#VISUALIZE SAMPLE DATA ####
pimage(toMatrix(numbers[2,]))

#IMAGE BLURR with a Gaussian Kernel ####
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

#PATTERN EXTRACTION ####
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

source("thinning.R")
#IMAGE THINNING ####
#source("thinning.R")
#is <- toMatrix(numbers[154,])
#pimage(is)
#pimage(is>128)
#is_thin <- thinImage(is>128)
#pimage(is_thin)
#(v <- toVector(is_thin))

#IMAGE EDGE DETECT with a Gaussian Kernel ####
#edge <- rbind(
#  c(0.0,-2.0, 0.0),
# c(-2.0, 8.0,-2.0),
#  c(0.0,-2.0, 0.0)
#)
#
#c3 <- convolve_2d(x_bb40, edge)
#pimage(c3)
#pimage(c3<0)

#END POINT DETECTION ####
#source("thinning.R")
#pimage(m>100)
#m_thin <- thinImage(m>100)
#pimage(m_thin)
#mc4 <- convolve_2d(m_thin, edge)
#pimage(mc4)
#pimage(mc4>4)
#sum(mc4>4)

#spatialfil() easily generates different Gaussian kernel types ####
cG <- convKernel(sigma = 1.5, k = "gaussian")   #Gaussian Filter
G <- cG$matrix
cS <- convKernel(sigma = 1.5, k = "sobel")      #Sobel Filter
S <- cS$matrix
cL <- convKernel(sigma = 1.5, k = "laplacian")  #Laplacian Filter
L <- cL$matrix
cK <- convKernel(sigma = 1.5, k = "LoG")        #Laplacian of Gaussian Filter
K <- cK$matrix
pimage(G)

#REINITIALIZE KEY DATASETS AS NEEDED TO REUSE CODE INSTEAD OF ADDING BLOCKS
#SENSITIVITY ANALYSIS ####
#A) Does Thinning help?
#B) Does PCA help?
#C) Does Outlier Removal help?

#THINNING DATASET LOOP ####
#load("numbers.rda")
#numbers_orig <- numbers
#numbers[numbers != 0] <- 0
#numbers <- numbers[1:5000,]
#x <- 0
#xt <- 0
#tv <- 0
#sample_orig_nbrs <- numbers_orig[sample(1:nrow(numbers_orig), 5000),]
#rm(numbers_orig)
#for (i in 1:5000){
#  x <- toMatrix(sample_orig_nbrs[i,])
#  xt <- thinImage(x>128)
#  tv <- toVector(xt)
#  numbers[i,] <- tv
#}

#PRINCIPAL COMPONENTS ANALYSIS ####
numbers <- sample_orig_nbrs
# Run all the subsequent code to get the k-Means analysis then run PCA on the gfeatures DF
pc_nbrs <- prcomp(features) #prcomp(numbers)
plot(pc_nbrs)
str(pc_nbrs)
plot(pc_nbrs$x[,1:2])
data_subspace <- pc_nbrs$x[,1:6] #use the first 6 PC
kc <- 16
kms_pc <- kmeans(data_subspace, centers = kc) #on this particular sample k = 16 looked promising
plot(data_subspace, col = kms_pc$cluster)

for(i in 1:kc)
  pimage(toMatrix(colMeans(numbers[kms_pc$cluster == i,])))
# There was no appreciative difference in running PCA against raw numbers or feature matrix

#SAMPLING WITH OUTLIER REMOVAL AND LABEL STORAGE ####
load("numbers.rda")
numbers$X <- as.vector(seq(1:nrow(numbers)), mode = "integer")
sample_orig_nbrs <- numbers[sample(1:nrow(numbers), 5500),]
lof_nbrs <- lof(sample_orig_nbrs)
sample_orig_nbrs$LOF <- lof_nbrs
tbl_df(sample_orig_nbrs)
sample_orig_nbrs <- sample_orig_nbrs %>% arrange(LOF)
sample_orig_nbrs <- slice(sample_orig_nbrs, 1:5000)
numbers[numbers != 0] <- 0
numbers <- numbers[1:5000,]
numbers <- sample_orig_nbrs
numbers$X   <- NULL
numbers$LOF <- NULL

#INITIAL DATA OPERATIONS with dplyr() ####
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
#save(gnumbers, file="gnumbers.rda")
gfeatures <- select(gnumbers, pixelSUM_ALL:pixelAVG_Q4)
#save(gfeatures, file="gfeatures.rda")

#CONVOLVE THE DATA SET FOR HORIZONTAL AND VERTICAL LINE PATTERNS ####
#Initialize Variables
x <- 0
obj = matrix(nrow = 5000, ncol = 2, byrow = TRUE) #create empty matrix to dump results into
#Execute For Loop ##a) 42000 (original images) ##b) 5000 (thinned images)
for(i in 1:5000){
x <- toMatrix(numbers[i,])
x_25 <- x>=quantile(x[x>0], 1-0) #Keep only X% darkest ##1-0.25 ##1-0
xh <- convolve_2d(x_25, conv_vc) #this attempts to find vertical lines 
xv <- convolve_2d(x_25, conv_hc) #this attempts to find horizontal lines 
xhq <- (xh>quantile(xh, .01)) # use (1-x)% of the highest values ##.97 ##.01
xvq <- (xv>quantile(xv, .01)) # use (1-x)% of the highest values ##.97 ##.01
obj[i,1] <- sum(toVector(xhq))
obj[i,2] <- sum(toVector(xvq))
}
colnames(obj) <- c("pixelSUM_Horiz","pixelSUM_Vert")

#pimage(x)
#pimage(x_25)
#pimage(xhq) 
#pimage(xvq) 

#CONVOLVE THE DATA SET FOR SOBEL FILTER from spatialfil() AND CENTROID (x,y) ####
#Initialize Variables
z <- 0
obj2 = matrix(nrow = 5000, ncol = 3, byrow = TRUE) #create empty matrix to dump results into
#Execute For Loop ##a) 42000 (original images) ##b) 5000 (thinned images)
for(i in 1:5000){
  z <- toMatrix(numbers[i,])
  z_25 <- z>=quantile(z[z>0], 1-.01) #Keep only 25% darkest ##1-.25 ##1-0
  zS <- convolve_2d(z, S) #this attempts to convolve with Sobel
  zSq <- (zS>quantile(zS, .01)) # use (1-x)% of the highest values ##.97 #.01
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

#UPDATE SELECTED FEATURES DATA FRAME WITH CONVOLUTION NUMBERS ####
gfeatures$pixelSUM_Horiz = obj[,"pixelSUM_Horiz"]
gfeatures$pixelSUM_Vert = obj[,"pixelSUM_Vert"]
gfeatures$pixelSUM_Sobel = obj2[,"pixelSUM_Sobel"]
gfeatures$pixelAVG_CentX = obj2[,"pixelAVG_CentX"]
gfeatures$pixelAVG_CentY = obj2[,"pixelAVG_CentY"]

gfeatures <- gfeatures %>% select(pixelAVG_ALL, 
                                  pixelPCT_H1,
                                  pixelPCT_Q1,
                                  pixelPCT_Q2,
                                  pixelPCT_Q3,
                                  pixelPCT_Q4,
                                  pixelSUM_Horiz,
                                  pixelSUM_Vert,
                                  pixelSUM_Sobel,
                                  pixelAVG_CentX,
                                  pixelAVG_CentY)

#SCALE SOME SELECTED FEATURES ####
#scaled_features_orig <- scaled_features
scaled_features <- gfeatures
scaled_features$pixelAVG_ALL <- scale(scaled_features$pixelAVG_ALL)
scaled_features$pixelSUM_Horiz <- scale(scaled_features$pixelSUM_Horiz)
scaled_features$pixelSUM_Vert  <- scale(scaled_features$pixelSUM_Vert)
scaled_features$pixelSUM_Sobel <- scale(scaled_features$pixelSUM_Sobel)
scaled_features$pixelAVG_CentX <- scale(scaled_features$pixelAVG_CentX)
scaled_features$pixelAVG_CentY <- scale(scaled_features$pixelAVG_CentY)
#scaled_features <- gfeatures %>% scale
#save(gfeatures, file="gfeatures.rda")

#CONFIRM CLUSTERING TENDENCY ####
#Distance Matrix
#sample_cviz <- scaled_features[sample(1:nrow(scaled_features), 5000),] # for Viz need smaller sample
#sample_cviz_df <- as.data.frame(sample_cviz) #next plot needs DF
#plot(sample_cviz_df$pixelSUM_ALL)
d_sample_cviz <- dist(scaled_features)    #For PCA round didn't need to resample already 5000 from Thin
#VAT(d_sample_cviz)
#iVAT(d_sample_cviz)

#CLUSTERING WITH SAMPLES OF 5000 AND DIFFERENT # OF CENTERS ####
sample_kms <- scaled_features #[sample(1:nrow(scaled_features), 5000),]
d_sample_kms <- dist(scaled_features)      #For PCA round didn't need to resample already 5000 from Thin
kms10 <- kmeans(sample_kms, centers = 10, nstart = 5)
kms15 <- kmeans(sample_kms, centers = 15, nstart = 5)
kms20 <- kmeans(sample_kms, centers = 20, nstart = 5)
kms30 <- kmeans(sample_kms, centers = 30, nstart = 5)

str(kms10)
str(kms20)
str(kms30)

kms10$size
kms20$size
kms30$size

100*(kms10$betweenss / kms10$totss)
100*(kms20$betweenss / kms20$totss)
100*(kms30$betweenss / kms30$totss)

100*(kms10$tot.withinss / kms10$totss)
100*(kms20$tot.withinss / kms20$totss)
100*(kms30$tot.withinss / kms30$totss)

plot(hclust(dist(kms10$centers)))
plot(hclust(dist(kms20$centers)))
plot(hclust(dist(kms30$centers)))

plot(hclust(as.dist(1-cor(t(kms10$centers)))))
plot(hclust(as.dist(1-cor(t(kms20$centers)))))
plot(hclust(as.dist(1-cor(t(kms30$centers)))))

fpc::cluster.stats(d_sample_kms, kms10$cluster, aggregateonly = TRUE) 
fpc::cluster.stats(d_sample_kms, kms20$cluster, aggregateonly = TRUE) 
fpc::cluster.stats(d_sample_kms, kms30$cluster, aggregateonly = TRUE) 

#CLUSTER OPTIMIZATION ####
#Total Within Sum of Squares (WSS) - Cohesion
ks <- seq(from = 10, to = 30, by = 2)
WSS <- sapply(ks, FUN=function(k) {
  kmeans(scaled_features, centers=k, nstart=5)$tot.withinss
})
plot(ks, WSS, type="l")
abline(v=c(10, 15, 20, 25), col="red", lty=2)

ks2 <- seq(from = 10, to = 40, by = 2)
WSS2 <- sapply(ks2, FUN=function(k) {
  kmeans(scaled_features, #[sample(1:nrow(scaled_features), 5000),], 
         centers=k, nstart=5)$tot.withinss
})
plot(ks2, WSS2, type="l")
abline(v=c(10, 15, 20, 25), col="red", lty=2)

ks3 <- seq(from = 10, to = 50, by = 5)
WSS3 <- sapply(ks3, FUN=function(k) {
  kmeans(scaled_features, #[sample(1:nrow(scaled_features), 5000),], 
         centers=k, nstart=5)$tot.withinss
})
plot(ks3, WSS3, type="l")
abline(v=c(10, 15, 20, 25), col="red", lty=2)

#Average Silhouette Width (ASW) - Cohesion and Separation
ASW <- sapply(ks, FUN=function(k) {
  fpc::cluster.stats(d_sample_kms, kmeans(sample_kms,
                                          centers=k,
                                          nstart=5)$cluster)$avg.silwidth
})

plot(ks, ASW, type="l")
#ks[which.max(ASW)] #10
abline(v=c(10, 15, 20, 25), col="red", lty=2)
#


#HIERARCHICAL CLUSTERING ####
sample_hc <- scaled_features #[sample(1:nrow(scaled_features), 5000),]
d_sample_hc <- dist(sample_hc)
hcl <- hclust(d_sample_hc, method="complete")
plot(as.dendrogram(hcl), leaflab = "none")
#rect.hclust(hcl, k=4)

#https://rpubs.com/gaston/dendrograms

#DBSCAN ####
sample_db <- scaled_features[sample(1:nrow(scaled_features), 1000),]
kNNdistplot(sample_db, k = 3)
abline(h=2, col="red")
(db <- dbscan(sample_db, eps=0.5, minPts=3))

#EXTERNAL CLUSTER VALIDATION ####
truth <- dplyr::inner_join(num_labels, sample_orig_nbrs, by = "X")
truth <- truth %>% select(label)
random10 <- sample(1:10, nrow(scaled_features), replace = TRUE)
random20 <- sample(1:20, nrow(scaled_features), replace = TRUE)
random30 <- sample(1:30, nrow(scaled_features), replace = TRUE)

val <- rbind(
  kms10 = c(
  unlist(fpc::cluster.stats(d_sample_kms, kms10$cluster, truth, compareonly = TRUE)),
  entropy = entropy(kms10$cluster, truth),
  purity = purity(kms10$cluster, truth)
  ),
  kms20 = c(
  unlist(fpc::cluster.stats(d_sample_kms, kms20$cluster, truth, compareonly = TRUE)),
  entropy = entropy(kms20$cluster, truth),
  purity = purity(kms20$cluster, truth)
  ),
  kms30 = c(
  unlist(fpc::cluster.stats(d_sample_kms, kms30$cluster, truth, compareonly = TRUE)),
  entropy = entropy(kms30$cluster, truth),
  purity = purity(kms30$cluster, truth)
  ),
  random10 = c(
  unlist(fpc::cluster.stats(d_sample_kms, random10, truth, compareonly = TRUE)), #recursive = FALSE),
  entropy = entropy(random10, truth),
  purity = purity(random10, truth)
  ),
  random20 = c(
  unlist(fpc::cluster.stats(d_sample_kms, random20, truth, compareonly = TRUE)),
  entropy = entropy(random20, truth),
  purity = purity(random20, truth)
  ),
  random30 = c(
  unlist(fpc::cluster.stats(d_sample_kms, random30, truth, compareonly = TRUE)),
  entropy = entropy(random30, truth),
  purity = purity(random30, truth)
  )
)
val


#CLUSTERING ON THE RAW DATA - REFERENCE CODE FROM CLASS ####
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
