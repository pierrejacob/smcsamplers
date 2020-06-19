## obtain data from
## https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary.html
## https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/binary/covtype.libsvm.binary.scale.bz2
## unzip to obtain file covtype.libsvm.binary.scale
xx <- read.csv("~/Downloads/covtype.libsvm.binary.scale", header = FALSE, sep = " ")
head(xx)
dim(xx)
xx <- xx[,1:11]
head(xx)
xxnumeric <- matrix(nrow = dim(xx)[1], ncol = dim(xx)[2])
xxnumeric[,1] <- as.numeric(xx[,1])
for (j in 2:dim(xx)[2]){
  xxnumeric[,j] <- as.numeric(sapply(strsplit(as.character(xx[,j]), ":"), function(v) v[2]))
}

head(xx)
head(xxnumeric)
summary(xxnumeric)
which(is.na(xxnumeric[,11]))
## remove NA
xxnumeric <- xxnumeric[-which(is.na(xxnumeric[,11])),]
summary(xxnumeric)

## random permutation
permut <- sample(1:dim(xxnumeric)[1], size = dim(xxnumeric)[1], replace = FALSE)
xxnumeric <- xxnumeric[permut,]

## split into training and test data
trainingset <- xxnumeric[1:4e5,]
testset <- xxnumeric[(4e5+1):dim(xxnumeric)[1],]

save(trainingset, testset, file = "~/Downloads/covtype.processed.RData")
