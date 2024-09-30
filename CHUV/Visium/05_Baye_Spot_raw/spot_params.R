if(disease == "breast"){
  nclus <- c(16, 16, 16, 14, 16)
  spcl.initmethod = c("mclust", "mclust", "kmeans", "kmeans", "mclust")
}else if(disease == "lung"){
  nclus <- c(11, 10, 12, 15, 15)
  spcl.initmethod = c("mclust", "mclust", "mclust", "kmeans", "kmeans")
}else if(disease == "dlbcl"){
  # nclus <- c(15, 16, 15, 15, 16, 17) # D1 D3 hair problem
  # spcl.initmethod = c("mclust", "kmeans", "mclust", "mclust", "mclust", "mclust")
  nclus <- c(15, 16, 19, 15, 16, 17)
  spcl.initmethod = c("kmeans", "kmeans", "kmeans", "mclust", "mclust", "mclust")
}
