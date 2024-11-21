rm(list = ls())
gc()
setwd("~/test")
# setwd("/mnt/f/07-CAUS/01-Linux-service/project_18DH-heterosis/Analysis/31-gwas_mph/test")

library(rbenchmark) # benchmark
library(bench)

# install.packages("combinat")                   # Install combinat package
# library("combinat") 
library(Rcpp)
# library(data.table)
# library(foreach)
# library(doParallel)
# library(BiocParallel)
library(getopt)
# library(biglm)
library(bigmemory)

# 测试参数
trait = "PH"
load(file = "Adata.RData")
load(file = "Ddata.RData")
load(paste0("trAdata_",trait,".RData"))
load(paste0("trDdata_",trait,".RData"))
load(paste0("TRAN_",trait,".RData"))
load(paste0("Y_",trait,".RData"))
load(file = "./total_ma_names.rdata")
nmar=500
TRAN = TRAN
eff_type = "ad"
epi1 = Adata
epi2 = Ddata
trdata1 = trAdata
trdata2 = trDdata
# total_ma_names = total_ma_names
no_cores = 4

# combinat::combn(1:5,2)
# utils::combn(1:5,2)
# 
# combos <- combinat::combn(total_ma_names,2)
# combos2 <- utils::combn(total_ma_names,2)# 这个计数有问题
# 
# 
# utils::combn()
# count <- choose(2100000000,2)
# count2 <-   nCm(2100000000,2)# 这个计数有问题
# rMVP
remove_bigmatrix <- function(x, desc_suffix=".desc", bin_suffix=".bin") {
  name <- basename(x)
  path <- dirname(x)
  
  descfile <- paste0(x, desc_suffix)
  binfile  <- paste0(x, bin_suffix)
  
  remove_var <- function(binfile, envir) {
    for (v in ls(envir = envir)) {
      if (is(get(v, envir = envir), "big.matrix")) {
        desc <- describe(get(v, envir = envir))@description
        if (desc$filename == binfile) {
          rm(list = v, envir = envir)
          gc()
        }
      }
    }
  }
  
  # remove_var(binfile, globalenv())
  remove_var(binfile, as.environment(-1L))
  
  if (file.exists(descfile)) {
    file.remove(descfile)
  }
  if (file.exists(binfile)) {
    file.remove(binfile)
  }
}

# utils::combn()
# x=length(total_ma_names)
# m=2
big.combn <- function (x, m, FUN = NULL, simplify = TRUE, ...) 
{
  stopifnot(length(m) == 1L, is.numeric(m))
  if (m < 0) 
    stop("m < 0", domain = NA)
  if (is.numeric(x) && length(x) == 1L && x > 0 && trunc(x) == 
      x) 
    x <- seq_len(x)
  n <- length(x)
  if (n < m) 
    stop("n < m", domain = NA)
  x0 <- x
  if (simplify) {
    if (is.factor(x)) 
      x <- as.integer(x)
  }
  m <- as.integer(m)
  e <- 0
  h <- m
  a <- seq_len(m)
  nofun <- is.null(FUN)
  if (!nofun && !is.function(FUN)) 
    stop("'FUN' must be a function or NULL")
  len.r <- length(r <- if (nofun) x[a] else FUN(x[a], ...))
  count <- as.numeric(round(choose(n, m)))
  if (simplify) {
    dim.use <- if (nofun) 
      c(m, count)
    else {
      d <- dim(r)
      if (length(d) > 1L) 
        c(d, count)
      else if (len.r != 1L) 
        c(len.r, count)
      else c(d, count)
    }
  }
  if (simplify) {
    # out <- matrix(r, nrow = len.r, ncol = count)
    remove_bigmatrix("combn")
    out <- filebacked.big.matrix(
      nrow = len.r,
      ncol = count,
      type = "integer", # char是c++的单个字符
      backingfile = "combn.bin", 
      backingpath = dirname("combn"), 
      descriptorfile = "combn.des",
      dimnames = c(NULL, NULL)
    )
    out[,1] <- r
  } else {
    out <- vector("list", count)
    out[[1L]] <- r
  }
  if (m > 0) {
    i <- 2
    nmmp1 <- n - m + 1
    while (a[1] != nmmp1) {
      if (e < n - h) {
        h <- 1
        e <- a[m]
        j <- 1
      } else {
        e <- a[m - h]
        h <- h + 1
        j <- 1:h
      }
      a[m - h + j] <- e + j
      r <- if (nofun) {
        x[a]
      } else FUN(x[a], ...)
      if (simplify) {
        out[, i] <- r
      } else out[[i]] <- r
      i <- i + 1
    }
  }
  if (simplify) {
    if (is.factor(x0)) {
      levels(out) <- levels(x0)
      # class(out) <- class(x0)
    }
    # dim(out) <- dim.use
  }
  out
}

combos <- big.combn(length(total_ma_names),2)
nmar <- length(total_ma_names)
big.combn.epi <- function (x, m=2, nmar,FUN = NULL, simplify = TRUE, ...) 
{
  stopifnot(length(m) == 1L, is.numeric(m))
  if (m < 0) 
    stop("m < 0", domain = NA)
  if (is.numeric(x) && length(x) == 1L && x > 0 && trunc(x) == 
      x) 
    x <- seq_len(x)
  n <- length(x)
  if (n < m) 
    stop("n < m", domain = NA)
  x0 <- x
  if (simplify) {
    if (is.factor(x)) 
      x <- as.integer(x)
  }
  m <- as.integer(m)
  e <- 0
  h <- m
  a <- seq_len(m)
  nofun <- is.null(FUN)
  if (!nofun && !is.function(FUN)) 
    stop("'FUN' must be a function or NULL")
  len.r <- length(r <- if (nofun) c(x[a][1]+1L,
                                    as.integer(x[a][2]+nmar),
                                    (x[a][1]+1L)*nmar-sum(1L:(x[a][1]-1L))+x[a][2]+1L)
                  else FUN(x[a], ...)
                  )
  count <- as.numeric(round(choose(n, m)))
  if (simplify) {
    dim.use <- if (nofun) 
      c(m, count)
    else {
      d <- dim(r)
      if (length(d) > 1L) 
        c(d, count)
      else if (len.r != 1L) 
        c(len.r, count)
      else c(d, count)
    }
  }
  if (simplify) {
    # out <- matrix(r, nrow = len.r, ncol = count)
    remove_bigmatrix("combn.epi")
    out <- filebacked.big.matrix(
      nrow = len.r,
      ncol = count,
      type = "integer", # char是c++的单个字符
      backingfile = "combn.epi.bin", 
      backingpath = dirname("combn.epi"), 
      descriptorfile = "combn.epi.des",
      dimnames = c(NULL, NULL)
    )
    out[,1] <- r
  } else {
    out <- vector("list", count)
    out[[1L]] <- r
  }
  if (m > 0) {
    i <- 2
    nmmp1 <- n - m + 1
    while (a[1] != nmmp1) {
      if (e < n - h) {
        h <- 1
        e <- a[m]
        j <- 1
      } else {
        e <- a[m - h]
        h <- h + 1
        j <- 1:h
      }
      a[m - h + j] <- e + j
      r <- if (nofun) {
        # x[a]
        c(x[a][1]+1L,
          as.integer(x[a][2]+nmar),
          (x[a][1]+1L)*nmar-sum(1L:(x[a][1]-1L))+x[a][2]+1L)
      } else FUN(x[a], ...)
      if (simplify) {
        out[, i] <- r
      } else out[[i]] <- r
      i <- i + 1
    }
  }
  if (simplify) {
    if (is.factor(x0)) {
      levels(out) <- levels(x0)
      # class(out) <- class(x0)
    }
    # dim(out) <- dim.use
  }
  out
}
epi_index <- big.combn.epi(length(total_ma_names),2,nmar = nmar)

remove_bigmatrix("epi_data")
epi_data <- filebacked.big.matrix(
  nrow = nrow(TRAN),
  ncol = ncol(combos)+2*nmar+ncol(Y),
  type = "double", # char是c++的单个字符
  backingfile = "epi_data.bin", 
  backingpath = dirname("epi_datai"), 
  descriptorfile = "epi_data.des",
  dimnames = c(NULL, NULL)
)

# 加载C++代码
sourceCpp(file = "./code_for_epi/source/epi.cpp")

# 调用 Rcpp 函数将 R 矩阵转换为 arma::mat
epi1 <- convertRMatrixToArmaMat(epi1)
epi2 <- convertRMatrixToArmaMat(epi2)
TRAN <- convertRMatrixToArmaMat(TRAN)

# 预计算lm需要的数据
tryCatch({
  ca_epi_data(epi_data@address, combos@address, epi1, epi2, TRAN, 
              Y, trdata1,trdata2, threads = 4)
}, error = function(e) {
  print(e)
})
# rm(epi1,epi2,Y,trdata1,trdata2,trAdata,trDdata,TRAN,combos)
# gc()
