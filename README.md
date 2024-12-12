高速、低内存占用的上位性分析

## 为什么做这个？
- 随着测序技术的发展，分子标记现在已经至少百万的水平，这对GWAS分析的计算挑战很大；
- 我在GWAS分析中，需要对两两SNP标记的上位性效应进行分析，模型目前已经简化成了回归模型；但是SNP两两组合的规模很大C(n,2)，当标记数n达到7万时，这个组合数已经超过R 32-bit 整数的范围；
- 对于上万行*上万列的矩阵，需要的内存很大，再加上还要进行矩阵运算，内存难以满足，即使1T的内存；
- 目前，我提取了一部分数据（大约3万标记）用R做lm分析，需要8-10天才能完成计算。偶尔服务器还会因为不明原因导致的内存不足半路失败。
## 要优化的for代码
```r
for_epi <- function(TRAN,eff_type,epi1,epi2,trdata1,trdata2,total_ma_names,trait,Y){
  nmar <- length(total_ma_names)
  res_eff <- matrix(0.0,nmar,nmar) #上三角存eff，下三角存pval
  rownames(res_eff) <- total_ma_names
  colnames(res_eff) <- total_ma_names

  for (m in 1:(nmar-1)) {
    for (n in (m+1):nmar) {
      epidesign <- TRAN%*%(epi1[,m]*epi2[,n])
      subdata <- data.frame(y=Y,trdata1[,m],trdata2[,n],epi=epidesign)
      colnames(subdata)[2:3] <- total_ma_names[c(m,n)
      fit.lm <- lm(as.formula(paste("y ~ -1 +",total_ma_names[m],"+",total_ma_names[n],"+ epi")),
                   data=subdata)
      infomat <- summary(fit.lm)$coefficient
      if ("epi"%in%rownames(infomat))
      {
        res_eff[m,n] <- summary(fit.lm)$coefficient["epi",1]
        res_eff[n,m] <- summary(fit.lm)$coefficient["epi",4]
      }else{
        res_eff[n,m] <- 1
      }
    }
    print(paste(trait,eff_type,"Epistasis scan invloving Marker",m,"completed."))
  }
  print(paste(trait,eff_type,"Epistasis completed."))
  write.table(res_eff,paste0("MapQTL_MPH_effect_pval_",trait,"_for_epi_",eff_type,".txt"),quote=FALSE)
  # write.table(res_pval,paste0("MapQTL_MPH_pval_",trait,"_for_epi_",eff_type,".txt"),quote=FALSE)
}
```
## 如何优化这个计算
- 整个分析流程，有两处计算耗时、耗内存的地方。 一个是构建两两标记的上位性效应，我把它单独拆成一个任务，另一个是lm回归分析；
- `bigmemory`可以通过内存映射的方式，将数据存到硬盘，减少对内存的占用；
- 用`C++`写代码，退而求其次，`Rcpp`；
- `Rcpp`中提供了`RcppArmadillo`包实现矩阵计算、特征分解做lm；
- `omp`可以加速`for`循环，充分利用多核优势。（我在`R`中尝试了`foreach+dopar`的方案【因为回归每次数据都不一样，别的并行都不合适，反正我没想出来】，可惜slurm集群中内存不会共享，除非自己有个很大内存的服务器才勉强能跑）

## 代码解析
### R中的代码
- 加载已经经过R各种处理过了的数据，用500个SNP做测试； https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/mpi_epi.r#L21
- 当需要重新跑代码时，这个函数用来删除bigmemory的backing文件； https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/mpi_epi.r#L51
- 修改了 `R utils` 中的`combn`函数，让它可以生成超过32-bit整数的组合，矩阵存到硬盘上，每列是对应的组合，如(1,2),(1,3),...,(1,500),(2,3),...,(499,500)；https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/mpi_epi.r#L84
- 将模型中Y,SNP1,SNP2,SNP1*SNP2这四部分矩阵组合成一个大矩阵，lm分析时，直接从中提取对应列；
- 生成了lm分析用到的那3列X的索引，因为Y始终在第一列，如(2,503,1003),(2,504,1004),...,(2,1001,1500),(3,504,1501),...,(500,1001,125751)； https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/mpi_epi.r#L173
- 转化`R`矩阵为`arma`的矩阵；https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/mpi_epi.r#L288-L290
- 预生成上边说的lm需要用的所有数据构成的矩阵；https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/mpi_epi.r#L291-L297
- 用来保存lm分析的结果主要是p值和effect值；https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/mpi_epi.r#L327
- 计算p值和eff值，先提取回归分析用到的数据，然后做lm分析；https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/mpi_epi.r#L359

### Rcpp中的代码
- 构造lm分析用到的数据，合并成一个大矩阵；https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/source/epi.cpp#L39
- 回归分析，不知道为啥arma计算的秩比实际小1，加阈值1e-7也不行；Eigen的计算没弄明白该怎么改，废弃掉了；https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/source/epi.cpp#L300
- 把上边这两步包起来，做for循环，并行；https://github.com/Sherry520/omp_epi/blob/5250a7bc2badc6706f4c8fd757cdd491dff0726f/source/epi.cpp#L538

## 可行性（有希望）
- 单线程的情况下，内存降低上百万倍，速度提升7倍
![image](https://github.com/user-attachments/assets/bfd1b9dd-438d-4372-bdbe-ee4229b0e21d)
![image](https://github.com/user-attachments/assets/b9f83834-b9cb-40f2-bad3-ff419f565dee)
- 之前好像就成功了这么一次，用的WSL的restudio-server，后来再也不行了，好像是把标记缩减到100跑出的
  ![image](https://github.com/user-attachments/assets/a9a2be78-d992-45db-a4d8-7aa428a7bb95)
