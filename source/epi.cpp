// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory)]]

#include <RcppArmadillo.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <omp.h>


using namespace std;
using namespace Rcpp;

// 函数声明，接受一个 R 矩阵并返回一个 arma::mat 对象
arma::mat rMatrixToArmaMat(const Rcpp::NumericMatrix& rMat) {
  return as<arma::mat>(rMat); // 使用 as<> 函数进行转换
}

// [[Rcpp::export]]
arma::mat convertRMatrixToArmaMat(Rcpp::NumericMatrix rMat) {
  return rMatrixToArmaMat(rMat);
}

// epi_data: big.matrix to save result
// combos: all combn of epi1 and epi2

// Logic for ca_epi_data
template <typename T>
void ca_epi_data(Rcpp::XPtr<BigMatrix> epi_data,
                                  const Rcpp::XPtr<BigMatrix> combos,
                                  const arma::mat& epi1,
                                  const arma::mat& epi2,
                                  const arma::mat& TRAN,
                                  int threads=0) {
  // omp_set_num_threads(threads);
  
  MatrixAccessor<T> epi_col = MatrixAccessor<T>(*epi_data);
  MatrixAccessor<int> combos_col = MatrixAccessor<int>(*combos);
  
  // Check dimensions
  if (epi_data->nrow() != TRAN.n_rows) {
    Rcpp::stop("Mismatch between TRAN rows and epi_data rows.");
  }
  if (epi_data->ncol() != combos->ncol()) {
    Rcpp::stop("Mismatch between combos columns and epi_data columns.");
  }
  
  // Rcpp::Rcout << "combos_col[1] is  " << combos_col[1][0] << std::endl;
  Rcpp::Rcout << "epi1 index 0 0 is  " << combos_col[0][0] << std::endl;
  Rcpp::Rcout << "epi1 index 0 1 is  " << combos_col[0][1] << std::endl;
  Rcpp::Rcout << "epi1 index 1 0 is  " << combos_col[1][0] << std::endl;
  Rcpp::Rcout << "epi1 index 1 1 is  " << combos_col[1][1] << std::endl;

  // #pragma omp parallel for
  for(size_t i=0; i < combos->ncol();i++) {
    // 测试哈达玛积功能
    // // 创建两个 arma::vec 类型的向量
    // arma::colvec u = {1, 2, 3};
    // arma::colvec v = {4, 5, 6};
    // 
    // // 计算哈达玛积（逐元素乘积）
    // arma::colvec result2 = u % v; // result 将会是 {4, 10, 18}
    // 
    // // 打印结果
    // std::cout << "Hadamard product: " << result2 << std::endl;
    
    arma::colvec col_epi1 = epi1.col(combos_col[i][0]-1);
    arma::colvec col_epi2 = epi2.col(combos_col[i][1]-1);
    // Rcpp::Rcout << "col_epi1 index " << i  << " is " << epi2.col(combos_col[i-1][1]) << std::endl;
    // Rcpp::Rcout << "epi2 index is  " << combos_col[i][1] << std::endl;
    
      // calculate hadamard product
      arma::colvec hadamard_product = col_epi1 % col_epi2;
      // std::cout << "Hadamard product: " << hadamard_product << std::endl;
      // matrix *
      arma::colvec result = TRAN * hadamard_product;
      
      // apply the result to big.matrix
      // epi_col[i] = result;
      // 将 result 写入到 big.matrix
      if (result.n_elem != epi_data->nrow()) {
        Rcpp::stop("Result vector length does not match epi_data row count.");
      }
      for (size_t j = 0; j < result.n_elem; j++) {
        if (j >= epi_data->nrow()) {
          Rcpp::stop("Row index exceeds epi_data rows.");
        }
        epi_col[i][j] = result[j];  // 两者数据类型不一样，得将 result[j] 写入 epi_data 的第 i 列第 j 行
    }
      // #pragma omp critical
      // {
      //   Rcpp::Rcout << "n_elem is  " << result.n_elem << std::endl;
      //   Rcpp::Rcout << "Processed column " << i << " of " << m << std::endl;
      // }
  }
  // return epi_data;
}
  
// [[Rcpp::export]]
void ca_epi_data(SEXP pEpi_data,
                 SEXP pCombos,
                 arma::mat& epi1,
                 arma::mat& epi2,
                 arma::mat& TRAN,
                 int threads=0){
  XPtr<BigMatrix> xpEpi_data(pEpi_data);
  XPtr<BigMatrix> xpCombos(pCombos);
  if(4 != xpCombos->matrix_type()){
    Rcpp::stop("big.matrix object of pCombos should be int type");
  }
  
  switch(xpEpi_data->matrix_type()) {
  case 1:
    return ca_epi_data<char>(xpEpi_data, xpCombos, epi1, epi2, TRAN, threads);
  case 2:
    return ca_epi_data<short>(xpEpi_data, xpCombos, epi1, epi2, TRAN, threads);
  case 4:
    return ca_epi_data<int>(xpEpi_data, xpCombos, epi1, epi2, TRAN, threads);
  case 8:
    return ca_epi_data<double>(xpEpi_data, xpCombos, epi1, epi2, TRAN, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}
