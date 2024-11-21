// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(bigmemory)]]

#include <RcppArmadillo.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <omp.h>

using namespace std;
using namespace Rcpp;

// function for convert R matrix class to arma::mat
arma::mat rMatrixToArmaMat(const Rcpp::NumericMatrix& rMat) {
  return as<arma::mat>(rMat);
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
  omp_set_num_threads(threads);
  
  MatrixAccessor<T> epi_col = MatrixAccessor<T>(*epi_data);
  MatrixAccessor<int> combos_col = MatrixAccessor<int>(*combos);
  
  // Check dimensions
  if (epi_data->nrow() != TRAN.n_rows) {
    Rcpp::stop("Mismatch between TRAN rows and epi_data rows.");
  }
  if (epi_data->ncol() != combos->ncol()) {
    Rcpp::stop("Mismatch between combos columns and epi_data columns.");
  }

  #pragma omp parallel for
  for(size_t i=0; i < combos->ncol();i++) {
    arma::colvec col_epi1 = epi1.col(combos_col[i][0]-1);
    arma::colvec col_epi2 = epi2.col(combos_col[i][1]-1);
    
    // calculate hadamard product
    arma::colvec hadamard_product = col_epi1 % col_epi2;
    
    // matrix *
    arma::colvec result = TRAN * hadamard_product;
      
    // apply the result to big.matrix
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
