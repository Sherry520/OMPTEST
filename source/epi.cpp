#include <RcppArmadillo.h>
// #include <RcppEigen.h>
#include <Rcpp.h>
#include <bigmemory/MatrixAccessor.hpp>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(bigmemory)]]

using namespace std;
using namespace Rcpp;
using namespace arma;
// using namespace Eigen;

// using namespace Eigen::Map;
// using namespace Eigen::MatrixXd;

// // function for convert R matrix class to arma::mat
// arma::mat rMatrixToArmaMat(const Rcpp::NumericMatrix& rMat) {
//   return Rcpp::as(rMat);
// }

// [[Rcpp::export]]
arma::mat TransferMatArma(Rcpp::NumericMatrix x) {
  arma::mat tx(x.begin(), x.nrow(), x.ncol(), false);
  return tx;
}

// arma::mat convertRMatrixToArmaMat(Rcpp::NumericMatrix rMat) {
//   return rMatrixToArmaMat(rMat);
// }

// epi_data: big.matrix to save result
// combos: all combn of epi1 and epi2

// Logic for ca_epi_data
template <typename T>
void ca_epi_data(Rcpp::XPtr<BigMatrix> epi_data,
                                  const Rcpp::XPtr<BigMatrix> combos,
                                  const arma::mat& epi1,
                                  const arma::mat& epi2,
                                  const arma::mat& TRAN,
                                  const Rcpp::NumericMatrix& Y,
                                  const Rcpp::NumericMatrix& trdata1,
                                  const Rcpp::NumericMatrix& trdata2,
                                  int threads=0) {
  omp_set_num_threads(threads);
  
  MatrixAccessor<T> epi_col = MatrixAccessor<T>(*epi_data);
  MatrixAccessor<int> combos_col = MatrixAccessor<int>(*combos);
  
  // Check dimensions
  if (epi_data->nrow() != TRAN.n_rows) {
    Rcpp::stop("Mismatch between TRAN rows and epi_data rows.");
  }
  if (epi_data->ncol() != combos->ncol()+1+2*trdata1.ncol()) {
    Rcpp::stop("Mismatch between combos columns and epi_data columns.");
  }
  
  // copy Y, trdata1 and trdata2 to epi_data
  #pragma omp parallel for
  for (size_t j = 0; j < Y.ncol(); j++) {
    if (j >= Y.ncol()) {
      Rcpp::stop("Col index exceeds Y cols.");
    }
    for (size_t k = 0; k < Y.nrow(); k++){
      if (k >= Y.nrow()) {
        Rcpp::stop("Row index exceeds trdata1 rows.");
      }
      // epi_col[i] is pointer,epi_col[i][j] is value.
      epi_col[j][k] = Y(k,j);
    }
  }
  #pragma omp parallel for
  for (size_t j = 0; j < trdata1.ncol(); j++) {
    if (j >= trdata1.ncol()) {
      Rcpp::stop("Col index exceeds trdata1 cols.");
    }
    for (size_t k = 0; k < trdata1.nrow(); k++){
      if (k >= trdata1.nrow()) {
        Rcpp::stop("Row index exceeds trdata1 rows.");
      }
      // epi_col[i] is pointer,epi_col[i][j] is value.
      epi_col[j+1][k] = trdata1(k,j);
    }
  }
  #pragma omp parallel for
  for (size_t j = 0; j < trdata2.ncol(); j++) {
    if (j >= trdata2.ncol()) {
      Rcpp::stop("Col index exceeds trdata1 cols.");
    }
    for (size_t k = 0; k < trdata2.nrow(); k++){
      if (k >= trdata2.nrow()) {
        Rcpp::stop("Row index exceeds trdata2 rows.");
      }
      // epi_col[i] is pointer,epi_col[i][j] is value.
      epi_col[j+1+trdata1.ncol()][k] = trdata2(k,j);
    }
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
       // epi_col[i] is pointer,epi_col[i][j] is value.
       epi_col[i+Y.ncol()+trdata1.ncol()+trdata2.ncol()][j] = result[j];
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
                 Rcpp::NumericMatrix& Y,
                 Rcpp::NumericMatrix& trdata1,
                 Rcpp::NumericMatrix& trdata2,
                 int threads=0){
  Rcpp::XPtr<BigMatrix> xpEpi_data(pEpi_data);
  Rcpp::XPtr<BigMatrix> xpCombos(pCombos);
  if(4 != xpCombos->matrix_type()){
    Rcpp::stop("big.matrix object of pCombos should be int type");
  }
  
  switch(xpEpi_data->matrix_type()) {
  case 1:
    return ca_epi_data<char>(xpEpi_data, xpCombos, epi1, epi2, TRAN,
                             Y, trdata1, trdata2, threads);
  case 2:
    return ca_epi_data<short>(xpEpi_data, xpCombos, epi1, epi2, TRAN,
                              Y, trdata1, trdata2, threads);
  case 4:
    return ca_epi_data<int>(xpEpi_data, xpCombos, epi1, epi2, TRAN,
                            Y, trdata1, trdata2, threads);
  case 8:
    return ca_epi_data<double>(xpEpi_data, xpCombos, epi1, epi2, TRAN,
                               Y, trdata1, trdata2, threads);
  default:
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}

// subset matrix from epi_data big.matrix with index from epi_index big.matrix
arma::mat extract_by_column_index(SEXP source_bigmat, int column_index, SEXP target_bigmat) {
  // 获取 source_bigmat 对象（用于提取索引）
  Rcpp::XPtr<BigMatrix> pSourceMat(source_bigmat);
  MatrixAccessor<int> source_col = MatrixAccessor<int>(*pSourceMat);
  
  size_t source_nrows = pSourceMat->nrow();
  size_t source_ncols = pSourceMat->ncol();
  
  // 检查列索引范围
  if (column_index < 0 || column_index >= source_ncols) {
    Rcpp::stop("Column index out of bounds.");
  }
  
  // 提取索引列
  Rcpp::IntegerVector indices(source_nrows + 1); // 加上首列Y
  indices[0] = 1;
  for (size_t i = 0; i < source_nrows; i++) {
    indices[i+1] = static_cast<int>(source_col[column_index][i]);
  }
  
  // 获取 target_bigmat 对象
  Rcpp::XPtr<BigMatrix> pTargetMat(target_bigmat);
  MatrixAccessor<double> target_col = MatrixAccessor<double>(*pTargetMat);
  
  size_t target_nrows = pTargetMat->nrow();
  size_t target_ncols = pTargetMat->ncol();
  
  // 创建返回矩阵：行数与目标矩阵一致，列数为索引大小
  arma::mat result(target_nrows, indices.size());
  
  // 填充返回矩阵
  for (size_t i = 0; i < indices.size(); i++) {
    int col_index = indices[i] - 1; // R 索引到 C++ 索引
    // Rcpp::Rcout << "col_index in C++ is " << col_index << std::endl;
    
    if (col_index < 0 || static_cast<size_t>(col_index) >= target_ncols) {
      Rcpp::Rcout << "col_index is " << col_index << std::endl;
      Rcpp::Rcout << "i is " << i << std::endl;
      Rcpp::stop("Index out of bounds in souce matrix.");
    }
    
    for (size_t j = 0; j < target_nrows; j++) {
      result(j, i) = target_col[col_index][j];
    }
  }
  // Rcpp::Rcout << "submatrix head 10 is " << result.rows(0,9) << std::endl;
  return result;
}

Rcpp::List lm_armadillo_with_p_values(arma::mat mat) {
  // // Eigen计算的se、t统计量、p值不对
  // Eigen::MatrixXd eigen_mat = Eigen::Map<Eigen::MatrixXd>(mat.memptr(),
  //                                                       mat.n_rows,
  //                                                       mat.n_cols);
  // Eigen::MatrixXd X = eigen_mat.block(0,1,eigen_mat.rows(),eigen_mat.cols()-1);
  // Eigen::VectorXd y = eigen_mat.col(0);
  // const int n(X.rows()), p(X.cols());
  // typedef Eigen::ColPivHouseholderQR<Eigen::MatrixXd> CPivQR;
  // typedef CPivQR::PermutationType Permutation;
  // const CPivQR PQR(X); // 使用 Householder 变换来进行 QR 分解
  // const Permutation Pmat(PQR.colsPermutation()); // 提取 QR 分解中的列置换矩阵
  // const int r(PQR.rank());
  // int df;
  // Rcpp::Rcerr << "X rank is: " << r << std::endl;
  // // Eigen::VectorXd X_sum12 = X.col(1) + X.col(2);
  // Eigen::VectorXd betahat(p), fitted, se(p),t_stats, p_values(p);
  // if (r == X.cols()) { // full rank case
  //   betahat = PQR.solve(y);
  //   fitted = X * betahat;
  //   se = Pmat * PQR.matrixQR().topRows(p).triangularView<Eigen::Upper>().solve(Eigen::MatrixXd::Identity(p, p)).rowwise().norm();
  //   // Calculate t-statistics
  //   t_stats = betahat.array() / se.array();
  //   df = n - p;
  //   // Calculate p-values
  //   // 遍历 t_stat 计算对应的 p 值
  //   for (size_t i = 0; i < t_stats.size(); ++i) {
  //     // 使用 abs(t_stat[i]) 计算单尾 p 值，并乘以 2 计算双尾 p 值
  //     p_values[i] = 2 * R::pt(-std::abs(t_stats[i]), df, true, false);
  //   }
  // } else {
  //   try{
  //     Eigen::MatrixXd Rinv = PQR.matrixQR().topLeftCorner(r, r).triangularView<Eigen::Upper>().solve(Eigen::MatrixXd::Identity(r, r));
  //     Rcpp::Rcerr << "Rinv is: " << Rinv << std::endl;
  //     Eigen::VectorXd effects(PQR.householderQ().adjoint() * y);
  //     Rcpp::Rcerr << "effects.head(r) is: " << effects.head(r) << std::endl;
  //     betahat.fill(NA_REAL);
  //     Rcpp::Rcout << "betahat.fill is: " << betahat << std::endl;
  //     Eigen::MatrixXd effect_head = effects.head(r);
  //     Rcpp::Rcout << "effects.head(r) =" << effect_head << std::endl;
  //     Rcpp::Rcout << "Rinv * effects.head(r) =" << Rinv * effect_head << std::endl;
  //     betahat.head(r) = Rinv * effects.head(r);
  //     Rcpp::Rcerr << "betahat.head is: " << betahat.head(r) << std::endl;
  //     Rcpp::Rcout << "dim of Pmat.cols()" << Pmat.cols() << std::endl;
  //     Rcpp::Rcout << "dim of Pmat.rows()" << Pmat.rows() << std::endl;
  //     betahat = Pmat * betahat;
  //     se.fill(::NA_REAL);
  //     se.head(r) = Rinv.rowwise().norm();
  //     se = Pmat * se;
  //     // create fitted values from effects
  //     effects.tail(X.rows() - r).setZero();
  //     fitted = PQR.householderQ() * effects;
  //     
  //     // Calculate t-statistics
  //     t_stats = betahat.array() / se.array();
  //     df = n-r;
  //     // Calculate p-values
  //     // Eigen::VectorXd p_values = 2 * Rcpp::pnorm(-(t_stats.array().abs()), r - 1); // Two-tailed test
  //     
  //     p_values.fill(NA_REAL); // 双尾 p 值
  //     // 遍历 t_stat 计算对应的 p 值
  //     for (size_t i = 0; i < t_stats.size(); ++i) {
  //       // 使用 abs(t_stat[i]) 计算单尾 p 值，并乘以 2 计算双尾 p 值
  //       p_values[i] = 2 * R::pt(-std::abs(t_stats[i]), df, true, false);
  //     }
  //   }catch (const std::exception& e) {
  //     Rcpp::Rcout << "Caught exception: " << e.what() << std::endl;
  //     Rcpp::stop("error occured in Rinv");
  //   }
  // 
  // 
  //   // Eigen::MatrixXd Rinv(PQR.matrixQR().topLeftCorner(r, r)
  //   //                  .triangularView<Eigen::Upper>().solve(Eigen::MatrixXd::Identity(r, r)));
  //   // Eigen::VectorXd effects(PQR.householderQ().adjoint() * y);
  //   // betahat.fill(::NA_REAL);
  //   // betahat.head(r) = Rinv * effects.head(r);
  //   // betahat = Pmat * betahat;
  //   // se.fill(::NA_REAL);
  //   // se.head(r) = Rinv.rowwise().norm();
  //   // se = Pmat * se;
  //   // // create fitted values from effects
  //   // effects.tail(X.rows() - r).setZero();
  //   // fitted = PQR.householderQ() * effects;
  // }
  // 
  // return Rcpp::List::create(
  //   // Rcpp::Named("input X is") = X,
  //   // Rcpp::Named("X[,1] + X[,2] is") = X_sum12,
  //   Rcpp::Named("coefficients") = betahat,
  //   Rcpp::Named("t value") = t_stats,
  //   // Rcpp::Named("fitted.values") = fitted,
  //   // Rcpp::Named("residuals") = resid,
  //   // Rcpp::Named("s") = s,
  //   Rcpp::Named("df.residual") = df,
  //   Rcpp::Named("rank") = r,
  //   Rcpp::Named("Std. Error") = se,
  //   Rcpp::Named("Pr(>|t|)") = p_values
  //   );
  

  // seprate matrix to y and x
  arma::mat X_mat = mat.cols(1, (mat.n_cols - 1));  // X 矩阵 (自变量)
  if(arma::rank(X_mat) < X_mat.n_cols-1){ // arma::rank计算的秩比实际小1
    Rcpp::Rcout << "rank(X_mat) is " << arma::rank(X_mat) << std::endl;
    Rcpp::Rcout << "head(X_mat) is " << X_mat.rows(0,5) << std::endl;
    X_mat = mat.cols(1, (mat.n_cols - 2));
    Rcpp::Rcout << "head(X_mat) is " << X_mat.rows(0,5) << std::endl;
  }
  arma::colvec y_vec = mat.col(0);  // y 向量 (因变量)
  // Rcpp::Rcout << "lm data of X_mat is " << X_mat << std::endl;
  // Rcpp::Rcout << "lm data of y_mat is " << y_vec << std::endl;

  // 检查矩阵维度
  if (X_mat.n_rows != y_vec.n_rows) {
    stop("Number of rows in X does not match length of y.");
  }


  // // 添加截距项
  // arma::mat X_aug = join_rows(ones(X_mat.n_rows), X_mat);
  // Rcpp::Rcout << "X_aug is " << X_aug.rows(0,10) << std::endl;
  // Rcpp::Rcout << "rank X_aug is " << arma::rank(X_aug) << std::endl;
  // // QR分解
  // arma::mat Q, R;
  // arma::qr(Q, R, X_aug);
  //
  // // 计算广义逆矩阵
  // arma::mat R_pinv = pinv(R);
  //
  // // 计算回归系数
  // arma::vec beta = R_pinv * Q.t() * y_vec;
  //
  // // 计算残差
  // arma::vec residuals = y_vec - X_aug * beta;
  //
  // // 计算残差平方和
  // double rss = dot(residuals, residuals);
  //
  // // 返回结果
  // Rcpp::List res = List::create(Named("coefficients") = beta,
  //                         // Named("residuals") = residuals,
  //                         Named("rss") = rss);
  // return res;


  // // 计算回归系数 (最小二乘法)
  // arma::colvec beta = solve(X_mat, y_vec);  // beta = (X'X)^(-1) X'y


  // 计算回归系数 (QR分解)
  arma::mat Q,R;
  arma::qr_econ(Q,R,X_mat); // 将 X 分解为 Q 和 R
  // arma::qr(Q,R,X_mat); // 将 X 分解为 Q 和 R
  // Rcpp::Rcout << "Q is " << Q << std::endl;
  Rcpp::Rcout << "R " << R << std::endl;

  // 解回归系数：beta = R^(-1) * Q^T * y

  // // 计算广义逆矩阵
  // arma::mat R_pinv = inv(R);

  // // 计算回归系数
  // arma::vec beta = R_pinv * Q.t() * y_vec;
  // 计算回归系数
  arma::vec beta = inv(R) * Q.t() * y_vec;

  // 计算预测值
  arma::colvec y_hat = X_mat * beta;

  // 计算残差
  arma::colvec residuals = y_vec - y_hat;

  // // 计算标准误差
  // arma::mat X_transpose = X_mat.t();
  // arma::mat var_beta = pinv(X_transpose * X_mat) * var(residuals);  //用广义逆 (X'X)^(-1) * Var(ε)
  // arma::colvec se_beta = sqrt(var_beta.diag());  // 提取标准误差
  // 
  // 计算残差方差
  double sigma_squared = dot(residuals, residuals) / (y_vec.n_elem - X_mat.n_cols);
  
  // 计算方差-协方差矩阵
  arma::mat X_transpose = X_mat.t();
  arma::mat var_beta = sigma_squared * inv(X_transpose * X_mat);
  
  // 计算标准误
  colvec se_beta = sqrt(var_beta.diag());

  // 计算 p 值 (双尾检验)
  int n = X_mat.n_rows;  // 样本数
  int p = X_mat.n_cols;  // 预测变量数
  int df = n - p;        // 自由度
  // 计算标准误差
  // double s2 = arma::dot(residuals, residuals) / (n -p); // std.errors of coefficients
  // arma::colvec se_beta = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X_mat)*X_mat)));

  // 计算 t 统计量
  arma::colvec t_stat = beta / se_beta;

  NumericVector p_values(t_stat.n_elem);  // 双尾 p 值
  // 遍历 t_stat 计算对应的 p 值
  for (size_t i = 0; i < t_stat.n_elem; ++i) {
    // 使用 abs(t_stat[i]) 计算单尾 p 值，并乘以 2 计算双尾 p 值
    p_values[i] = 2 * (1 - R::pt(std::abs(t_stat[i]), df,true,false));
  }

  // // 计算TSS
  // double mean_y = mean(y_vec);
  // vec centered_y = y_vec - mean_y;
  // double TSS = dot(centered_y, centered_y);
  // 
  // // 计算RSS
  // double RSS = dot(residuals, residuals);
  // 
  // // 计算SSR
  // double SSR = TSS - RSS;
  // 
  // // 计算R-squared
  // double R_squared = SSR / TSS;
  // 
  // // 计算调整后的R-square
  // double adjusted_R_squared = 1 - (1 - R_squared) * (n - 1) / (n - p - 1);

  // 返回结果
  return List::create(
    Named("coefficients") = beta,
    // Named("residuals") = residuals,
    Named("standard_errors") = se_beta,
    Named("t_stat") = t_stat,
    Named("p_values") = p_values
    // Named("r_squared") = R_squared,
    // Named("Adjusted r_squared") = adjusted_R_squared
  );
}

// Logic for epi_eff_pval

// [[Rcpp::export]]
Rcpp::List ca_epi_pval_eff(Rcpp::XPtr<BigMatrix> epi_pval,
                 Rcpp::XPtr<BigMatrix> epi_eff,
                 const Rcpp::XPtr<BigMatrix> epi_data,
                 const Rcpp::XPtr<BigMatrix> epi_index,
                 int threads=0) {
  omp_set_num_threads(threads);
  
  MatrixAccessor<double> pval_col = MatrixAccessor<double>(*epi_pval);
  MatrixAccessor<double> eff_col = MatrixAccessor<double>(*epi_eff);
  MatrixAccessor<double> epi_col = MatrixAccessor<double>(*epi_data);
  MatrixAccessor<int> epi_index_col = MatrixAccessor<int>(*epi_index);
  
  arma::icolvec subdata_epi_index;
  Rcpp::List result;
  size_t j=0;
  
  // Rcpp::Environment stats("package:stats");
  // Rcpp::Function lm = stats["lm"];
  // Rcpp::Function as_formula = stats["as.formula"];
  // Rcpp::Environment base("package:base");
  // Rcpp::Function colnames = base["colnames"];
  // Rcpp::Function c = base["c"];
  
  // for(size_t i=0; i < epi_index->ncol();i++) {
    
  for(size_t i=2658; i < 2659;i++) {
    try {
      arma::mat subdata = extract_by_column_index(epi_index,i,epi_data);
      result = lm_armadillo_with_p_values(subdata);
      // // NumericMatrix rsubdata =Rcpp::wrap(subdata);
      // Rcpp::NumericMatrix mat(subdata.n_rows, subdata.n_cols-1);
      // std::copy(subdata.begin(), subdata.end(), mat.begin());
      // 
      // int nrow = mat.nrow();
      // int ncol = mat.ncol();
      // 
      // NumericVector first_col(nrow);
      // NumericMatrix rest_mat(nrow, ncol - 1);
      // 
      // // 提取第一列
      // for (int i = 0; i < nrow; ++i) {
      //   first_col[i] = mat(i, 0);
      // }
      // 
      // // 提取其余列
      // for (int i = 0; i < nrow; ++i) {
      //   for (int j = 1; j < ncol; ++j) {
      //     rest_mat(i, j - 1) = mat(i, j);
      //   }
      // }
      // 
      // // 获取矩阵的属性列表
      // List attr_list = rest_mat.attr("attributes");
      // 
      // // 添加新的列名属性
      // attr_list["colnames"] = CharacterVector::creat("M1","M2","epi");
      // 
      // // 将属性列表赋给矩阵
      // rest_mat.attr("attributes") = attr_list;
      // 
      // // 将C++数据转换为R对象
      // Environment base("package:stats");
      // Function lm_fun = base["lm"];
      // Function summary_fun = base["summary"];
      // 
      // // 拟合线性模型
      // Rcpp::List model = lm_fun(Named("formula", y ~ M1+M2+epi), Named("data", DataFrame::create(Named("y", first_col), Named("X", rest_mat))));
      
      // // 获取系数
      // Rcpp::NumericVector coefficients = summary_fun(model)["coefficients"];
  
      
      // double a = 1/3;
      // double b =-1/3;
      // Rcpp::Rcout << " a 1/3 + b -1/3 = " << a+b  << std::endl;
      // // 下边这样写会报错<<重载
      // Rcpp::Rcout << "lm result for epi_index " << i << " is " << result << std::endl;
    // 正常处理代码
    } catch (const std::exception& e) {
      j++;
      Rcpp::Rcout << "Processed " << i << "st epi_index error"<< std::endl;
      Rcpp::Rcout << "Caught exception: " << e.what() << std::endl;
      // 可以选择在这里处理异常，比如记录日志或跳过
      continue; // 跳过当前迭代，继续下一个
    }
    
    
    
    // // 打开文件
    // std::ofstream file("lm_result.list.txt");
    // if (!file.is_open()) {
    //   Rcpp::stop("Unable to open file!");
    // }
    // 
    // // 遍历列表并写入文件
    // for (int i = 0; i < result.size(); ++i) {
    //   file << "Element " << i << ": " << result[i] << std::endl;
    // }
    // 
    // // 关闭文件
    // file.close();
  }
  Rcpp::Rcout << j << "caught exception times" << std::endl;
  // return model;
  return result;
  // return NULL;
    
    // // 获取big.matrix的内存地址
    // double* mat = (double*)xpMat->matrix();
    // 
    // // 创建arma::mat对象，不复制内存
    // arma::mat M(mat, xpMat->nrow(), xpMat->ncol(), false);
    // 
    // // 提取指定列
    // return M.cols(colidx);
    // 
    // 
    // // matrix *
    // arma::colvec result = TRAN * hadamard_product;
    // 
    // // apply the result to big.matrix
    // if (result.n_elem != epi_data->nrow()) {
    //   Rcpp::stop("Result vector length does not match epi_data row count.");
    // }
    // for (size_t j = 0; j < result.n_elem; j++) {
    //   if (j >= epi_data->nrow()) {
    //     Rcpp::stop("Row index exceeds epi_data rows.");
    //   }
    //   // epi_col[i] is pointer,epi_col[i][j] is value.
    //   epi_col[i+Y.ncol()+trdata1.ncol()+trdata2.ncol()][j] = result[j];
    // }
    
    // #pragma omp critical
    // {
    //   Rcpp::Rcout << "Processed column " << i << " of " << m << std::endl;
    // }
  // }
  // return epi_data;
}
