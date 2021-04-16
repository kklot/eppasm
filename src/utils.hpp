#include <iostream>
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/CXX11/Tensor>
#pragma once

namespace epp {
  template<typename T> using scalar  = Eigen::Tensor<T, 0>;
  template<typename T> using vector  = Eigen::Tensor<T, 1>;
  template<typename T> using matrix  = Eigen::Tensor<T, 2>;
  template<typename T> using cube    = Eigen::Tensor<T, 3>;
  template<typename T> using array4D = Eigen::Tensor<T, 4>;
  template<typename T> using array5D = Eigen::Tensor<T, 5>;
  template<size_t x  > using index   = Eigen::array<Eigen::Index, x>;
  // Mapping
  template<typename T> using map_of_vector = Eigen::TensorMap<vector<T>>;
  template<typename T> using map_of_matrix = Eigen::TensorMap<matrix<T>>;
  template<typename T> using map_of_cube   = Eigen::TensorMap<cube<T>>;
  template<typename T> using map_of_4D     = Eigen::TensorMap<array4D<T>>;
  template<typename T> using map_of_5D     = Eigen::TensorMap<array5D<T>>;
  // tensor contractions and manipulations
  index<1> by_rows({0});
  index<1> by_cols({1});
  index<1> by_rack({2});
  index<2> by_bay({0, 1}); // row and col
  index<2> transpose({1, 0});

  // power method
  template<typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1> 
  power_method(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& Jac, 
               T e = 1e-6, T R0 = 1., int max_iter = 1e3)
  {
    Eigen::SparseMatrix<T> Jac_sparse = Jac.sparseView();
    Eigen::Matrix<T, Eigen::Dynamic, 1> v(Jac.rows()); 
    v.setRandom();
    T r = R0, u;
    for (int i = 0; i < max_iter; i++) 
    {
      v = Jac_sparse * v;
      u = v.norm();
      v = v / u;
      if (std::abs( r - u ) > e) {
        if (i==max_iter)
          Rf_warning("stop at iter:", max_iter, "with error: ", std::abs(r - u));
        r = u;
      } else {
        break;
      }
    }
    return v;
  }
} 

SEXP get_value(SEXP list, const char *str) {
  SEXP out = R_NilValue, names = GET_NAMES(list);
  int l_len = GET_LENGTH(list);
  for (int i = 0; i < l_len; i++ ) {
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) {
      out = VECTOR_ELT(list, i);
        break;
    }
  }
  if ( out == R_NilValue )
    Rf_warning("%s is NULL, check ?prepare_fp_for_Cpp", str);
  return out;
}

bool has_value(SEXP list, const char *str) {
  SEXP names = GET_NAMES(list);
  int l_len = GET_LENGTH(list);
  for (int i = 0; i < l_len; i++ )
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) {
      if (VECTOR_ELT(list, i) == R_NilValue) {
        Rf_warning("%s is NULL", str);
        return false;
      }
      else
        return true;
    }
  Rf_warning("%s does not exist.", str);
  return false;
}

template<typename Type>
epp::cube<Type> expand_2to3 (const epp::matrix<Type>& b, const int replicates)
{
  epp::index<3> cube_dims = {replicates, b.dimension(0), b.dimension(1)};
  epp::cube<Type> out = b
    .reshape(epp::index<2>({1, (int) b.size()})) // transpose
    .broadcast(epp::index<2>({replicates, 1}))
    .reshape(cube_dims); // replicate columns
  return out; // use for mulitplying tensor e.g. DS X AGE X SEX * AGE x SEX
}

template<typename Type>
epp::cube<Type> expand_right_2to3 (const epp::matrix<Type>& b, const int replicates)
{
  epp::index<3> cube_dims = {b.dimension(0), b.dimension(1), replicates};
  epp::cube<Type> out = b
    .reshape(epp::index<1>({(int) b.size()})) // prolong
    .broadcast(epp::index<1>({replicates}))
    .reshape(cube_dims); // replicate
  return out; // use for mulitplying tensor e.g. AGE X SEX x DS * AGE x SEX
}

template<typename Type>
epp::matrix<Type> expand_2to4 (const epp::matrix<Type>& b, const int replicates)
{
  epp::matrix<Type> out = b
    .reshape(Eigen::array<int, 2>({1, (int) b.size()})) // to one long row
    .broadcast(Eigen::array<int, 2>({replicates, 1})); // replicate each
  return out; // use for mulitplying tensor TS X DS X AGE X SEX * AGE x SEX
} 

typedef std::vector<double>    dvec;
typedef std::vector<int>       ivec;

Eigen::array<long, 1> get_dim_1D(SEXP r_in, const char *str) {
  int len = GET_LENGTH(get_value(r_in, str));
  Eigen::array<long, 1> out = { len };
  return out;
}

Eigen::array<long, 2> get_dim_2D(SEXP r_in, const char *str) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(get_value(r_in, str), R_DimSymbol));
  if (r_dims == R_NilValue)
    Rf_error("%s dimension is wrong! expecting a 2-D", str);
  int * rdims = INTEGER(r_dims);
  UNPROTECT(1);
  Eigen::array<long, 2> out = { rdims[0], rdims[1] };
  return out;
}

Eigen::array<long, 3> get_dim_3D(SEXP r_in, const char *str) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(get_value(r_in, str), R_DimSymbol));
  if (r_dims == R_NilValue || Rf_length(r_dims) != 3)
    Rf_error("%s dimension is wrong! expecting a 3-D", str);
  int *rdims = INTEGER(r_dims);
  UNPROTECT(1);
  Eigen::array<long, 3> out = { rdims[0], rdims[1], rdims[2] };
  return out;
}

Eigen::array<long, 4> get_dim_4D(SEXP r_in, const char *str) {
  SEXP r_dims = Rf_protect(Rf_getAttrib(get_value(r_in, str), R_DimSymbol));
  if (r_dims == R_NilValue || Rf_length(r_dims) != 4)
    Rf_error("%s dimension is wrong! expecting a 4-D", str);
  int *rdims = INTEGER(r_dims);
  UNPROTECT(1);
  Eigen::array<long, 4> out = {rdims[0], rdims[1], rdims[2], rdims[3]};
  return out;
}

// Tensor array NA/INF to zero: num/0.0 or 0.0/0.0
template <class K>
void replace_na_with (K& A, double B = 0.0) {
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
    if (std::isnan(*i) || std::isinf(*i))
      *i = B;
}

// Tensor array all to zero
template <class K>
void zeroing (K& A) {
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
      *i = 0.0;
}

// Tensor sum
template <class K>
double sumArray (K& A) {
  double sum = 0;
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
    sum += *i;
  return sum;
}

// stand vector sum
template <class K>
double sum_vector (K& A) {
  double sum = 0;
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
    sum += *i;
  return sum;
}

// Tensor array multiply each cell
template <class K>
K multiply_with (K A, double X) { // not by reference
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
      *i *= X;
  return A;
}

// Tensor array multiply each cell
template <class K>
void multiply_with_inplace (K& A, double X) { // by reference
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
      *i *= X;
}

// Tensor array multiply each cell
template <class K>
K add_to_each (K A, double X) { // not by reference
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
      *i += X;
  return A;
}

// Tensor array multiply each cell
template <class K>
void add_to_each_inplace (K& A, double X) { // by reference
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
      *i += X;
}

// Tensor array substract each cell
template <class K>
K substract_from_each (K A, double X) { // not by reference
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
      *i -= X;
  return A;
}

// Tensor array fill
template <class K>
K fill_with (K A, double B) {
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
      *i = B;
  return A;
}

// Tensor array fill inplace
template <class K>
void fill_with_inplace (K& A, double B) {
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
      *i = B;
}

template <class K>
double min (K& A) {
  double * mini = std::min_element( A.origin(), A.origin() + A.size());
  return * mini;
}

template <class K>
double max (K& A) {
  double * maxi = std::max_element( A.origin(), A.origin() + A.size());
  return * maxi;
}

template <class K>
void replace_elem_with (K& A, double B = 1, double C = 1) {
  for (auto i = A.data(); i < (A.data() + A.size()); ++i)
    if (*i > B) *i = C;
}

template<typename Type>
Eigen::Tensor<Type, 2> sumByAG2D (
  const Eigen::Tensor<Type, 2>& B, 
  const Eigen::Tensor<Type, 1>& pace_v, 
  int new_row) 
{
  Eigen::Tensor<Type, 2> A(new_row, B.dimension(1));
  A.setZero();
  for (int i = 0; i < B.dimension(1); ++i) {
    int r_i = pace_v(0); // first age group
    for (int j = 0; j < B.dimension(0); ++j) {
      if (pace_v(j) != r_i) 
        ++r_i;
      A(r_i - 1, i) += B(j, i);
    } // end age-groups
  } // end columns
  return A;
}

template<typename Type>
Eigen::Tensor<Type, 1> sumByAG1D (
  const Eigen::Tensor<Type, 1>& B, 
  const Eigen::Tensor<Type, 1>& pace_v, 
  int new_row)
{
  Eigen::Tensor<Type, 1> A(new_row);
  A.setZero();
  int r_i = pace_v(0); // first age group
  for (int i = 0; i < B.dimension(0); ++i) {
    if (pace_v(i) != r_i)
      ++r_i;
    A(r_i - 1) += B(i);
  } // end age-groups
  return A;
}

template<typename Type>
Eigen::Tensor<Type, 2> expand_age_group(
  const Eigen::Tensor<Type, 2>& B, 
  const Eigen::Tensor<Type, 1>& pace_v) 
{
  int new_row = sumArray(pace_v);
  Eigen::Tensor<Type, 2> A(new_row, B.dimension(1));
  for (int c = 0; c < B.dimension(1); c++) 
  {
    int cusu = 0;
    for (int r = 0; r < B.dimension(0); r++) 
    {
      for (int id = 0; id < pace_v(r); id++)
        A(cusu + id, c) = B(r, c);
      cusu += pace_v(r);
    }
  }
  return A;
}
