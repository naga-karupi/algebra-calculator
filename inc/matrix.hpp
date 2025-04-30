/**
 * @file matrix.cpp
 * @author Nishinaga Rikuto (you@domain.com)
 * @brief definition of matrix class
 * @version 0.1
 * @date 2025-04-24
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#pragma once

#include <complex>
#include <vector>
#include <tuple>

namespace my_mt
{

/// todo: i wanna use concept but c++20 is not supported,
/// so i will use sfinae...?
/// if initialized ROW = 0, COLUMNS = 0, use vector(heap memory)
template <typename TYPE, size_t ROW = 0, size_t COL = 0>
class Matrix 
{
    TYPE m_data[ROW][COL];
    ::std::vector<::std::vector<TYPE>> m_data_vector;

    inline TYPE& get(size_t row, size_t col);
public:
    Matrix() = default;
    Matrix(const TYPE& value);
    Matrix(const TYPE (&arr)[ROW][COL]);
    Matrix(const size_t row, size_t col);
    Matrix<TYPE, ROW, COL>& operator = (const Matrix<TYPE, ROW, COL>& other);
    TYPE& operator () (size_t row, size_t col);

    inline size_t rows() const noexcept;
    inline size_t cols() const noexcept;

    Matrix<TYPE, ROW, COL> operator + (const Matrix<TYPE, ROW, COL>& other) const;
    Matrix<TYPE, ROW, COL> operator - (const Matrix<TYPE, ROW, COL>& other) const;
    Matrix<TYPE, ROW, COL> operator * (const TYPE& scalar) const;
    Matrix<TYPE, ROW, COL> operator / (const TYPE& scalar) const;
    Matrix<TYPE, ROW, COL> operator * (const Matrix<TYPE, COL, ROW>& other) const;

    inline bool isZero() const;
    inline bool isIdentity() const;
    inline bool isSymmetric() const;

    inline static Matrix<TYPE, ROW, COL> Identity() noexcept;
    inline static Matrix<TYPE, ROW, COL> Zero() noexcept;

    TYPE determinant() const;
    TYPE trace() const;
    TYPE norm() const;
    TYPE squareNorm() const;

    Matrix<TYPE, ROW, COL>  inverse() const;
    Matrix<TYPE, ROW, COL>  transpose() const;
    Matrix<TYPE, ROW, COL>  adjoint() const;

    Matrix<TYPE, ROW, COL>  block(size_t row, size_t col, size_t rows, size_t cols) const;
};

#include "LU_decomposition.hpp"
#include "detail/matrix.ipp"


} // namespace my_mt
