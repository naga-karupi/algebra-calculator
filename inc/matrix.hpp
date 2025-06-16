/**
 * @file matrix.cpp
 * @author Nishinaga Rikuto
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
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <execution>


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

    
    inline TYPE& set(size_t row, size_t col);
public:
    Matrix() = default;
    Matrix(const TYPE& value);
    Matrix(const TYPE (&arr)[ROW][COL]);
    Matrix(const size_t row, size_t col);
    Matrix(const Matrix<TYPE, ROW, COL>& other) = default;
    Matrix(const Matrix<TYPE, ROW, COL>&& other);
    Matrix<TYPE, ROW, COL>& operator = (const Matrix<TYPE, ROW, COL>& other);
    template <size_t ROW2, size_t COL2>
    Matrix<TYPE, 0, 0>& operator = (const Matrix<TYPE, ROW2, COL2>& other);

    TYPE& operator () (size_t row, size_t col);
    TYPE& operator () (int row, int col);

    TYPE& operator () (size_t n);
    TYPE& operator () (int n);

    inline TYPE get(size_t row, size_t col) const;

    inline size_t rows() const noexcept;
    inline size_t cols() const noexcept;

    // Matrix<TYPE, ROW, COL> operator + (const Matrix<TYPE, ROW, COL>& other) const;
    // Matrix<TYPE, ROW, COL> operator + (const Matrix<TYPE, 0, 0>& other) const;
    template <size_t ROW2, size_t COL2>
    Matrix<TYPE, ROW, COL> operator + (const Matrix<TYPE, ROW2, COL2>& other) const;

    // Matrix<TYPE, ROW, COL> operator - (const Matrix<TYPE, ROW, COL>& other) const;
    // Matrix<TYPE, ROW, COL> operator - (const Matrix<TYPE, 0, 0>& other) const;
    template <size_t ROW2, size_t COL2>
    Matrix<TYPE, ROW, COL> operator - (const Matrix<TYPE, ROW2, COL2>& other) const;

    Matrix<TYPE, ROW, COL> operator * (const TYPE& scalar) const;
    Matrix<TYPE, ROW, COL> operator / (const TYPE& scalar) const;

    inline bool isZero() const;
    inline bool isIdentity() const;
    inline bool isSymmetric() const;

    inline static Matrix<TYPE, ROW, COL> Identity() noexcept;
    inline static Matrix<TYPE, 0, 0> Identity(size_t row, size_t col) noexcept;
    inline static Matrix<TYPE, ROW, COL> Zero() noexcept;
    inline static Matrix<TYPE, 0, 0> Zero(size_t row, size_t col) noexcept;

    TYPE determinant() const;
    TYPE trace() const;
    TYPE norm() const;
    TYPE squaredNorm() const;

    Matrix<TYPE, ROW, COL> inverse() const;
    Matrix<TYPE, ROW, COL> transpose() const;
    Matrix<TYPE, ROW, COL> adjoint() const;

    Matrix<TYPE, 0, 0> block(size_t row, size_t col, size_t rows, size_t cols) const;

    void resize(size_t row, size_t col);

    void fill(const TYPE& value);
};

#include "helper_class.hpp"

template <typename TYPE, size_t ROW, size_t COL>
void print(const Matrix<TYPE, ROW, COL>& matrix);

template <typename TYPE, size_t ROW, size_t COL,size_t ROW2, size_t COL2>
typename ResultMatrixType<TYPE, ROW, COL, ROW2, COL2>::type operator * (const Matrix<TYPE, ROW2, COL2>& other);

template <typename TYPE, size_t SIZE>
using Vector = Matrix<TYPE, SIZE, 1>;
template <typename TYPE, size_t SIZE>
using VectorA = Matrix<TYPE, SIZE, 0>;




#include "LU_decomposition.hpp"
#if MULTI == 0
#include "multi.hpp"
#endif
#include "detail/matrix.ipp"




} // namespace my_mt
