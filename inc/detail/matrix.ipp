/**
 * @file matrix.ipp
 * @author Nishinaga Rikuto (you@domain.com)
 * @brief implementation of matrix class
 * @version 0.1
 * @date 2025-04-24
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <complex>
#include "../LU_decomposition.hpp"

 template <typename TYPE, size_t ROW, size_t COL>
 class Matrix;

template <typename TYPE, size_t ROW, size_t COL>
inline TYPE Matrix<TYPE, ROW, COL>::get(size_t row, size_t col) const
{
    if constexpr (ROW == 0 && COL == 0) 
    {
        if(row >= this->m_data_vector.size() || col >= this->m_data_vector[row].size()) 
        {
            throw ::std::out_of_range("Index out of range");
        }
        return this->m_data_vector[row][col];
    } 
    else 
    {
        if(row >= ROW || col >= COL) 
        {
            throw ::std::out_of_range("Index out of range");
        }
        return this->m_data[row][col];
    }
}

template <typename TYPE, size_t ROW, size_t COL>
inline TYPE& Matrix<TYPE, ROW, COL>::set(size_t row, size_t col) 
{
    if constexpr (ROW == 0 && COL == 0) 
    {
        if(row >= this->m_data_vector.size() || col >= this->m_data_vector[row].size()) 
        {
            throw ::std::out_of_range("Index out of range");
        }
        return this->m_data_vector[row][col];
    } 
    else 
    {
        if(row >= ROW || col >= COL) 
        {
            throw ::std::out_of_range("Index out of range");
        }
        return this->m_data[row][col];
    }
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL>::Matrix(const TYPE& value) 
{
    if constexpr (ROW == 0 && COL == 0) 
    {
        this->m_data_vector.resize(ROW, ::std::vector<TYPE>(COL, value));
    } 
    else 
    {
        for (size_t i = 0; i < ROW; ++i) 
        {
            for (size_t j = 0; j < COL; ++j) 
            {
                this->set(i, j) = value;
            }
        }
    }
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            this->set(i, j) = value;
        }
    }
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL>::Matrix(const TYPE (&arr)[ROW][COL]) 
{
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            this->get(i, j) = arr[i][j];
        }
    }
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL>::Matrix(const size_t row, const size_t col)
{
    static_assert(ROW == 0 && COL == 0, "This constructor is allowed only allocated matrix");

    this->m_data_vector.resize(row);

    for(auto& v: m_data_vector)
    {
        v.resize(col);
    }
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL>& Matrix<TYPE, ROW, COL>::operator = (const Matrix<TYPE, ROW, COL>& other) 
{
    if (this != &other) {
        for (size_t i = 0; i < ROW; ++i) 
        {
            for (size_t j = 0; j < COL; ++j) 
            {
                this->set(i, j) = other(i, j);
            }
        }
    }
    return *this;
}

template <typename TYPE, size_t ROW, size_t COL>
TYPE& Matrix<TYPE, ROW, COL>::operator () (size_t row, size_t col) 
{
    return this->set(row, col);
}

template <typename TYPE, size_t ROW, size_t COL>
inline size_t Matrix<TYPE, ROW, COL>::rows() const noexcept 
{
    if constexpr(    ROW == 0 && COL == 0) 
    {
        return this->m_data_vector.size();
    } 
    else 
    {
        return ROW;
    }
}

template <typename TYPE, size_t ROW, size_t COL>
inline size_t Matrix<TYPE, ROW, COL>::cols() const noexcept 
{
    if constexpr(ROW == 0 && COL == 0) 
    {
        return this->m_data_vector[0].size();
    } 
    else 
    {
        return COL;
    }
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator + (const Matrix<TYPE, ROW, COL>& other) const 
{
    Matrix<TYPE, ROW, COL> result;
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            result(i, j) = this->get(i, j) + other.get(i, j);
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator - (const Matrix<TYPE, ROW, COL>& other) const 
{
    Matrix<TYPE, ROW, COL> result;
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            result(i, j) = this->get(i, j) - other.get(i, j);
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator * (const TYPE& scalar) const 
{
    Matrix<TYPE, ROW, COL> result;
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            result(i, j) = this->get(i, j) * scalar;
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator / (const TYPE& scalar) const 
{
    Matrix<TYPE, ROW, COL> result;
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            result(i, j) = this->get(i, j) / scalar;
        }
    }
    return result;
}

// TODO:  implement matrix multiplication by Hashimoto
template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator*(const Matrix<TYPE, COL, ROW>& other) const 
{
    Matrix<TYPE, ROW, ROW> result;
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < ROW; ++j) 
        {
            result(i, j) = 0;
            for (size_t k = 0; k < COL; ++k) 
            {
                result(i, j) += this->get(i, k) * other.get(k, j);
            }
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
bool Matrix<TYPE, ROW, COL>::isZero() const 
{
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            if (::std::fabs(this->get(i, j)) > 1e-8) 
            {
                return false;
            }
        }
    }
    return true;
}

template <typename TYPE, size_t ROW, size_t COL>
bool Matrix<TYPE, ROW, COL>::isIdentity() const 
{
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            if (i == j && ::std::fabs(this->get(i, j) - TYPE(1)) > 1e-8) 
            {
                return false;
            } 
            else if (i != j && ::std::fabs(this->get(i, j)) > 1e-8) 
            {
                return false;
            }
        }
    }
    return true;
}

template <typename TYPE, size_t ROW, size_t COL>
bool Matrix<TYPE, ROW, COL>::isSymmetric() const 
{
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            if (::std::fabs(this->get(i, j) - this->get(j, i)) > 1e-8) 
            {
                return false;
            }
        }
    }
    return true;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::Identity() noexcept
{
    Matrix<TYPE, ROW, COL> result;
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            result(i, j) = (i == j) ? TYPE(1) : TYPE(0);
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::Zero() noexcept
{
    return Matrix<TYPE, ROW, COL>(TYPE(0));
}

template <typename TYPE, size_t ROW, size_t COL>
TYPE Matrix<TYPE, ROW, COL>::determinant() const 
{
    TYPE ret = calculate_determinant(*this);
    return ret;
}

template <typename TYPE, size_t ROW, size_t COL>
TYPE Matrix<TYPE, ROW, COL>::trace() const 
{
    TYPE sum = TYPE(0);
    for (size_t i = 0; i < ::std::min(ROW, COL); ++i) 
    {
        sum += this->get(i, i);
    }
    return sum;
}

template <typename TYPE, size_t ROW, size_t COL>
TYPE Matrix<TYPE, ROW, COL>::norm() const 
{
    TYPE sum = TYPE(0);
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            sum += this->get(i, j) * this->get(i, j);
        }
    }
    return ::std::sqrt(sum);
}

template <typename TYPE, size_t ROW, size_t COL>
TYPE Matrix<TYPE, ROW, COL>::squareNorm() const 
{
    return this->norm();
}


template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::inverse() const 
{
    // Implement matrix inversion here
    // For simplicity, we will just return an identity matrix for now
    return Matrix<TYPE, ROW, COL>::Identity();
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::transpose() const 
{
    static_assert(ROW == COL, "Transpose is only defined for square matrices");

    Matrix<TYPE, ROW, COL> result;
    for (size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            result(j, i) = this->get(i, j);
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::adjoint() const 
{
    auto t_mt = this->transpose();
    Matrix<TYPE, ROW, COL> result;

    for(size_t i = 0; i < ROW; ++i) 
    {
        for (size_t j = 0; j < COL; ++j) 
        {
            result(i, j) = ::std::conj(t_mt(i, j));
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::block(size_t row, size_t col, size_t rows, size_t cols) const 
{
    Matrix<TYPE, 0, 0> result;
    for (size_t i = 0; i < rows; ++i) 
    {
        for (size_t j = 0; j < cols; ++j) 
        {
            result(i, j) = this->get(row + i, col + j);
        }
    }
    return result;
}



