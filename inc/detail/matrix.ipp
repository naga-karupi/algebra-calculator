/**
 * @file matrix.ipp
 * @author Nishinaga Rikuto
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
            throw ::std::out_of_range("Index out of range 1");
        }
        return this->m_data_vector[row][col];
    } 
    else 
    {
        if(row >= ROW || col >= COL) 
        {
            throw ::std::out_of_range("Index out of range 2");
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
            throw ::std::out_of_range("Index out of range 3");
        }
        return this->m_data_vector[row][col];
    } 
    else 
    {
        if(row >= ROW || col >= COL) 
        {
            throw ::std::out_of_range("Index out of range 4");
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

template<typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL>::Matrix(const Matrix<TYPE, ROW, COL>&& other)
{
    if constexpr (ROW == 0 && COL == 0) 
    {
        this->m_data_vector = std::move(other.m_data_vector);
    } 
    else 
    {
        for (size_t i = 0; i < ROW; ++i) 
        {
            for (size_t j = 0; j < COL; ++j) 
            {
                this->set(i, j) = other.get(i, j);
            }
        }
    }
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL>& Matrix<TYPE, ROW, COL>::operator = (const Matrix<TYPE, ROW, COL>& other) 
{
    if constexpr (ROW == 0 && COL == 0) 
    {
        this->m_data_vector = other.m_data_vector;
    } 
    else 
    {
        for (size_t i = 0; i < ROW; ++i) 
        {
            for (size_t j = 0; j < COL; ++j) 
            {
                this->set(i, j) = other.get(i, j);
            }
        }
    }

    return *this;
}

template <typename TYPE, size_t ROW, size_t COL>
template <size_t ROW2, size_t COL2>
Matrix<TYPE, 0, 0>& Matrix<TYPE, ROW, COL>::operator = (const Matrix<TYPE, ROW2, COL2>& other)
{
    static_assert(ROW == 0 && COL == 0, "This operator is allowed only allocated matrix");

    if (this->rows() != other.rows() || this->cols() != other.cols()) 
    {
        throw ::std::invalid_argument("Matrix dimensions do not match");
    }
    this->m_data_vector.resize(other.rows());
    for(auto& v: m_data_vector)
    {
        v.resize(other.cols());
    }
    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < this->cols(); ++j) 
        {
            this->set(i, j) = other.get(i, j);
        }
    }

    return *this;
}


template <typename TYPE, size_t ROW, size_t COL>
TYPE& Matrix<TYPE, ROW, COL>::operator () (size_t row, size_t col) 
{
    if (row >= this->rows() || col >= this->cols()) 
    {
        throw ::std::out_of_range("Index out of range 5");
    }
    return this->set(row, col);
}

template <typename TYPE, size_t ROW, size_t COL>
TYPE& Matrix<TYPE, ROW, COL>::operator () (int row, int col) 
{
    if (size_t(row) >= this->rows() || size_t(col) >= this->cols()) 
    {
        throw ::std::out_of_range("Index out of range 6");
    }
    return this->set(row, col);
}

template <typename TYPE, size_t ROW, size_t COL>
TYPE& Matrix<TYPE, ROW, COL>::operator () (size_t n)
{
    if(this->cols() != 1)
    {
        throw ::std::invalid_argument("Matrix is not a vector");
    }
    if (n < 0 || n >= this->rows()) 
    {
        throw ::std::out_of_range("Index out of range 7");
    }
    return this->set(n, 0);
}

template <typename TYPE, size_t ROW, size_t COL>
TYPE& Matrix<TYPE, ROW, COL>::operator () (int n)
{
    if(this->cols() != 1)
    {
        throw ::std::invalid_argument("Matrix is not a vector");
    }
    if (n < 0 || n >= this->rows()) 
    {
        throw ::std::out_of_range("Index out of range 8");
    }
    return this->set(n, 0);
}

template <typename TYPE, size_t ROW, size_t COL>
inline size_t Matrix<TYPE, ROW, COL>::rows() const noexcept 
{
    if constexpr(ROW == 0 && COL == 0) 
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
    if (this->rows() == 0) 
    {
        throw ::std::invalid_argument("Matrix is not allocated");
    }

    if constexpr(ROW == 0 && COL == 0) 
    {
        return this->m_data_vector[0].size();
    } 
    else 
    {
        return COL;
    }
}

// template <typename TYPE, size_t ROW, size_t COL>
// Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator + (const Matrix<TYPE, ROW, COL>& other) const 
// {
//     if(this->rows() != other.rows() || this->cols() != other.cols()) 
//     {
//         throw ::std::invalid_argument("Matrix dimensions do not match");
//     }
//     Matrix<TYPE, ROW, COL> result;
//     if constexpr(ROW == 0 && COL == 0) 
//     {
//         result.resize(this->rows(), this->cols());
//     }

//     for (size_t i = 0; i < this->rows(); ++i) 
//     {
//         for (size_t j = 0; j < this->cols(); ++j) 
//         {
//             result(i, j) = this->get(i, j) + other.get(i, j);
//         }
//     }
//     return result;
// }

// template <typename TYPE, size_t ROW, size_t COL>
// Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator + (const Matrix<TYPE, 0, 0>& other) const
// {
//     if(this->rows() != other.rows() || this->cols() != other.cols()) 
//     {
//         throw ::std::invalid_argument("Matrix dimensions do not match");
//     }
//     Matrix<TYPE, ROW, COL> result;
//     if constexpr(ROW == 0 && COL == 0) 
//     {
//         result.resize(this->rows(), this->cols());
//     }

//     for (size_t i = 0; i < this->rows(); ++i) 
//     {
//         for (size_t j = 0; j < this->cols(); ++j) 
//         {
//             result(i, j) = this->get(i, j) + other.get(i, j);
//         }
//     }
//     return result;
// }

template <typename TYPE, size_t ROW, size_t COL>
template <size_t ROW2, size_t COL2>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator + (const Matrix<TYPE, ROW2, COL2>& other) const
{
    if(this->rows() != other.rows() || this->cols() != other.cols()) 
    {
        throw ::std::invalid_argument("Matrix dimensions do not match");
    }
    Matrix<TYPE, ROW, COL> result;
    if constexpr(ROW == 0 && COL == 0) 
    {
        result.resize(this->rows(), this->cols());
    }

    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < this->cols(); ++j) 
        {
            result(i, j) = this->get(i, j) + other.get(i, j);
        }
    }
    return result;
}

// template <typename TYPE, size_t ROW, size_t COL>
// Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator - (const Matrix<TYPE, ROW, COL>& other) const 
// {
//     if(this->rows() != other.rows() || this->cols() != other.cols()) 
//     {
//         throw ::std::invalid_argument("Matrix dimensions do not match");
//     }
//     Matrix<TYPE, ROW, COL> result;
//     for (size_t i = 0; i < this->rows(); ++i) 
//     {
//         for (size_t j = 0; j < this->cols(); ++j) 
//         {
//             result(i, j) = this->get(i, j) - other.get(i, j);
//         }
//     }
//     return result;
// }

// template <typename TYPE, size_t ROW, size_t COL>
// Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator - (const Matrix<TYPE, 0, 0>& other) const
// {
//     Matrix<TYPE, ROW, COL> result;
//     for (size_t i = 0; i < this->rows(); ++i) 
//     {
//         for (size_t j = 0; j < this->cols(); ++j) 
//         {
//             result(i, j) = this->get(i, j) - other.get(i, j);
//         }
//     }
//     return result;
// }

template <typename TYPE, size_t ROW, size_t COL>
template <size_t ROW2, size_t COL2>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator - (const Matrix<TYPE, ROW2, COL2>& other) const
{
    if(this->rows() != other.rows() || this->cols() != other.cols()) 
    {
        throw ::std::invalid_argument("Matrix dimensions do not match");
    }
    Matrix<TYPE, ROW, COL> result;
    if constexpr(ROW == 0 && COL == 0) 
    {
        result.resize(this->rows(), this->cols());
    }
    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < this->cols(); ++j) 
        {
            result.set(i, j) = this->get(i, j) - other.get(i, j);
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator * (const TYPE& scalar) const 
{
    Matrix<TYPE, ROW, COL> result;
    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < cols(); ++j) 
        {
            result.set(i, j) = this->get(i, j) * scalar;
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::operator / (const TYPE& scalar) const 
{
    Matrix<TYPE, ROW, COL> result;

    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < this->cols(); ++j) 
        {
            result.set(i, j) = this->get(i, j) / scalar;
        }
    }
    return result;
}

// TODO:  implement matrix multiplication by Hashimoto
template <typename TYPE, size_t ROW, size_t COL,size_t ROW2, size_t COL2>
typename ResultMatrixType<TYPE, ROW, COL, ROW2, COL2>::type operator * (const Matrix<TYPE, ROW, COL>& lhs, const Matrix<TYPE, ROW2, COL2>& rhs)
{
    if (lhs.cols() != rhs.rows())
    {
        throw ::std::invalid_argument("Matrix dimensions do not match");
    }

    using ret_mat = ResultMatrixType<TYPE, ROW, COL, ROW2, COL2>;

    Matrix<typename ret_mat::type, ret_mat::ret_row, ret_mat::ret_col> result;
    if constexpr (ret_mat::ret_row == 0 && ret_mat::ret_col  == 0) 
    {
        result.resize(lhs.rows(), rhs.cols());
    }

    for (size_t i = 0; i < lhs.rows(); ++i) 
    {
        for (size_t j = 0; j < rhs.cols(); ++j) 
        {
            result(i, j) = TYPE{};
            for (size_t k = 0; k < lhs.cols(); ++k) 
            {
                result(i, j) += lhs.get(i, k) * rhs.get(k, j);
            }
        }
    }

    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
inline bool Matrix<TYPE, ROW, COL>::isZero() const 
{
    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < this->cols(); ++j) 
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
    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < this->cols(); ++j) 
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
    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < this->cols(); ++j) 
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
Matrix<TYPE, 0, 0> Matrix<TYPE, ROW, COL>::Identity(size_t row, size_t col) noexcept
{
    Matrix<TYPE, 0, 0> result;
    result.resize(row, col);
    if (row != col) 
    {
        throw ::std::invalid_argument("Matrix must be square");
    }
    for (size_t i = 0; i < row; ++i) 
    {
        for (size_t j = 0; j < col; ++j) 
        {
            result.set(i, j) = (i == j) ? TYPE(1) : TYPE(0);
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
Matrix<TYPE, 0, 0> Matrix<TYPE, ROW, COL>::Zero(size_t row, size_t col) noexcept
{
    Matrix<TYPE, 0, 0> result;
    result.resize(row, col);
    for (size_t i = 0; i < row; ++i) 
    {
        for (size_t j = 0; j < col; ++j) 
        {
            result(i, j) = TYPE(0);
        }
    }
    return result;
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
    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < this->cols(); ++j) 
        {
            sum += this->get(i, j) * this->get(i, j);
        }
    }
    return ::std::sqrt(sum);
}

template <typename TYPE, size_t ROW, size_t COL>
TYPE Matrix<TYPE, ROW, COL>::squaredNorm() const 
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
    Matrix<TYPE, COL, ROW> result;
    if constexpr(ROW == 0 && COL == 0) 
    {
        result.resize(this->cols(), this->rows());
    }
    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < this->cols(); ++j) 
        {
            result.set(j, i) = this->get(i, j);
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> Matrix<TYPE, ROW, COL>::adjoint() const 
{
    auto t_mt = this->transpose();
    Matrix<TYPE, COL, ROW> result;
    if constexpr(ROW == 0 && COL == 0) 
    {
        result.resize(this->cols(), this->rows());
    }

    for(size_t i = 0; i < t_mt.rows(); ++i) 
    {
        for (size_t j = 0; j < t_mt.cols(); ++j) 
        {
            result.set(i, j) = ::std::conj(t_mt.set(i, j));
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, 0, 0> Matrix<TYPE, ROW, COL>::block(size_t row, size_t col, size_t rows, size_t cols) const 
{
    Matrix<TYPE, 0, 0> result;
    if (row + rows > this->rows() || col + cols > this->cols()) 
    {
        throw ::std::out_of_range("Block size exceeds matrix dimensions");
    }
    result.resize(rows, cols);
    for (size_t i = 0; i < rows; ++i) 
    {
        for (size_t j = 0; j < cols; ++j) 
        {
            result(i, j) = this->get(row + i, col + j);
        }
    }
    return result;
}

template <typename TYPE, size_t ROW, size_t COL>
void Matrix<TYPE, ROW, COL>::resize(size_t row, size_t col) 
{
    static_assert(ROW == 0 && COL == 0, "This function is only available for dynamic matrix");

    this->m_data_vector.resize(row);
    for(auto& v: m_data_vector)
    {
        v.resize(col);
    }
}

template <typename TYPE, size_t ROW, size_t COL>
Matrix<TYPE, ROW, COL> operator * (const TYPE& scalar, const Matrix<TYPE, ROW, COL>& matrix) 
{
    return matrix * scalar;
}

template <typename TYPE, size_t ROW, size_t COL>
void Matrix<TYPE, ROW, COL>::fill(const TYPE& value) 
{
    for (size_t i = 0; i < this->rows(); ++i) 
    {
        for (size_t j = 0; j < this->cols(); ++j) 
        {
            this->set(i, j) = value;
        }
    }
}
template <typename TYPE, size_t ROW, size_t COL>
void print(const Matrix<TYPE, ROW, COL>& matrix) 
{
    for (size_t i = 0; i < matrix.rows(); ++i) 
    {
        for (size_t j = 0; j < matrix.cols(); ++j) 
        {
            std::cout << matrix.get(i, j) << " ";
        }
        std::cout << "\n";
    }
}

