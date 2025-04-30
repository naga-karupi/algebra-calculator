/**
 * @file LU_decomposition.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2025-04-28
 * 
 * @copyright Copyright (c) 2025
 * 
 */

 #pragma once
 
#include <complex>
#include <tuple>
#include "matrix.hpp"

template <typename TYPE, size_t COLUMNS, size_t ROWS>
::std::pair<::my_mt::Matrix<double, COLUMNS, ROWS>, ::my_mt::Matrix<double, COLUMNS, ROWS>> 
LU_decompose(const my_mt::Matrix<double, COLUMNS, ROWS>& A)
{
    static_assert(ROWS == COLUMNS, "Matrix must be square for LU decomposition");

    // Initialize L and U matrices
    // L is lower triangular, U is upper triangular
    my_mt::Matrix<TYPE, COLUMNS, ROWS> L;
    my_mt::Matrix<TYPE, COLUMNS, ROWS> U;

    for(size_t j = 0; j < COLUMNS; j++)
    {
        U(0, j) = A(0, j);
    }

    for(size_t i = 0; i < ROWS; i++)
    {
        L(i, 0) = A(i, 0) / U(0, 0);
    }

    // LU decomposition algorithm
    // Inner product Gaussian method
    for(size_t i = 0; i < ROWS; i++)
    {
        for (size_t j = 0; j < COLUMNS; j++)
        {
            bool pass_L = false;
            bool pass_U = false;

            if(i <= j)
            {
                if (i == j)
                    L(i, j) = 1;
                else
                    L(i, j) = 0;
                
                pass_L = true;
            }
            else 
            {
                L(i, j) = A(i, j);
            }

            if (i > j)
            {
                U(i, j) = 0;
                pass_U = true;
            }
            else
            {
                U(i, j) = A(i, j);
            }

            for(size_t k = 0; k < i && k < j; k++)
            {
                if (not pass_U)
                {
                    U(i, j) -= L(i, k) * U(k, j);
                }

                if (not pass_L)
                {
                    L(i, j) -= L(i, k) * U(k, j);
                }
            }

            L(i, j) /= U(j, j);

        }
    }
    
    return {L, U};
}

template <typename TYPE, size_t COLUMNS, size_t ROWS>
TYPE calculate_determinant(const my_mt::Matrix<TYPE, ROWS, COLUMNS>& matrix)
{
    TYPE det = 1;
    auto L_and_U = LU_decompose<ROWS, COLUMNS>(matrix);

    auto L = L_and_U[0];
    auto U = L_and_U[1];

    for (size_t i = 0; i < COLUMNS; i++)
    {
        det *= U(i, i);
        // std::cout << "U(i, i) = " << U(i, i) << std::endl;
    }

    return det;
}
