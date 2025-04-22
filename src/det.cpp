/**
 * @file det.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2025-04-16
 * 
 * @copyright Copyright (c) 2025
 * 
 */

/// USE LU DECOMPOSITION
#include <iostream>
#include <array>
#include <vector>
#include <random>

#include <Eigen/Dense>

void show_22matrix(Eigen::Matrix<double, 2, 2>mat)
{
    for(auto i = 0lu; i < 2; i++){
        for (auto j = 0lu; j < 2; j++)
        {
            std::cout << mat(i, j) << ' ';
        }
        std::cout << "\n";
    }
}

template <size_t COLUMNS, size_t ROWS>
std::array<Eigen::Matrix<double, COLUMNS, ROWS>, 2> LU_decompose(const Eigen::Matrix<double, COLUMNS, ROWS> A)
{
    static_assert(ROWS == COLUMNS, "Matrix must be square for LU decomposition");

    // Initialize L and U matrices
    // L is lower triangular, U is upper triangular
    Eigen::Matrix<double, COLUMNS, ROWS> L;
    Eigen::Matrix<double, COLUMNS, ROWS> U;

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

template <size_t COLUMNS, size_t ROWS>
double calculate_determinant(const Eigen::Matrix<double, ROWS, COLUMNS>& matrix)
{


    double det = 1;
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



int main()
{
    constexpr size_t N = 100;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);
    
    Eigen::Matrix<double, N, N> mat;
    for(auto i = 0; i < N; i++)
    {
        for(auto j = 0; j < N; j++)
        {
            mat(i, j) = dis(gen);
        }
    }

    // calc determine by eigen lib
    auto det_eigen = mat.determinant();
    std::cout << "det_eigen = " << det_eigen << std::endl;

    auto LU = LU_decompose<N, N>(mat);

    // std::cout << "mat" << std::endl;
    // show_22matrix(mat);
    // std::cout << "L" << std::endl;
    // show_22matrix(LU[0]);
    // std::cout << "U" << std::endl;
    // show_22matrix(LU[1]);

    auto det = calculate_determinant<N, N>(mat);

    std::cout << "det = " << det << std::endl;

    return 0;
}

