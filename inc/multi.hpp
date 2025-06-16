/**
 * @file multi.hpp
 * @author Yuri Hashimoto
 * @brief implementation of multi matrix
 * @version 0.1
 * @date 2025-05-07
 * 
 * @copyright Copyright (c) 2025
 * 
 */

 #pragma once
#include "matrix.hpp"

 template<typename T, size_t ROW, size_t COL>
::my_mt::Matrix<T, ROW, COL> add(const ::my_mt::Matrix<T, ROW, COL>& A, const ::my_mt::Matrix<T, ROW, COL>& B) {
    ::my_mt::Matrix<T, ROW, COL> C;
    for (size_t i = 0; i < ROW; ++i)
        for (size_t j = 0; j < COL; ++j)
            C(i, j) = A.get(i, j) + B.get(i, j);
    return C;
}

template<typename T, size_t ROW, size_t COL>
::my_mt::Matrix<T, ROW, COL> subtract(const ::my_mt::Matrix<T, ROW, COL>& A, const ::my_mt::Matrix<T, ROW, COL>& B) {
    ::my_mt::Matrix<T, ROW, COL> C;
    for (size_t i = 0; i < ROW; ++i)
        for (size_t j = 0; j < COL; ++j)
            C(i, j) = A.get(i, j) - B.get(i, j);
    return C;
}

constexpr size_t nextPowerOfTwo(size_t n) {
    size_t power = 1;
    while (power < n) power <<= 1;
    return power;
}

template<typename T, size_t SIZE>
::my_mt::Matrix<T, SIZE, SIZE> strassen(const ::my_mt::Matrix<T, SIZE, SIZE>& A, const ::my_mt::Matrix<T, SIZE, SIZE>& B) {
    if constexpr (SIZE == 1) {
        ::my_mt::Matrix<T, 1, 1> C;
        C(0, 0) = A.get(0, 0) * B.get(0, 0);
        return C;
    } else {
        constexpr size_t K = SIZE / 2;
        using Mat = ::my_mt::Matrix<T, K, K>;

        Mat A11, A12, A21, A22;
        Mat B11, B12, B21, B22;

        for (size_t i = 0; i < K; ++i) {
            for (size_t j = 0; j < K; ++j) {
                A11(i,j) = A.get(i,j);             A12(i,j) = A.get(i,j+K);
                A21(i,j) = A.get(i+K,j);           A22(i,j) = A.get(i+K,j+K);
                B11(i,j) = B.get(i,j);             B12(i,j) = B.get(i,j+K);
                B21(i,j) = B.get(i+K,j);           B22(i,j) = B.get(i+K,j+K);
            }
        }

        auto M1 = strassen(add(A11, A22), add(B11, B22));
        auto M2 = strassen(add(A21, A22), B11);
        auto M3 = strassen(A11, subtract(B12, B22));
        auto M4 = strassen(A22, subtract(B21, B11));
        auto M5 = strassen(add(A11, A12), B22);
        auto M6 = strassen(subtract(A21, A11), add(B11, B12));
        auto M7 = strassen(subtract(A12, A22), add(B21, B22));

        Mat C11 = add(subtract(add(M1, M4), M5), M7);
        Mat C12 = add(M3, M5);
        Mat C21 = add(M2, M4);
        Mat C22 = add(subtract(add(M1, M3), M2), M6);

        ::my_mt::Matrix<T, SIZE, SIZE> C;
        for (size_t i = 0; i < K; ++i) {
            for (size_t j = 0; j < K; ++j) {
                C(i,j)     = C11(i,j);
                C(i,j+K)   = C12(i,j);
                C(i+K,j)   = C21(i,j);
                C(i+K,j+K) = C22(i,j);
            }
        }

        return C;
    }
}

// 任意サイズの行列に対応したStrassen乗算
template<typename T, size_t M, size_t K, size_t N>
::my_mt::Matrix<T, M, N> multiply(const ::my_mt::Matrix<T, M, K>& A, const ::my_mt::Matrix<T, K, N>& B) {
    constexpr size_t maxDim = nextPowerOfTwo(std::max({M, K, N}));
    ::my_mt::Matrix<T, maxDim, maxDim> A_pad, B_pad;
    A_pad.fill(0);
    B_pad.fill(0);

    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < K; ++j)
            A_pad(i, j) = A.get(i, j);
    for (size_t i = 0; i < K; ++i)
        for (size_t j = 0; j < N; ++j)
            B_pad(i, j) = B.get(i, j);

    ::my_mt::Matrix<T, maxDim, maxDim> C_pad = strassen(A_pad, B_pad);

    ::my_mt::Matrix<T, M, N> C;
    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < N; ++j)
            C(i, j) = C_pad(i, j);

    return C;
}