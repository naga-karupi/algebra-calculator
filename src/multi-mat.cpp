/**
 * @file multi-mat.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2025-04-17
 * 
 * @copyright Copyright (c) 2025
 * 
 */

 
#include <array>
#include <iostream>
#include <algorithm>
#include <cassert>

namespace my_mt {

/*
template<typename T, size_t ROW, size_t COL>
class Matrix {
public:
    std::array<std::array<T, COL>, ROW> data{};

    T& operator()(size_t i, size_t j) { return data[i][j]; }
    const T& operator()(size_t i, size_t j) const { return data[i][j]; }

    void fill(T value) {
        for (auto& row : data)
            row.fill(value);
    }

    void print() const {
        for (size_t i = 0; i < ROW; ++i) {
            for (size_t j = 0; j < COL; ++j)
                std::cout << data[i][j] << " ";
            std::cout << "\n";
        }
    }
};
*/

template<typename T, size_t ROW, size_t COL>
Matrix<T, ROW, COL> add(const Matrix<T, ROW, COL>& A, const Matrix<T, ROW, COL>& B) {
    Matrix<T, ROW, COL> C;
    for (size_t i = 0; i < ROW; ++i)
        for (size_t j = 0; j < COL; ++j)
            C(i, j) = A(i, j) + B(i, j);
    return C;
}

template<typename T, size_t ROW, size_t COL>
Matrix<T, ROW, COL> subtract(const Matrix<T, ROW, COL>& A, const Matrix<T, ROW, COL>& B) {
    Matrix<T, ROW, COL> C;
    for (size_t i = 0; i < ROW; ++i)
        for (size_t j = 0; j < COL; ++j)
            C(i, j) = A(i, j) - B(i, j);
    return C;
}

constexpr size_t nextPowerOfTwo(size_t n) {
    size_t power = 1;
    while (power < n) power <<= 1;
    return power;
}

template<typename T, size_t SIZE>
Matrix<T, SIZE, SIZE> strassen(const Matrix<T, SIZE, SIZE>& A, const Matrix<T, SIZE, SIZE>& B) {
    if constexpr (SIZE == 1) {
        Matrix<T, 1, 1> C;
        C(0, 0) = A(0, 0) * B(0, 0);
        return C;
    } else {
        constexpr size_t K = SIZE / 2;
        using Mat = Matrix<T, K, K>;

        Mat A11, A12, A21, A22;
        Mat B11, B12, B21, B22;

        for (size_t i = 0; i < K; ++i) {
            for (size_t j = 0; j < K; ++j) {
                A11(i,j) = A(i,j);             A12(i,j) = A(i,j+K);
                A21(i,j) = A(i+K,j);           A22(i,j) = A(i+K,j+K);
                B11(i,j) = B(i,j);             B12(i,j) = B(i,j+K);
                B21(i,j) = B(i+K,j);           B22(i,j) = B(i+K,j+K);
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

        Matrix<T, SIZE, SIZE> C;
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
Matrix<T, M, N> multiply(const Matrix<T, M, K>& A, const Matrix<T, K, N>& B) {
    constexpr size_t maxDim = nextPowerOfTwo(std::max({M, K, N}));
    Matrix<T, maxDim, maxDim> A_pad, B_pad;
    A_pad.fill(0);
    B_pad.fill(0);

    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < K; ++j)
            A_pad(i, j) = A(i, j);
    for (size_t i = 0; i < K; ++i)
        for (size_t j = 0; j < N; ++j)
            B_pad(i, j) = B(i, j);

    Matrix<T, maxDim, maxDim> C_pad = strassen(A_pad, B_pad);

    Matrix<T, M, N> C;
    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < N; ++j)
            C(i, j) = C_pad(i, j);

    return C;
}

} // namespace my_mt

/*
// Example usage
int main() {
    using namespace my_mt;

    Matrix<int, 22, 10> A;
    Matrix<int, 10, 8> B;

    A(0,0) = 12; A(0,1) = 20; A(0,2)= 16; A(0,3) = 13; A(0,4) =  5; A(0,5) = 11; A(0,6) = 10; A(0,7) = 20; A(0,8) =  8; A(0,9) =  4;
    A(1,0) =  4; A(1,1) = 11; A(1,2)= 14; A(1,3) = 10; A(1,4) =  6; A(1,5) =  2; A(1,6) =  4; A(1,7) = 18; A(1,8) =  9; A(1,9) = 15;
    A(2,0) =  2; A(2,1) = 20; A(2,2)=  1; A(2,3) = 17; A(2,4) = 14; A(2,5) = 17; A(2,6) = 19; A(2,7) =  7; A(2,8) = 10; A(2,9) =  9;
    A(3,0) =  7; A(3,1) =  6; A(3,2)= 13; A(3,3) =  2; A(3,4) = 20; A(3,5) =  5; A(3,6) = 11; A(3,7) = 18; A(3,8) =  2; A(3,9) = 11;
    A(4,0) = 19; A(4,1) = 19; A(4,2)=  8; A(4,3) = 17; A(4,4) = 12; A(4,5) = 18; A(4,6) = 10; A(4,7) =  9; A(4,8) =  3; A(4,9) = 17;
    A(5,0) =  8; A(5,1) =  3; A(5,2)= 13; A(5,3) = 10; A(5,4) = 21; A(5,5) = 17; A(5,6) =  3; A(5,7) = 18; A(5,8) = 14; A(5,9) =  3;
    A(6,0) =  6; A(6,1) = 13; A(6,2)=  4; A(6,3) =  2; A(6,4) = 15; A(6,5) =  7; A(6,6) =  2; A(6,7) = 12; A(6,8) =  5; A(6,9) =  2;
    A(7,0) = 17; A(7,1) = 15; A(7,2)=  7; A(7,3) =  5; A(7,4) =  4; A(7,5) =  6; A(7,6) = 15; A(7,7) =  6; A(7,8) =  6; A(7,9) = 20;
    A(8,0) =  2; A(8,1) =  3; A(8,2)=  4; A(8,3) = 15; A(8,4) = 10; A(8,5) =  9; A(8,6) =  2; A(8,7) = 13; A(8,8) =  8; A(8,9) = 15;
    A(9,0) = 15; A(9,1) = 15; A(9,2)= 15; A(9,3) = 16; A(9,4) =  8; A(9,5) =  4; A(9,6) =  9; A(9,7) =  7; A(9,8) = 14; A(9,9) =  3;
    A(10,0) = 17; A(10,1) = 12; A(10,2)=  8; A(10,3) =  6; A(10,4) = 13; A(10,5) = 18; A(10,6) =  6; A(10,7) = 14; A(10,8) = 13; A(10,9) =  9;
    A(11,0) =  2; A(11,1) = 20; A(11,2)= 10; A(11,3) =  8; A(11,4) = 17; A(11,5) = 14; A(11,6) = 18; A(11,7) =  7; A(11,8) =  6; A(11,9) =  2;
    A(12,0) = 10; A(12,1) =  6; A(12,2)= 18; A(12,3) = 11; A(12,4) =  2; A(12,5) = 20; A(12,6) = 16; A(12,7) = 21; A(12,8) = 12; A(12,9) = 11;
    A(13,0) = 11; A(13,1) = 18; A(13,2)=  3; A(13,3) = 17; A(13,4) =  5; A(13,5) =  9; A(13,6) = 14; A(13,7) = 11; A(13,8) = 11; A(13,9) =  6;
    A(14,0) =  4; A(14,1) =  7; A(14,2)= 21; A(14,3) =  2; A(14,4) = 21; A(14,5) =  5; A(14,6) = 16; A(14,7) = 13; A(14,8) =  4; A(14,9) = 15;
    A(15,0) =  8; A(15,1) = 14; A(15,2)= 20; A(15,3) = 19; A(15,4) =  7; A(15,5) = 19; A(15,6) =  7; A(15,7) = 14; A(15,8) = 20; A(15,9) =  2;
    A(16,0) = 20; A(16,1) = 11; A(16,2)=  9; A(16,3) =  8; A(16,4) = 11; A(16,5) =  6; A(16,6) =  2; A(16,7) = 10; A(16,8) = 15; A(16,9) = 19;
    A(17,0) =  8; A(17,1) =  8; A(17,2)= 18; A(17,3) = 20; A(17,4) = 17; A(17,5) =  6; A(17,6) = 14; A(17,7) =  9; A(17,8) = 16; A(17,9) = 17;
    A(18,0) = 12; A(18,1) = 19; A(18,2)= 15; A(18,3) = 17; A(18,4) =  9; A(18,5) =  5; A(18,6) =  3; A(18,7) =  5; A(18,8) =  3; A(18,9) =  9;
    A(19,0) =  5; A(19,1) =  6; A(19,2)=  2; A(19,3) = 14; A(19,4) = 20; A(19,5) =  9; A(19,6) =  3; A(19,7) = 20; A(19,8) = 15; A(19,9) =  5;
    A(20,0) =  2; A(20,1) =  5; A(20,2)=  5; A(20,3) = 20; A(20,4) =  5; A(20,5) = 11; A(20,6) = 10; A(20,7) = 16; A(20,8) = 21; A(20,9) = 16;
    A(21,0) =  5; A(21,1) = 14; A(21,2)=  4; A(21,3) = 19; A(21,4) = 19; A(21,5) = 11; A(21,6) = 10; A(21,7) =  5; A(21,8) =  1; A(21,9) =  7;

    B(0,0) =  7; B(0,1) =  1; B(0,2) =  8; B(0,3) = 18; B(0,4) =  5; B(0,5) = 19; B(0,6) =  5; B(0,7) =  8;
    B(1,0) =  1; B(1,1) = 17; B(1,2) =  7; B(1,3) =  7; B(1,4) = 14; B(1,5) = 14; B(1,6) = 15; B(1,7) = 14;
    B(2,0) = 14; B(2,1) = 16; B(2,2) = 17; B(2,3) = 20; B(2,4) =  4; B(2,5) = 21; B(2,6) =  5; B(2,7) =  2;
    B(3,0) =  3; B(3,1) =  3; B(3,2) =  4; B(3,3) =  6; B(3,4) = 17; B(3,5) = 21; B(3,6) = 16; B(3,7) =  2;
    B(4,0) = 21; B(4,1) = 12; B(4,2) =  9; B(4,3) = 15; B(4,4) =  2; B(4,5) = 17; B(4,6) =  4; B(4,7) = 16;
    B(5,0) =  4; B(5,1) =  8; B(5,2) = 12; B(5,3) = 15; B(5,4) = 20; B(5,5) = 17; B(5,6) = 20; B(5,7) =  3;
    B(6,0) =  9; B(6,1) =  4; B(6,2) = 14; B(6,3) =  9; B(6,4) = 12; B(6,5) = 19; B(6,6) = 12; B(6,7) =  5;
    B(7,0) =  1; B(7,1) =  5; B(7,2) = 19; B(7,3) =  6; B(7,4) = 13; B(7,5) = 15; B(7,6) = 15; B(7,7) = 14;
    B(8,0) = 12; B(8,1) = 16; B(8,2) = 19; B(8,3) = 15; B(8,4) = 13; B(8,5) = 16; B(8,6) =  4; B(8,7) = 19;
    B(9,0) =  6; B(9,1) = 21; B(9,2) = 13; B(9,3) = 18; B(9,4) = 16; B(9,5) = 10; B(9,6) = 20; B(9,7) = 17;

   auto C = multiply(A, B);

    C.print(); // Output result

    return 0;
}

*/
