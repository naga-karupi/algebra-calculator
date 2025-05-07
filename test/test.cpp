/**
 * @file test.cpp
 * @author your name (you@domain.com)
 * @brief test file for matrix class
 * @version 0.1
 * @date 2025-04-30
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "../inc/matrix.hpp"

#include <random>
#include <complex>

#include <gtest/gtest.h>

TEST(MatrixTest, MatrixConstructor)
{
    my_mt::Matrix<double, 3, 3> mat1;
    EXPECT_EQ(mat1.rows(), 3);
    EXPECT_EQ(mat1.cols(), 3);

    my_mt::Matrix<double, 0, 0> mat2(3, 4);
    EXPECT_EQ(mat2.rows(), 3);
    EXPECT_EQ(mat2.cols(), 4);

    my_mt::Matrix<double, 2, 2> mat3(5.0);
    EXPECT_EQ(mat3(0, 0), 5.0);
    EXPECT_EQ(mat3(1, 1), 5.0);
}

TEST(MatrixTest, MatrixAddition)
{
    my_mt::Matrix<double, 2, 2> mat1;
    mat1(0, 0) = 1.0;
    mat1(0, 1) = 2.0;
    mat1(1, 0) = 3.0;
    mat1(1, 1) = 4.0;

    my_mt::Matrix<double, 2, 2> mat2;
    mat2(0, 0) = 5.0;
    mat2(0, 1) = 6.0;
    mat2(1, 0) = 7.0;
    mat2(1, 1) = 8.0;

    auto result = mat1 + mat2;

    EXPECT_EQ(result(0, 0), 6.0);
    EXPECT_EQ(result(0, 1), 8.0);
    EXPECT_EQ(result(1, 0), 10.0);
    EXPECT_EQ(result(1, 1), 12.0);
}
TEST(MatrixTest, MatrixSubtraction)
{
    my_mt::Matrix<double, 2, 2> mat1;
    mat1(0, 0) = 5.0;
    mat1(0, 1) = 6.0;
    mat1(1, 0) = 7.0;
    mat1(1, 1) = 8.0;

    my_mt::Matrix<double, 2, 2> mat2;
    mat2(0, 0) = 1.0;
    mat2(0, 1) = 2.0;
    mat2(1, 0) = 3.0;
    mat2(1, 1) = 4.0;

    auto result = mat1 - mat2;

    EXPECT_EQ(result(0, 0), 4.0);
    EXPECT_EQ(result(0, 1), 4.0);
    EXPECT_EQ(result(1, 0), 4.0);
    EXPECT_EQ(result(1, 1), 4.0);
}

TEST(MatrixTest, MatrixScalarMultiplication)
{
    my_mt::Matrix<double, 2, 2> mat;
    mat(0, 0) = 1.0;
    mat(0, 1) = 2.0;
    mat(1, 0) = 3.0;
    mat(1, 1) = 4.0;

    auto result = mat * 2.0;

    EXPECT_EQ(result(0, 0), 2.0);
    EXPECT_EQ(result(0, 1), 4.0);
    EXPECT_EQ(result(1, 0), 6.0);
    EXPECT_EQ(result(1, 1), 8.0);
}

TEST(MatrixTest, MatrixScalarDivision)
{
    my_mt::Matrix<double, 2, 2> mat;
    mat(0, 0) = 2.0;
    mat(0, 1) = 4.0;
    mat(1, 0) = 6.0;
    mat(1, 1) = 8.0;

    auto result = mat / 2.0;

    EXPECT_EQ(result(0, 0), 1.0);
    EXPECT_EQ(result(0, 1), 2.0);
    EXPECT_EQ(result(1, 0), 3.0);
    EXPECT_EQ(result(1, 1), 4.0);
}

TEST(MatrixTest, MatrixMultiplication)
{
    my_mt::Matrix<double, 2, 2> mat1;
    mat1(0, 0) = 1.0;
    mat1(0, 1) = 2.0;
    mat1(1, 0) = 3.0;
    mat1(1, 1) = 4.0;

    my_mt::Matrix<double, 2, 2> mat2;
    mat2(0, 0) = 5.0;
    mat2(0, 1) = 6.0;
    mat2(1, 0) = 7.0;
    mat2(1, 1) = 8.0;

    auto result = mat1 * mat2;

    EXPECT_EQ(result(0, 0), 19.0);
    EXPECT_EQ(result(0, 1), 22.0);
    EXPECT_EQ(result(1, 0), 43.0);
    EXPECT_EQ(result(1, 1), 50.0);
}

TEST(MatrixTest, MatrixZero)
{
    my_mt::Matrix<double, 2, 2> zeroMat = my_mt::Matrix<double, 2, 2>::Zero();

    EXPECT_TRUE(zeroMat.isZero());
    EXPECT_FALSE(zeroMat.isIdentity());
}

TEST(MatrixTest, MatrixIdentity)
{
    my_mt::Matrix<double, 2, 2> identityMat = my_mt::Matrix<double, 2, 2>::Identity();

    EXPECT_FALSE(identityMat.isZero());
    EXPECT_TRUE(identityMat.isIdentity());
}

TEST(MatrixTest, MatrixSymmetric)
{
    my_mt::Matrix<double, 2, 2> symMat;
    symMat(0, 0) = 1.0;
    symMat(0, 1) = 2.0;
    symMat(1, 0) = 2.0;
    symMat(1, 1) = 3.0;

    EXPECT_TRUE(symMat.isSymmetric());

    my_mt::Matrix<double, 2, 2> nonSymMat;
    nonSymMat(0, 0) = 1.0;
    nonSymMat(0, 1) = 2.0;
    nonSymMat(1, 0) = 3.0;
    nonSymMat(1, 1) = 4.0;

    EXPECT_FALSE(nonSymMat.isSymmetric());
}

TEST(MatrixTest, MatrixDeterminant)
{
    my_mt::Matrix<double, 2, 2> mat;
    mat(0, 0) = 1.0;
    mat(0, 1) = 2.0;
    mat(1, 0) = 3.0;
    mat(1, 1) = 4.0;

    EXPECT_LE(std::fabs(mat.determinant() - (-2.0)), 1e-8);
}

TEST(MatrixTest, MatrixTrace)
{
    my_mt::Matrix<double, 2, 2> mat;
    mat(0, 0) = 1.0;
    mat(0, 1) = 2.0;
    mat(1, 0) = 3.0;
    mat(1, 1) = 4.0;

    EXPECT_EQ(mat.trace(), 5.0);
}

TEST(MatrixTest, MatrixNorm)
{
    my_mt::Matrix<double, 2, 2> mat;
    mat(0, 0) = 1.0;
    mat(0, 1) = 2.0;
    mat(1, 0) = 3.0;
    mat(1, 1) = 4.0;

    EXPECT_DOUBLE_EQ(mat.norm(), ::std::sqrt(30.0));
}

TEST(MatrixTest, MatrixTranspose)
{
    my_mt::Matrix<double, 2, 2> mat;
    mat(0, 0) = 1.0;
    mat(0, 1) = 2.0;
    mat(1, 0) = 3.0;
    mat(1, 1) = 4.0;

    auto result = mat.transpose();

    EXPECT_EQ(result(0, 0), 1.0);
    EXPECT_EQ(result(0, 1), 3.0);
    EXPECT_EQ(result(1, 0), 2.0);
    EXPECT_EQ(result(1, 1), 4.0);
}

TEST(MatrixTest, MatrixAdjoint)
{
    my_mt::Matrix<std::complex<double>, 2, 2> mat;
    mat(0, 0) = std::complex<double>(1.0, 1.0);
    mat(0, 1) = std::complex<double>(2.0, 2.0);
    mat(1, 0) = std::complex<double>(3.0, 3.0);
    mat(1, 1) = std::complex<double>(4.0, 4.0);
    auto result = mat.adjoint();
    EXPECT_EQ(result(0, 0), std::complex<double>(1.0, -1.0));
    EXPECT_EQ(result(0, 1), std::complex<double>(3.0, -3.0));
    EXPECT_EQ(result(1, 0), std::complex<double>(2.0, -2.0));
    EXPECT_EQ(result(1, 1), std::complex<double>(4.0, -4.0));
}
