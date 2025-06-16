/**
 * @file helper_class.hpp
 * @author your name (you@domain.com)
 * @brief helper class to determine the result matrix type of matrix operations
 * @version 0.1
 * @date 2025-05-28
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#pragma once

#include <cstddef> 

template <typename TYPE, size_t ROW, size_t COL>
class Matrix;

/**
 * @brief helper struct to determine the result matrix type of matrix operations
 * 
 * @tparam TYPE 
 * @tparam R1 
 * @tparam C1 
 * @tparam R2 
 * @tparam C2 
 */
template <typename TYPE, size_t R1, size_t C1, size_t R2, size_t C2>
struct ResultMatrixType {
    static constexpr size_t res_row = (R1 == 0 || C1 == 0 || R2 == 0 || C2 == 0) ? 0 : R1;
    static constexpr size_t res_col = (R1 == 0 || C1 == 0 || R2 == 0 || C2 == 0) ? 0 : C2;

    using type = Matrix<TYPE, res_row, res_col>;
};
